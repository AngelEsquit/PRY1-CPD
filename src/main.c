// Fluid screensaver, sequential version (comparison baseline for the
// OpenMP parallel version). Solves incompressible Navier-Stokes with Jos
// Stam's "Stable Fluids" method: see solver.c for the actual math.
//
// Per-frame flow: move sources (n-body), inject ink/force, solve velocity,
// advect each ink channel through it, fade the ink a little, render.
//
// Build: make (see Makefile). Usage: ./ss -n 128 -f 6

#include <SDL2/SDL.h>
#include <stdio.h>
#include <stdlib.h>

#include "background.h"
#include "common.h"
#include "config.h"
#include "fields.h"
#include "nbody.h"
#include "render.h"
#include "solver.h"
#include "sources.h"

// Grid resolution, global so the IX() macro stays readable. Declared
// extern in common.h so every other module can use it.
int grid_n = GRID_DEFAULT;

int main(int argc, char *argv[]) {
  Config config;
  FluidFields fields;
  InkSource *sources = NULL;
  Star *stars = NULL;
  SDL_Window *window = NULL;
  SDL_Renderer *renderer = NULL;
  SDL_Texture *texture = NULL;
  SDL_Event event;

  int running = 1;
  int exit_code = EXIT_SUCCESS;
  int accumulated_frames = 0;
  Uint64 prev_ticks, current_ticks, clock_frequency;
  Uint64 last_report_ticks;
  double elapsed_seconds, current_fps = 0.0;
  char window_title[128];

  // Parse command-line arguments.
  int args_result = parse_arguments(argc, argv, &config);
  if (args_result == -1)
    return EXIT_SUCCESS; // Help was requested.
  if (args_result == 0)
    return EXIT_FAILURE; // Invalid arguments.

  grid_n = config.grid_n; // Sets the resolution used by the IX() macro.
  srand(config.seed);

  // Initialize SDL and detect the screen (for fullscreen).
  if (SDL_Init(SDL_INIT_VIDEO) != 0) {
    fprintf(stderr, "Error: could not initialize SDL: %s\n", SDL_GetError());
    return EXIT_FAILURE;
  }

  if (config.fullscreen) {
    SDL_DisplayMode dm;
    if (SDL_GetCurrentDisplayMode(0, &dm) == 0 ||
        SDL_GetDesktopDisplayMode(0, &dm) == 0) {
      config.window_width = dm.w;
      config.window_height = dm.h;
    }
  }

  printf("Fluid simulation (Navier-Stokes) - SEQUENTIAL version\n");
  printf("  Grid         : %d x %d cells (%d interior cells)\n", config.grid_n,
         config.grid_n, config.grid_n * config.grid_n);
  printf("  Sources      : %d\n", config.num_sources);
  printf("  Window       : %d x %d px%s\n", config.window_width,
         config.window_height, config.fullscreen ? " (fullscreen)" : "");
  printf("  Seed         : %u\n", config.seed);
  printf("  visc / diff  : %.5f / %.5f\n", (double)config.viscosity,
         (double)config.diffusion);

  // Allocate the simulation field memory.
  if (!allocate_fields(&fields, config.grid_n)) {
    SDL_Quit();
    return EXIT_FAILURE;
  }

  // Allocate and initialize the ink sources.
  sources = (InkSource *)malloc((size_t)config.num_sources * sizeof(InkSource));
  if (sources == NULL) {
    fprintf(stderr, "Error: could not allocate memory for the sources.\n");
    SDL_Quit();
    free_fields(&fields);
    return EXIT_FAILURE;
  }
  init_sources(sources, config.num_sources, config.grid_n);

  // Allocate and initialize the background stars.
  stars = (Star *)malloc((size_t)BACKGROUND_STARS_COUNT * sizeof(Star));
  if (stars == NULL) {
    fprintf(stderr, "Error: could not allocate memory for the background.\n");
    SDL_Quit();
    free(sources);
    free_fields(&fields);
    return EXIT_FAILURE;
  }
  init_background(stars, BACKGROUND_STARS_COUNT, config.window_width,
                  config.window_height);

  // Create the window and renderer.
  SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "linear");

  window =
      SDL_CreateWindow("Fluid screensaver - sequential", SDL_WINDOWPOS_CENTERED,
                       SDL_WINDOWPOS_CENTERED, config.window_width,
                       config.window_height, SDL_WINDOW_SHOWN);
  if (window == NULL) {
    fprintf(stderr, "Error: could not create the window: %s\n", SDL_GetError());
    exit_code = EXIT_FAILURE;
    goto cleanup;
  }

  // Fullscreen is requested on the window after creation, not via a flag
  // to SDL_CreateWindow, because tiling window managers can drop a
  // fullscreen request made before the window is actually mapped.
  if (config.fullscreen) {
    if (SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN_DESKTOP) != 0) {
      fprintf(stderr, "Warning: could not enter fullscreen: %s\n",
              SDL_GetError());
    }
    SDL_ShowCursor(SDL_DISABLE);
  }

  // Gets the actual size SDL assigned (post-fullscreen, if requested).
  SDL_GetWindowSize(window, &config.window_width, &config.window_height);

  renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
  if (renderer == NULL) {
    // Retry with software rendering if acceleration isn't available.
    fprintf(stderr,
            "Notice: no accelerated renderer (%s). "
            "Trying software.\n",
            SDL_GetError());
    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
  }
  if (renderer == NULL) {
    fprintf(stderr, "Error: could not create the renderer: %s\n",
            SDL_GetError());
    exit_code = EXIT_FAILURE;
    goto cleanup;
  }

  texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888,
                              SDL_TEXTUREACCESS_STREAMING, config.window_width,
                              config.window_height);
  if (texture == NULL) {
    fprintf(stderr, "Error: could not create the texture: %s\n",
            SDL_GetError());
    exit_code = EXIT_FAILURE;
    goto cleanup;
  }
  // The ink texture carries alpha proportional to brightness (see
  // render.c). Without this, SDL_RenderCopy would treat it as opaque and
  // fully cover the background regardless of pixel alpha.
  SDL_SetTextureBlendMode(texture, SDL_BLENDMODE_BLEND);

  // Main loop. Performance measurement setup first.
  clock_frequency = SDL_GetPerformanceFrequency();
  prev_ticks = SDL_GetPerformanceCounter();
  last_report_ticks = prev_ticks;

  while (running) {
    // Window and keyboard events.
    while (SDL_PollEvent(&event)) {
      if (event.type == SDL_QUIT) {
        running = 0;
      } else if (event.type == SDL_KEYDOWN) {
        SDL_Keycode key = event.key.keysym.sym;
        if (key == SDLK_ESCAPE || key == SDLK_q) {
          running = 0;
        } else if (key == SDLK_r) {
          clear_fields(&fields);
          init_sources(sources, config.num_sources, config.grid_n);
        }
      }
    }

    // One simulation step. Move the sources first.
    update_nbody_sources(sources, config.num_sources, config.grid_n);

    // Inject ink and force from the sources.
    inject_sources(sources, config.num_sources, &fields,
                   (float)config.window_width / (float)config.window_height);

    // Solve the velocity field for this frame.
    velocity_step(&fields, config.viscosity, DT_DEFAULT);

    // Advect each ink color channel through that velocity.
    ink_step(&fields.ink_r, &fields.ink_r_p, fields.vel_x, fields.vel_y,
             config.diffusion, DT_DEFAULT, fields.total_cells);
    ink_step(&fields.ink_g, &fields.ink_g_p, fields.vel_x, fields.vel_y,
             config.diffusion, DT_DEFAULT, fields.total_cells);
    ink_step(&fields.ink_b, &fields.ink_b_p, fields.vel_x, fields.vel_y,
             config.diffusion, DT_DEFAULT, fields.total_cells);

    // Fade the ink so it doesn't saturate the screen.
    dissipate_ink(&fields);

    // Advance the background stars' twinkle.
    update_background(stars, BACKGROUND_STARS_COUNT);

    // Render this frame.
    render_ink(texture, &fields, config.window_width, config.window_height);

    SDL_SetRenderDrawColor(renderer, BACKGROUND_COLOR_R, BACKGROUND_COLOR_G,
                           BACKGROUND_COLOR_B, 255);
    SDL_RenderClear(renderer);
    draw_background(renderer, stars, BACKGROUND_STARS_COUNT);
    SDL_RenderCopy(renderer, texture, NULL, NULL);
    SDL_RenderPresent(renderer);

    // Measure FPS, averaged over roughly half-second windows.
    accumulated_frames++;
    current_ticks = SDL_GetPerformanceCounter();
    elapsed_seconds =
        (double)(current_ticks - last_report_ticks) / (double)clock_frequency;

    if (elapsed_seconds >= 0.5) {
      current_fps = (double)accumulated_frames / elapsed_seconds;
      accumulated_frames = 0;
      last_report_ticks = current_ticks;

      snprintf(window_title, sizeof(window_title),
               "Sequential fluids | N=%d | sources=%d | FPS= %.2f",
               config.grid_n, config.num_sources, current_fps);
      SDL_SetWindowTitle(window, window_title);

      printf("FPS= %.2f\n", current_fps);
      fflush(stdout);
    }

    prev_ticks = current_ticks;
  }

  (void)prev_ticks; // Unused, kept for clarity. Silences the compiler.

// Release resources.
cleanup:
  if (texture != NULL)
    SDL_DestroyTexture(texture);
  if (renderer != NULL)
    SDL_DestroyRenderer(renderer);
  if (window != NULL)
    SDL_DestroyWindow(window);
  SDL_Quit();

  free(sources);
  free(stars);
  free_fields(&fields);

  return exit_code;
}
