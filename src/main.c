/* ===========================================================================
 * main.c
 * ---------------------------------------------------------------------------
 * Fluid screensaver based on the Navier-Stokes equations for incompressible
 * fluids, solved with Jos Stam's "Stable Fluids" method.
 *
 * SEQUENTIAL VERSION (comparison baseline for the OpenMP parallel version).
 *
 * Physical model
 * --------------
 * The two coupled equations are solved over a regular 2D grid:
 *
 *   Velocity:  du/dt = -(u . grad)u + visc * lap(u) + f      with div(u) = 0
 *   Ink:       dp/dt = -(u . grad)p + diff * lap(p) + s
 *
 * Each time step is broken down into four operators (Stam, 1999):
 *   1. add_source  -> applies the external forces f and the ink sources s
 *   2. diffuse     -> solves viscous diffusion implicitly
 *   3. advect      -> transports the field via semi-Lagrangian tracing
 *   4. project     -> Hodge projection: removes divergence (incompressible)
 *
 * Steps 2 and 4 require solving a sparse linear system; here iterative
 * Gauss-Seidel relaxation is used (see the note in solver.c on its impact
 * on the parallel version).
 *
 * Trigonometry: each ink source injects velocity in a direction that
 * rotates over time, computed with sin() and cos() over its own phase
 * (see sources.c).
 *
 * Project modules (see include/ and src/):
 *   common.h        - constants, macros and grid indexing
 *   utils.{h,c}     - general-purpose functions
 *   config.{h,c}    - command-line arguments
 *   fields.{h,c}    - simulation field memory
 *   sources.{h,c}   - ink sources
 *   nbody.{h,c}     - source movement (n-body system)
 *   solver.{h,c}    - numerical core (Stable Fluids)
 *   render.{h,c}    - dump to an SDL texture
 *   main.c          - main program (this file)
 *
 * Build:  make            (see Makefile)
 * Usage:  ./Screensaver -n 128 -f 6
 * ===========================================================================
 */

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

/* Grid resolution; global so the IX() macro stays readable. The project's
 * only definition (declared extern in common.h so every other module can
 * use it). */
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

  // Args
  int args_result = parse_arguments(argc, argv, &config);
  if (args_result == -1)
    return EXIT_SUCCESS; /* help was requested       */
  if (args_result == 0)
    return EXIT_FAILURE; /* invalid arguments        */

  grid_n = config.grid_n; /* sets the resolution used by the IX() macro  */
  srand(config.seed);

  // SDL init / screen detection
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

  // Allocate field memory
  if (!allocate_fields(&fields, config.grid_n)) {
    SDL_Quit();
    return EXIT_FAILURE;
  }

  // Init ink sources
  sources = (InkSource *)malloc((size_t)config.num_sources * sizeof(InkSource));
  if (sources == NULL) {
    fprintf(stderr, "Error: could not allocate memory for the sources.\n");
    SDL_Quit();
    free_fields(&fields);
    return EXIT_FAILURE;
  }
  init_sources(sources, config.num_sources, config.grid_n);

  // Init stars / bg
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

  // Window rendering / creation
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

  /* Fullscreen is requested on the window after creation, not via a
   * flag to SDL_CreateWindow due to issues w/ tiling WM's. */
  if (config.fullscreen) {
    if (SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN_DESKTOP) != 0) {
      fprintf(stderr, "Warning: could not enter fullscreen: %s\n",
              SDL_GetError());
    }
    SDL_ShowCursor(SDL_DISABLE);
  }

  /* Gets the actual size SDL assigned (post-fullscreen, if requested) */
  SDL_GetWindowSize(window, &config.window_width, &config.window_height);

  renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
  if (renderer == NULL) {
    /* Retry with software rendering if acceleration isn't available */
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
  /* The ink texture now carries alpha proportional to brightness (see
   * render.c); without this SDL_RenderCopy would treat it as opaque and
   * fully cover the background (stars) regardless of pixel alpha. */
  SDL_SetTextureBlendMode(texture, SDL_BLENDMODE_BLEND);

  // Main loop

  // Performance  measurements
  clock_frequency = SDL_GetPerformanceFrequency();
  prev_ticks = SDL_GetPerformanceCounter();
  last_report_ticks = prev_ticks;

  while (running) {
    // Window / kb events
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

    // Simulation per step:

    // Move the n bodies / source emitters
    update_nbody_sources(sources, config.num_sources, config.grid_n);

    // Inject ink from the sources
    inject_sources(sources, config.num_sources, &fields,
                   (float)config.window_width / (float)config.window_height);

    // Update velocity field (stable fluids)
    velocity_step(&fields, config.viscosity, DT_DEFAULT);

    // Update each rgb channel's scalar field
    ink_step(&fields.ink_r, &fields.ink_r_p, fields.vel_x, fields.vel_y,
             config.diffusion, DT_DEFAULT, fields.total_cells);
    ink_step(&fields.ink_g, &fields.ink_g_p, fields.vel_x, fields.vel_y,
             config.diffusion, DT_DEFAULT, fields.total_cells);
    ink_step(&fields.ink_b, &fields.ink_b_p, fields.vel_x, fields.vel_y,
             config.diffusion, DT_DEFAULT, fields.total_cells);

    // Dissipate rgb channels (avoid saturating screen w/ ink)
    dissipate_ink(&fields);

    // Update stars in bg
    update_background(stars, BACKGROUND_STARS_COUNT);

    // Draw / render
    render_ink(texture, &fields, config.window_width, config.window_height);

    SDL_SetRenderDrawColor(renderer, BACKGROUND_COLOR_R, BACKGROUND_COLOR_G,
                           BACKGROUND_COLOR_B, 255);
    SDL_RenderClear(renderer);
    draw_background(renderer, stars, BACKGROUND_STARS_COUNT);
    SDL_RenderCopy(renderer, texture, NULL, NULL);
    SDL_RenderPresent(renderer);

    // FPS measurement (0.5s timeframe avg)
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

  (void)prev_ticks; // make compiler happy :)

// Resource release
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
