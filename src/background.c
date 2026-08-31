#include "background.h"
#include "common.h"
#include "utils.h"

#include <math.h>

// Places every star at a random position with random brightness, phase and
// twinkle speed.
void init_background(Star *stars, int count, int width, int height)
{
    int e;

    for (e = 0; e < count; e++) {
        stars[e].x = random_range(0.0f, (float)(width - 1));
        stars[e].y = random_range(0.0f, (float)(height - 1));

        // Multiplying two random values skews most stars dim, with only a
        // few bright ones. Matches a real starry sky better than a flat
        // random range would.
        stars[e].base_brightness = random_range(0.15f, 1.0f) *
                                    random_range(0.15f, 1.0f);

        stars[e].phase       = random_range(0.0f, 2.0f * PI);
        stars[e].phase_speed = random_range(0.01f, 0.04f);

        // Most stars are 1 pixel, a few are 2, as if slightly closer.
        stars[e].radius = (random_range(0.0f, 1.0f) > 0.92f) ? 2 : 1;
    }
}

// Advances every star's twinkle phase by one frame.
void update_background(Star *stars, int count)
{
    int e;

    for (e = 0; e < count; e++) {
        stars[e].phase += stars[e].phase_speed;
        if (stars[e].phase > 2.0f * PI) {
            stars[e].phase -= 2.0f * PI;
        }
    }
}

// Draws every star with the renderer directly (not a texture). Brightness
// oscillates between roughly 30% and 100% of its base value so it flickers
// without ever going fully dark.
void draw_background(SDL_Renderer *renderer, const Star *stars, int count)
{
    int e;

    for (e = 0; e < count; e++) {
        float twinkle    = 0.65f + 0.35f * sinf(stars[e].phase);
        float brightness = stars[e].base_brightness * twinkle;
        Uint8 value       = (Uint8)clamp(brightness * 255.0f, 0.0f, 255.0f);

        // Slightly bluish white instead of pure white, to match the
        // background's cool base color.
        SDL_SetRenderDrawColor(renderer,
                               (Uint8)(value * 0.85f),
                               (Uint8)(value * 0.92f),
                               value, 255);

        if (stars[e].radius <= 1) {
            SDL_RenderDrawPoint(renderer, (int)stars[e].x,
                                (int)stars[e].y);
        } else {
            SDL_Rect r;
            r.x = (int)stars[e].x;
            r.y = (int)stars[e].y;
            r.w = stars[e].radius;
            r.h = stars[e].radius;
            SDL_RenderFillRect(renderer, &r);
        }
    }
}
