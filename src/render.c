#include "render.h"
#include "common.h"
#include "utils.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>

/* Exponent of alpha's gamma curve (see render_ink): the smaller it is, the
 * faster alpha reaches opaque with little ink brightness. */
#define ALPHA_GAMMA 0.7f

/*
 * Bilinearly interpolates "field" (an N x N grid with a ghost ring) at the
 * continuous position (fx, fy), expressed in the same space as a pixel
 * index of the original N x N texture: fx=0 is the center of interior cell
 * 1, fx=grid_n-1 is the center of interior cell grid_n.
 */
static float sample_bilinear(const float *field, float fx, float fy)
{
    int i0, i1, j0, j1;
    float weight_x1, weight_x0, weight_y1, weight_y0;

    fx = clamp(fx, 0.0f, (float)(grid_n - 1));
    fy = clamp(fy, 0.0f, (float)(grid_n - 1));

    i0 = (int)fx;  i1 = (i0 < grid_n - 1) ? i0 + 1 : i0;
    j0 = (int)fy;  j1 = (j0 < grid_n - 1) ? j0 + 1 : j0;

    weight_x1 = fx - (float)i0;  weight_x0 = 1.0f - weight_x1;
    weight_y1 = fy - (float)j0;  weight_y0 = 1.0f - weight_y1;

    /* +1 shifts from pixel space (0..N-1) to the interior cell index
     * (1..N) inside the grid with its ghost ring. */
    return weight_x0 * (weight_y0 * field[IX(i0 + 1, j0 + 1)] +
                        weight_y1 * field[IX(i0 + 1, j1 + 1)]) +
           weight_x1 * (weight_y0 * field[IX(i1 + 1, j0 + 1)] +
                        weight_y1 * field[IX(i1 + 1, j1 + 1)]);
}

/*
 * Dumps the ink density into the SDL texture, one texture pixel per
 * destination pixel (see render.h for why: SDL2's software renderer
 * doesn't do linear scaling).
 */
void render_ink(SDL_Texture *texture, const FluidFields *fields,
                int texture_width, int texture_height)
{
    void  *raw_pixels = NULL;
    int    pitch_bytes = 0;
    Uint32 *pixels;
    int    i, j;
    float  fx, fy;
    float  value_r, value_g, value_b, average;
    const float scale_x = (float)grid_n / (float)texture_width;
    const float scale_y = (float)grid_n / (float)texture_height;

    if (SDL_LockTexture(texture, NULL, &raw_pixels, &pitch_bytes) != 0) {
        fprintf(stderr, "Warning: could not lock the texture: %s\n",
                SDL_GetError());
        return;
    }

    pixels = (Uint32 *)raw_pixels;

    for (j = 0; j < texture_height; j++) {
        Uint32 *row = pixels + (size_t)j * ((size_t)pitch_bytes / sizeof(Uint32));
        fy = ((float)j + 0.5f) * scale_y - 0.5f;

        for (i = 0; i < texture_width; i++) {
            fx = ((float)i + 0.5f) * scale_x - 0.5f;

            /* A Reinhard-style mapping (x / (x + k)) is applied instead of
             * a hard clip, so the buildup of many sources compresses
             * smoothly toward white instead of clipping abruptly. */
            value_r = sample_bilinear(fields->ink_r, fx, fy) * BRIGHTNESS_FACTOR;
            value_g = sample_bilinear(fields->ink_g, fx, fy) * BRIGHTNESS_FACTOR;
            value_b = sample_bilinear(fields->ink_b, fx, fy) * BRIGHTNESS_FACTOR;

            value_r = value_r / (value_r + CONTRAST_FACTOR);
            value_g = value_g / (value_g + CONTRAST_FACTOR);
            value_b = value_b / (value_b + CONTRAST_FACTOR);

            /* Pushes each channel away from the pixel's average gray to
             * recover saturation without changing the mean brightness,
             * keeping the result from looking washed-out/white. */
            average = (value_r + value_g + value_b) / 3.0f;
            value_r = clamp(average + (value_r - average) * SATURATION_FACTOR, 0.0f, 1.0f);
            value_g = clamp(average + (value_g - average) * SATURATION_FACTOR, 0.0f, 1.0f);
            value_b = clamp(average + (value_b - average) * SATURATION_FACTOR, 0.0f, 1.0f);

            /* Alpha proportional to the ink's brightness (instead of a
             * fixed 0xFF): where there's no ink the pixel is transparent
             * and the background (stars/vignette) drawn earlier on the
             * renderer shows through; where the ink saturates, the pixel
             * is opaque and covers it completely. The highest of the three
             * channels is used as a brightness approximation (keeps a
             * color heavily saturated in a single channel, with a low
             * average, from looking more transparent than it should). */
            {
                float max_brightness, alpha;

                max_brightness = value_r;
                if (value_g > max_brightness) max_brightness = value_g;
                if (value_b > max_brightness) max_brightness = value_b;

                /* Using max_brightness directly as alpha leaves most of
                 * the blob (anything not already saturated to full
                 * brightness) semi-transparent, and SDL_BLENDMODE_BLEND
                 * mixes it with the dark background behind it: the result
                 * looks washed out, as if the ink's color lost intensity.
                 * A gamma curve (<1) pushes alpha to opaque much faster
                 * than the real brightness, so the blob looks as vivid as
                 * it did before the background was added, and only the
                 * edge where the ink drops to nearly zero stays
                 * transparent (so the stars show through there). */
                alpha = powf(max_brightness, ALPHA_GAMMA);

                row[i] = ((Uint32)(alpha * 255.0f) << 24) |
                         ((Uint32)(value_r * 255.0f) << 16) |
                         ((Uint32)(value_g * 255.0f) <<  8) |
                         ((Uint32)(value_b * 255.0f));
            }
        }
    }

    SDL_UnlockTexture(texture);
}
