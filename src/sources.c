#include "sources.h"
#include "common.h"
#include "utils.h"

#include <string.h>
#include <math.h>

/*
 * Initializes the sources with pseudo-random position, velocity, mass,
 * color, phase and force. Positions are kept away from the border so the
 * jet can develop before hitting the walls.
 */
void init_sources(InkSource *sources, int count, int resolution)
{
    int f;
    const int margin = (resolution / 8 > 2) ? resolution / 8 : 2;

    for (f = 0; f < count; f++) {
        sources[f].pos_x = (float)margin +
            random_range(0.0f, (float)(resolution - 2 * margin));
        sources[f].pos_y = (float)margin +
            random_range(0.0f, (float)(resolution - 2 * margin));

        /* Small initial velocity: the main motion comes from mutual
         * gravitational attraction once the simulation gets going. */
        sources[f].vel_x = random_range(-0.3f, 0.3f);
        sources[f].vel_y = random_range(-0.3f, 0.3f);
        sources[f].mass  = random_range(25.0f, 70.0f);

        /* Saturated color: one dominant channel and two partial ones */
        sources[f].color_r = random_range(0.15f, 1.0f);
        sources[f].color_g = random_range(0.15f, 1.0f);
        sources[f].color_b = random_range(0.15f, 1.0f);

        sources[f].phase       = random_range(0.0f, 2.0f * PI);
        sources[f].angular_vel = random_range(-0.045f, 0.045f);
        /* With force in [25,60] the advection displacement right next to
         * the source (dt*grid_n*|v|) reaches ~15-19 cells per frame:
         * several times the diameter of the injection neighborhood, so the
         * ink "jumps" instead of flowing, looking like falling droplets
         * instead of a continuous jet. At [6,14] the jump stays around
         * ~5-6 cells, on the order of the injection neighborhood itself,
         * and consecutive frames overlap into a continuous stroke. */
        sources[f].force     = random_range(6.0f, 14.0f);
        sources[f].flow_rate = random_range(60.0f, 110.0f);
    }
}

/*
 * Deposits into the "previous" buffers (which act as the source term) the
 * ink and momentum of each source for the current frame.
 *
 * TRIGONOMETRIC ELEMENT: each source's jet direction rotates over time as
 * (cos(phase), sin(phase)), which generates spiral vortices. Ink is
 * deposited over a neighborhood of cells with Gaussian falloff (instead of
 * a uniform value) so the blob is born already round: with few cells the
 * shape reads as a diamond/pixelated blob as soon as it's created, because
 * there aren't enough samples to read as a circle, and since ink diffusion
 * is 0 by default nothing smooths it afterward except motion. Each cell's
 * weight is computed against the source's continuous position (not the
 * nearest cell), so the blob slides continuously as it moves instead of
 * jumping cell to cell.
 */
/* How narrow the hot core is relative to the color blob (fraction of
 * sigma); smaller = a more concentrated, pointlike spot. */
#define CORE_SIGMA_FACTOR      0.32f
/* How much white brightness is added at the exact center, as a fraction of
 * the source's flow rate (so a stronger jet also has a brighter core,
 * instead of a fixed white that would stand out differently per source). */
#define CORE_BRIGHTNESS_FACTOR 0.9f

void inject_sources(InkSource *sources, int count, FluidFields *fields,
                    float aspect)
{
    int f, dx, dy, cell_i, cell_j, center_x, center_y;
    float direction_x, direction_y;
    const int   radius = source_radius(grid_n);
    const float sigma  = source_sigma(grid_n);
    const float scale  = source_scale(grid_n);
    /* render.c uses a different scale on x and y to fill a window that
     * isn't square (see the comment in sources.h), so an equal sigma on
     * both axes would look like an ellipse on screen. Shrinking sigma_x by
     * "aspect" compensates for that stretch so the blob looks circular. */
    const float sigma_x = sigma / aspect;
    /* The neighborhood's area (~sigma^2) grows with the scale; dividing by
     * scale^2 keeps the total ink deposited per source invariant,
     * regardless of how many cells the circle uses to draw itself.
     * Shrinking sigma_x reduces that area by a factor of "aspect", so it's
     * multiplied back in to avoid losing ink just from correcting the
     * shape. */
    const float normalizer = aspect / (scale * scale);
    /* "Hot" core: a second Gaussian, much narrower and centered at the same
     * point, added equally to all three color channels (instead of the
     * source's own color). Adding pure white pushes the exact center of
     * the source toward white/bright, distinguishing it from its trail's
     * color; moving away from the center the core fades out much faster
     * than the color blob (via the CORE_SIGMA_FACTOR factor), so the
     * effect is a bright point that quickly blends into the source's own
     * color. Uses the same "normalizer" as the color blob: since both
     * sigmas scale the same way with grid resolution, the same
     * resolution-correction factor still applies. */
    const float core_sigma_x = sigma_x * CORE_SIGMA_FACTOR;
    const float core_sigma   = sigma   * CORE_SIGMA_FACTOR;

    /* The source buffers are cleared every frame */
    const size_t bytes = (size_t)fields->total_cells * sizeof(float);
    memset(fields->vel_x_p, 0, bytes);
    memset(fields->vel_y_p, 0, bytes);
    memset(fields->ink_r_p, 0, bytes);
    memset(fields->ink_g_p, 0, bytes);
    memset(fields->ink_b_p, 0, bytes);

    for (f = 0; f < count; f++) {
        /* Advance the jet's angular phase */
        sources[f].phase += sources[f].angular_vel;
        if (sources[f].phase > 2.0f * PI) {
            sources[f].phase -= 2.0f * PI;
        }

        direction_x = cosf(sources[f].phase) * sources[f].force;
        direction_y = sinf(sources[f].phase) * sources[f].force;

        /* The position is continuous (moved by the n-body system); only
         * floor() is used to locate the neighborhood's base cell, each
         * cell's weight uses the real continuous position. */
        center_x = (int)floorf(sources[f].pos_x);
        center_y = (int)floorf(sources[f].pos_y);

        for (dy = -radius; dy <= radius; dy++) {
            for (dx = -radius; dx <= radius; dx++) {
                cell_i = center_x + dx;
                cell_j = center_y + dy;

                /* Only write within the interior cells */
                if (cell_i < 1 || cell_i > grid_n ||
                    cell_j < 1 || cell_j > grid_n) {
                    continue;
                }

                /* Gaussian kernel against the source's real continuous
                 * position (not the neighborhood's integer center), so the
                 * blob slides smoothly between cells. */
                float diff_x = (float)cell_i - sources[f].pos_x;
                float diff_y = (float)cell_j - sources[f].pos_y;
                float dist2 = (diff_x * diff_x) / (sigma_x * sigma_x) +
                              (diff_y * diff_y) / (sigma   * sigma);
                float weight = expf(-dist2 / 2.0f) * normalizer;

                /* Same calculation as "weight" but with the hot core's
                 * narrow sigma, so it falls to zero much faster with
                 * distance and only contributes brightness right at the
                 * center. */
                float dist2_core = (diff_x * diff_x) / (core_sigma_x * core_sigma_x) +
                                   (diff_y * diff_y) / (core_sigma   * core_sigma);
                float weight_core  = expf(-dist2_core / 2.0f) * normalizer;
                float brightness   = CORE_BRIGHTNESS_FACTOR * sources[f].flow_rate *
                                      weight_core;

                fields->vel_x_p[IX(cell_i, cell_j)] += direction_x * weight;
                fields->vel_y_p[IX(cell_i, cell_j)] += direction_y * weight;

                fields->ink_r_p[IX(cell_i, cell_j)] +=
                    sources[f].flow_rate * sources[f].color_r * weight + brightness;
                fields->ink_g_p[IX(cell_i, cell_j)] +=
                    sources[f].flow_rate * sources[f].color_g * weight + brightness;
                fields->ink_b_p[IX(cell_i, cell_j)] +=
                    sources[f].flow_rate * sources[f].color_b * weight + brightness;
            }
        }
    }
}

/*
 * Multiplies the ink by a factor below 1 so it fades out gradually; without
 * this the screen would end up saturated white.
 */
void dissipate_ink(FluidFields *fields)
{
    int i;
    for (i = 0; i < fields->total_cells; i++) {
        fields->ink_r[i] *= DISSIPATION;
        fields->ink_g[i] *= DISSIPATION;
        fields->ink_b[i] *= DISSIPATION;
    }
}
