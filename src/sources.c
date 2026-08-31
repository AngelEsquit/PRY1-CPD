#include "sources.h"
#include "common.h"
#include "utils.h"

#include <string.h>
#include <math.h>

// Sets each source to a random position, velocity, mass, color, phase and
// force. Positions stay away from the border so the jet can develop before
// hitting a wall.
void init_sources(InkSource *sources, int count, int resolution)
{
    int f;
    const int margin = (resolution / 8 > 2) ? resolution / 8 : 2;

    for (f = 0; f < count; f++) {
        sources[f].pos_x = (float)margin +
            random_range(0.0f, (float)(resolution - 2 * margin));
        sources[f].pos_y = (float)margin +
            random_range(0.0f, (float)(resolution - 2 * margin));

        // Small initial velocity. Real motion comes from mutual gravity
        // once the sim gets going.
        sources[f].vel_x = random_range(-0.3f, 0.3f);
        sources[f].vel_y = random_range(-0.3f, 0.3f);
        sources[f].mass  = random_range(25.0f, 70.0f);

        // Saturated color: one dominant channel, two partial ones.
        sources[f].color_r = random_range(0.15f, 1.0f);
        sources[f].color_g = random_range(0.15f, 1.0f);
        sources[f].color_b = random_range(0.15f, 1.0f);

        sources[f].phase       = random_range(0.0f, 2.0f * PI);
        sources[f].angular_vel = random_range(-0.045f, 0.045f);
        // Force stays in [6,14]. Any higher and the ink jumps cell to cell
        // instead of flowing (looks like droplets, not a jet), since each
        // frame's displacement would outrun the injection blob's size.
        sources[f].force     = random_range(6.0f, 14.0f);
        sources[f].flow_rate = random_range(60.0f, 110.0f);
    }
}

// How narrow the hot core is relative to the color blob. Smaller means a
// more concentrated, pointlike spot.
#define CORE_SIGMA_FACTOR      0.32f
// How much white brightness gets added at the exact center, as a fraction
// of the source's flow rate.
#define CORE_BRIGHTNESS_FACTOR 0.9f

// Deposits ink and momentum for every source into the "previous" buffers
// (the source term for this frame's solve). Jet direction spins over time
// via (cos(phase), sin(phase)), producing the spiral shape. Ink is spread
// over a small Gaussian-weighted neighborhood instead of a single cell so
// the blob looks round instead of pixelated, and the weight is computed
// against each source's real (fractional) position so it slides smoothly
// as the source moves instead of jumping cell to cell.
void inject_sources(InkSource *sources, int count, FluidFields *fields,
                    float aspect)
{
    int f, dx, dy, cell_i, cell_j, center_x, center_y;
    float direction_x, direction_y;
    const int   radius = source_radius(grid_n);
    const float sigma  = source_sigma(grid_n);
    const float scale  = source_scale(grid_n);
    // The window usually isn't square, so render.c scales x and y
    // differently. Shrinking sigma_x by "aspect" here cancels that out so
    // the blob still looks circular on screen.
    const float sigma_x = sigma / aspect;
    // Shrinking sigma_x also shrinks the blob's area, which would reduce
    // the total ink deposited. This multiplies that back in so total ink
    // per source stays the same regardless of resolution or aspect ratio.
    const float normalizer = aspect / (scale * scale);
    // A second, much narrower Gaussian centered on the same point, added
    // equally to all 3 color channels. Pushes the exact center toward
    // white/bright and fades out fast, so each source has a bright core
    // that blends into its own color a short distance out.
    const float core_sigma_x = sigma_x * CORE_SIGMA_FACTOR;
    const float core_sigma   = sigma   * CORE_SIGMA_FACTOR;

    // Source buffers get cleared and rebuilt every frame.
    const size_t bytes = (size_t)fields->total_cells * sizeof(float);
    memset(fields->vel_x_p, 0, bytes);
    memset(fields->vel_y_p, 0, bytes);
    memset(fields->ink_r_p, 0, bytes);
    memset(fields->ink_g_p, 0, bytes);
    memset(fields->ink_b_p, 0, bytes);

    for (f = 0; f < count; f++) {
        // Spin the jet's direction forward one step.
        sources[f].phase += sources[f].angular_vel;
        if (sources[f].phase > 2.0f * PI) {
            sources[f].phase -= 2.0f * PI;
        }

        direction_x = cosf(sources[f].phase) * sources[f].force;
        direction_y = sinf(sources[f].phase) * sources[f].force;

        // floor() just picks which cell to start scanning the neighborhood
        // from. Each cell's actual weight below still uses the source's
        // real, continuous position.
        center_x = (int)floorf(sources[f].pos_x);
        center_y = (int)floorf(sources[f].pos_y);

        for (dy = -radius; dy <= radius; dy++) {
            for (dx = -radius; dx <= radius; dx++) {
                cell_i = center_x + dx;
                cell_j = center_y + dy;

                // Skip cells outside the real (non-ghost) grid.
                if (cell_i < 1 || cell_i > grid_n ||
                    cell_j < 1 || cell_j > grid_n) {
                    continue;
                }

                // Gaussian weight against the source's continuous
                // position, not the integer cell center.
                float diff_x = (float)cell_i - sources[f].pos_x;
                float diff_y = (float)cell_j - sources[f].pos_y;
                float dist2 = (diff_x * diff_x) / (sigma_x * sigma_x) +
                              (diff_y * diff_y) / (sigma   * sigma);
                float weight = expf(-dist2 / 2.0f) * normalizer;

                // Same idea, but with the narrow "hot core" sigma, so it
                // only contributes right at the center.
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

// Fades every ink channel by a factor below 1 each frame. Without this the
// screen would just keep accumulating ink and saturate to white forever.
void dissipate_ink(FluidFields *fields)
{
    int i;
    #pragma omp parallel for schedule(dynamic, 16)
    for (i = 0; i < fields->total_cells; i++) {
        fields->ink_r[i] *= DISSIPATION;
        fields->ink_g[i] *= DISSIPATION;
        fields->ink_b[i] *= DISSIPATION;
    }
}
