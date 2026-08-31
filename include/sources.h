#ifndef SOURCES_H
#define SOURCES_H

#include "fields.h"

/* ===========================================================================
 * sources.h
 * Ink sources: initialization, per-frame injection and dissipation.
 * ======================================================================== */

/*
 * Radius (in cells) of the injection neighborhood, calibrated for a
 * reference grid of SOURCE_GRID_REF cells and scaled with the actual
 * resolution (see source_radius/source_sigma) so the blob always covers the
 * same fraction of the screen regardless of grid resolution: with a fixed
 * radius, a finer grid shrinks how many screen pixels that same handful of
 * cells occupies, and the ink is born looking like "a handful of pixels"
 * instead of a smooth blob.
 *
 * Each cell's Gaussian weight is divided by the square of the scale factor
 * (see inject_sources): the neighborhood's area (~radius*sigma) grows with
 * the scale, and without that correction the total ink deposited per source
 * would grow with it and saturate the screen white (already happened once:
 * see history). Dividing by scale^2 keeps the total invariant and only
 * improves how many cells draw the circle.
 */
#define SOURCE_RADIUS_BASE   2
#define SOURCE_SIGMA_BASE    1.1f
#define SOURCE_GRID_REF      128

static inline float source_scale(int resolution)
{
    float scale = (float)resolution / (float)SOURCE_GRID_REF;
    return (scale < 1.0f) ? 1.0f : scale;
}

static inline int source_radius(int resolution)
{
    float scale  = source_scale(resolution);
    int   radius = (int)(SOURCE_RADIUS_BASE * scale + 0.5f);
    return (radius < SOURCE_RADIUS_BASE) ? SOURCE_RADIUS_BASE : radius;
}

static inline float source_sigma(int resolution)
{
    return SOURCE_SIGMA_BASE * source_scale(resolution);
}

/*
 * Struct: InkSource
 * A point emitter that injects color and momentum into the grid. Sources
 * move across the grid as an n-body system: each one gravitationally
 * attracts the others (see nbody.h).
 */
typedef struct {
    float pos_x, pos_y;          /* continuous position on the grid (cells)   */
    float vel_x, vel_y;          /* drift velocity (cells/frame)              */
    float mass;                  /* gravitational mass (n-body)               */
    float color_r, color_g, color_b; /* ink color in [0,1]                    */
    float phase;                 /* current angular phase (radians)           */
    float angular_vel;           /* jet spin rate (rad/frame)                 */
    float force;                 /* magnitude of the injected velocity        */
    float flow_rate;             /* amount of ink injected per frame          */
} InkSource;

/*
 * Initializes the sources with pseudo-random position, velocity, mass,
 * color, phase and force. Positions are kept away from the border so the
 * jet can develop before hitting the walls.
 */
void init_sources(InkSource *sources, int count, int resolution);

/*
 * Deposits into the "previous" buffers (which act as the source term) the
 * ink and momentum of each source for the current frame.
 *
 * TRIGONOMETRIC ELEMENT: each source's jet direction rotates over time as
 * (cos(phase), sin(phase)), which generates spiral vortices. Ink is
 * deposited over a neighborhood of cells with Gaussian falloff to smooth
 * out the injection and avoid harsh pointwise values.
 *
 * "aspect" is window_width/window_height. The grid is always square
 * (grid_n x grid_n) but the window doesn't have to be: render.c uses a
 * different scale on x and y to fill the window, so a grid cell doesn't map
 * to a square on screen but to a rectangle "aspect" times wider than tall.
 * Without correcting for this, an isotropic Gaussian kernel (same sigma on
 * x and y, in cell units) would look on screen like an ellipse wider than
 * tall. Here sigma on x is shrunk by that same factor so that, after
 * render.c's anisotropic scaling, the blob looks circular again.
 */
void inject_sources(InkSource *sources, int count, FluidFields *fields,
                    float aspect);

/*
 * Multiplies the ink by a factor below 1 so it fades out gradually; without
 * this the screen would end up saturated white.
 */
void dissipate_ink(FluidFields *fields);

#endif /* SOURCES_H */
