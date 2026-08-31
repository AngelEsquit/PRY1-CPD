#include "nbody.h"
#include "common.h"
#include "sources.h"
#include "utils.h"

#include <math.h>

/*
 * NBODY_G:           gravitational constant (arbitrary units, hand-tuned
 *                     together with the masses in init_sources() so the
 *                     motion is noticeable but not chaotic).
 * NBODY_SOFTENING:   effective minimum distance between two sources; keeps
 *                     the force (proportional to 1/dist^2) from blowing up
 *                     to infinity on a close encounter.
 * NBODY_VEL_MAX:     maximum velocity per axis (cells/frame); limits the
 *                     "slingshot effect" of a close encounter so explicit
 *                     integration stays stable.
 * NBODY_RESTITUTION: fraction of velocity kept when bouncing off a wall
 *                     (<1 = inelastic bounce, damped).
 */
#define NBODY_G            0.3f
#define NBODY_SOFTENING    6.0f
#define NBODY_VEL_MAX      5.0f
#define NBODY_RESTITUTION  0.9f

/*
 * Advances one frame of the n-body system formed by the sources: each
 * source attracts the others according to the law of universal gravitation
 * (F = G*m1*m2/d^2, softened) and bounces elastically off the grid's edges
 * to stay always visible on screen.
 *
 * The "frame" is used as the unit of time (same as angular_vel, already
 * expressed in rad/frame), so it doesn't depend on the fluid simulation's
 * dt: the sources' motion is a separate phenomenon that only shares the
 * grid with the solver.
 */
void update_nbody_sources(InkSource *sources, int count, int resolution)
{
    float accel_x[SOURCES_MAX] = { 0 };
    float accel_y[SOURCES_MAX] = { 0 };
    int   i, j;
    const int radius = source_radius(resolution);
    /* Sources can't leave the injection neighborhood (see source_radius(),
     * defined in sources.h) without losing valid interior cells. */
    const float lower_bound = (float)(radius + 1);
    const float upper_bound = (float)(resolution - radius);

    for (i = 0; i < count; i++) {
        for (j = i + 1; j < count; j++) {
            float diff_x   = sources[j].pos_x - sources[i].pos_x;
            float diff_y   = sources[j].pos_y - sources[i].pos_y;
            float dist2    = diff_x * diff_x + diff_y * diff_y +
                              NBODY_SOFTENING * NBODY_SOFTENING;
            float inv_dist = 1.0f / sqrtf(dist2);
            /* common factor G/d^3 (softened); multiplied by the other
             * body's mass gives the acceleration per Newton's 2nd law, and
             * the 3rd law (action-reaction) avoids recomputing the
             * symmetric pair. */
            float factor   = NBODY_G * inv_dist * inv_dist * inv_dist;

            accel_x[i] += factor * sources[j].mass * diff_x;
            accel_y[i] += factor * sources[j].mass * diff_y;
            accel_x[j] -= factor * sources[i].mass * diff_x;
            accel_y[j] -= factor * sources[i].mass * diff_y;
        }
    }

    for (i = 0; i < count; i++) {
        sources[i].vel_x = clamp(sources[i].vel_x + accel_x[i],
                                 -NBODY_VEL_MAX, NBODY_VEL_MAX);
        sources[i].vel_y = clamp(sources[i].vel_y + accel_y[i],
                                 -NBODY_VEL_MAX, NBODY_VEL_MAX);

        sources[i].pos_x += sources[i].vel_x;
        sources[i].pos_y += sources[i].vel_y;

        /* Elastic (damped) bounce off the grid's edges: the source stays
         * confined to the box instead of escaping the screen. */
        if (sources[i].pos_x < lower_bound) {
            sources[i].pos_x = lower_bound;
            sources[i].vel_x = -sources[i].vel_x * NBODY_RESTITUTION;
        } else if (sources[i].pos_x > upper_bound) {
            sources[i].pos_x = upper_bound;
            sources[i].vel_x = -sources[i].vel_x * NBODY_RESTITUTION;
        }

        if (sources[i].pos_y < lower_bound) {
            sources[i].pos_y = lower_bound;
            sources[i].vel_y = -sources[i].vel_y * NBODY_RESTITUTION;
        } else if (sources[i].pos_y > upper_bound) {
            sources[i].pos_y = upper_bound;
            sources[i].vel_y = -sources[i].vel_y * NBODY_RESTITUTION;
        }
    }
}
