#include "nbody.h"
#include "common.h"
#include "sources.h"
#include "utils.h"

#include <math.h>

// NBODY_G: gravitational constant, hand-tuned with the masses in
// init_sources() so the motion is noticeable but not chaotic.
// NBODY_SOFTENING: minimum effective distance between two sources. Keeps
// the force (proportional to 1/dist^2) from blowing up on a close
// encounter.
// NBODY_VEL_MAX: max velocity per axis (cells/frame). Caps the slingshot
// effect from a close encounter so the integration stays stable.
// NBODY_RESTITUTION: fraction of velocity kept when bouncing off a wall
// (below 1 means a damped, inelastic bounce).
#define NBODY_G            0.3f
#define NBODY_SOFTENING    6.0f
#define NBODY_VEL_MAX      5.0f
#define NBODY_RESTITUTION  0.9f

// Advances one frame of the n-body system. Each source attracts every
// other source (F = G*m1*m2/d^2, softened) and bounces off the grid's
// edges to stay on screen.
//
// This runs on a per-frame time unit, not the fluid solver's dt. Source
// movement is a separate system that only happens to share the grid.
void update_nbody_sources(InkSource *sources, int count, int resolution)
{
    float accel_x[SOURCES_MAX] = { 0 };
    float accel_y[SOURCES_MAX] = { 0 };
    int   i, j;
    const int radius = source_radius(resolution);
    // Sources can't leave the injection neighborhood (source_radius() in
    // sources.h) without losing valid interior cells.
    const float lower_bound = (float)(radius + 1);
    const float upper_bound = (float)(resolution - radius);

    // O(count^2): every source pulls on every other source once.
    for (i = 0; i < count; i++) {
        for (j = i + 1; j < count; j++) {
            float diff_x   = sources[j].pos_x - sources[i].pos_x;
            float diff_y   = sources[j].pos_y - sources[i].pos_y;
            float dist2    = diff_x * diff_x + diff_y * diff_y +
                              NBODY_SOFTENING * NBODY_SOFTENING;
            float inv_dist = 1.0f / sqrtf(dist2);
            // Common factor G/d^3. Multiply by the other body's mass to get
            // acceleration. Action-reaction lets both sources get updated
            // from a single pair check instead of computing it twice.
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

        // Damped bounce off the grid's edges, keeps the source on screen.
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
