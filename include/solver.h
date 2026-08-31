#ifndef SOLVER_H
#define SOLVER_H

#include "fields.h"

/* ===========================================================================
 * solver.h
 * ---------------------------------------------------------------------------
 * Numerical core of the "Stable Fluids" method: diffusion, advection,
 * projection and the full velocity and ink steps.
 *
 * This is the module that concentrates most of the computational cost and,
 * therefore, the main candidate for OpenMP parallelization in the next
 * phase of the project (see the note on Gauss-Seidel in solver.c).
 * ======================================================================== */

/*
 * Full step of one channel's ink density:
 * source -> diffusion -> advection.
 */
void ink_step(float **ink, float **ink_prev,
             const float *vel_x, const float *vel_y,
             float diffusion, float dt, int total_cells);

/*
 * Full step of the velocity field:
 * forces -> viscous diffusion -> projection -> self-advection -> projection.
 * It's projected twice because both diffusion and advection reintroduce
 * divergence into the field.
 */
void velocity_step(FluidFields *fields, float viscosity, float dt);

#endif /* SOLVER_H */
