#include "solver.h"
#include "common.h"
#include "utils.h"

/* ===========================================================================
 * solver.c -- the numerical core (Jos Stam's "Stable Fluids", 1999)
 * ---------------------------------------------------------------------------
 * Each frame, velocity_step() and ink_step() run the same four operators in
 * sequence on their respective fields:
 *
 *   add_source  -> apply_boundary is NOT called here; sources.c already
 *                  wrote directly into the *_p ("previous") buffers, this
 *                  just adds dt * that into the real field.
 *   diffuse     -> solve_linear(): an implicit (unconditionally stable)
 *                  diffusion solve via Gauss-Seidel relaxation. This is the
 *                  hotspot -- GAUSS_SEIDEL_ITERS sweeps of the whole grid,
 *                  called multiple times per frame.
 *   advect      -> semi-Lagrangian backward trace + bilinear interpolation.
 *                  Also unconditionally stable, at the cost of some
 *                  numerical smoothing.
 *   project     -> ONLY applied to velocity, never to ink. Removes
 *                  divergence (enforces incompressibility) via a Poisson
 *                  solve, which is also what produces the swirling motion.
 *
 * apply_boundary() is the plumbing all of the above shares: it fills the
 * grid's 1-cell ghost ring (see its own comment below) so that diffuse/
 * advect/project never need special-case logic for edge cells -- they just
 * read whatever is in the ghost ring as if it were a real neighbor.
 *
 * velocity_step() is the only one that calls project(), and calls it twice
 * (once after diffusion, once after advection) because both of those steps
 * reintroduce a little divergence into the field, undoing the previous
 * projection.
 * ======================================================================== */

/*
 * Applies the boundary conditions over the ring of ghost cells.
 *   bnd_type = BND_SCALAR -> the border copies the interior neighbor (Neumann)
 *   bnd_type = BND_VEL_X  -> the x component is inverted at vertical walls
 *   bnd_type = BND_VEL_Y  -> the y component is inverted at horizontal walls
 * Inverting the sign is equivalent to a solid wall: the fluid bounces
 * instead of passing through it.
 */
static void apply_boundary(int bnd_type, float *field)
{
    int i;

    for (i = 1; i <= grid_n; i++) {
        field[IX(0, i)]         = (bnd_type == BND_VEL_X)
                                  ? -field[IX(1, i)]      : field[IX(1, i)];
        field[IX(grid_n + 1, i)] = (bnd_type == BND_VEL_X)
                                  ? -field[IX(grid_n, i)] : field[IX(grid_n, i)];
        field[IX(i, 0)]         = (bnd_type == BND_VEL_Y)
                                  ? -field[IX(i, 1)]      : field[IX(i, 1)];
        field[IX(i, grid_n + 1)] = (bnd_type == BND_VEL_Y)
                                  ? -field[IX(i, grid_n)] : field[IX(i, grid_n)];
    }

    /* The four corners are averaged from their two neighbors */
    field[IX(0, 0)] = 0.5f * (field[IX(1, 0)] + field[IX(0, 1)]);
    field[IX(0, grid_n + 1)] =
        0.5f * (field[IX(1, grid_n + 1)] + field[IX(0, grid_n)]);
    field[IX(grid_n + 1, 0)] =
        0.5f * (field[IX(grid_n, 0)] + field[IX(grid_n + 1, 1)]);
    field[IX(grid_n + 1, grid_n + 1)] =
        0.5f * (field[IX(grid_n, grid_n + 1)] + field[IX(grid_n + 1, grid_n)]);
}

/*
 * Adds the source field's contributions, scaled by dt, into the destination
 * field. This is the external-forces / ink-injection operator.
 */
static void add_source(float *dest, const float *source, float dt,
                       int total_cells)
{
    int i;
    for (i = 0; i < total_cells; i++) {
        dest[i] += dt * source[i];
    }
}

/*
 * Solves the sparse linear system   x = (x0 + a * neighbor_sum) / c
 * via Gauss-Seidel relaxation.
 *
 * IMPORTANT NOTE FOR THE PARALLEL VERSION:
 * Gauss-Seidel uses values already updated within the same iteration,
 * which creates a data dependency between neighboring cells and makes this
 * loop NOT directly parallelizable. The OpenMP version replaces it with
 * Jacobi (which reads from a separate buffer and therefore has no
 * dependencies) or with a red-black sweep. Gauss-Seidel is kept here as the
 * sequential reference method, and the one that converges fastest per
 * iteration.
 */
static void solve_linear(int bnd_type, float *field, const float *field_prev,
                         float a, float c)
{
    int iter, i, j;
    const float inv_c = 1.0f / c;

    for (iter = 0; iter < GAUSS_SEIDEL_ITERS; iter++) {
        for (j = 1; j <= grid_n; j++) {
            for (i = 1; i <= grid_n; i++) {
                field[IX(i, j)] = (field_prev[IX(i, j)] +
                                   a * (field[IX(i - 1, j)] + field[IX(i + 1, j)] +
                                        field[IX(i, j - 1)] + field[IX(i, j + 1)]))
                                  * inv_c;
            }
        }
        apply_boundary(bnd_type, field);
    }
}

/*
 * Implicit diffusion: solves  x - a*lap(x) = x0,  with a = dt*coef*N*N.
 * The implicit scheme is unconditionally stable (that's the key to Stam's
 * method, as opposed to explicit diffusion).
 */
static void diffuse(int bnd_type, float *field, const float *field_prev,
                    float coefficient, float dt)
{
    const float a = dt * coefficient * (float)grid_n * (float)grid_n;
    solve_linear(bnd_type, field, field_prev, a, 1.0f + 4.0f * a);
}

/*
 * Semi-Lagrangian advection: for each cell, trace backward in time along
 * the velocity field and bilinearly interpolate the field's value at the
 * origin position. This scheme is also unconditionally stable, regardless
 * of velocity magnitude.
 */
static void advect(int bnd_type, float *field, const float *field_prev,
                   const float *vel_x, const float *vel_y, float dt)
{
    int   i, j, i0, i1, j0, j1;
    float origin_x, origin_y;
    float weight_i1, weight_i0, weight_j1, weight_j0;
    const float dt_grid = dt * (float)grid_n;
    const float limit   = (float)grid_n + 0.5f;

    for (j = 1; j <= grid_n; j++) {
        for (i = 1; i <= grid_n; i++) {
            /* Backward trace of the particle arriving at (i,j) */
            origin_x = (float)i - dt_grid * vel_x[IX(i, j)];
            origin_y = (float)j - dt_grid * vel_y[IX(i, j)];

            /* Clamp the origin to the domain to avoid reading out of bounds */
            origin_x = clamp(origin_x, 0.5f, limit);
            origin_y = clamp(origin_y, 0.5f, limit);

            i0 = (int)origin_x;  i1 = i0 + 1;
            j0 = (int)origin_y;  j1 = j0 + 1;

            /* Bilinear interpolation weights */
            weight_i1 = origin_x - (float)i0;  weight_i0 = 1.0f - weight_i1;
            weight_j1 = origin_y - (float)j0;  weight_j0 = 1.0f - weight_j1;

            field[IX(i, j)] =
                weight_i0 * (weight_j0 * field_prev[IX(i0, j0)] +
                             weight_j1 * field_prev[IX(i0, j1)]) +
                weight_i1 * (weight_j0 * field_prev[IX(i1, j0)] +
                             weight_j1 * field_prev[IX(i1, j1)]);
        }
    }

    apply_boundary(bnd_type, field);
}

/*
 * Hodge projection: decomposes the velocity field into a divergence-free
 * part plus the gradient of a pressure field, and keeps only the former.
 * This enforces the incompressibility condition div(u) = 0, and is what
 * produces the fluid's characteristic swirls.
 */
static void project(float *vel_x, float *vel_y, float *pressure,
                    float *divergence)
{
    int i, j;
    const float h = 1.0f / (float)grid_n;

    /* 1. Compute the divergence of the velocity field */
    for (j = 1; j <= grid_n; j++) {
        for (i = 1; i <= grid_n; i++) {
            divergence[IX(i, j)] = -0.5f * h *
                (vel_x[IX(i + 1, j)] - vel_x[IX(i - 1, j)] +
                 vel_y[IX(i, j + 1)] - vel_y[IX(i, j - 1)]);
            pressure[IX(i, j)] = 0.0f;
        }
    }
    apply_boundary(BND_SCALAR, divergence);
    apply_boundary(BND_SCALAR, pressure);

    /* 2. Solve the Poisson equation  lap(p) = div(u) */
    solve_linear(BND_SCALAR, pressure, divergence, 1.0f, 4.0f);

    /* 3. Subtract the pressure gradient from the velocity */
    for (j = 1; j <= grid_n; j++) {
        for (i = 1; i <= grid_n; i++) {
            vel_x[IX(i, j)] -= 0.5f *
                (pressure[IX(i + 1, j)] - pressure[IX(i - 1, j)]) / h;
            vel_y[IX(i, j)] -= 0.5f *
                (pressure[IX(i, j + 1)] - pressure[IX(i, j - 1)]) / h;
        }
    }
    apply_boundary(BND_VEL_X, vel_x);
    apply_boundary(BND_VEL_Y, vel_y);
}

/* Swaps two field pointers (avoids copying whole arrays). */
static void swap_fields(float **field_a, float **field_b)
{
    float *tmp = *field_a;
    *field_a = *field_b;
    *field_b = tmp;
}

/*
 * Full step of one channel's ink density:
 * source -> diffusion -> advection.
 */
void ink_step(float **ink, float **ink_prev,
             const float *vel_x, const float *vel_y,
             float diffusion, float dt, int total_cells)
{
    add_source(*ink, *ink_prev, dt, total_cells);
    swap_fields(ink_prev, ink);
    diffuse(BND_SCALAR, *ink, *ink_prev, diffusion, dt);
    swap_fields(ink_prev, ink);
    advect(BND_SCALAR, *ink, *ink_prev, vel_x, vel_y, dt);
}

/*
 * Full step of the velocity field:
 * forces -> viscous diffusion -> projection -> self-advection -> projection.
 * It's projected twice because both diffusion and advection reintroduce
 * divergence into the field.
 */
void velocity_step(FluidFields *fields, float viscosity, float dt)
{
    add_source(fields->vel_x, fields->vel_x_p, dt, fields->total_cells);
    add_source(fields->vel_y, fields->vel_y_p, dt, fields->total_cells);

    swap_fields(&fields->vel_x_p, &fields->vel_x);
    diffuse(BND_VEL_X, fields->vel_x, fields->vel_x_p, viscosity, dt);
    swap_fields(&fields->vel_y_p, &fields->vel_y);
    diffuse(BND_VEL_Y, fields->vel_y, fields->vel_y_p, viscosity, dt);

    project(fields->vel_x, fields->vel_y, fields->pressure, fields->divergence);

    swap_fields(&fields->vel_x_p, &fields->vel_x);
    swap_fields(&fields->vel_y_p, &fields->vel_y);
    advect(BND_VEL_X, fields->vel_x, fields->vel_x_p,
          fields->vel_x_p, fields->vel_y_p, dt);
    advect(BND_VEL_Y, fields->vel_y, fields->vel_y_p,
          fields->vel_x_p, fields->vel_y_p, dt);

    project(fields->vel_x, fields->vel_y, fields->pressure, fields->divergence);
}
