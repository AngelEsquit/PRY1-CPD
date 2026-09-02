#include "solver.h"
#include "common.h"
#include "utils.h"

// Stable Fluids solver (Stam 1999). Every field goes through the same 4
// steps each frame: add_source -> diffuse -> advect -> project (project is
// velocity-only). Grid has a 1-cell ghost border (real cells are 1..grid_n);
// IX(i,j) maps (col,row) to an index into the flat 1D array backing every
// field.

// Fills the ghost border so the interior update never special-cases edges.
// BND_SCALAR: border copies its neighbor. BND_VEL_X/Y: copies but flips
// sign across the wall that axis is perpendicular to, so velocity reflects
// off walls instead of leaking through.
static void apply_boundary(int bnd_type, float *field) {
  int i;

  for (i = 1; i <= grid_n; i++) {
    // Left/right walls: flip for BND_VEL_X, else copy.
    field[IX(0, i)] =
        (bnd_type == BND_VEL_X) ? -field[IX(1, i)] : field[IX(1, i)];
    field[IX(grid_n + 1, i)] =
        (bnd_type == BND_VEL_X) ? -field[IX(grid_n, i)] : field[IX(grid_n, i)];
    // Top/bottom walls: flip for BND_VEL_Y, else copy.
    field[IX(i, 0)] =
        (bnd_type == BND_VEL_Y) ? -field[IX(i, 1)] : field[IX(i, 1)];
    field[IX(i, grid_n + 1)] =
        (bnd_type == BND_VEL_Y) ? -field[IX(i, grid_n)] : field[IX(i, grid_n)];
  }

  // Corners are outside the loop above: average their 2 neighbors.
  field[IX(0, 0)] = 0.5f * (field[IX(1, 0)] + field[IX(0, 1)]);
  field[IX(0, grid_n + 1)] =
      0.5f * (field[IX(1, grid_n + 1)] + field[IX(0, grid_n)]);
  field[IX(grid_n + 1, 0)] =
      0.5f * (field[IX(grid_n, 0)] + field[IX(grid_n + 1, 1)]);
  field[IX(grid_n + 1, grid_n + 1)] =
      0.5f * (field[IX(grid_n, grid_n + 1)] + field[IX(grid_n + 1, grid_n)]);
}

// u += dt * f: merges this frame's injected force/ink (`source`, written by
// sources.c) into the field.
static void add_source(float *dest, const float *source, float dt,
                       int total_cells) {
  int i;
  // Each cell only touches itself, safe to split across threads as-is.
  #pragma omp parallel for schedule(static)
  for (i = 0; i < total_cells; i++) {
    dest[i] += dt * source[i];
  }
}

// One Gauss-Seidel color of the checkerboard (parity = (i + j) % 2), called
// from within an active OpenMP parallel region (see solve_linear). Every
// cell's 4 neighbors are the opposite color, so same-parity cells have no
// data dependencies and can update in parallel.
static void relax_color(float *field, const float *field_prev, float a,
                        float inv_c, int parity) {
  #pragma omp for schedule(static)
  for (int j = 1; j <= grid_n; j++) {
    // First column of this row matching `parity`; the loop then steps by
    // 2 to stay on that color.
    int start_i = (((1 + j) % 2) == parity) ? 1 : 2;

    for (int i = start_i; i <= grid_n; i += 2) {
      // field = field_prev + a * sum(neighbors), normalized by inv_c =
      // 1/c. One sweep of Gauss-Seidel relaxation for A*field = field_prev.
      field[IX(i, j)] = (field_prev[IX(i, j)] +
                         a * (field[IX(i - 1, j)] + field[IX(i + 1, j)] +
                              field[IX(i, j - 1)] + field[IX(i, j + 1)])) *
                        inv_c;
    }
  }
}

// Solves the sparse system A*field = field_prev (one equation per cell) by
// running relax_color GAUSS_SEIDEL_ITERS times, red cells then black.
//
// A single OpenMP parallel region encloses all iterations to avoid
// thread-team fork/join overhead per sub-step. The implicit barrier at the
// end of each `omp for` synchronizes red and black sweeps; `omp single`
// restricts the boundary update to one thread.
static void solve_linear(int bnd_type, float *field, const float *field_prev,
                         float a, float c) {
  const float inv_c = 1.0f / c;

  #pragma omp parallel
  {
    for (int iter = 0; iter < GAUSS_SEIDEL_ITERS; iter++) {
      relax_color(field, field_prev, a, inv_c, 0);
      relax_color(field, field_prev, a, inv_c, 1);

      #pragma omp single
      {
        apply_boundary(bnd_type, field);
      }
    }
  }
}

// Diffusion term nu*laplacian(u), solved implicitly via solve_linear.
// a = nu*dt*n^2 is the per-cell coupling strength to its neighbors.
static void diffuse(int bnd_type, float *field, const float *field_prev,
                    float coefficient, float dt) {
  const float a = dt * coefficient * (float)grid_n * (float)grid_n;
  solve_linear(bnd_type, field, field_prev, a, 1.0f + 4.0f * a);
}

// Advection term -(u.grad)u via semi-Lagrangian backtrace: for each cell,
// trace backward along velocity by dt to find the source point, then
// bilinearly sample field_prev there.
static void advect(int bnd_type, float *field, const float *field_prev,
                   const float *vel_x, const float *vel_y, float dt) {
  int i, j, i0, i1, j0, j1;
  float origin_x, origin_y;
  float weight_i1, weight_i0, weight_j1, weight_j0;
  const float dt_grid = dt * (float)grid_n;
  const float limit = (float)grid_n + 0.5f;

  // Each destination cell only reads field_prev, never writes it, so the
  // grid splits across threads freely. collapse(2) merges both loops so
  // narrow grids (fewer rows than threads) still spread across all of them.
  #pragma omp parallel for collapse(2) \
      private(i0, i1, j0, j1, origin_x, origin_y, \
              weight_i1, weight_i0, weight_j1, weight_j0) schedule(static)
  for (j = 1; j <= grid_n; j++) {
    for (i = 1; i <= grid_n; i++) {
      // Backtrace: (origin_x, origin_y) is cell (i,j)'s position moved
      // backward one time step along its own velocity. dt_grid is dt
      // scaled to grid units. Generally lands between cells.
      origin_x = (float)i - dt_grid * vel_x[IX(i, j)];
      origin_y = (float)j - dt_grid * vel_y[IX(i, j)];

      // Clamp to valid range in case a large velocity sent the origin
      // outside the grid.
      origin_x = clamp(origin_x, 0.5f, limit);
      origin_y = clamp(origin_y, 0.5f, limit);

      // Enclosing 2x2 cell box: i0/j0 floor, i1/j1 the next cell over.
      i0 = (int)origin_x;
      i1 = i0 + 1;
      j0 = (int)origin_y;
      j1 = j0 + 1;

      // Fractional offset of the origin into that box, in [0, 1] per axis.
      weight_i1 = origin_x - (float)i0;
      weight_i0 = 1.0f - weight_i1;
      weight_j1 = origin_y - (float)j0;
      weight_j0 = 1.0f - weight_j1;

      // Bilinear interpolation of the 4 surrounding cells' old values.
      field[IX(i, j)] = weight_i0 * (weight_j0 * field_prev[IX(i0, j0)] +
                                     weight_j1 * field_prev[IX(i0, j1)]) +
                        weight_i1 * (weight_j0 * field_prev[IX(i1, j0)] +
                                     weight_j1 * field_prev[IX(i1, j1)]);
    }
  }

  apply_boundary(bnd_type, field);
}

// Velocity-only pressure projection. Enforces incompressibility
// (div(u) = 0), which diffuse/advect can violate slightly. 3 stages:
// measure divergence, solve for a pressure field that cancels it, subtract
// its gradient from velocity.
static void project(float *vel_x, float *vel_y, float *pressure,
                    float *divergence) {
  int i, j;
  const float h = 1.0f / (float)grid_n;

  // Stage 1: divergence of the velocity field per cell, central
  // difference: -0.5h * (dvel_x/dx + dvel_y/dy).
  #pragma omp parallel for collapse(2) schedule(static)
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

  // Stage 2: Poisson equation laplacian(pressure) = divergence, solved with
  // the same Gauss-Seidel relaxation as diffuse() (relax_color, 20 iters).
  solve_linear(BND_SCALAR, pressure, divergence, 1.0f, 4.0f);

  // Stage 3: vel -= grad(pressure), central difference. Subtracting the
  // pressure gradient makes the resulting field divergence-free.
  #pragma omp parallel for collapse(2) schedule(static)
  for (j = 1; j <= grid_n; j++) {
    for (i = 1; i <= grid_n; i++) {
      vel_x[IX(i, j)] -=
          0.5f * (pressure[IX(i + 1, j)] - pressure[IX(i - 1, j)]) / h;
      vel_y[IX(i, j)] -=
          0.5f * (pressure[IX(i, j + 1)] - pressure[IX(i, j - 1)]) / h;
    }
  }
  apply_boundary(BND_VEL_X, vel_x);
  apply_boundary(BND_VEL_Y, vel_y);
}

// Swaps which buffer counts as "current" vs "previous" (pointer swap, no
// copy). Needed because diffuse/advect read every cell's old value while
// writing every cell's new value, so each field keeps 2 buffers and flips
// which is which after each step.
static void swap_fields(float **field_a, float **field_b) {
  float *tmp = *field_a;
  *field_a = *field_b;
  *field_b = tmp;
}

// One color channel, one frame: inject, diffuse, advect. No project():
// that's velocity-only, ink just rides along.
void ink_step(float **ink, float **ink_prev, const float *vel_x,
              const float *vel_y, float diffusion, float dt, int total_cells) {
  add_source(*ink, *ink_prev, dt, total_cells);
  swap_fields(ink_prev, ink);
  diffuse(BND_SCALAR, *ink, *ink_prev, diffusion, dt);
  swap_fields(ink_prev, ink);
  advect(BND_SCALAR, *ink, *ink_prev, vel_x, vel_y, dt);
}

// Velocity, one frame: force, diffuse, project, self-advect, project again
// (diffuse and advect each reintroduce a little divergence, so projection
// runs twice).
void velocity_step(FluidFields *fields, float viscosity, float dt) {
  add_source(fields->vel_x, fields->vel_x_p, dt, fields->total_cells);
  add_source(fields->vel_y, fields->vel_y_p, dt, fields->total_cells);

  swap_fields(&fields->vel_x_p, &fields->vel_x);
  diffuse(BND_VEL_X, fields->vel_x, fields->vel_x_p, viscosity, dt);
  swap_fields(&fields->vel_y_p, &fields->vel_y);
  diffuse(BND_VEL_Y, fields->vel_y, fields->vel_y_p, viscosity, dt);

  project(fields->vel_x, fields->vel_y, fields->pressure, fields->divergence);

  // Velocity carries itself: vel_x_p/vel_y_p serve as both the field being
  // advected and the flow doing the advecting.
  swap_fields(&fields->vel_x_p, &fields->vel_x);
  swap_fields(&fields->vel_y_p, &fields->vel_y);
  advect(BND_VEL_X, fields->vel_x, fields->vel_x_p, fields->vel_x_p,
         fields->vel_y_p, dt);
  advect(BND_VEL_Y, fields->vel_y, fields->vel_y_p, fields->vel_x_p,
         fields->vel_y_p, dt);

  project(fields->vel_x, fields->vel_y, fields->pressure, fields->divergence);
}
