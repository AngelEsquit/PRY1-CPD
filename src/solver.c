#include "solver.h"
#include "common.h"
#include "utils.h"

// The fluid math (Stable Fluids, Stam 1999). Every field goes through the
// same 4 steps each frame: add_source -> diffuse -> advect -> project
// (project is velocity-only). Grid has a 1-cell ghost border (real cells
// are 1..grid_n); IX(i,j) turns a (col,row) into an index into the flat 1D
// array backing every field.

// Fills the ghost border so the real math never special-cases edge cells.
// BND_SCALAR: border copies its neighbor. BND_VEL_X/Y: copies but flips
// sign on the walls that axis is perpendicular to. That is what makes
// velocity bounce off walls instead of leaking through.
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

  // Corners aren't covered by the loop above, average their 2 neighbors.
  field[IX(0, 0)] = 0.5f * (field[IX(1, 0)] + field[IX(0, 1)]);
  field[IX(0, grid_n + 1)] =
      0.5f * (field[IX(1, grid_n + 1)] + field[IX(0, grid_n)]);
  field[IX(grid_n + 1, 0)] =
      0.5f * (field[IX(grid_n, 0)] + field[IX(grid_n + 1, 1)]);
  field[IX(grid_n + 1, grid_n + 1)] =
      0.5f * (field[IX(grid_n, grid_n + 1)] + field[IX(grid_n + 1, grid_n)]);
}

// Merges this frame's injected force/ink (in `source`, already written by
// sources.c) into the real field, scaled by dt.
static void add_source(float *dest, const float *source, float dt,
                       int total_cells) {
  int i;
  // Each cell only touches itself, safe to split across threads as-is.
  #pragma omp parallel for schedule(dynamic, 16)
  for (i = 0; i < total_cells; i++) {
    dest[i] += dt * source[i]; // Just add it in, cell by cell.
  }
}

// Relaxes a single color of the checkerboard (parity = (i + j) % 2) from
// within an active OpenMP parallel region (see solve_linear). The 4 neighbors
// of any colored cell always have the opposite color, so cells of the same
// parity have no data dependencies between them and can be computed in
// parallel without race conditions.
static void relax_color(float *field, const float *field_prev, float a,
                        float inv_c, int parity) {
  #pragma omp for schedule(dynamic, 16)
  for (int j = 1; j <= grid_n; j++) {
    int start_i = (((1 + j) % 2) == parity) ? 1 : 2;

    for (int i = start_i; i <= grid_n; i += 2) {
      field[IX(i, j)] = (field_prev[IX(i, j)] +
                         a * (field[IX(i - 1, j)] + field[IX(i + 1, j)] +
                              field[IX(i, j - 1)] + field[IX(i, j + 1)])) *
                        inv_c;
    }
  }
}

// Iteratively solves new_value = (starting_value + a * sum_of_4_neighbors) / c
// using Red-Black Gauss-Seidel relaxation. Each iteration is split into two
// sub-steps: red cells (parity 0) then black cells (parity 1).
//
// A single OpenMP parallel region encloses all GAUSS_SEIDEL_ITERS iterations
// to avoid thread team fork/join overhead on every sub-step. Implicit
// barriers at the end of each omp for ensure synchronization between red
// and black sweeps, while omp single ensures only one thread updates the
// boundary conditions.
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

// Smooths the field toward its neighbors' average. Just works out the a/c
// constants solve_linear() needs from the diffusion rate, dt, and grid
// size, then hands off.
static void diffuse(int bnd_type, float *field, const float *field_prev,
                    float coefficient, float dt) {
  const float a = dt * coefficient * (float)grid_n * (float)grid_n;
  solve_linear(bnd_type, field, field_prev, a, 1.0f + 4.0f * a);
}

// Carries the field along the velocity field. For each cell: trace
// backward along velocity to find where its new value should come from,
// then bilinearly interpolate (blend the 4 nearest cells) since that point
// almost never lands exactly on a cell.
static void advect(int bnd_type, float *field, const float *field_prev,
                   const float *vel_x, const float *vel_y, float dt) {
  int i, j, i0, i1, j0, j1;
  float origin_x, origin_y;
  float weight_i1, weight_i0, weight_j1, weight_j0;
  const float dt_grid = dt * (float)grid_n;
  const float limit = (float)grid_n + 0.5f;

  // Each destination cell only reads field_prev (never writes it), so the
  // grid can be split across threads freely. collapse(2) merges both loops
  // into one range so narrow grids (fewer rows than threads) still spread
  // work across all of them.
  #pragma omp parallel for collapse(2) \
      private(i0, i1, j0, j1, origin_x, origin_y, \
              weight_i1, weight_i0, weight_j1, weight_j0) schedule(dynamic, 16)
  for (j = 1; j <= grid_n; j++) {
    for (i = 1; i <= grid_n; i++) {
      // Step backward along this cell's own velocity.
      origin_x = (float)i - dt_grid * vel_x[IX(i, j)];
      origin_y = (float)j - dt_grid * vel_y[IX(i, j)];

      // Keep it inside the grid so the lookup below stays in bounds.
      origin_x = clamp(origin_x, 0.5f, limit);
      origin_y = clamp(origin_y, 0.5f, limit);

      // The 4 real cells surrounding that (fractional) origin point.
      i0 = (int)origin_x;
      i1 = i0 + 1;
      j0 = (int)origin_y;
      j1 = j0 + 1;

      // How close the origin is to each side, used to weight the blend.
      weight_i1 = origin_x - (float)i0;
      weight_i0 = 1.0f - weight_i1;
      weight_j1 = origin_y - (float)j0;
      weight_j0 = 1.0f - weight_j1;

      // Blend those 4 cells' old values by distance (bilinear interpolation).
      field[IX(i, j)] = weight_i0 * (weight_j0 * field_prev[IX(i0, j0)] +
                                     weight_j1 * field_prev[IX(i0, j1)]) +
                        weight_i1 * (weight_j0 * field_prev[IX(i1, j0)] +
                                     weight_j1 * field_prev[IX(i1, j1)]);
    }
  }

  apply_boundary(bnd_type, field);
}

// Velocity-only. Real fluid can't locally bunch up or thin out, so this
// measures where that's happening, solves for a pressure field that would
// cancel it, and subtracts that pressure's push from velocity. This is
// also what makes the flow swirl instead of just spreading outward.
static void project(float *vel_x, float *vel_y, float *pressure,
                    float *divergence) {
  int i, j;
  const float h = 1.0f / (float)grid_n;

  // Measure how much more velocity leaves each cell than enters it.
  #pragma omp parallel for collapse(2) schedule(dynamic, 16)
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

  // Solve for the pressure field that would cancel that imbalance out.
  solve_linear(BND_SCALAR, pressure, divergence, 1.0f, 4.0f);

  // Push velocity away from high pressure, toward low.
  #pragma omp parallel for collapse(2) schedule(dynamic, 16)
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

// Swaps which buffer counts as "current" vs "previous". No data copied,
// just the pointers. Needed because diffuse/advect read every cell's old
// value while writing every cell's new value, so each field keeps 2
// buffers and flips which is which after each step.
static void swap_fields(float **field_a, float **field_b) {
  float *tmp = *field_a;
  *field_a = *field_b;
  *field_b = tmp;
}

// One color channel, one frame: inject, smooth, carry along velocity.
// No project(): that's only for velocity, ink just rides along.
// Takes float** so swap_fields() can actually repoint *ink/*ink_prev.
void ink_step(float **ink, float **ink_prev, const float *vel_x,
              const float *vel_y, float diffusion, float dt, int total_cells) {
  add_source(*ink, *ink_prev, dt, total_cells);
  swap_fields(ink_prev, ink);
  diffuse(BND_SCALAR, *ink, *ink_prev, diffusion, dt);
  swap_fields(ink_prev, ink);
  advect(BND_SCALAR, *ink, *ink_prev, vel_x, vel_y, dt);
}

// Velocity, one frame: force, smooth, fix imbalance, self-carry, fix
// imbalance again (smoothing and self-carrying each reintroduce a little
// imbalance, so it needs cleaning up twice).
void velocity_step(FluidFields *fields, float viscosity, float dt) {
  add_source(fields->vel_x, fields->vel_x_p, dt, fields->total_cells);
  add_source(fields->vel_y, fields->vel_y_p, dt, fields->total_cells);

  swap_fields(&fields->vel_x_p, &fields->vel_x);
  diffuse(BND_VEL_X, fields->vel_x, fields->vel_x_p, viscosity, dt);
  swap_fields(&fields->vel_y_p, &fields->vel_y);
  diffuse(BND_VEL_Y, fields->vel_y, fields->vel_y_p, viscosity, dt);

  project(fields->vel_x, fields->vel_y, fields->pressure, fields->divergence);

  // Velocity carries itself: vel_x_p/vel_y_p are both the field being
  // moved and the flow doing the moving.
  swap_fields(&fields->vel_x_p, &fields->vel_x);
  swap_fields(&fields->vel_y_p, &fields->vel_y);
  advect(BND_VEL_X, fields->vel_x, fields->vel_x_p, fields->vel_x_p,
         fields->vel_y_p, dt);
  advect(BND_VEL_Y, fields->vel_y, fields->vel_y_p, fields->vel_x_p,
         fields->vel_y_p, dt);

  project(fields->vel_x, fields->vel_y, fields->pressure, fields->divergence);
}
