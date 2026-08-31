#include "solver.h"
#include "common.h"
#include "utils.h"

// =============================================================================
// solver.c -- the actual fluid math (Jos Stam's "Stable Fluids" method, 1999)
// -----------------------------------------------------------------------------
// Every frame, each field (velocity, and each ink color channel separately)
// gets pushed through the same 4 steps, in this order:
//
//   1. add_source  -- dump in whatever was injected this frame (force / ink)
//   2. diffuse     -- blur/smooth the field a little (viscosity for velocity,
//                      diffusion rate for ink)
//   3. advect      -- carry the field along wherever the velocity is pointing
//   4. project     -- (velocity only) cancel out any "ballooning" so the
//                      fluid doesn't compress or expand anywhere -- this is
//                      also what makes it swirl instead of just spreading out
//
// The grid is one cell bigger than the simulated area on every side (a
// "ghost ring" -- see apply_boundary() below). Every real cell is indexed
// 1..grid_n; index 0 and grid_n+1 are the fake border cells.
//
// IX(i, j) (defined in common.h) converts a (column, row) position into a
// single index, because the grid is stored as one flat 1D array instead of
// a real 2D array. Moving one column over (i+1) moves 1 slot in memory;
// moving one row down (j+1) jumps forward by a whole row's width. It shows
// up on almost every line below because that's just how you read/write "the
// cell above me" or "the cell to my left" out of a flat array.
// =============================================================================

// Fills in the fake border ring so the real math never needs a special case
// for cells sitting at the edge of the grid.
//   bnd_type = BND_SCALAR -> border cell just copies its real neighbor
//                             (used for ink and pressure: nothing bounces,
//                             it just can't leak past the wall)
//   bnd_type = BND_VEL_X  -> border cell copies its neighbor but FLIPS THE
//                             SIGN, on the left/right walls (used for
//                             velocity's x component)
//   bnd_type = BND_VEL_Y  -> same flip, but on the top/bottom walls (used
//                             for velocity's y component)
// Flipping the sign is what makes velocity "bounce": if real fluid right
// next to a wall is flowing INTO it, the ghost cell holds the same speed
// flowing back OUT, so the solver reads it as an equal push coming back at
// it and the net result cancels out instead of leaking through.
static void apply_boundary(int bnd_type, float *field) {
  int i;

  // Walk every real edge cell and fix up its matching ghost cell on each of
  // the 4 sides at once (left, right, bottom, top).
  for (i = 1; i <= grid_n; i++) {
    // Left wall ghost cell = copy (or flip) the real cell just inside it
    field[IX(0, i)] =
        (bnd_type == BND_VEL_X) ? -field[IX(1, i)] : field[IX(1, i)];
    // Right wall ghost cell, same idea
    field[IX(grid_n + 1, i)] =
        (bnd_type == BND_VEL_X) ? -field[IX(grid_n, i)] : field[IX(grid_n, i)];
    // Bottom wall ghost cell
    field[IX(i, 0)] =
        (bnd_type == BND_VEL_Y) ? -field[IX(i, 1)] : field[IX(i, 1)];
    // Top wall ghost cell
    field[IX(i, grid_n + 1)] =
        (bnd_type == BND_VEL_Y) ? -field[IX(i, grid_n)] : field[IX(i, grid_n)];
  }

  // The 4 corner ghost cells don't belong to any single wall above (the
  // loop only fills in edges, not corners), so just average each corner's
  // two neighboring ghost cells to give it a sane value too.
  field[IX(0, 0)] = 0.5f * (field[IX(1, 0)] + field[IX(0, 1)]);
  field[IX(0, grid_n + 1)] =
      0.5f * (field[IX(1, grid_n + 1)] + field[IX(0, grid_n)]);
  field[IX(grid_n + 1, 0)] =
      0.5f * (field[IX(grid_n, 0)] + field[IX(grid_n + 1, 1)]);
  field[IX(grid_n + 1, grid_n + 1)] =
      0.5f * (field[IX(grid_n, grid_n + 1)] + field[IX(grid_n + 1, grid_n)]);
}

// Adds whatever was injected this frame into the real field. "source" here
// is really the *_p buffer that sources.c already wrote the injected
// force/ink into earlier this frame -- this just merges it in, scaled by
// how much time (dt) has passed.
static void add_source(float *dest, const float *source, float dt,
                       int total_cells) {
  int i;
  // Every cell in the whole grid (ghost ring included -- doesn't matter
  // here, boundaries get fixed up separately by whoever calls this).
  for (i = 0; i < total_cells; i++) {
    dest[i] += dt * source[i];
  }
}

// Iteratively solves for a new grid where every cell's value is a blend of
// where it started (field_prev) and its neighbors' CURRENT values, using
// the formula:
//
//     new_value = (starting_value + a * sum_of_4_neighbors) / c
//
// This is used for two different things depending on what a/c are set to:
// smoothing a field toward its neighbors (diffusion), or solving for a
// pressure field that cancels out unwanted expansion/contraction (used
// inside project()). Both are "every cell depends on its neighbors, which
// depend on their neighbors" problems -- there's no way to compute the
// exact answer in one pass, so instead this nudges every cell a little
// closer to correct, over and over, until it's converged enough.
//
// IMPORTANT FOR THE PARALLEL VERSION: this writes each cell's new value
// directly back into `field`, then immediately reads that fresh value again
// while computing its right-hand neighbor. That's what makes it converge
// fast, but it also means cell (i+1, j) literally cannot be computed at the
// same time as cell (i, j) without a race -- this loop is NOT safely
// parallelizable as written. The parallel version needs a different
// approach here (Jacobi, which only reads last pass's values, or a
// red/black checkerboard split).
static void solve_linear(int bnd_type, float *field, const float *field_prev,
                         float a, float c) {
  int iter, i, j;
  const float inv_c = 1.0f / c;

  // Repeat the whole-grid nudge GAUSS_SEIDEL_ITERS times -- each pass gets
  // closer to the true answer; 20 passes is "good enough for real-time
  // graphics," not an exact solve.
  for (iter = 0; iter < GAUSS_SEIDEL_ITERS; iter++) {
    // Sweep every real cell, top-left to bottom-right.
    for (j = 1; j <= grid_n; j++) {
      for (i = 1; i <= grid_n; i++) {
        // Blend this cell's starting value with its 4 neighbors' current
        // values, weighted by `a`, then normalize by `c`.
        field[IX(i, j)] = (field_prev[IX(i, j)] +
                           a * (field[IX(i - 1, j)] + field[IX(i + 1, j)] +
                                field[IX(i, j - 1)] + field[IX(i, j + 1)])) *
                          inv_c;
      }
    }
    // Refresh the ghost ring after every pass -- the interior cells right
    // next to the border need a correct "fake neighbor" to read on the
    // NEXT pass, otherwise edge cells would slowly drift wrong.
    apply_boundary(bnd_type, field);
  }
}

// Smooths the field toward its own neighbors' average, at a rate controlled
// by `coefficient` (viscosity for velocity, the diffusion rate for ink).
// This doesn't do the smoothing itself -- it just works out the two numbers
// (`a`, `c`) that solve_linear() needs to produce that smoothing rate, then
// hands off to it. Bigger `coefficient`, bigger grid, or bigger `dt` all
// push `a` up, which makes each cell pull harder toward its neighbors.
static void diffuse(int bnd_type, float *field, const float *field_prev,
                    float coefficient, float dt) {
  const float a = dt * coefficient * (float)grid_n * (float)grid_n;
  solve_linear(bnd_type, field, field_prev, a, 1.0f + 4.0f * a);
}

// Moves the field along wherever the velocity is pointing. For every cell,
// instead of asking "where does my stuff go," it asks the reverse: "walk
// backward along the velocity -- where did the stuff that just arrived here
// come from?" Then it copies (interpolated) whatever was at that earlier
// spot into this cell now.
static void advect(int bnd_type, float *field, const float *field_prev,
                   const float *vel_x, const float *vel_y, float dt) {
  int i, j, i0, i1, j0, j1;
  float origin_x, origin_y;
  float weight_i1, weight_i0, weight_j1, weight_j0;
  const float dt_grid = dt * (float)grid_n;
  const float limit = (float)grid_n + 0.5f;

  for (j = 1; j <= grid_n; j++) {
    for (i = 1; i <= grid_n; i++) {
      // Step backward from this cell along its own velocity to find where
      // its new value should be copied from.
      origin_x = (float)i - dt_grid * vel_x[IX(i, j)];
      origin_y = (float)j - dt_grid * vel_y[IX(i, j)];

      // If that backward step lands outside the grid, pull it back to the
      // nearest valid spot -- otherwise the lookup below would read
      // out-of-bounds memory.
      origin_x = clamp(origin_x, 0.5f, limit);
      origin_y = clamp(origin_y, 0.5f, limit);

      // The origin point almost never lands exactly on a cell (e.g. it
      // might be (4.3, 5.7)), so grab the 4 real cells surrounding it.
      i0 = (int)origin_x;
      i1 = i0 + 1;
      j0 = (int)origin_y;
      j1 = j0 + 1;

      // Work out how close the origin is to each side, so the 4
      // surrounding cells can be blended by distance in the next step
      // (this is bilinear interpolation: closer cells count for more).
      weight_i1 = origin_x - (float)i0;
      weight_i0 = 1.0f - weight_i1;
      weight_j1 = origin_y - (float)j0;
      weight_j0 = 1.0f - weight_j1;

      // Blend the 4 surrounding cells' old values by those weights to get
      // this cell's new value.
      field[IX(i, j)] = weight_i0 * (weight_j0 * field_prev[IX(i0, j0)] +
                                     weight_j1 * field_prev[IX(i0, j1)]) +
                        weight_i1 * (weight_j0 * field_prev[IX(i1, j0)] +
                                     weight_j1 * field_prev[IX(i1, j1)]);
    }
  }

  // Fix up the ghost ring now that every real cell has a new value.
  apply_boundary(bnd_type, field);
}

// Only ever called on velocity, never on ink. Real fluid can't locally
// bunch up or thin out anywhere -- whatever flows into a spot has to flow
// back out somewhere. This function measures where that rule is currently
// being broken, figures out a "pressure" field that would push things back
// into balance, and subtracts that push from velocity. What's left over is
// velocity that can no longer pile up or empty out anywhere -- it can only
// curl around itself, which is exactly what produces the swirls you see on
// screen.
static void project(float *vel_x, float *vel_y, float *pressure,
                    float *divergence) {
  int i, j;
  const float h = 1.0f / (float)grid_n;

  // Step 1: for every cell, measure how much more velocity is leaving than
  // entering (a positive number = this cell is acting like a source,
  // pushing fluid outward; negative = acting like a drain). This should be
  // zero everywhere for a real incompressible fluid, but after diffuse/
  // advect it usually isn't.
  for (j = 1; j <= grid_n; j++) {
    for (i = 1; i <= grid_n; i++) {
      divergence[IX(i, j)] = -0.5f * h *
                             (vel_x[IX(i + 1, j)] - vel_x[IX(i - 1, j)] +
                              vel_y[IX(i, j + 1)] - vel_y[IX(i, j - 1)]);
      pressure[IX(i, j)] = 0.0f; // start the solve below from a blank slate
    }
  }
  apply_boundary(BND_SCALAR, divergence);
  apply_boundary(BND_SCALAR, pressure);

  // Step 2: solve for the pressure field that would need to exist to
  // exactly cancel out the imbalance measured above. Reuses the same
  // neighbor-blending solver as diffuse() -- different equation, same
  // "every cell depends on its neighbors" iterative trick.
  solve_linear(BND_SCALAR, pressure, divergence, 1.0f, 4.0f);

  // Step 3: push velocity away from high pressure and toward low pressure --
  // this is the actual correction that removes the imbalance.
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

// Swaps which buffer is "current" and which is "previous" without copying
// any actual data -- just trades the two pointers. Needed because several
// steps above (diffuse, advect) need to read every cell's OLD value while
// writing every cell's NEW value, so the code keeps two buffers per field
// and flips which one is "the one to write into" after each step, instead
// of copying the whole grid around.
static void swap_fields(float **field_a, float **field_b) {
  float *tmp = *field_a;
  *field_a = *field_b;
  *field_b = tmp;
}

// Runs one color channel's ink through the full pipeline for this frame:
// add whatever was injected -> let it blur a little -> carry it along by
// the velocity field. No project() here -- that step only matters for
// velocity (it's the incompressibility fix); ink is just along for the
// ride and has nothing to stay "balanced."
//
// Takes float** (pointer to the field pointer) instead of float* because
// swap_fields() above needs to actually change which buffer `*ink` and
// `*ink_prev` point at -- in C, a function can only modify a caller's
// variable if it's handed the address of that variable.
void ink_step(float **ink, float **ink_prev, const float *vel_x,
              const float *vel_y, float diffusion, float dt, int total_cells) {
  add_source(*ink, *ink_prev, dt, total_cells);
  swap_fields(ink_prev, ink); // *ink_prev now holds what was just injected
  diffuse(BND_SCALAR, *ink, *ink_prev, diffusion, dt); // write smoothed result into *ink
  swap_fields(ink_prev, ink); // flip again so *ink_prev holds the smoothed result
  advect(BND_SCALAR, *ink, *ink_prev, vel_x, vel_y, dt); // write final carried result into *ink
}

// Runs the velocity field through the full frame pipeline: add forces,
// smooth (viscosity), cancel out any ballooning, let velocity carry
// itself along (yes -- velocity moves itself, since it's also a flow),
// then cancel out ballooning again. It has to be done twice because both
// the smoothing step and the self-carrying step above each nudge some
// ballooning back in as a side effect, undoing the first cleanup.
void velocity_step(FluidFields *fields, float viscosity, float dt) {
  // Add this frame's injected force into both velocity components.
  add_source(fields->vel_x, fields->vel_x_p, dt, fields->total_cells);
  add_source(fields->vel_y, fields->vel_y_p, dt, fields->total_cells);

  // Smooth (viscosity) each velocity component.
  swap_fields(&fields->vel_x_p, &fields->vel_x);
  diffuse(BND_VEL_X, fields->vel_x, fields->vel_x_p, viscosity, dt);
  swap_fields(&fields->vel_y_p, &fields->vel_y);
  diffuse(BND_VEL_Y, fields->vel_y, fields->vel_y_p, viscosity, dt);

  // First cleanup: remove whatever ballooning the smoothing step above
  // just introduced.
  project(fields->vel_x, fields->vel_y, fields->pressure, fields->divergence);

  // Let velocity carry itself along its own flow (self-advection). Note
  // vel_x_p/vel_y_p are passed as BOTH the "old field being moved" AND the
  // "flow doing the moving" -- velocity pushes itself around, unlike ink
  // which gets pushed around by a separate field.
  swap_fields(&fields->vel_x_p, &fields->vel_x);
  swap_fields(&fields->vel_y_p, &fields->vel_y);
  advect(BND_VEL_X, fields->vel_x, fields->vel_x_p, fields->vel_x_p,
         fields->vel_y_p, dt);
  advect(BND_VEL_Y, fields->vel_y, fields->vel_y_p, fields->vel_x_p,
         fields->vel_y_p, dt);

  // Second cleanup: remove whatever ballooning self-advection just
  // introduced. After this, velocity is done for the frame.
  project(fields->vel_x, fields->vel_y, fields->pressure, fields->divergence);
}
