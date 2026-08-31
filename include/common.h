#ifndef COMMON_H
#define COMMON_H

/* ===========================================================================
 * common.h
 * ---------------------------------------------------------------------------
 * Constants, macros and values shared by every module of the fluid
 * screensaver (Navier-Stokes / Stable Fluids).
 * ======================================================================== */

/* ---------------------------------------------------------------------------
 * Configuration constants and validation limits (defensive programming)
 * ------------------------------------------------------------------------- */
#define GRID_MIN              16    /* minimum grid resolution                  */
#define GRID_MAX             1024    /* maximum grid resolution                  */
#define GRID_DEFAULT          256
#define SOURCES_MIN             1    /* at least one ink source                  */
#define SOURCES_MAX           256
#define SOURCES_DEFAULT          6
#define WINDOW_WIDTH_MIN       640    /* required by the project spec             */
#define WINDOW_HEIGHT_MIN      480
#define WINDOW_WIDTH_DEFAULT  1920
#define WINDOW_HEIGHT_DEFAULT 1080
#define GAUSS_SEIDEL_ITERS       20    /* iterations of the linear solver          */

/* M_PI is not part of the C11 standard (POSIX only), so it's defined here
 * for the program to compile portably with -std=c11.                        */
#define PI 3.14159265358979323846f

/* Default physical parameters (adjustable via command line) */
#define DT_DEFAULT          0.07f   /* time step                                */
#define BRIGHTNESS_FACTOR    0.6f   /* brightness factor applied to the render  */
#define CONTRAST_FACTOR      1.0f   /* saturates slower the higher it is (keeps
                                        ink from looking white where many
                                        sources overlap)                       */
#define SATURATION_FACTOR    2.4f   /* pushes each channel away from the
                                        pixel's average gray (>1 = more vivid);
                                        being symmetric around the average, it
                                        doesn't push the color toward white   */
#define VISC_DEFAULT        0.0000f /* kinematic viscosity of the fluid         */
#define DIFF_DEFAULT        0.0001f /* ink diffusion coefficient: without this
                                        (0.0) ink is only smoothed by the fixed
                                        shape of the injection kernel, and shows
                                        a hard/abrupt edge right as it's born;
                                        with diffusion on, each source is born
                                        already blurring, like real ink in water */
#define DISSIPATION         0.995f /* ink retained per frame (multiplicative
                                       factor); at 0.3 the ink lost 70% of its
                                       value every frame and faded out before it
                                       could be advected out of its origin cell */

/* ---------------------------------------------------------------------------
 * Grid indexing.
 * The grid has (N+2)x(N+2) cells: N x N interior cells (indices 1..N) plus a
 * ring of ghost cells at the border for the boundary conditions. Stored as a
 * 1D array in row-major order.
 *
 * "grid_n" is defined once in main.c and declared here as extern so the
 * IX() macro is available (and readable) in every module that operates on
 * the grid.
 * ------------------------------------------------------------------------- */
extern int grid_n;

#define IX(i, j) ((i) + (grid_n + 2) * (j))

/* Boundary codes used by apply_boundary() (see solver.c) */
#define BND_SCALAR   0  /* scalar field: copies the interior neighbor        */
#define BND_VEL_X    1  /* x velocity component: reflected at walls         */
#define BND_VEL_Y    2  /* y velocity component: reflected at walls         */

#endif /* COMMON_H */
