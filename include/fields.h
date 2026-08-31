#ifndef FIELDS_H
#define FIELDS_H

/* ===========================================================================
 * fields.h
 * Allocation, freeing and clearing of the simulation fields.
 * ======================================================================== */

/*
 * Struct: FluidFields
 * All the scalar and vector fields of the simulation. Each field has a
 * "current" buffer and a "previous" one, swapped between operators. Ink is
 * carried in three independent channels (RGB) so colors can mix on screen.
 */
typedef struct {
    float *vel_x,   *vel_y;      /* current velocity field                    */
    float *vel_x_p, *vel_y_p;    /* previous velocity / force accumulator     */
    float *ink_r, *ink_g, *ink_b;        /* ink density per channel           */
    float *ink_r_p, *ink_g_p, *ink_b_p;  /* previous ink buffers              */
    float *pressure, *divergence;/* projection step scratch fields            */
    int    total_cells;          /* (N+2)*(N+2)                               */
} FluidFields;

/*
 * Allocates and zero-initializes every simulation field.
 * Returns 1 if all memory was obtained successfully, 0 if some allocation
 * failed (in which case whatever was already allocated is freed before
 * returning).
 */
int allocate_fields(FluidFields *fields, int resolution);

/* Frees every field and sets the pointers to NULL. */
void free_fields(FluidFields *fields);

/* Zeroes out every field (used when resetting with the R key). */
void clear_fields(FluidFields *fields);

#endif /* FIELDS_H */
