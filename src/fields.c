#include "fields.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Allocates and zero-initializes every simulation field.
 * Returns 1 if all memory was obtained successfully, 0 if some allocation
 * failed (in which case whatever was already allocated is freed before
 * returning).
 */
int allocate_fields(FluidFields *fields, int resolution)
{
    int i;
    /* Pointers grouped together to allocate and validate in a single loop */
    float **all_fields[] = {
        &fields->vel_x,   &fields->vel_y,
        &fields->vel_x_p, &fields->vel_y_p,
        &fields->ink_r,   &fields->ink_g,   &fields->ink_b,
        &fields->ink_r_p, &fields->ink_g_p, &fields->ink_b_p,
        &fields->pressure, &fields->divergence
    };
    const int field_count = (int)(sizeof(all_fields) / sizeof(all_fields[0]));

    fields->total_cells = (resolution + 2) * (resolution + 2);

    /* Initialize to NULL so free_fields() is safe on failure */
    for (i = 0; i < field_count; i++) {
        *all_fields[i] = NULL;
    }

    for (i = 0; i < field_count; i++) {
        *all_fields[i] = (float *)calloc((size_t)fields->total_cells,
                                         sizeof(float));
        if (*all_fields[i] == NULL) {
            fprintf(stderr,
                    "Error: not enough memory for a %dx%d grid.\n"
                    "Try a smaller value for -n.\n",
                    resolution, resolution);
            /* Free whatever was successfully allocated */
            while (--i >= 0) {
                free(*all_fields[i]);
                *all_fields[i] = NULL;
            }
            return 0;
        }
    }

    return 1;
}

/* Frees every field and sets the pointers to NULL. */
void free_fields(FluidFields *fields)
{
    free(fields->vel_x);      fields->vel_x      = NULL;
    free(fields->vel_y);      fields->vel_y      = NULL;
    free(fields->vel_x_p);    fields->vel_x_p    = NULL;
    free(fields->vel_y_p);    fields->vel_y_p    = NULL;
    free(fields->ink_r);      fields->ink_r      = NULL;
    free(fields->ink_g);      fields->ink_g      = NULL;
    free(fields->ink_b);      fields->ink_b      = NULL;
    free(fields->ink_r_p);    fields->ink_r_p    = NULL;
    free(fields->ink_g_p);    fields->ink_g_p    = NULL;
    free(fields->ink_b_p);    fields->ink_b_p    = NULL;
    free(fields->pressure);   fields->pressure   = NULL;
    free(fields->divergence); fields->divergence = NULL;
}

/* Zeroes out every field (used when resetting with the R key). */
void clear_fields(FluidFields *fields)
{
    const size_t bytes = (size_t)fields->total_cells * sizeof(float);
    memset(fields->vel_x,      0, bytes);
    memset(fields->vel_y,      0, bytes);
    memset(fields->vel_x_p,    0, bytes);
    memset(fields->vel_y_p,    0, bytes);
    memset(fields->ink_r,      0, bytes);
    memset(fields->ink_g,      0, bytes);
    memset(fields->ink_b,      0, bytes);
    memset(fields->ink_r_p,    0, bytes);
    memset(fields->ink_g_p,    0, bytes);
    memset(fields->ink_b_p,    0, bytes);
    memset(fields->pressure,   0, bytes);
    memset(fields->divergence, 0, bytes);
}
