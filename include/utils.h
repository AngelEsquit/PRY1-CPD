#ifndef UTILS_H
#define UTILS_H

/* ===========================================================================
 * utils.h
 * General-purpose functions used by several modules.
 * ======================================================================== */

/* Returns a uniform pseudo-random float in [min, max]. */
float random_range(float min, float max);

/* Clamps a float value to the interval [min, max]. */
float clamp(float value, float min, float max);

#endif /* UTILS_H */
