#include "utils.h"

#include <stdlib.h>

// Returns a uniform pseudo-random float in [min, max].
float random_range(float min, float max)
{
    return min + ((float)rand() / (float)RAND_MAX) * (max - min);
}

// Clamps a float value to the interval [min, max].
float clamp(float value, float min, float max)
{
    if (value < min) return min;
    if (value > max) return max;
    return value;
}
