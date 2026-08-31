#include "config.h"
#include "common.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <time.h>

/* Prints the program's usage screen. */
static void print_help(const char *program_name)
{
    printf(
        "Fluid screensaver (Navier-Stokes) - sequential version\n\n"
        "Usage: %s [options]\n\n"
        "Options:\n"
        "  -n <N>       Grid resolution: N x N cells            [%d..%d] (def. %d)\n"
        "  -f <F>       Number of ink sources                   [%d..%d] (def. %d)\n"
        "  -W <width>   Window width in pixels                  (min. %d, def. %d)\n"
        "  -H <height>  Window height in pixels                 (min. %d, def. %d)\n"
        "  -s <seed>    Pseudo-random seed                      (def. system clock)\n"
        "  -v <visc>    Kinematic viscosity of the fluid        (def. %.5f)\n"
        "  -d <diff>    Ink diffusion coefficient                (def. %.5f)\n"
        "  -p, -F       Fullscreen (exact size of the current screen)\n"
        "               (def. off)\n"
        "  -h           Shows this help and exits\n\n"
        "Runtime controls:\n"
        "  ESC / Q      Quit\n"
        "  R            Reset the simulation\n\n",
        program_name,
        GRID_MIN, GRID_MAX, GRID_DEFAULT,
        SOURCES_MIN, SOURCES_MAX, SOURCES_DEFAULT,
        WINDOW_WIDTH_MIN, WINDOW_WIDTH_DEFAULT,
        WINDOW_HEIGHT_MIN, WINDOW_HEIGHT_DEFAULT,
        VISC_DEFAULT, DIFF_DEFAULT);
}

/*
 * Converts a string to an int, validating that it's a complete numeric
 * value and that it falls within [min, max].
 * Returns 1 if the conversion succeeded, 0 on error (and prints the reason
 * to stderr).
 */
static int read_int(const char *text, const char *option_name,
                    int min, int max, int *out)
{
    char *leftover = NULL;
    long  value;

    if (text == NULL || *text == '\0') {
        fprintf(stderr, "Error: option %s requires a value.\n", option_name);
        return 0;
    }

    errno = 0;
    value = strtol(text, &leftover, 10);

    if (*leftover != '\0') {
        fprintf(stderr, "Error: '%s' is not a valid integer for %s.\n",
                text, option_name);
        return 0;
    }
    if (errno == ERANGE || value < (long)INT_MIN || value > (long)INT_MAX) {
        fprintf(stderr, "Error: the value of %s is out of range.\n", option_name);
        return 0;
    }
    if (value < (long)min || value > (long)max) {
        fprintf(stderr, "Error: %s must be between %d and %d (got %ld).\n",
                option_name, min, max, value);
        return 0;
    }

    *out = (int)value;
    return 1;
}

/*
 * Converts a string to a float, validating format and range.
 * Returns 1 if the conversion succeeded, 0 on error.
 */
static int read_float(const char *text, const char *option_name,
                      float min, float max, float *out)
{
    char *leftover = NULL;
    double value;

    if (text == NULL || *text == '\0') {
        fprintf(stderr, "Error: option %s requires a value.\n", option_name);
        return 0;
    }

    errno = 0;
    value = strtod(text, &leftover);

    if (*leftover != '\0') {
        fprintf(stderr, "Error: '%s' is not a valid number for %s.\n",
                text, option_name);
        return 0;
    }
    if (errno == ERANGE || value < (double)min || value > (double)max) {
        fprintf(stderr, "Error: %s must be between %g and %g (got %g).\n",
                option_name, (double)min, (double)max, value);
        return 0;
    }

    *out = (float)value;
    return 1;
}

/*
 * Walks argv and fills in the config struct.
 * Returns  1 if everything is correct,
 *          0 if there was an argument error,
 *         -1 if the user asked for help (-h): the program should exit cleanly.
 */
int parse_arguments(int argc, char *argv[], Config *config)
{
    int index;

    /* Default values */
    config->grid_n         = GRID_DEFAULT;
    config->num_sources    = SOURCES_DEFAULT;
    config->window_width   = WINDOW_WIDTH_DEFAULT;
    config->window_height  = WINDOW_HEIGHT_DEFAULT;
    config->seed           = (unsigned int)time(NULL);
    config->viscosity      = VISC_DEFAULT;
    config->diffusion      = DIFF_DEFAULT;
    config->fullscreen     = 0;

    for (index = 1; index < argc; index++) {
        const char *option = argv[index];

        if (strcmp(option, "-h") == 0 || strcmp(option, "--help") == 0) {
            print_help(argv[0]);
            return -1;
        }

        if (strcmp(option, "-p") == 0 || strcmp(option, "-F") == 0 ||
            strcmp(option, "--fullscreen") == 0) {
            config->fullscreen = 1;
            continue;
        }

        /* Every other option requires an associated value */
        if (index + 1 >= argc) {
            fprintf(stderr, "Error: missing value for option '%s'.\n", option);
            return 0;
        }

        if (strcmp(option, "-n") == 0) {
            if (!read_int(argv[++index], "-n", GRID_MIN, GRID_MAX,
                          &config->grid_n)) return 0;
        } else if (strcmp(option, "-f") == 0) {
            if (!read_int(argv[++index], "-f", SOURCES_MIN, SOURCES_MAX,
                          &config->num_sources)) return 0;
        } else if (strcmp(option, "-W") == 0) {
            if (!read_int(argv[++index], "-W", WINDOW_WIDTH_MIN, 7680,
                          &config->window_width)) return 0;
        } else if (strcmp(option, "-H") == 0) {
            if (!read_int(argv[++index], "-H", WINDOW_HEIGHT_MIN, 4320,
                          &config->window_height)) return 0;
        } else if (strcmp(option, "-s") == 0) {
            int seed_read;
            if (!read_int(argv[++index], "-s", 0, INT_MAX,
                          &seed_read)) return 0;
            config->seed = (unsigned int)seed_read;
        } else if (strcmp(option, "-v") == 0) {
            if (!read_float(argv[++index], "-v", 0.0f, 1.0f,
                            &config->viscosity)) return 0;
        } else if (strcmp(option, "-d") == 0) {
            if (!read_float(argv[++index], "-d", 0.0f, 1.0f,
                            &config->diffusion)) return 0;
        } else {
            fprintf(stderr, "Error: unknown option '%s'.\n", option);
            fprintf(stderr, "Run '%s -h' to see usage.\n", argv[0]);
            return 0;
        }
    }

    return 1;
}
