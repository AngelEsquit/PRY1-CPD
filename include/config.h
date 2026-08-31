#ifndef CONFIG_H
#define CONFIG_H

/* ===========================================================================
 * config.h
 * Reading and validation of command-line arguments.
 * ======================================================================== */

/*
 * Struct: Config
 * Groups every parameter read from the command line.
 */
typedef struct {
    int          grid_n;          /* interior cells per side (N)               */
    int          num_sources;     /* number of ink sources                     */
    int          window_width;
    int          window_height;
    int          fullscreen;      /* 1 = exact size of the current screen      */
    unsigned int seed;            /* PRNG seed                                 */
    float        viscosity;
    float        diffusion;
} Config;

/*
 * Walks argv and fills in the config struct.
 * Returns  1 if everything is correct,
 *          0 if there was an argument error,
 *         -1 if the user asked for help (-h): the program should exit cleanly.
 */
int parse_arguments(int argc, char *argv[], Config *config);

#endif /* CONFIG_H */
