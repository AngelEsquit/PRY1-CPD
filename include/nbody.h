#ifndef NBODY_H
#define NBODY_H

#include "sources.h"

/* ===========================================================================
 * nbody.h
 * N-body system that drives the movement of the ink sources.
 * ======================================================================== */

/*
 * Advances one frame of the n-body system formed by the sources: each
 * source attracts the others according to the law of universal gravitation
 * (softened to avoid infinite forces on close encounters), and bounces
 * elastically off the grid's edges to stay always visible on screen.
 */
void update_nbody_sources(InkSource *sources, int count, int resolution);

#endif /* NBODY_H */
