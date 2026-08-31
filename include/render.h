#ifndef RENDER_H
#define RENDER_H

#include <SDL2/SDL.h>
#include "fields.h"

/* ===========================================================================
 * render.h
 * Dumps the ink density into an SDL texture to draw it.
 * ======================================================================== */

/*
 * Dumps the ink density into the SDL texture, sampling the N x N field with
 * its own bilinear interpolation at the destination texture's resolution
 * (texture_width x texture_height, normally the window size).
 *
 * Scaling isn't delegated to SDL_RenderCopy because SDL2's software
 * renderer ignores SDL_HINT_RENDER_SCALE_QUALITY and always scales with
 * nearest-neighbor: that looks fine far from a source (the ink is already
 * spread over many cells with small differences between them) but looks
 * pixelated right where it's born (density drops sharply within 2-3 cells).
 */
void render_ink(SDL_Texture *texture, const FluidFields *fields,
                int texture_width, int texture_height);

#endif /* RENDER_H */
