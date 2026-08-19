#ifndef RENDER_H
#define RENDER_H

#include <SDL2/SDL.h>
#include "campos.h"

/* ===========================================================================
 * render.h
 * Volcado de la densidad de tinta a una textura SDL para dibujarla.
 * ======================================================================== */

/*
 * Vuelca la densidad de tinta a la textura de SDL. Cada celda interior de la
 * malla se convierte en un pixel de una textura de N x N que luego SDL escala
 * al tamano de la ventana con filtrado lineal.
 */
void renderizar_tinta(SDL_Texture *textura, const CamposFluido *campos);

#endif /* RENDER_H */
