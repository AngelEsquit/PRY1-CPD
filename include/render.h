#ifndef RENDER_H
#define RENDER_H

#include <SDL2/SDL.h>
#include "campos.h"

/* ===========================================================================
 * render.h
 * Volcado de la densidad de tinta a una textura SDL para dibujarla.
 * ======================================================================== */

/*
 * Vuelca la densidad de tinta a la textura de SDL, muestreando el campo N x N
 * con interpolacion bilineal propia a la resolucion de la textura destino
 * (ancho_textura x alto_textura, normalmente el tamano de la ventana).
 *
 * No se delega el escalado a SDL_RenderCopy porque el renderizador por
 * software de SDL2 ignora SDL_HINT_RENDER_SCALE_QUALITY y siempre escala con
 * el vecino mas cercano: eso se ve bien lejos de una fuente (la tinta ya esta
 * repartida en muchas celdas con diferencias pequenas entre si) pero se ve
 * pixelado justo donde nace (la densidad cae bruscamente en 2-3 celdas).
 */
void renderizar_tinta(SDL_Texture *textura, const CamposFluido *campos,
                      int ancho_textura, int alto_textura);

#endif /* RENDER_H */
