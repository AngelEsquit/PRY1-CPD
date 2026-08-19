#include "render.h"
#include "comun.h"
#include "utilidades.h"

#include <stddef.h>
#include <stdio.h>

/*
 * Vuelca la densidad de tinta a la textura de SDL. Cada celda interior de la
 * malla se convierte en un pixel de una textura de N x N que luego SDL escala
 * al tamano de la ventana con filtrado lineal.
 */
void renderizar_tinta(SDL_Texture *textura, const CamposFluido *campos)
{
    void  *pixeles_crudos = NULL;
    int    pitch_bytes    = 0;
    Uint32 *pixeles;
    int    i, j;
    float  valor_r, valor_g, valor_b, promedio;

    if (SDL_LockTexture(textura, NULL, &pixeles_crudos, &pitch_bytes) != 0) {
        fprintf(stderr, "Advertencia: no se pudo bloquear la textura: %s\n",
                SDL_GetError());
        return;
    }

    pixeles = (Uint32 *)pixeles_crudos;

    for (j = 0; j < malla_n; j++) {
        Uint32 *fila = pixeles + (size_t)j * ((size_t)pitch_bytes / sizeof(Uint32));
        for (i = 0; i < malla_n; i++) {
            /* La celda interior (i+1, j+1) corresponde al pixel (i, j).
             * Se aplica un mapeo tipo Reinhard (x / (x + k)) en lugar de un
             * recorte duro, para que la acumulacion de muchas fuentes se
             * comprima suavemente hacia blanco en vez de saturar de golpe. */
            valor_r = campos->tinta_r[IX(i + 1, j + 1)] * BRILLO_FACTOR;
            valor_g = campos->tinta_g[IX(i + 1, j + 1)] * BRILLO_FACTOR;
            valor_b = campos->tinta_b[IX(i + 1, j + 1)] * BRILLO_FACTOR;

            valor_r = valor_r / (valor_r + CONTRASTE_FACTOR);
            valor_g = valor_g / (valor_g + CONTRASTE_FACTOR);
            valor_b = valor_b / (valor_b + CONTRASTE_FACTOR);

            /* Empuja cada canal lejos del gris promedio del pixel para
             * recuperar saturacion sin alterar el brillo medio, evitando
             * que el resultado se vea lavado/blanco. */
            promedio = (valor_r + valor_g + valor_b) / 3.0f;
            valor_r = acotar(promedio + (valor_r - promedio) * SATURACION_FACTOR, 0.0f, 1.0f);
            valor_g = acotar(promedio + (valor_g - promedio) * SATURACION_FACTOR, 0.0f, 1.0f);
            valor_b = acotar(promedio + (valor_b - promedio) * SATURACION_FACTOR, 0.0f, 1.0f);

            fila[i] = ((Uint32)0xFF << 24) |
                      ((Uint32)(valor_r * 255.0f) << 16) |
                      ((Uint32)(valor_g * 255.0f) <<  8) |
                      ((Uint32)(valor_b * 255.0f));
        }
    }

    SDL_UnlockTexture(textura);
}
