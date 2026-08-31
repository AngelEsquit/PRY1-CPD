#include "render.h"
#include "comun.h"
#include "utilidades.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>

/* Exponente de la curva gamma del alfa (ver renderizar_tinta): mientras mas
 * chico, mas rapido llega el alfa a opaco con poco brillo de tinta. */
#define ALFA_GAMMA 0.35f

/*
 * Interpola bilinealmente "campo" (malla N x N con anillo fantasma) en la
 * posicion continua (fx, fy), expresada en el mismo espacio que un indice de
 * pixel de la textura N x N original: fx=0 es el centro de la celda interior
 * 1, fx=malla_n-1 es el centro de la celda interior malla_n.
 */
static float muestrear_bilineal(const float *campo, float fx, float fy)
{
    int i0, i1, j0, j1;
    float peso_x1, peso_x0, peso_y1, peso_y0;

    fx = acotar(fx, 0.0f, (float)(malla_n - 1));
    fy = acotar(fy, 0.0f, (float)(malla_n - 1));

    i0 = (int)fx;  i1 = (i0 < malla_n - 1) ? i0 + 1 : i0;
    j0 = (int)fy;  j1 = (j0 < malla_n - 1) ? j0 + 1 : j0;

    peso_x1 = fx - (float)i0;  peso_x0 = 1.0f - peso_x1;
    peso_y1 = fy - (float)j0;  peso_y0 = 1.0f - peso_y1;

    /* +1 desplaza del espacio de pixel (0..N-1) al indice de celda interior
     * (1..N) dentro de la malla con anillo fantasma. */
    return peso_x0 * (peso_y0 * campo[IX(i0 + 1, j0 + 1)] +
                       peso_y1 * campo[IX(i0 + 1, j1 + 1)]) +
           peso_x1 * (peso_y0 * campo[IX(i1 + 1, j0 + 1)] +
                       peso_y1 * campo[IX(i1 + 1, j1 + 1)]);
}

/*
 * Vuelca la densidad de tinta a la textura de SDL, con un pixel de la
 * textura por cada pixel de destino (ver render.h para el porque de esto:
 * el renderizador por software de SDL2 no hace escalado lineal).
 */
void renderizar_tinta(SDL_Texture *textura, const CamposFluido *campos,
                      int ancho_textura, int alto_textura)
{
    void  *pixeles_crudos = NULL;
    int    pitch_bytes    = 0;
    Uint32 *pixeles;
    int    i, j;
    float  fx, fy;
    float  valor_r, valor_g, valor_b, promedio;
    const float escala_x = (float)malla_n / (float)ancho_textura;
    const float escala_y = (float)malla_n / (float)alto_textura;

    if (SDL_LockTexture(textura, NULL, &pixeles_crudos, &pitch_bytes) != 0) {
        fprintf(stderr, "Advertencia: no se pudo bloquear la textura: %s\n",
                SDL_GetError());
        return;
    }

    pixeles = (Uint32 *)pixeles_crudos;

    /* Cada pixel se calcula de forma independiente (solo lee los campos de
     * tinta, no los modifica). Nota: fy depende solo de j y no de i, asi que
     * con collapse(2) se recalcula por cada (i,j) en vez de una vez por
     * fila; es un recalculo barato (una resta y una multiplicacion) y a
     * cambio incluso mallas muy angostas reparten trabajo entre todos los
     * hilos. */
    for (j = 0; j < alto_textura; j++) {
        for (i = 0; i < ancho_textura; i++) {
            Uint32 *fila = pixeles + (size_t)j * ((size_t)pitch_bytes / sizeof(Uint32));
            fy = ((float)j + 0.5f) * escala_y - 0.5f;
            fx = ((float)i + 0.5f) * escala_x - 0.5f;

            /* Se aplica un mapeo tipo Reinhard (x / (x + k)) en lugar de un
             * recorte duro, para que la acumulacion de muchas fuentes se
             * comprima suavemente hacia blanco en vez de saturar de golpe. */
            valor_r = muestrear_bilineal(campos->tinta_r, fx, fy) * BRILLO_FACTOR;
            valor_g = muestrear_bilineal(campos->tinta_g, fx, fy) * BRILLO_FACTOR;
            valor_b = muestrear_bilineal(campos->tinta_b, fx, fy) * BRILLO_FACTOR;

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

            /* Alfa proporcional al brillo de la tinta (en vez de 0xFF fijo):
             * donde no hay tinta el pixel es transparente y deja ver el
             * fondo (estrellas/vinieta) dibujado antes en el renderer; donde
             * la tinta satura, el pixel es opaco y la tapa por completo. El
             * canal mas alto de los tres se usa como aproximacion de
             * brillo (evita que un color muy saturado en un solo canal, con
             * promedio bajo, se vea mas transparente de lo que deberia). */
            {
                float brillo_max, alfa;

                brillo_max = valor_r;
                if (valor_g > brillo_max) brillo_max = valor_g;
                if (valor_b > brillo_max) brillo_max = valor_b;

                /* Usar brillo_max directo como alfa deja gran parte de la
                 * mancha (todo lo que no este ya saturado a full brillo)
                 * semi-transparente, y SDL_BLENDMODE_BLEND la mezcla con el
                 * fondo oscuro de mas atras: el resultado se ve deslavado,
                 * como si el color de la tinta perdiera intensidad. Una
                 * curva gamma (<1) empuja el alfa a opaco mucho mas rapido
                 * que el brillo real, asi que la mancha se ve tan vivida
                 * como antes de agregar el fondo, y solo el borde donde la
                 * tinta cae casi a cero se queda transparente (para que ahi
                 * se vean las estrellas). */
                alfa = powf(brillo_max, ALFA_GAMMA);

                fila[i] = ((Uint32)(alfa * 255.0f) << 24) |
                          ((Uint32)(valor_r * 255.0f) << 16) |
                          ((Uint32)(valor_g * 255.0f) <<  8) |
                          ((Uint32)(valor_b * 255.0f));
            }
        }
    }

    SDL_UnlockTexture(textura);
}
