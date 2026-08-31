#include "fondo.h"
#include "comun.h"
#include "utilidades.h"

#include <math.h>

void inicializar_fondo(Estrella *estrellas, int cantidad, int ancho, int alto)
{
    int e;

    for (e = 0; e < cantidad; e++) {
        estrellas[e].x = aleatorio_rango(0.0f, (float)(ancho - 1));
        estrellas[e].y = aleatorio_rango(0.0f, (float)(alto - 1));

        /* La mayoria quedan tenues (fondo) y solo unas pocas brillan fuerte
         * (primer plano), imitando la distribucion real de un cielo
         * estrellado en vez de un brillo uniforme. */
        estrellas[e].brillo_base = aleatorio_rango(0.15f, 1.0f) *
                                    aleatorio_rango(0.15f, 1.0f);

        estrellas[e].fase     = aleatorio_rango(0.0f, 2.0f * PI);
        estrellas[e].vel_fase = aleatorio_rango(0.01f, 0.04f);

        /* La gran mayoria son de 1 pixel; muy pocas de 2, como si estuvieran
         * un poco mas cerca. */
        estrellas[e].radio = (aleatorio_rango(0.0f, 1.0f) > 0.92f) ? 2 : 1;
    }
}

void actualizar_fondo(Estrella *estrellas, int cantidad)
{
    int e;

    for (e = 0; e < cantidad; e++) {
        estrellas[e].fase += estrellas[e].vel_fase;
        if (estrellas[e].fase > 2.0f * PI) {
            estrellas[e].fase -= 2.0f * PI;
        }
    }
}

void dibujar_fondo(SDL_Renderer *renderizador, const Estrella *estrellas,
                   int cantidad)
{
    int e;

    /* Titileo suave: oscila entre ~30% y 100% del brillo base en vez de
     * apagarse del todo, para que parpadee sin llegar a desaparecer. */
    for (e = 0; e < cantidad; e++) {
        float titileo = 0.65f + 0.35f * sinf(estrellas[e].fase);
        float brillo  = estrellas[e].brillo_base * titileo;
        Uint8 valor   = (Uint8)acotar(brillo * 255.0f, 0.0f, 255.0f);

        /* Blanco ligeramente azulado en vez de blanco puro, para que
         * combine con el tono frio del color base del fondo. */
        SDL_SetRenderDrawColor(renderizador,
                               (Uint8)(valor * 0.85f),
                               (Uint8)(valor * 0.92f),
                               valor, 255);

        if (estrellas[e].radio <= 1) {
            SDL_RenderDrawPoint(renderizador, (int)estrellas[e].x,
                                (int)estrellas[e].y);
        } else {
            SDL_Rect r;
            r.x = (int)estrellas[e].x;
            r.y = (int)estrellas[e].y;
            r.w = estrellas[e].radio;
            r.h = estrellas[e].radio;
            SDL_RenderFillRect(renderizador, &r);
        }
    }
}
