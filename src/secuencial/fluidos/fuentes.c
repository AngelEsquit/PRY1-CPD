#include "fuentes.h"
#include "comun.h"
#include "utilidades.h"

#include <string.h>
#include <math.h>

/*
 * Inicializa las fuentes con posicion, velocidad, masa, color, fase y fuerza
 * pseudoaleatorias. Las posiciones se mantienen alejadas del borde para que
 * el chorro se desarrolle antes de chocar con los muros.
 */
void inicializar_fuentes(FuenteTinta *fuentes, int cantidad, int resolucion)
{
    int f;
    const int margen = (resolucion / 8 > 2) ? resolucion / 8 : 2;

    for (f = 0; f < cantidad; f++) {
        fuentes[f].pos_x = (float)margen +
            aleatorio_rango(0.0f, (float)(resolucion - 2 * margen));
        fuentes[f].pos_y = (float)margen +
            aleatorio_rango(0.0f, (float)(resolucion - 2 * margen));

        /* Velocidad inicial pequena: el movimiento principal lo genera la
         * atraccion gravitacional mutua una vez arranca la simulacion. */
        fuentes[f].vel_x = aleatorio_rango(-0.3f, 0.3f);
        fuentes[f].vel_y = aleatorio_rango(-0.3f, 0.3f);
        fuentes[f].masa  = aleatorio_rango(25.0f, 70.0f);

        /* Color saturado: un canal dominante y los otros dos parciales */
        fuentes[f].color_r = aleatorio_rango(0.15f, 1.0f);
        fuentes[f].color_g = aleatorio_rango(0.15f, 1.0f);
        fuentes[f].color_b = aleatorio_rango(0.15f, 1.0f);

        fuentes[f].fase        = aleatorio_rango(0.0f, 2.0f * PI);
        fuentes[f].vel_angular = aleatorio_rango(-0.045f, 0.045f);
        /* Con fuerza en [25,60] el desplazamiento por adveccion junto a la
         * fuente (dt*malla_n*|v|) llega a ~15-19 celdas por frame: varias
         * veces el diametro del vecindario de inyeccion, asi que la tinta
         * "salta" en vez de fluir, y se ve como gotas cayendo en vez de un
         * chorro continuo. En [6,14] el salto queda en ~5-6 celdas, del
         * orden del propio vecindario de inyeccion, y frames consecutivos
         * se solapan en un trazo continuo. */
        fuentes[f].fuerza      = aleatorio_rango(6.0f, 14.0f);
        fuentes[f].caudal      = aleatorio_rango(60.0f, 110.0f);
    }
}

/*
 * Deposita en los buffers "previos" (que actuan como termino fuente) la tinta
 * y la cantidad de movimiento de cada fuente para el frame actual.
 *
 * ELEMENTO TRIGONOMETRICO: la direccion del chorro de cada fuente rota en el
 * tiempo segun (cos(fase), sin(fase)), lo que genera vortices en espiral.
 * La tinta se deposita en un vecindario de 5x5 celdas con caida gaussiana
 * (en vez de un valor uniforme) para que la mancha nazca ya redondeada: con
 * solo 3x3 celdas la forma se percibe como un rombo/pixelado apenas se crea,
 * porque no hay suficientes muestras para leerse como un circulo, y como la
 * difusion de tinta esta en 0 por defecto nada la suaviza despues salvo el
 * movimiento.
 */
#define FUENTE_SIGMA  1.1f

void inyectar_fuentes(FuenteTinta *fuentes, int cantidad, CamposFluido *campos)
{
    int f, dx, dy, celda_i, celda_j, centro_x, centro_y;
    float direccion_x, direccion_y;

    /* Los buffers fuente se limpian cada frame */
    const size_t bytes = (size_t)campos->celdas_total * sizeof(float);
    memset(campos->vel_x_p,   0, bytes);
    memset(campos->vel_y_p,   0, bytes);
    memset(campos->tinta_r_p, 0, bytes);
    memset(campos->tinta_g_p, 0, bytes);
    memset(campos->tinta_b_p, 0, bytes);

    for (f = 0; f < cantidad; f++) {
        /* Avance de la fase angular del chorro */
        fuentes[f].fase += fuentes[f].vel_angular;
        if (fuentes[f].fase > 2.0f * PI) {
            fuentes[f].fase -= 2.0f * PI;
        }

        direccion_x = cosf(fuentes[f].fase) * fuentes[f].fuerza;
        direccion_y = sinf(fuentes[f].fase) * fuentes[f].fuerza;

        /* La posicion es continua (la mueve el sistema de n-cuerpos); se
         * redondea a la celda mas cercana solo para ubicar el vecindario de
         * inyeccion. */
        centro_x = (int)(fuentes[f].pos_x + 0.5f);
        centro_y = (int)(fuentes[f].pos_y + 0.5f);

        for (dy = -FUENTE_RADIO; dy <= FUENTE_RADIO; dy++) {
            for (dx = -FUENTE_RADIO; dx <= FUENTE_RADIO; dx++) {
                /* Nucleo gaussiano normalizado al centro (peso=1 en dx=dy=0),
                 * cae suavemente hasta los bordes del vecindario 5x5. */
                float dist2 = (float)(dx * dx + dy * dy);
                float peso  = expf(-dist2 / (2.0f * FUENTE_SIGMA * FUENTE_SIGMA));

                celda_i = centro_x + dx;
                celda_j = centro_y + dy;

                /* Solo se escribe dentro de las celdas interiores */
                if (celda_i < 1 || celda_i > malla_n ||
                    celda_j < 1 || celda_j > malla_n) {
                    continue;
                }

                campos->vel_x_p[IX(celda_i, celda_j)] += direccion_x * peso;
                campos->vel_y_p[IX(celda_i, celda_j)] += direccion_y * peso;

                campos->tinta_r_p[IX(celda_i, celda_j)] +=
                    fuentes[f].caudal * fuentes[f].color_r * peso;
                campos->tinta_g_p[IX(celda_i, celda_j)] +=
                    fuentes[f].caudal * fuentes[f].color_g * peso;
                campos->tinta_b_p[IX(celda_i, celda_j)] +=
                    fuentes[f].caudal * fuentes[f].color_b * peso;
            }
        }
    }
}

/*
 * Multiplica la tinta por un factor menor que 1 para que se desvanezca poco a
 * poco; sin esto la pantalla terminaria saturada de blanco.
 */
void disipar_tinta(CamposFluido *campos)
{
    int i;
    for (i = 0; i < campos->celdas_total; i++) {
        campos->tinta_r[i] *= DISIPACION;
        campos->tinta_g[i] *= DISIPACION;
        campos->tinta_b[i] *= DISIPACION;
    }
}
