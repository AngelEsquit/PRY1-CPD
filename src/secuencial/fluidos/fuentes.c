#include "fuentes.h"
#include "comun.h"
#include "utilidades.h"

#include <string.h>
#include <math.h>

/*
 * Inicializa las fuentes con posicion, color, fase y fuerza pseudoaleatorias.
 * Las posiciones se mantienen alejadas del borde para que el chorro se
 * desarrolle antes de chocar con los muros.
 */
void inicializar_fuentes(FuenteTinta *fuentes, int cantidad, int resolucion)
{
    int f;
    const int margen = (resolucion / 8 > 2) ? resolucion / 8 : 2;

    for (f = 0; f < cantidad; f++) {
        fuentes[f].celda_x = margen +
            (int)aleatorio_rango(0.0f, (float)(resolucion - 2 * margen));
        fuentes[f].celda_y = margen +
            (int)aleatorio_rango(0.0f, (float)(resolucion - 2 * margen));

        /* Color saturado: un canal dominante y los otros dos parciales */
        fuentes[f].color_r = aleatorio_rango(0.15f, 1.0f);
        fuentes[f].color_g = aleatorio_rango(0.15f, 1.0f);
        fuentes[f].color_b = aleatorio_rango(0.15f, 1.0f);

        fuentes[f].fase        = aleatorio_rango(0.0f, 2.0f * PI);
        fuentes[f].vel_angular = aleatorio_rango(-0.045f, 0.045f);
        fuentes[f].fuerza      = aleatorio_rango(25.0f, 60.0f);
        fuentes[f].caudal      = aleatorio_rango(60.0f, 110.0f);
    }
}

/*
 * Deposita en los buffers "previos" (que actuan como termino fuente) la tinta
 * y la cantidad de movimiento de cada fuente para el frame actual.
 *
 * ELEMENTO TRIGONOMETRICO: la direccion del chorro de cada fuente rota en el
 * tiempo segun (cos(fase), sin(fase)), lo que genera vortices en espiral.
 * La tinta se deposita en un vecindario de 3x3 celdas para suavizar la
 * inyeccion y evitar valores puntuales muy bruscos.
 */
void inyectar_fuentes(FuenteTinta *fuentes, int cantidad, CamposFluido *campos)
{
    int f, dx, dy, celda_i, celda_j;
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

        for (dy = -1; dy <= 1; dy++) {
            for (dx = -1; dx <= 1; dx++) {
                celda_i = fuentes[f].celda_x + dx;
                celda_j = fuentes[f].celda_y + dy;

                /* Solo se escribe dentro de las celdas interiores */
                if (celda_i < 1 || celda_i > malla_n ||
                    celda_j < 1 || celda_j > malla_n) {
                    continue;
                }

                campos->vel_x_p[IX(celda_i, celda_j)] += direccion_x;
                campos->vel_y_p[IX(celda_i, celda_j)] += direccion_y;

                campos->tinta_r_p[IX(celda_i, celda_j)] +=
                    fuentes[f].caudal * fuentes[f].color_r;
                campos->tinta_g_p[IX(celda_i, celda_j)] +=
                    fuentes[f].caudal * fuentes[f].color_g;
                campos->tinta_b_p[IX(celda_i, celda_j)] +=
                    fuentes[f].caudal * fuentes[f].color_b;
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
