#ifndef FUENTES_H
#define FUENTES_H

#include "campos.h"

/* ===========================================================================
 * fuentes.h
 * Fuentes de tinta: inicializacion, inyeccion por frame y disipacion.
 * ======================================================================== */

/*
 * Estructura: FuenteTinta
 * Emisor puntual que inyecta color y cantidad de movimiento en la malla.
 */
typedef struct {
    int   celda_x, celda_y;      /* posicion de la fuente en la malla         */
    float color_r, color_g, color_b; /* color de la tinta en [0,1]            */
    float fase;                  /* fase angular actual (radianes)            */
    float vel_angular;           /* rapidez de giro del chorro (rad/frame)    */
    float fuerza;                /* magnitud de la velocidad inyectada        */
    float caudal;                /* cantidad de tinta inyectada por frame     */
} FuenteTinta;

/*
 * Inicializa las fuentes con posicion, color, fase y fuerza pseudoaleatorias.
 * Las posiciones se mantienen alejadas del borde para que el chorro se
 * desarrolle antes de chocar con los muros.
 */
void inicializar_fuentes(FuenteTinta *fuentes, int cantidad, int resolucion);

/*
 * Deposita en los buffers "previos" (que actuan como termino fuente) la tinta
 * y la cantidad de movimiento de cada fuente para el frame actual.
 *
 * ELEMENTO TRIGONOMETRICO: la direccion del chorro de cada fuente rota en el
 * tiempo segun (cos(fase), sin(fase)), lo que genera vortices en espiral.
 * La tinta se deposita en un vecindario de 3x3 celdas para suavizar la
 * inyeccion y evitar valores puntuales muy bruscos.
 */
void inyectar_fuentes(FuenteTinta *fuentes, int cantidad, CamposFluido *campos);

/*
 * Multiplica la tinta por un factor menor que 1 para que se desvanezca poco a
 * poco; sin esto la pantalla terminaria saturada de blanco.
 */
void disipar_tinta(CamposFluido *campos);

#endif /* FUENTES_H */
