#ifndef FUENTES_H
#define FUENTES_H

#include "campos.h"

/* ===========================================================================
 * fuentes.h
 * Fuentes de tinta: inicializacion, inyeccion por frame y disipacion.
 * ======================================================================== */

/*
 * Radio (en celdas) del vecindario de inyeccion de cada fuente (ver
 * inyectar_fuentes). Es publico porque tambien delimita cuanto puede
 * acercarse una fuente al borde de la malla sin perder celdas interiores
 * validas (ver ncuerpos.c).
 */
#define FUENTE_RADIO  2

/*
 * Estructura: FuenteTinta
 * Emisor puntual que inyecta color y cantidad de movimiento en la malla.
 * Las fuentes se desplazan por la malla como un sistema de n-cuerpos: cada
 * una atrae gravitacionalmente a las demas (ver ncuerpos.h).
 */
typedef struct {
    float pos_x, pos_y;          /* posicion continua en la malla (celdas)    */
    float vel_x, vel_y;          /* velocidad de desplazamiento (celdas/frame)*/
    float masa;                  /* masa gravitacional (n-cuerpos)            */
    float color_r, color_g, color_b; /* color de la tinta en [0,1]            */
    float fase;                  /* fase angular actual (radianes)            */
    float vel_angular;           /* rapidez de giro del chorro (rad/frame)    */
    float fuerza;                /* magnitud de la velocidad inyectada        */
    float caudal;                /* cantidad de tinta inyectada por frame     */
} FuenteTinta;

/*
 * Inicializa las fuentes con posicion, velocidad, masa, color, fase y fuerza
 * pseudoaleatorias. Las posiciones se mantienen alejadas del borde para que
 * el chorro se desarrolle antes de chocar con los muros.
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
