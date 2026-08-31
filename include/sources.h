#ifndef SOURCES_H
#define SOURCES_H

#include "fields.h"

/* ===========================================================================
 * fuentes.h
 * Fuentes de tinta: inicializacion, inyeccion por frame y disipacion.
 * ======================================================================== */

/*
 * Radio (en celdas) del vecindario de inyeccion, calibrado para una malla
 * de referencia de FUENTE_MALLA_REF celdas y escalado con la resolucion
 * real (ver fuente_radio/fuente_sigma) para que la mancha ocupe siempre la
 * misma fraccion de pantalla sin importar cuantas celdas tenga la malla:
 * con el radio fijo, una malla mas fina reduce cuantos pixeles de pantalla
 * ocupa ese mismo puñado de celdas, y el nacimiento de la tinta se ve como
 * "un puñado de pixeles" en vez de una mancha suave.
 *
 * El peso gaussiano de cada celda se divide entre el cuadrado del factor
 * de escala (ver inyectar_fuentes): el area del vecindario (~radio*sigma)
 * crece con la escala, y sin esa correccion la tinta total depositada por
 * fuente crecería con ella y saturaría la pantalla de blanco (ya paso una
 * vez: ver historial). Dividir entre escala^2 mantiene el total invariante
 * y solo mejora cuantas celdas dibujan el circulo.
 */
#define FUENTE_RADIO_BASE   2
#define FUENTE_SIGMA_BASE   1.1f
#define FUENTE_MALLA_REF    128

static inline float fuente_escala(int resolucion)
{
    float escala = (float)resolucion / (float)FUENTE_MALLA_REF;
    return (escala < 1.0f) ? 1.0f : escala;
}

static inline int fuente_radio(int resolucion)
{
    float escala = fuente_escala(resolucion);
    int   radio  = (int)(FUENTE_RADIO_BASE * escala + 0.5f);
    return (radio < FUENTE_RADIO_BASE) ? FUENTE_RADIO_BASE : radio;
}

static inline float fuente_sigma(int resolucion)
{
    return FUENTE_SIGMA_BASE * fuente_escala(resolucion);
}

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
 * La tinta se deposita en un vecindario de celdas con caida gaussiana para
 * suavizar la inyeccion y evitar valores puntuales muy bruscos.
 *
 * "aspecto" es ancho_ventana/alto_ventana. La malla es siempre cuadrada
 * (malla_n x malla_n) pero la ventana no tiene por que serlo: render.c usa
 * una escala distinta en x y en y para llenar la ventana, asi que una celda
 * de la malla no ocupa un cuadrado de pantalla sino un rectangulo de
 * "aspecto" veces mas ancho que alto. Sin corregir esto, un nucleo gaussiano
 * isotropico (mismo sigma en x y en y, en unidades de celda) se ve en
 * pantalla como una elipse mas ancha que alta. Aqui se encoge el sigma en x
 * por ese mismo factor para que, tras el escalado anisotropico de render.c,
 * la mancha vuelva a verse circular.
 */
void inyectar_fuentes(FuenteTinta *fuentes, int cantidad, CamposFluido *campos,
                      float aspecto);

/*
 * Multiplica la tinta por un factor menor que 1 para que se desvanezca poco a
 * poco; sin esto la pantalla terminaria saturada de blanco.
 */
void disipar_tinta(CamposFluido *campos);

#endif /* SOURCES_H */
