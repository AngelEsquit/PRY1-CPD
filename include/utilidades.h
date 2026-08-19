#ifndef UTILIDADES_H
#define UTILIDADES_H

/* ===========================================================================
 * utilidades.h
 * Funciones de proposito general usadas por varios modulos.
 * ======================================================================== */

/* Devuelve un flotante pseudoaleatorio uniforme en [minimo, maximo]. */
float aleatorio_rango(float minimo, float maximo);

/* Limita un valor flotante al intervalo [minimo, maximo]. */
float acotar(float valor, float minimo, float maximo);

#endif /* UTILIDADES_H */
