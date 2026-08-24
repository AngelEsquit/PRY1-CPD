#ifndef NCUERPOS_H
#define NCUERPOS_H

#include "fuentes.h"

/* ===========================================================================
 * ncuerpos.h
 * Sistema de n-cuerpos que gobierna el movimiento de las fuentes de tinta.
 * ======================================================================== */

/*
 * Avanza un frame el sistema de n-cuerpos formado por las fuentes: cada
 * fuente atrae a las demas segun la ley de gravitacion universal (suavizada
 * para evitar fuerzas infinitas en encuentros cercanos), y rebota
 * elasticamente contra los bordes de la malla para permanecer siempre
 * visible en pantalla.
 */
void actualizar_fuentes_nbody(FuenteTinta *fuentes, int cantidad, int resolucion);

#endif /* NCUERPOS_H */
