#ifndef SOLVER_H
#define SOLVER_H

#include "fields.h"

/* ===========================================================================
 * solver.h
 * ---------------------------------------------------------------------------
 * Nucleo numerico del metodo "Stable Fluids": difusion, adveccion,
 * proyeccion y los pasos completos de velocidad y tinta.
 *
 * Este es el modulo que concentra la mayor parte del costo computacional y,
 * por lo tanto, el candidato principal para paralelizar con OpenMP en la
 * siguiente fase del proyecto (ver nota sobre Gauss-Seidel en solver.c).
 * ======================================================================== */

/*
 * Paso completo de la densidad de tinta de un canal:
 * fuente -> difusion -> adveccion.
 */
void paso_tinta(float **tinta, float **tinta_previa,
                const float *vel_x, const float *vel_y,
                float difusion, float dt, int celdas_total);

/*
 * Paso completo del campo de velocidad:
 * fuerzas -> difusion viscosa -> proyeccion -> autoadveccion -> proyeccion.
 * Se proyecta dos veces porque tanto la difusion como la adveccion vuelven a
 * introducir divergencia en el campo.
 */
void paso_velocidad(CamposFluido *campos, float viscosidad, float dt);

#endif /* SOLVER_H */
