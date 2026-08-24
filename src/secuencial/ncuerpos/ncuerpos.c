#include "ncuerpos.h"
#include "comun.h"
#include "utilidades.h"

#include <math.h>

/*
 * NBODY_G:           constante de gravitacion (unidades arbitrarias, ajustada
 *                     a mano junto con las masas de inicializar_fuentes() para
 *                     que el movimiento se note pero no sea caotico).
 * NBODY_SUAVIZADO:   distancia minima efectiva entre dos fuentes; evita que
 *                     la fuerza (proporcional a 1/dist^2) se dispare a
 *                     infinito en un encuentro cercano.
 * NBODY_VEL_MAX:     velocidad maxima por eje (celdas/frame); limita el
 *                     "efecto honda" de un encuentro cercano para que la
 *                     integracion explicita se mantenga estable.
 * NBODY_RESTITUCION: fraccion de la velocidad que se conserva al rebotar
 *                     contra un muro (<1 = rebote inelastico, se amortigua).
 */
#define NBODY_G            0.3f
#define NBODY_SUAVIZADO    6.0f
#define NBODY_VEL_MAX      5.0f
#define NBODY_RESTITUCION  0.9f

/*
 * Avanza un frame el sistema de n-cuerpos formado por las fuentes: cada
 * fuente atrae a las demas segun la ley de gravitacion universal
 * (F = G*m1*m2/d^2, suavizada) y rebota elasticamente contra los bordes de
 * la malla para permanecer siempre visible en pantalla.
 *
 * El "frame" se usa como unidad de tiempo (igual que vel_angular, que ya se
 * expresa en rad/frame), asi que no depende del dt de la simulacion de
 * fluidos: el movimiento de las fuentes es un fenomeno aparte que solo
 * comparte la malla con el solver.
 */
void actualizar_fuentes_nbody(FuenteTinta *fuentes, int cantidad, int resolucion)
{
    float acel_x[FUENTES_MAX] = { 0 };
    float acel_y[FUENTES_MAX] = { 0 };
    int   i, j;
    /* Las fuentes no pueden salir del vecindario de inyeccion (ver
     * FUENTE_RADIO, definido en fuentes.h) sin perder celdas interiores
     * validas. */
    const float limite_inf = (float)(FUENTE_RADIO + 1);
    const float limite_sup = (float)(resolucion - FUENTE_RADIO);

    for (i = 0; i < cantidad; i++) {
        for (j = i + 1; j < cantidad; j++) {
            float dif_x    = fuentes[j].pos_x - fuentes[i].pos_x;
            float dif_y    = fuentes[j].pos_y - fuentes[i].pos_y;
            float dist2    = dif_x * dif_x + dif_y * dif_y +
                              NBODY_SUAVIZADO * NBODY_SUAVIZADO;
            float inv_dist = 1.0f / sqrtf(dist2);
            /* factor comun G/d^3 (suavizada); multiplicado por la masa del
             * otro cuerpo da la aceleracion segun la 2a ley de Newton, y la
             * 3a ley (accion-reaccion) evita recalcular el par simetrico. */
            float factor   = NBODY_G * inv_dist * inv_dist * inv_dist;

            acel_x[i] += factor * fuentes[j].masa * dif_x;
            acel_y[i] += factor * fuentes[j].masa * dif_y;
            acel_x[j] -= factor * fuentes[i].masa * dif_x;
            acel_y[j] -= factor * fuentes[i].masa * dif_y;
        }
    }

    for (i = 0; i < cantidad; i++) {
        fuentes[i].vel_x = acotar(fuentes[i].vel_x + acel_x[i],
                                  -NBODY_VEL_MAX, NBODY_VEL_MAX);
        fuentes[i].vel_y = acotar(fuentes[i].vel_y + acel_y[i],
                                  -NBODY_VEL_MAX, NBODY_VEL_MAX);

        fuentes[i].pos_x += fuentes[i].vel_x;
        fuentes[i].pos_y += fuentes[i].vel_y;

        /* Rebote elastico (amortiguado) contra los bordes de la malla: la
         * fuente queda confinada a la caja en vez de escapar de la pantalla. */
        if (fuentes[i].pos_x < limite_inf) {
            fuentes[i].pos_x = limite_inf;
            fuentes[i].vel_x = -fuentes[i].vel_x * NBODY_RESTITUCION;
        } else if (fuentes[i].pos_x > limite_sup) {
            fuentes[i].pos_x = limite_sup;
            fuentes[i].vel_x = -fuentes[i].vel_x * NBODY_RESTITUCION;
        }

        if (fuentes[i].pos_y < limite_inf) {
            fuentes[i].pos_y = limite_inf;
            fuentes[i].vel_y = -fuentes[i].vel_y * NBODY_RESTITUCION;
        } else if (fuentes[i].pos_y > limite_sup) {
            fuentes[i].pos_y = limite_sup;
            fuentes[i].vel_y = -fuentes[i].vel_y * NBODY_RESTITUCION;
        }
    }
}
