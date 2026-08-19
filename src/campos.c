#include "campos.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Reserva e inicializa en cero todos los campos de la simulacion.
 * Retorna 1 si toda la memoria se obtuvo correctamente, 0 si alguna reserva
 * fallo (en cuyo caso libera lo ya reservado antes de retornar).
 */
int reservar_campos(CamposFluido *campos, int resolucion)
{
    int i;
    /* Punteros agrupados para reservar y validar en un solo ciclo */
    float **todos_los_campos[] = {
        &campos->vel_x,     &campos->vel_y,
        &campos->vel_x_p,   &campos->vel_y_p,
        &campos->tinta_r,   &campos->tinta_g,   &campos->tinta_b,
        &campos->tinta_r_p, &campos->tinta_g_p, &campos->tinta_b_p,
        &campos->presion,   &campos->divergencia
    };
    const int cantidad_campos = (int)(sizeof(todos_los_campos) /
                                      sizeof(todos_los_campos[0]));

    campos->celdas_total = (resolucion + 2) * (resolucion + 2);

    /* Inicializar en NULL para que liberar_campos() sea seguro ante fallos */
    for (i = 0; i < cantidad_campos; i++) {
        *todos_los_campos[i] = NULL;
    }

    for (i = 0; i < cantidad_campos; i++) {
        *todos_los_campos[i] = (float *)calloc((size_t)campos->celdas_total,
                                               sizeof(float));
        if (*todos_los_campos[i] == NULL) {
            fprintf(stderr,
                    "Error: memoria insuficiente para una malla de %dx%d.\n"
                    "Intente con un valor menor de -n.\n",
                    resolucion, resolucion);
            /* Liberar lo que si se pudo reservar */
            while (--i >= 0) {
                free(*todos_los_campos[i]);
                *todos_los_campos[i] = NULL;
            }
            return 0;
        }
    }

    return 1;
}

/* Libera todos los campos y deja los punteros en NULL. */
void liberar_campos(CamposFluido *campos)
{
    free(campos->vel_x);       campos->vel_x       = NULL;
    free(campos->vel_y);       campos->vel_y       = NULL;
    free(campos->vel_x_p);     campos->vel_x_p     = NULL;
    free(campos->vel_y_p);     campos->vel_y_p     = NULL;
    free(campos->tinta_r);     campos->tinta_r     = NULL;
    free(campos->tinta_g);     campos->tinta_g     = NULL;
    free(campos->tinta_b);     campos->tinta_b     = NULL;
    free(campos->tinta_r_p);   campos->tinta_r_p   = NULL;
    free(campos->tinta_g_p);   campos->tinta_g_p   = NULL;
    free(campos->tinta_b_p);   campos->tinta_b_p   = NULL;
    free(campos->presion);     campos->presion     = NULL;
    free(campos->divergencia); campos->divergencia = NULL;
}

/* Pone todos los campos en cero (usado al reiniciar con la tecla R). */
void limpiar_campos(CamposFluido *campos)
{
    const size_t bytes = (size_t)campos->celdas_total * sizeof(float);
    memset(campos->vel_x,       0, bytes);
    memset(campos->vel_y,       0, bytes);
    memset(campos->vel_x_p,     0, bytes);
    memset(campos->vel_y_p,     0, bytes);
    memset(campos->tinta_r,     0, bytes);
    memset(campos->tinta_g,     0, bytes);
    memset(campos->tinta_b,     0, bytes);
    memset(campos->tinta_r_p,   0, bytes);
    memset(campos->tinta_g_p,   0, bytes);
    memset(campos->tinta_b_p,   0, bytes);
    memset(campos->presion,     0, bytes);
    memset(campos->divergencia, 0, bytes);
}
