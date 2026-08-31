#ifndef FIELDS_H
#define FIELDS_H

/* ===========================================================================
 * campos.h
 * Reserva, liberacion y limpieza de los campos de la simulacion.
 * ======================================================================== */

/*
 * Estructura: CamposFluido
 * Todos los campos escalares y vectoriales de la simulacion. Cada campo tiene
 * un buffer "actual" y uno "previo" que se intercambian entre operadores.
 * La tinta se lleva en tres canales independientes (RGB) para poder mezclar
 * colores en pantalla.
 */
typedef struct {
    float *vel_x,   *vel_y;      /* campo de velocidad actual                 */
    float *vel_x_p, *vel_y_p;    /* velocidad previa / acumulador de fuerzas  */
    float *tinta_r, *tinta_g, *tinta_b;      /* densidad de tinta por canal   */
    float *tinta_r_p, *tinta_g_p, *tinta_b_p;/* buffers previos de tinta      */
    float *presion, *divergencia;/* auxiliares del paso de proyeccion         */
    int    celdas_total;         /* (N+2)*(N+2)                               */
} CamposFluido;

/*
 * Reserva e inicializa en cero todos los campos de la simulacion.
 * Retorna 1 si toda la memoria se obtuvo correctamente, 0 si alguna reserva
 * fallo (en cuyo caso libera lo ya reservado antes de retornar).
 */
int reservar_campos(CamposFluido *campos, int resolucion);

/* Libera todos los campos y deja los punteros en NULL. */
void liberar_campos(CamposFluido *campos);

/* Pone todos los campos en cero (usado al reiniciar con la tecla R). */
void limpiar_campos(CamposFluido *campos);

#endif /* FIELDS_H */
