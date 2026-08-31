#ifndef CONFIGURACION_H
#define CONFIGURACION_H

/* ===========================================================================
 * configuracion.h
 * Lectura y validacion de los argumentos de linea de comandos.
 * ======================================================================== */

/*
 * Estructura: Configuracion
 * Agrupa todos los parametros leidos de la linea de comandos.
 */
typedef struct {
    int          malla_n;        /* celdas interiores por lado (N)            */
    int          num_fuentes;    /* cantidad de fuentes de tinta              */
    int          ventana_ancho;
    int          ventana_alto;
    unsigned int semilla;        /* semilla del generador pseudoaleatorio     */
    float        dt;
    float        viscosidad;
    float        difusion;
    int          nbody;          /* 1 = las fuentes se mueven por gravitacion */
} Configuracion;

/*
 * Recorre argv y llena la estructura de configuracion.
 * Retorna  1 si todo esta correcto,
 *          0 si hubo un error de argumentos,
 *         -1 si el usuario pidio la ayuda (-h): el programa debe terminar bien.
 */
int procesar_argumentos(int argc, char *argv[], Configuracion *config);

#endif /* CONFIGURACION_H */
