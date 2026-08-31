#include "configuracion.h"
#include "comun.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <time.h>

/* Imprime el modo de uso del programa. */
static void imprimir_ayuda(const char *nombre_programa)
{
    printf(
        "Screensaver de fluidos (Navier-Stokes) - version secuencial\n\n"
        "Uso: %s [opciones]\n\n"
        "Opciones:\n"
        "  -n <N>       Resolucion de la malla: N x N celdas   [%d..%d] (def. %d)\n"
        "  -f <F>       Cantidad de fuentes de tinta           [%d..%d] (def. %d)\n"
        "  -W <ancho>   Ancho de la ventana en pixeles         (min. %d, def. %d)\n"
        "  -H <alto>    Alto de la ventana en pixeles          (min. %d, def. %d)\n"
        "  -s <semilla> Semilla pseudoaleatoria                (def. reloj del sistema)\n"
        "  -t <dt>      Paso de tiempo de la simulacion        (def. %.3f)\n"
        "  -v <visc>    Viscosidad cinematica del fluido       (def. %.5f)\n"
        "  -d <diff>    Coeficiente de difusion de la tinta    (def. %.5f)\n"
        "  -b           Activa el sistema de n-cuerpos: las fuentes se mueven\n"
        "               por atraccion gravitacional mutua  (def. desactivado)\n"
        "  -p, -F       Pantalla completa (tamano exacto de la pantalla actual)\n"
        "               (def. desactivado)\n"
        "  -h           Muestra esta ayuda y termina\n\n"
        "Controles en ejecucion:\n"
        "  ESC / Q      Salir\n"
        "  R            Reiniciar la simulacion\n\n",
        nombre_programa,
        MALLA_MIN, MALLA_MAX, MALLA_DEFAULT,
        FUENTES_MIN, FUENTES_MAX, FUENTES_DEFAULT,
        VENTANA_ANCHO_MIN, VENTANA_ANCHO_DEF,
        VENTANA_ALTO_MIN, VENTANA_ALTO_DEF,
        DT_DEFAULT, VISC_DEFAULT, DIFF_DEFAULT);
}

/*
 * Convierte una cadena a entero validando que sea numerica completa y que
 * caiga dentro de [minimo, maximo].
 * Retorna 1 si la conversion fue exitosa, 0 en caso de error (e imprime el
 * motivo en stderr).
 */
static int leer_entero(const char *texto, const char *nombre_opcion,
                       int minimo, int maximo, int *destino)
{
    char *sobrante = NULL;
    long  valor;

    if (texto == NULL || *texto == '\0') {
        fprintf(stderr, "Error: la opcion %s requiere un valor.\n", nombre_opcion);
        return 0;
    }

    errno = 0;
    valor = strtol(texto, &sobrante, 10);

    if (*sobrante != '\0') {
        fprintf(stderr, "Error: '%s' no es un entero valido para %s.\n",
                texto, nombre_opcion);
        return 0;
    }
    if (errno == ERANGE || valor < (long)INT_MIN || valor > (long)INT_MAX) {
        fprintf(stderr, "Error: el valor de %s esta fuera de rango.\n", nombre_opcion);
        return 0;
    }
    if (valor < (long)minimo || valor > (long)maximo) {
        fprintf(stderr, "Error: %s debe estar entre %d y %d (recibido %ld).\n",
                nombre_opcion, minimo, maximo, valor);
        return 0;
    }

    *destino = (int)valor;
    return 1;
}

/*
 * Convierte una cadena a flotante validando formato y rango.
 * Retorna 1 si la conversion fue exitosa, 0 en caso de error.
 */
static int leer_flotante(const char *texto, const char *nombre_opcion,
                         float minimo, float maximo, float *destino)
{
    char *sobrante = NULL;
    double valor;

    if (texto == NULL || *texto == '\0') {
        fprintf(stderr, "Error: la opcion %s requiere un valor.\n", nombre_opcion);
        return 0;
    }

    errno = 0;
    valor = strtod(texto, &sobrante);

    if (*sobrante != '\0') {
        fprintf(stderr, "Error: '%s' no es un numero valido para %s.\n",
                texto, nombre_opcion);
        return 0;
    }
    if (errno == ERANGE || valor < (double)minimo || valor > (double)maximo) {
        fprintf(stderr, "Error: %s debe estar entre %g y %g (recibido %g).\n",
                nombre_opcion, (double)minimo, (double)maximo, valor);
        return 0;
    }

    *destino = (float)valor;
    return 1;
}

/*
 * Recorre argv y llena la estructura de configuracion.
 * Retorna  1 si todo esta correcto,
 *          0 si hubo un error de argumentos,
 *         -1 si el usuario pidio la ayuda (-h): el programa debe terminar bien.
 */
int procesar_argumentos(int argc, char *argv[], Configuracion *config)
{
    int indice;

    /* Valores por defecto */
    config->malla_n       = MALLA_DEFAULT;
    config->num_fuentes   = FUENTES_DEFAULT;
    config->ventana_ancho = VENTANA_ANCHO_DEF;
    config->ventana_alto  = VENTANA_ALTO_DEF;
    config->semilla       = (unsigned int)time(NULL);
    config->dt            = DT_DEFAULT;
    config->viscosidad    = VISC_DEFAULT;
    config->difusion      = DIFF_DEFAULT;
    config->nbody         = 0;
    config->pantalla_completa = 0;

    for (indice = 1; indice < argc; indice++) {
        const char *opcion = argv[indice];

        if (strcmp(opcion, "-h") == 0 || strcmp(opcion, "--help") == 0) {
            imprimir_ayuda(argv[0]);
            return -1;
        }

        if (strcmp(opcion, "-b") == 0 || strcmp(opcion, "--nbody") == 0) {
            config->nbody = 1;
            continue;
        }

        if (strcmp(opcion, "-p") == 0 || strcmp(opcion, "-F") == 0 ||
            strcmp(opcion, "--fullscreen") == 0 || strcmp(opcion, "--pantalla-completa") == 0) {
            config->pantalla_completa = 1;
            continue;
        }

        /* Todas las demas opciones requieren un valor asociado */
        if (indice + 1 >= argc) {
            fprintf(stderr, "Error: falta el valor para la opcion '%s'.\n", opcion);
            return 0;
        }

        if (strcmp(opcion, "-n") == 0) {
            if (!leer_entero(argv[++indice], "-n", MALLA_MIN, MALLA_MAX,
                             &config->malla_n)) return 0;
        } else if (strcmp(opcion, "-f") == 0) {
            if (!leer_entero(argv[++indice], "-f", FUENTES_MIN, FUENTES_MAX,
                             &config->num_fuentes)) return 0;
        } else if (strcmp(opcion, "-W") == 0) {
            if (!leer_entero(argv[++indice], "-W", VENTANA_ANCHO_MIN, 7680,
                             &config->ventana_ancho)) return 0;
        } else if (strcmp(opcion, "-H") == 0) {
            if (!leer_entero(argv[++indice], "-H", VENTANA_ALTO_MIN, 4320,
                             &config->ventana_alto)) return 0;
        } else if (strcmp(opcion, "-s") == 0) {
            int semilla_leida;
            if (!leer_entero(argv[++indice], "-s", 0, INT_MAX,
                             &semilla_leida)) return 0;
            config->semilla = (unsigned int)semilla_leida;
        } else if (strcmp(opcion, "-t") == 0) {
            if (!leer_flotante(argv[++indice], "-t", 0.001f, 1.0f,
                               &config->dt)) return 0;
        } else if (strcmp(opcion, "-v") == 0) {
            if (!leer_flotante(argv[++indice], "-v", 0.0f, 1.0f,
                               &config->viscosidad)) return 0;
        } else if (strcmp(opcion, "-d") == 0) {
            if (!leer_flotante(argv[++indice], "-d", 0.0f, 1.0f,
                               &config->difusion)) return 0;
        } else {
            fprintf(stderr, "Error: opcion desconocida '%s'.\n", opcion);
            fprintf(stderr, "Ejecute '%s -h' para ver el modo de uso.\n", argv[0]);
            return 0;
        }
    }

    return 1;
}
