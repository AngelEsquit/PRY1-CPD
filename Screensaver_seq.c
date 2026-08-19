/* ===========================================================================
 * screensaver_seq.c
 * ---------------------------------------------------------------------------
 * Screensaver de fluidos basado en las ecuaciones de Navier-Stokes para
 * fluidos incompresibles, resueltas con el metodo "Stable Fluids" de Jos Stam.
 *
 * VERSION SECUENCIAL (base de comparacion para la version paralela con OpenMP).
 *
 * Modelo fisico
 * -------------
 * Se resuelven las dos ecuaciones acopladas sobre una malla regular 2D:
 *
 *   Velocidad:  du/dt = -(u . grad)u + visc * lap(u) + f      con div(u) = 0
 *   Tinta:      dp/dt = -(u . grad)p + diff * lap(p) + s
 *
 * Cada paso de tiempo se descompone en cuatro operadores (Stam, 1999):
 *   1. add_source  -> aplica las fuerzas externas f y las fuentes de tinta s
 *   2. diffuse     -> resuelve la difusion viscosa de forma implicita
 *   3. advect      -> transporta el campo con trazado semi-Lagrangiano
 *   4. project     -> proyeccion de Hodge: elimina la divergencia (incompresible)
 *
 * Los pasos 2 y 4 requieren resolver un sistema lineal disperso; aqui se usa
 * relajacion iterativa de Gauss-Seidel (ver nota en lin_solve sobre su impacto
 * en la version paralela).
 *
 * Trigonometria: cada fuente de tinta inyecta velocidad en una direccion que
 * rota en el tiempo, calculada con sin() y cos() sobre una fase propia.
 *
 * Compilacion:  make            (ver Makefile)
 * Uso:          ./screensaver_seq -n 128 -f 6
 * ===========================================================================
 */

#include <SDL2/SDL.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <limits.h>

/* ---------------------------------------------------------------------------
 * Constantes de configuracion y limites de validacion (programacion defensiva)
 * ------------------------------------------------------------------------- */
#define MALLA_MIN            16    /* resolucion minima de la malla            */
#define MALLA_MAX          1024    /* resolucion maxima de la malla            */
#define MALLA_DEFAULT       128
#define FUENTES_MIN           1    /* al menos una fuente de tinta             */
#define FUENTES_MAX         256
#define FUENTES_DEFAULT       6
#define VENTANA_ANCHO_MIN   640    /* exigido por el enunciado del proyecto    */
#define VENTANA_ALTO_MIN    480
#define VENTANA_ANCHO_DEF   800
#define VENTANA_ALTO_DEF    600
#define ITER_GAUSS_SEIDEL    20    /* iteraciones del solver lineal            */

/* M_PI no forma parte del estandar C11 (solo de POSIX), se define aqui para
 * que el programa compile de forma portable con -std=c11.                    */
#define PI 3.14159265358979323846f

/* Parametros fisicos por defecto (ajustables por linea de comandos) */
#define DT_DEFAULT         0.07f   /* paso de tiempo                           */
#define BRILLO_FACTOR       0.6f   /* factor de brillo aplicado al render      */
#define CONTRASTE_FACTOR     1.0f   /* satura mas lento cuanto mas alto (evita
                                        que la acumulacion de tinta se vea
                                        blanca al superponerse muchas fuentes) */
#define SATURACION_FACTOR    1.6f   /* aleja cada canal del gris promedio del
                                        pixel (>1 = mas vivido); al ser
                                        simetrico respecto al promedio no
                                        empuja el color hacia blanco        */
#define VISC_DEFAULT       0.0000f /* viscosidad cinematica del fluido         */
#define DIFF_DEFAULT       0.0000f /* difusion de la tinta                     */
#define DISIPACION         0.995f  /* desvanecimiento de la tinta por frame    */

/* ---------------------------------------------------------------------------
 * Indexacion de la malla.
 * La malla tiene (N+2)x(N+2) celdas: N x N celdas interiores (indices 1..N)
 * mas un anillo de celdas fantasma en el borde para las condiciones de
 * frontera. Se almacena en un arreglo 1D en orden row-major.
 * ------------------------------------------------------------------------- */
#define IX(i, j) ((i) + (malla_n + 2) * (j))

/* Codigos de frontera usados por aplicar_frontera() */
#define BND_ESCALAR   0  /* campo escalar: se copia el vecino interior        */
#define BND_VEL_X     1  /* componente x de velocidad: se refleja en muros    */
#define BND_VEL_Y     2  /* componente y de velocidad: se refleja en muros    */

/* ---------------------------------------------------------------------------
 * Estructura: FuenteTinta
 * Emisor puntual que inyecta color y cantidad de movimiento en la malla.
 * ------------------------------------------------------------------------- */
typedef struct {
    int   celda_x, celda_y;      /* posicion de la fuente en la malla         */
    float color_r, color_g, color_b; /* color de la tinta en [0,1]            */
    float fase;                  /* fase angular actual (radianes)            */
    float vel_angular;           /* rapidez de giro del chorro (rad/frame)    */
    float fuerza;                /* magnitud de la velocidad inyectada        */
    float caudal;                /* cantidad de tinta inyectada por frame     */
} FuenteTinta;

/* ---------------------------------------------------------------------------
 * Estructura: Configuracion
 * Agrupa todos los parametros leidos de la linea de comandos.
 * ------------------------------------------------------------------------- */
typedef struct {
    int          malla_n;        /* celdas interiores por lado (N)            */
    int          num_fuentes;    /* cantidad de fuentes de tinta              */
    int          ventana_ancho;
    int          ventana_alto;
    unsigned int semilla;        /* semilla del generador pseudoaleatorio     */
    float        dt;
    float        viscosidad;
    float        difusion;
} Configuracion;

/* ---------------------------------------------------------------------------
 * Estructura: CamposFluido
 * Todos los campos escalares y vectoriales de la simulacion. Cada campo tiene
 * un buffer "actual" y uno "previo" que se intercambian entre operadores.
 * La tinta se lleva en tres canales independientes (RGB) para poder mezclar
 * colores en pantalla.
 * ------------------------------------------------------------------------- */
typedef struct {
    float *vel_x,   *vel_y;      /* campo de velocidad actual                 */
    float *vel_x_p, *vel_y_p;    /* velocidad previa / acumulador de fuerzas  */
    float *tinta_r, *tinta_g, *tinta_b;      /* densidad de tinta por canal   */
    float *tinta_r_p, *tinta_g_p, *tinta_b_p;/* buffers previos de tinta      */
    float *presion, *divergencia;/* auxiliares del paso de proyeccion         */
    int    celdas_total;         /* (N+2)*(N+2)                               */
} CamposFluido;

/* Resolucion de la malla; global para que la macro IX() sea legible. */
static int malla_n = MALLA_DEFAULT;

/* ===========================================================================
 * SECCION 1: utilidades generales
 * ======================================================================== */

/* Devuelve un flotante pseudoaleatorio uniforme en [minimo, maximo]. */
static float aleatorio_rango(float minimo, float maximo)
{
    return minimo + ((float)rand() / (float)RAND_MAX) * (maximo - minimo);
}

/* Limita un valor flotante al intervalo [minimo, maximo]. */
static float acotar(float valor, float minimo, float maximo)
{
    if (valor < minimo) return minimo;
    if (valor > maximo) return maximo;
    return valor;
}

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

/* ===========================================================================
 * SECCION 2: lectura y validacion de argumentos (programacion defensiva)
 * ======================================================================== */

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
static int procesar_argumentos(int argc, char *argv[], Configuracion *config)
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

    for (indice = 1; indice < argc; indice++) {
        const char *opcion = argv[indice];

        if (strcmp(opcion, "-h") == 0 || strcmp(opcion, "--help") == 0) {
            imprimir_ayuda(argv[0]);
            return -1;
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

/* ===========================================================================
 * SECCION 3: manejo de memoria de los campos
 * ======================================================================== */

/*
 * Reserva e inicializa en cero todos los campos de la simulacion.
 * Retorna 1 si toda la memoria se obtuvo correctamente, 0 si alguna reserva
 * fallo (en cuyo caso libera lo ya reservado antes de retornar).
 */
static int reservar_campos(CamposFluido *campos, int resolucion)
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
static void liberar_campos(CamposFluido *campos)
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
static void limpiar_campos(CamposFluido *campos)
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

/* ===========================================================================
 * SECCION 4: nucleo numerico (Stable Fluids)
 * ======================================================================== */

/*
 * Aplica las condiciones de frontera sobre el anillo de celdas fantasma.
 *   tipo_borde = BND_ESCALAR -> el borde copia al vecino interior (Neumann)
 *   tipo_borde = BND_VEL_X   -> se invierte la componente x en muros verticales
 *   tipo_borde = BND_VEL_Y   -> se invierte la componente y en muros horizontales
 * Invertir el signo equivale a un muro solido: el fluido rebota en lugar de
 * atravesar la pared.
 */
static void aplicar_frontera(int tipo_borde, float *campo)
{
    int i;

    for (i = 1; i <= malla_n; i++) {
        campo[IX(0, i)]           = (tipo_borde == BND_VEL_X)
                                    ? -campo[IX(1, i)]       : campo[IX(1, i)];
        campo[IX(malla_n + 1, i)] = (tipo_borde == BND_VEL_X)
                                    ? -campo[IX(malla_n, i)] : campo[IX(malla_n, i)];
        campo[IX(i, 0)]           = (tipo_borde == BND_VEL_Y)
                                    ? -campo[IX(i, 1)]       : campo[IX(i, 1)];
        campo[IX(i, malla_n + 1)] = (tipo_borde == BND_VEL_Y)
                                    ? -campo[IX(i, malla_n)] : campo[IX(i, malla_n)];
    }

    /* Las cuatro esquinas se promedian a partir de sus dos vecinos */
    campo[IX(0, 0)] = 0.5f * (campo[IX(1, 0)] + campo[IX(0, 1)]);
    campo[IX(0, malla_n + 1)] =
        0.5f * (campo[IX(1, malla_n + 1)] + campo[IX(0, malla_n)]);
    campo[IX(malla_n + 1, 0)] =
        0.5f * (campo[IX(malla_n, 0)] + campo[IX(malla_n + 1, 1)]);
    campo[IX(malla_n + 1, malla_n + 1)] =
        0.5f * (campo[IX(malla_n, malla_n + 1)] + campo[IX(malla_n + 1, malla_n)]);
}

/*
 * Suma al campo destino las contribuciones del campo fuente escaladas por dt.
 * Es el operador de fuerzas externas / inyeccion de tinta.
 */
static void agregar_fuente(float *destino, const float *fuente, float dt,
                           int celdas_total)
{
    int i;
    for (i = 0; i < celdas_total; i++) {
        destino[i] += dt * fuente[i];
    }
}

/*
 * Resuelve el sistema lineal disperso   x = (x0 + a * suma_vecinos) / c
 * mediante relajacion de Gauss-Seidel.
 *
 * NOTA IMPORTANTE PARA LA VERSION PARALELA:
 * Gauss-Seidel usa los valores ya actualizados de la misma iteracion, lo que
 * crea una dependencia de datos entre celdas vecinas y hace que este ciclo NO
 * sea directamente paralelizable. En la version con OpenMP se sustituye por
 * Jacobi (que lee de un buffer separado y por tanto no tiene dependencias) o
 * por un recorrido red-black. Aqui se conserva Gauss-Seidel por ser el metodo
 * secuencial de referencia y el que converge mas rapido por iteracion.
 */
static void resolver_lineal(int tipo_borde, float *campo, const float *campo_previo,
                            float a, float c)
{
    int iteracion, i, j;
    const float inverso_c = 1.0f / c;

    for (iteracion = 0; iteracion < ITER_GAUSS_SEIDEL; iteracion++) {
        for (j = 1; j <= malla_n; j++) {
            for (i = 1; i <= malla_n; i++) {
                campo[IX(i, j)] = (campo_previo[IX(i, j)] +
                                   a * (campo[IX(i - 1, j)] + campo[IX(i + 1, j)] +
                                        campo[IX(i, j - 1)] + campo[IX(i, j + 1)]))
                                  * inverso_c;
            }
        }
        aplicar_frontera(tipo_borde, campo);
    }
}

/*
 * Difusion implicita: resuelve  x - a*lap(x) = x0,  con a = dt*coef*N*N.
 * El esquema implicito es incondicionalmente estable (esa es la clave del
 * metodo de Stam frente a una difusion explicita).
 */
static void difundir(int tipo_borde, float *campo, const float *campo_previo,
                     float coeficiente, float dt)
{
    const float a = dt * coeficiente * (float)malla_n * (float)malla_n;
    resolver_lineal(tipo_borde, campo, campo_previo, a, 1.0f + 4.0f * a);
}

/*
 * Adveccion semi-Lagrangiana: para cada celda se retrocede en el tiempo
 * siguiendo el campo de velocidad y se interpola bilinealmente el valor del
 * campo en la posicion de origen. Este esquema tambien es incondicionalmente
 * estable, sin importar la magnitud de la velocidad.
 */
static void advectar(int tipo_borde, float *campo, const float *campo_previo,
                     const float *vel_x, const float *vel_y, float dt)
{
    int   i, j, i0, i1, j0, j1;
    float origen_x, origen_y;
    float peso_i1, peso_i0, peso_j1, peso_j0;
    const float dt_malla = dt * (float)malla_n;
    const float limite   = (float)malla_n + 0.5f;

    for (j = 1; j <= malla_n; j++) {
        for (i = 1; i <= malla_n; i++) {
            /* Trazado hacia atras de la particula que llega a (i,j) */
            origen_x = (float)i - dt_malla * vel_x[IX(i, j)];
            origen_y = (float)j - dt_malla * vel_y[IX(i, j)];

            /* Se acota el origen al dominio para no leer fuera del arreglo */
            origen_x = acotar(origen_x, 0.5f, limite);
            origen_y = acotar(origen_y, 0.5f, limite);

            i0 = (int)origen_x;  i1 = i0 + 1;
            j0 = (int)origen_y;  j1 = j0 + 1;

            /* Pesos de la interpolacion bilineal */
            peso_i1 = origen_x - (float)i0;  peso_i0 = 1.0f - peso_i1;
            peso_j1 = origen_y - (float)j0;  peso_j0 = 1.0f - peso_j1;

            campo[IX(i, j)] =
                peso_i0 * (peso_j0 * campo_previo[IX(i0, j0)] +
                           peso_j1 * campo_previo[IX(i0, j1)]) +
                peso_i1 * (peso_j0 * campo_previo[IX(i1, j0)] +
                           peso_j1 * campo_previo[IX(i1, j1)]);
        }
    }

    aplicar_frontera(tipo_borde, campo);
}

/*
 * Proyeccion de Hodge: descompone el campo de velocidad en una parte sin
 * divergencia mas el gradiente de un campo de presion, y se queda solo con la
 * primera. Esto impone la condicion de incompresibilidad div(u) = 0 y es lo
 * que produce los remolinos caracteristicos del fluido.
 */
static void proyectar(float *vel_x, float *vel_y, float *presion,
                      float *divergencia)
{
    int i, j;
    const float h = 1.0f / (float)malla_n;

    /* 1. Calcular la divergencia del campo de velocidad */
    for (j = 1; j <= malla_n; j++) {
        for (i = 1; i <= malla_n; i++) {
            divergencia[IX(i, j)] = -0.5f * h *
                (vel_x[IX(i + 1, j)] - vel_x[IX(i - 1, j)] +
                 vel_y[IX(i, j + 1)] - vel_y[IX(i, j - 1)]);
            presion[IX(i, j)] = 0.0f;
        }
    }
    aplicar_frontera(BND_ESCALAR, divergencia);
    aplicar_frontera(BND_ESCALAR, presion);

    /* 2. Resolver la ecuacion de Poisson  lap(p) = div(u) */
    resolver_lineal(BND_ESCALAR, presion, divergencia, 1.0f, 4.0f);

    /* 3. Restar el gradiente de presion a la velocidad */
    for (j = 1; j <= malla_n; j++) {
        for (i = 1; i <= malla_n; i++) {
            vel_x[IX(i, j)] -= 0.5f *
                (presion[IX(i + 1, j)] - presion[IX(i - 1, j)]) / h;
            vel_y[IX(i, j)] -= 0.5f *
                (presion[IX(i, j + 1)] - presion[IX(i, j - 1)]) / h;
        }
    }
    aplicar_frontera(BND_VEL_X, vel_x);
    aplicar_frontera(BND_VEL_Y, vel_y);
}

/* Intercambia dos punteros a campo (evita copiar arreglos completos). */
static void intercambiar(float **campo_a, float **campo_b)
{
    float *temporal = *campo_a;
    *campo_a = *campo_b;
    *campo_b = temporal;
}

/*
 * Paso completo de la densidad de tinta de un canal:
 * fuente -> difusion -> adveccion.
 */
static void paso_tinta(float **tinta, float **tinta_previa,
                       const float *vel_x, const float *vel_y,
                       float difusion, float dt, int celdas_total)
{
    agregar_fuente(*tinta, *tinta_previa, dt, celdas_total);
    intercambiar(tinta_previa, tinta);
    difundir(BND_ESCALAR, *tinta, *tinta_previa, difusion, dt);
    intercambiar(tinta_previa, tinta);
    advectar(BND_ESCALAR, *tinta, *tinta_previa, vel_x, vel_y, dt);
}

/*
 * Paso completo del campo de velocidad:
 * fuerzas -> difusion viscosa -> proyeccion -> autoadveccion -> proyeccion.
 * Se proyecta dos veces porque tanto la difusion como la adveccion vuelven a
 * introducir divergencia en el campo.
 */
static void paso_velocidad(CamposFluido *campos, float viscosidad, float dt)
{
    agregar_fuente(campos->vel_x, campos->vel_x_p, dt, campos->celdas_total);
    agregar_fuente(campos->vel_y, campos->vel_y_p, dt, campos->celdas_total);

    intercambiar(&campos->vel_x_p, &campos->vel_x);
    difundir(BND_VEL_X, campos->vel_x, campos->vel_x_p, viscosidad, dt);
    intercambiar(&campos->vel_y_p, &campos->vel_y);
    difundir(BND_VEL_Y, campos->vel_y, campos->vel_y_p, viscosidad, dt);

    proyectar(campos->vel_x, campos->vel_y, campos->presion, campos->divergencia);

    intercambiar(&campos->vel_x_p, &campos->vel_x);
    intercambiar(&campos->vel_y_p, &campos->vel_y);
    advectar(BND_VEL_X, campos->vel_x, campos->vel_x_p,
             campos->vel_x_p, campos->vel_y_p, dt);
    advectar(BND_VEL_Y, campos->vel_y, campos->vel_y_p,
             campos->vel_x_p, campos->vel_y_p, dt);

    proyectar(campos->vel_x, campos->vel_y, campos->presion, campos->divergencia);
}

/* ===========================================================================
 * SECCION 5: fuentes de tinta
 * ======================================================================== */

/*
 * Inicializa las fuentes con posicion, color, fase y fuerza pseudoaleatorias.
 * Las posiciones se mantienen alejadas del borde para que el chorro se
 * desarrolle antes de chocar con los muros.
 */
static void inicializar_fuentes(FuenteTinta *fuentes, int cantidad, int resolucion)
{
    int f;
    const int margen = (resolucion / 8 > 2) ? resolucion / 8 : 2;

    for (f = 0; f < cantidad; f++) {
        fuentes[f].celda_x = margen +
            (int)aleatorio_rango(0.0f, (float)(resolucion - 2 * margen));
        fuentes[f].celda_y = margen +
            (int)aleatorio_rango(0.0f, (float)(resolucion - 2 * margen));

        /* Color saturado: un canal dominante y los otros dos parciales */
        fuentes[f].color_r = aleatorio_rango(0.15f, 1.0f);
        fuentes[f].color_g = aleatorio_rango(0.15f, 1.0f);
        fuentes[f].color_b = aleatorio_rango(0.15f, 1.0f);

        fuentes[f].fase        = aleatorio_rango(0.0f, 2.0f * PI);
        fuentes[f].vel_angular = aleatorio_rango(-0.045f, 0.045f);
        fuentes[f].fuerza      = aleatorio_rango(25.0f, 60.0f);
        fuentes[f].caudal      = aleatorio_rango(60.0f, 110.0f);
    }
}

/*
 * Deposita en los buffers "previos" (que actuan como termino fuente) la tinta
 * y la cantidad de movimiento de cada fuente para el frame actual.
 *
 * ELEMENTO TRIGONOMETRICO: la direccion del chorro de cada fuente rota en el
 * tiempo segun (cos(fase), sin(fase)), lo que genera vortices en espiral.
 * La tinta se deposita en un vecindario de 3x3 celdas para suavizar la
 * inyeccion y evitar valores puntuales muy bruscos.
 */
static void inyectar_fuentes(FuenteTinta *fuentes, int cantidad,
                             CamposFluido *campos)
{
    int f, dx, dy, celda_i, celda_j;
    float direccion_x, direccion_y;

    /* Los buffers fuente se limpian cada frame */
    const size_t bytes = (size_t)campos->celdas_total * sizeof(float);
    memset(campos->vel_x_p,   0, bytes);
    memset(campos->vel_y_p,   0, bytes);
    memset(campos->tinta_r_p, 0, bytes);
    memset(campos->tinta_g_p, 0, bytes);
    memset(campos->tinta_b_p, 0, bytes);

    for (f = 0; f < cantidad; f++) {
        /* Avance de la fase angular del chorro */
        fuentes[f].fase += fuentes[f].vel_angular;
        if (fuentes[f].fase > 2.0f * PI) {
            fuentes[f].fase -= 2.0f * PI;
        }

        direccion_x = cosf(fuentes[f].fase) * fuentes[f].fuerza;
        direccion_y = sinf(fuentes[f].fase) * fuentes[f].fuerza;

        for (dy = -1; dy <= 1; dy++) {
            for (dx = -1; dx <= 1; dx++) {
                celda_i = fuentes[f].celda_x + dx;
                celda_j = fuentes[f].celda_y + dy;

                /* Solo se escribe dentro de las celdas interiores */
                if (celda_i < 1 || celda_i > malla_n ||
                    celda_j < 1 || celda_j > malla_n) {
                    continue;
                }

                campos->vel_x_p[IX(celda_i, celda_j)] += direccion_x;
                campos->vel_y_p[IX(celda_i, celda_j)] += direccion_y;

                campos->tinta_r_p[IX(celda_i, celda_j)] +=
                    fuentes[f].caudal * fuentes[f].color_r;
                campos->tinta_g_p[IX(celda_i, celda_j)] +=
                    fuentes[f].caudal * fuentes[f].color_g;
                campos->tinta_b_p[IX(celda_i, celda_j)] +=
                    fuentes[f].caudal * fuentes[f].color_b;
            }
        }
    }
}

/*
 * Multiplica la tinta por un factor menor que 1 para que se desvanezca poco a
 * poco; sin esto la pantalla terminaria saturada de blanco.
 */
static void disipar_tinta(CamposFluido *campos)
{
    int i;
    for (i = 0; i < campos->celdas_total; i++) {
        campos->tinta_r[i] *= DISIPACION;
        campos->tinta_g[i] *= DISIPACION;
        campos->tinta_b[i] *= DISIPACION;
    }
}

/* ===========================================================================
 * SECCION 6: renderizado
 * ======================================================================== */

/*
 * Vuelca la densidad de tinta a la textura de SDL. Cada celda interior de la
 * malla se convierte en un pixel de una textura de N x N que luego SDL escala
 * al tamano de la ventana con filtrado lineal.
 */
static void renderizar_tinta(SDL_Texture *textura, const CamposFluido *campos)
{
    void  *pixeles_crudos = NULL;
    int    pitch_bytes    = 0;
    Uint32 *pixeles;
    int    i, j;
    float  valor_r, valor_g, valor_b, promedio;

    if (SDL_LockTexture(textura, NULL, &pixeles_crudos, &pitch_bytes) != 0) {
        fprintf(stderr, "Advertencia: no se pudo bloquear la textura: %s\n",
                SDL_GetError());
        return;
    }

    pixeles = (Uint32 *)pixeles_crudos;

    for (j = 0; j < malla_n; j++) {
        Uint32 *fila = pixeles + (size_t)j * ((size_t)pitch_bytes / sizeof(Uint32));
        for (i = 0; i < malla_n; i++) {
            /* La celda interior (i+1, j+1) corresponde al pixel (i, j).
             * Se aplica un mapeo tipo Reinhard (x / (x + k)) en lugar de un
             * recorte duro, para que la acumulacion de muchas fuentes se
             * comprima suavemente hacia blanco en vez de saturar de golpe. */
            valor_r = campos->tinta_r[IX(i + 1, j + 1)] * BRILLO_FACTOR;
            valor_g = campos->tinta_g[IX(i + 1, j + 1)] * BRILLO_FACTOR;
            valor_b = campos->tinta_b[IX(i + 1, j + 1)] * BRILLO_FACTOR;

            valor_r = valor_r / (valor_r + CONTRASTE_FACTOR);
            valor_g = valor_g / (valor_g + CONTRASTE_FACTOR);
            valor_b = valor_b / (valor_b + CONTRASTE_FACTOR);

            /* Empuja cada canal lejos del gris promedio del pixel para
             * recuperar saturacion sin alterar el brillo medio, evitando
             * que el resultado se vea lavado/blanco. */
            promedio = (valor_r + valor_g + valor_b) / 3.0f;
            valor_r = acotar(promedio + (valor_r - promedio) * SATURACION_FACTOR, 0.0f, 1.0f);
            valor_g = acotar(promedio + (valor_g - promedio) * SATURACION_FACTOR, 0.0f, 1.0f);
            valor_b = acotar(promedio + (valor_b - promedio) * SATURACION_FACTOR, 0.0f, 1.0f);

            fila[i] = ((Uint32)0xFF << 24) |
                      ((Uint32)(valor_r * 255.0f) << 16) |
                      ((Uint32)(valor_g * 255.0f) <<  8) |
                      ((Uint32)(valor_b * 255.0f));
        }
    }

    SDL_UnlockTexture(textura);
}

/* ===========================================================================
 * SECCION 7: programa principal
 * ======================================================================== */

int main(int argc, char *argv[])
{
    Configuracion  config;
    CamposFluido   campos;
    FuenteTinta   *fuentes = NULL;
    SDL_Window    *ventana = NULL;
    SDL_Renderer  *renderizador = NULL;
    SDL_Texture   *textura = NULL;
    SDL_Event      evento;

    int    ejecutando = 1;
    int    codigo_salida = EXIT_SUCCESS;
    int    frames_acumulados = 0;
    Uint64 ticks_previos, ticks_actuales, frecuencia_reloj;
    Uint64 ticks_ultimo_reporte;
    double segundos_transcurridos, fps_actual = 0.0;
    char   titulo_ventana[128];

    /* --- 1. Argumentos --------------------------------------------------- */
    int resultado_args = procesar_argumentos(argc, argv, &config);
    if (resultado_args == -1) return EXIT_SUCCESS;   /* se pidio la ayuda    */
    if (resultado_args ==  0) return EXIT_FAILURE;   /* argumentos invalidos */

    malla_n = config.malla_n;   /* fija la resolucion usada por la macro IX() */
    srand(config.semilla);

    printf("Simulacion de fluidos (Navier-Stokes) - version SECUENCIAL\n");
    printf("  Malla        : %d x %d celdas (%d celdas interiores)\n",
           config.malla_n, config.malla_n, config.malla_n * config.malla_n);
    printf("  Fuentes      : %d\n", config.num_fuentes);
    printf("  Ventana      : %d x %d px\n",
           config.ventana_ancho, config.ventana_alto);
    printf("  Semilla      : %u\n", config.semilla);
    printf("  dt / visc / diff : %.3f / %.5f / %.5f\n\n",
           (double)config.dt, (double)config.viscosidad, (double)config.difusion);

    /* --- 2. Memoria de los campos ---------------------------------------- */
    if (!reservar_campos(&campos, config.malla_n)) {
        return EXIT_FAILURE;
    }

    fuentes = (FuenteTinta *)malloc((size_t)config.num_fuentes *
                                    sizeof(FuenteTinta));
    if (fuentes == NULL) {
        fprintf(stderr, "Error: no se pudo reservar memoria para las fuentes.\n");
        liberar_campos(&campos);
        return EXIT_FAILURE;
    }
    inicializar_fuentes(fuentes, config.num_fuentes, config.malla_n);

    /* --- 3. Inicializacion de SDL ---------------------------------------- */
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "Error: no se pudo inicializar SDL: %s\n", SDL_GetError());
        codigo_salida = EXIT_FAILURE;
        goto limpieza;
    }

    /* Filtrado lineal para que la malla escalada no se vea pixelada */
    SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "linear");

    ventana = SDL_CreateWindow("Screensaver de fluidos - secuencial",
                               SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                               config.ventana_ancho, config.ventana_alto,
                               SDL_WINDOW_SHOWN);
    if (ventana == NULL) {
        fprintf(stderr, "Error: no se pudo crear la ventana: %s\n", SDL_GetError());
        codigo_salida = EXIT_FAILURE;
        goto limpieza;
    }

    renderizador = SDL_CreateRenderer(ventana, -1, SDL_RENDERER_ACCELERATED);
    if (renderizador == NULL) {
        /* Reintento por software si no hay aceleracion disponible */
        fprintf(stderr, "Aviso: sin renderizador acelerado (%s). "
                        "Se intenta por software.\n", SDL_GetError());
        renderizador = SDL_CreateRenderer(ventana, -1, SDL_RENDERER_SOFTWARE);
    }
    if (renderizador == NULL) {
        fprintf(stderr, "Error: no se pudo crear el renderizador: %s\n",
                SDL_GetError());
        codigo_salida = EXIT_FAILURE;
        goto limpieza;
    }

    textura = SDL_CreateTexture(renderizador, SDL_PIXELFORMAT_ARGB8888,
                                SDL_TEXTUREACCESS_STREAMING,
                                config.malla_n, config.malla_n);
    if (textura == NULL) {
        fprintf(stderr, "Error: no se pudo crear la textura: %s\n", SDL_GetError());
        codigo_salida = EXIT_FAILURE;
        goto limpieza;
    }

    /* --- 4. Ciclo principal ---------------------------------------------- */
    frecuencia_reloj     = SDL_GetPerformanceFrequency();
    ticks_previos        = SDL_GetPerformanceCounter();
    ticks_ultimo_reporte = ticks_previos;

    while (ejecutando) {
        /* 4.1 Eventos de ventana y teclado */
        while (SDL_PollEvent(&evento)) {
            if (evento.type == SDL_QUIT) {
                ejecutando = 0;
            } else if (evento.type == SDL_KEYDOWN) {
                SDL_Keycode tecla = evento.key.keysym.sym;
                if (tecla == SDLK_ESCAPE || tecla == SDLK_q) {
                    ejecutando = 0;
                } else if (tecla == SDLK_r) {
                    limpiar_campos(&campos);
                    inicializar_fuentes(fuentes, config.num_fuentes,
                                        config.malla_n);
                }
            }
        }

        /* 4.2 Un paso de simulacion */
        inyectar_fuentes(fuentes, config.num_fuentes, &campos);

        paso_velocidad(&campos, config.viscosidad, config.dt);

        paso_tinta(&campos.tinta_r, &campos.tinta_r_p,
                   campos.vel_x, campos.vel_y,
                   config.difusion, config.dt, campos.celdas_total);
        paso_tinta(&campos.tinta_g, &campos.tinta_g_p,
                   campos.vel_x, campos.vel_y,
                   config.difusion, config.dt, campos.celdas_total);
        paso_tinta(&campos.tinta_b, &campos.tinta_b_p,
                   campos.vel_x, campos.vel_y,
                   config.difusion, config.dt, campos.celdas_total);

        disipar_tinta(&campos);

        /* 4.3 Dibujo */
        renderizar_tinta(textura, &campos);
        SDL_RenderClear(renderizador);
        SDL_RenderCopy(renderizador, textura, NULL, NULL);
        SDL_RenderPresent(renderizador);

        /* 4.4 Medicion de FPS (promedio sobre ventanas de ~0.5 s) */
        frames_acumulados++;
        ticks_actuales = SDL_GetPerformanceCounter();
        segundos_transcurridos = (double)(ticks_actuales - ticks_ultimo_reporte) /
                                 (double)frecuencia_reloj;

        if (segundos_transcurridos >= 0.5) {
            fps_actual = (double)frames_acumulados / segundos_transcurridos;
            frames_acumulados = 0;
            ticks_ultimo_reporte = ticks_actuales;

            snprintf(titulo_ventana, sizeof(titulo_ventana),
                     "Fluidos secuencial | N=%d | fuentes=%d | FPS= %.2f",
                     config.malla_n, config.num_fuentes, fps_actual);
            SDL_SetWindowTitle(ventana, titulo_ventana);

            printf("FPS= %.2f\n", fps_actual);
            fflush(stdout);
        }

        ticks_previos = ticks_actuales;
    }

    (void)ticks_previos;  /* se conserva por claridad del ciclo de tiempo */

    /* --- 5. Liberacion ordenada de recursos ------------------------------ */
limpieza:
    if (textura      != NULL) SDL_DestroyTexture(textura);
    if (renderizador != NULL) SDL_DestroyRenderer(renderizador);
    if (ventana      != NULL) SDL_DestroyWindow(ventana);
    SDL_Quit();

    free(fuentes);
    liberar_campos(&campos);

    return codigo_salida;
}