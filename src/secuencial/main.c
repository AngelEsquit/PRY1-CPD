/* ===========================================================================
 * main.c
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
 * relajacion iterativa de Gauss-Seidel (ver nota en solver.c sobre su impacto
 * en la version paralela).
 *
 * Trigonometria: cada fuente de tinta inyecta velocidad en una direccion que
 * rota en el tiempo, calculada con sin() y cos() sobre una fase propia
 * (ver fuentes.c).
 *
 * Modulos del proyecto (ver include/ y src/):
 *   comun.h              - constantes, macros e indexacion de la malla
 *   utilidades.{h,c}     - funciones de proposito general
 *   configuracion.{h,c}  - argumentos de linea de comandos
 *   campos.{h,c}         - memoria de los campos de la simulacion
 *   fuentes.{h,c}        - fuentes de tinta
 *   solver.{h,c}         - nucleo numerico (Stable Fluids)
 *   render.{h,c}         - volcado a textura SDL
 *   main.c               - programa principal (este archivo)
 *
 * Compilacion:  make            (ver Makefile)
 * Uso:          ./screensaver_seq -n 128 -f 6
 * ===========================================================================
 */

#include <SDL2/SDL.h>
#include <stdio.h>
#include <stdlib.h>

#include "comun.h"
#include "configuracion.h"
#include "campos.h"
#include "fuentes.h"
#include "solver.h"
#include "render.h"

/* Resolucion de la malla; global para que la macro IX() sea legible.
 * Unica definicion del proyecto (declarada extern en comun.h para que el
 * resto de los modulos la puedan usar). */
int malla_n = MALLA_DEFAULT;

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
