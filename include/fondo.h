#ifndef FONDO_H
#define FONDO_H

#include <SDL2/SDL.h>

/* ===========================================================================
 * fondo.h
 * Fondo del screensaver: un campo de estrellas fijas con titileo, dibujado
 * con el renderer de SDL (no con la textura de tinta) justo antes de volcar
 * la tinta encima. La textura de tinta ahora se dibuja con alfa proporcional
 * al brillo (ver render.c), asi que donde no hay tinta el fondo se ve.
 * ======================================================================== */

#define FONDO_ESTRELLAS_CANTIDAD 220

/* Color base del fondo: un azul/purpura muy oscuro en vez de negro plano,
 * para dar sensacion de profundidad sin competir con los colores de la
 * tinta. Se aplica como color de SDL_RenderClear (ver main.c). */
#define FONDO_COLOR_R 6
#define FONDO_COLOR_G 8
#define FONDO_COLOR_B 20

typedef struct {
    float x, y;         /* posicion fija en pixeles de pantalla            */
    float brillo_base;  /* brillo maximo del titileo, en [0,1]             */
    float fase;         /* fase actual del titileo (radianes)              */
    float vel_fase;      /* velocidad angular del titileo (rad/frame)       */
    int   radio;         /* tamano en pixeles (1 = un punto, 2 = 2x2, etc.) */
} Estrella;

/*
 * Coloca "cantidad" estrellas en posiciones aleatorias dentro de la ventana
 * (ancho x alto en pixeles), con brillo base, fase inicial y velocidad de
 * titileo tambien aleatorios.
 */
void inicializar_fondo(Estrella *estrellas, int cantidad, int ancho, int alto);

/* Avanza la fase de titileo de cada estrella un frame. */
void actualizar_fondo(Estrella *estrellas, int cantidad);

/*
 * Pinta las estrellas directamente con el renderer de SDL. Debe llamarse
 * DESPUES de SDL_RenderClear (que ya dejo el color base) y ANTES de copiar
 * la textura de tinta (que ahora es semi-transparente donde no hay tinta,
 * dejando ver lo que se dibuje aqui).
 */
void dibujar_fondo(SDL_Renderer *renderizador, const Estrella *estrellas,
                   int cantidad);

#endif /* FONDO_H */
