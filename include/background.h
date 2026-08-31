#ifndef BACKGROUND_H
#define BACKGROUND_H

#include <SDL2/SDL.h>

/* ===========================================================================
 * background.h
 * Screensaver background: a field of fixed, twinkling stars, drawn with
 * SDL's renderer (not the ink texture) right before the ink is drawn on
 * top. The ink texture is now drawn with alpha proportional to brightness
 * (see render.c), so the background shows through wherever there's no ink.
 * ======================================================================== */

#define BACKGROUND_STARS_COUNT 220

/* Base background color: a very dark blue/purple instead of flat black, to
 * give a sense of depth without competing with the ink colors. Applied as
 * the SDL_RenderClear color (see main.c). */
#define BACKGROUND_COLOR_R 6
#define BACKGROUND_COLOR_G 8
#define BACKGROUND_COLOR_B 20

typedef struct {
    float x, y;          /* fixed position in screen pixels                 */
    float base_brightness;/* maximum twinkle brightness, in [0,1]            */
    float phase;          /* current twinkle phase (radians)                 */
    float phase_speed;    /* twinkle angular speed (rad/frame)               */
    int   radius;         /* size in pixels (1 = one dot, 2 = 2x2, etc.)     */
} Star;

/*
 * Places "count" stars at random positions within the window (width x
 * height in pixels), with random base brightness, initial phase and
 * twinkle speed too.
 */
void init_background(Star *stars, int count, int width, int height);

/* Advances every star's twinkle phase by one frame. */
void update_background(Star *stars, int count);

/*
 * Paints the stars directly with SDL's renderer. Must be called AFTER
 * SDL_RenderClear (which already laid down the base color) and BEFORE
 * copying the ink texture (which is now semi-transparent wherever there's
 * no ink, letting whatever is drawn here show through).
 */
void draw_background(SDL_Renderer *renderer, const Star *stars, int count);

#endif /* BACKGROUND_H */
