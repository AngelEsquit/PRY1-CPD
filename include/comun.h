#ifndef COMUN_H
#define COMUN_H

/* ===========================================================================
 * comun.h
 * ---------------------------------------------------------------------------
 * Constantes, macros y valores compartidos por todos los modulos del
 * screensaver de fluidos (Navier-Stokes / Stable Fluids).
 * ======================================================================== */

/* ---------------------------------------------------------------------------
 * Constantes de configuracion y limites de validacion (programacion defensiva)
 * ------------------------------------------------------------------------- */
#define MALLA_MIN            16    /* resolucion minima de la malla            */
#define MALLA_MAX          1024    /* resolucion maxima de la malla            */
#define MALLA_DEFAULT       256
#define FUENTES_MIN           1    /* al menos una fuente de tinta             */
#define FUENTES_MAX         256
#define FUENTES_DEFAULT       6
#define VENTANA_ANCHO_MIN   640    /* exigido por el enunciado del proyecto    */
#define VENTANA_ALTO_MIN    480
#define VENTANA_ANCHO_DEF  1920 
#define VENTANA_ALTO_DEF   1080 
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
#define SATURACION_FACTOR    2.4f   /* aleja cada canal del gris promedio del
                                        pixel (>1 = mas vivido); al ser
                                        simetrico respecto al promedio no
                                        empuja el color hacia blanco        */
#define VISC_DEFAULT       0.0000f /* viscosidad cinematica del fluido         */
#define DIFF_DEFAULT       0.0001f /* difusion de la tinta: sin esto (0.0) la
                                       tinta solo se suaviza por la forma fija
                                       del nucleo de inyeccion, y se ve con un
                                       borde duro/brusco justo al nacer; con
                                       difusion activa cada fuente nace ya
                                       difuminandose, como tinta real en agua */
#define DISIPACION         0.995f /* retencion de tinta por frame (factor
                                      multiplicativo); con 0.3 la tinta perdia
                                      el 70% de su valor en cada frame y se
                                      apagaba antes de poder advectarse fuera
                                      de la celda de origen                   */

/* ---------------------------------------------------------------------------
 * Indexacion de la malla.
 * La malla tiene (N+2)x(N+2) celdas: N x N celdas interiores (indices 1..N)
 * mas un anillo de celdas fantasma en el borde para las condiciones de
 * frontera. Se almacena en un arreglo 1D en orden row-major.
 *
 * "malla_n" se define una sola vez en main.c y se declara aqui como extern
 * para que la macro IX() este disponible (y sea legible) en todos los
 * modulos que trabajan sobre la malla.
 * ------------------------------------------------------------------------- */
extern int malla_n;

#define IX(i, j) ((i) + (malla_n + 2) * (j))

/* Codigos de frontera usados por aplicar_frontera() (ver solver.c) */
#define BND_ESCALAR   0  /* campo escalar: se copia el vecino interior        */
#define BND_VEL_X     1  /* componente x de velocidad: se refleja en muros    */
#define BND_VEL_Y     2  /* componente y de velocidad: se refleja en muros    */

#endif /* COMUN_H */
