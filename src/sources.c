#include "sources.h"
#include "common.h"
#include "utils.h"

#include <string.h>
#include <math.h>

/*
 * Inicializa las fuentes con posicion, velocidad, masa, color, fase y fuerza
 * pseudoaleatorias. Las posiciones se mantienen alejadas del borde para que
 * el chorro se desarrolle antes de chocar con los muros.
 */
void inicializar_fuentes(FuenteTinta *fuentes, int cantidad, int resolucion)
{
    int f;
    const int margen = (resolucion / 8 > 2) ? resolucion / 8 : 2;

    for (f = 0; f < cantidad; f++) {
        fuentes[f].pos_x = (float)margen +
            aleatorio_rango(0.0f, (float)(resolucion - 2 * margen));
        fuentes[f].pos_y = (float)margen +
            aleatorio_rango(0.0f, (float)(resolucion - 2 * margen));

        /* Velocidad inicial pequena: el movimiento principal lo genera la
         * atraccion gravitacional mutua una vez arranca la simulacion. */
        fuentes[f].vel_x = aleatorio_rango(-0.3f, 0.3f);
        fuentes[f].vel_y = aleatorio_rango(-0.3f, 0.3f);
        fuentes[f].masa  = aleatorio_rango(25.0f, 70.0f);

        /* Color saturado: un canal dominante y los otros dos parciales */
        fuentes[f].color_r = aleatorio_rango(0.15f, 1.0f);
        fuentes[f].color_g = aleatorio_rango(0.15f, 1.0f);
        fuentes[f].color_b = aleatorio_rango(0.15f, 1.0f);

        fuentes[f].fase        = aleatorio_rango(0.0f, 2.0f * PI);
        fuentes[f].vel_angular = aleatorio_rango(-0.045f, 0.045f);
        /* Con fuerza en [25,60] el desplazamiento por adveccion junto a la
         * fuente (dt*malla_n*|v|) llega a ~15-19 celdas por frame: varias
         * veces el diametro del vecindario de inyeccion, asi que la tinta
         * "salta" en vez de fluir, y se ve como gotas cayendo en vez de un
         * chorro continuo. En [6,14] el salto queda en ~5-6 celdas, del
         * orden del propio vecindario de inyeccion, y frames consecutivos
         * se solapan en un trazo continuo. */
        fuentes[f].fuerza      = aleatorio_rango(6.0f, 14.0f);
        fuentes[f].caudal      = aleatorio_rango(60.0f, 110.0f);
    }
}

/*
 * Deposita en los buffers "previos" (que actuan como termino fuente) la tinta
 * y la cantidad de movimiento de cada fuente para el frame actual.
 *
 * ELEMENTO TRIGONOMETRICO: la direccion del chorro de cada fuente rota en el
 * tiempo segun (cos(fase), sin(fase)), lo que genera vortices en espiral.
 * La tinta se deposita en un vecindario de celdas con caida gaussiana (en
 * vez de un valor uniforme) para que la mancha nazca ya redondeada: con
 * pocas celdas la forma se percibe como un rombo/pixelado apenas se crea,
 * porque no hay suficientes muestras para leerse como un circulo, y como la
 * difusion de tinta esta en 0 por defecto nada la suaviza despues salvo el
 * movimiento. El peso de cada celda se calcula contra la posicion continua
 * de la fuente (no la celda mas cercana), para que al desplazarse la mancha
 * se deslice de forma continua en vez de saltar de celda en celda.
 */
/* Que tan angosto es el nucleo caliente respecto a la mancha de color
 * (fraccion de sigma); mas chico = punto mas concentrado y puntual. */
#define NUCLEO_SIGMA_FACTOR  0.32f
/* Cuanto brillo blanco se suma en el centro exacto, como fraccion del
 * caudal de la fuente (asi un chorro mas fuerte tiene un nucleo mas
 * brillante tambien, en vez de un blanco fijo que se notaria distinto
 * segun la fuente). */
#define NUCLEO_BRILLO_FACTOR 0.9f

void inyectar_fuentes(FuenteTinta *fuentes, int cantidad, CamposFluido *campos,
                      float aspecto)
{
    int f, dx, dy, celda_i, celda_j, centro_x, centro_y;
    float direccion_x, direccion_y;
    const int   radio  = fuente_radio(malla_n);
    const float sigma  = fuente_sigma(malla_n);
    const float escala = fuente_escala(malla_n);
    /* render.c usa una escala distinta en x y en y para llenar una ventana
     * que no es cuadrada (ver comentario en fuentes.h), asi que un sigma
     * igual en ambos ejes se veria en pantalla como una elipse. Encoger
     * sigma_x por "aspecto" compensa ese estiramiento para que la mancha
     * se vea circular. */
    const float sigma_x = sigma / aspecto;
    /* El area del vecindario (~sigma^2) crece con la escala; dividir entre
     * escala^2 mantiene la tinta total depositada por fuente invariante,
     * sin importar cuantas celdas use el circulo para dibujarse. Encoger
     * sigma_x reduce esa area en un factor "aspecto", asi que se multiplica
     * de vuelta para no perder tinta solo por corregir la forma. */
    const float normalizador = aspecto / (escala * escala);
    /* Nucleo "caliente": un segundo gaussiano, mucho mas angosto y centrado
     * en el mismo punto, que se suma por igual a los tres canales de color
     * (en vez de al color propio de la fuente). Sumar blanco puro empuja el
     * centro exacto de la fuente hacia blanco/brillante, distinguiendolo del
     * color de su estela; al alejarse del centro el nucleo se apaga mucho
     * mas rapido que la mancha de color (por el factor NUCLEO_SIGMA_FACTOR),
     * as que el efecto es un punto brillante que se funde rapidamente en el
     * color propio de la fuente. Usa el mismo "normalizador" que la mancha
     * de color: como ambos sigmas escalan igual con la resolucion de malla,
     * el mismo factor de correccion por resolucion sigue aplicando. */
    const float nucleo_sigma_x  = sigma_x * NUCLEO_SIGMA_FACTOR;
    const float nucleo_sigma    = sigma   * NUCLEO_SIGMA_FACTOR;

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

        /* La posicion es continua (la mueve el sistema de n-cuerpos); solo
         * se usa floor() para ubicar la celda base del vecindario a
         * recorrer, el peso de cada celda usa la posicion continua real. */
        centro_x = (int)floorf(fuentes[f].pos_x);
        centro_y = (int)floorf(fuentes[f].pos_y);

        for (dy = -radio; dy <= radio; dy++) {
            for (dx = -radio; dx <= radio; dx++) {
                celda_i = centro_x + dx;
                celda_j = centro_y + dy;

                /* Solo se escribe dentro de las celdas interiores */
                if (celda_i < 1 || celda_i > malla_n ||
                    celda_j < 1 || celda_j > malla_n) {
                    continue;
                }

                /* Nucleo gaussiano contra la posicion continua real de la
                 * fuente (no contra el centro entero del vecindario), para
                 * que la mancha se deslice suavemente entre celdas. */
                float dif_x = (float)celda_i - fuentes[f].pos_x;
                float dif_y = (float)celda_j - fuentes[f].pos_y;
                float dist2 = (dif_x * dif_x) / (sigma_x * sigma_x) +
                              (dif_y * dif_y) / (sigma  * sigma);
                float peso  = expf(-dist2 / 2.0f) * normalizador;

                /* Mismo calculo que "peso" pero con el sigma angosto del
                 * nucleo caliente, asi que cae a cero mucho mas rapido con
                 * la distancia y solo aporta brillo justo en el centro. */
                float dist2_nucleo = (dif_x * dif_x) / (nucleo_sigma_x * nucleo_sigma_x) +
                                     (dif_y * dif_y) / (nucleo_sigma   * nucleo_sigma);
                float peso_nucleo  = expf(-dist2_nucleo / 2.0f) * normalizador;
                float brillo       = NUCLEO_BRILLO_FACTOR * fuentes[f].caudal *
                                      peso_nucleo;

                campos->vel_x_p[IX(celda_i, celda_j)] += direccion_x * peso;
                campos->vel_y_p[IX(celda_i, celda_j)] += direccion_y * peso;

                campos->tinta_r_p[IX(celda_i, celda_j)] +=
                    fuentes[f].caudal * fuentes[f].color_r * peso + brillo;
                campos->tinta_g_p[IX(celda_i, celda_j)] +=
                    fuentes[f].caudal * fuentes[f].color_g * peso + brillo;
                campos->tinta_b_p[IX(celda_i, celda_j)] +=
                    fuentes[f].caudal * fuentes[f].color_b * peso + brillo;
            }
        }
    }
}

/*
 * Multiplica la tinta por un factor menor que 1 para que se desvanezca poco a
 * poco; sin esto la pantalla terminaria saturada de blanco.
 */
void disipar_tinta(CamposFluido *campos)
{
    int i;
    for (i = 0; i < campos->celdas_total; i++) {
        campos->tinta_r[i] *= DISIPACION;
        campos->tinta_g[i] *= DISIPACION;
        campos->tinta_b[i] *= DISIPACION;
    }
}
