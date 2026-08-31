#include "utils.h"

#include <stdlib.h>

/* Devuelve un flotante pseudoaleatorio uniforme en [minimo, maximo]. */
float aleatorio_rango(float minimo, float maximo)
{
    return minimo + ((float)rand() / (float)RAND_MAX) * (maximo - minimo);
}

/* Limita un valor flotante al intervalo [minimo, maximo]. */
float acotar(float valor, float minimo, float maximo)
{
    if (valor < minimo) return minimo;
    if (valor > maximo) return maximo;
    return valor;
}
