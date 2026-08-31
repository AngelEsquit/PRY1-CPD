#include "solver.h"
#include "comun.h"
#include "utilidades.h"

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
void paso_tinta(float **tinta, float **tinta_previa,
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
void paso_velocidad(CamposFluido *campos, float viscosidad, float dt)
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
