# Desarrollo del proyecto: estructura y flujo de trabajo

## 1. Objetivo del proyecto

Este repositorio implementa un screensaver de fluidos en 2D basado en las ecuaciones de Navier-Stokes para fluidos incompresibles, utilizando la aproximacion de Jos Stam conocida como Stable Fluids.

La version actual es la secuencial, y la estructura del proyecto esta preparada para que la version paralela con OpenMP pueda desarrollarse sin romper la organizacion ni la claridad del codigo.

---

## 2. Estado actual del repositorio

La base funcional del proyecto ya esta implementada y compilando correctamente en una version secuencial. En esta fase no se ha incorporado paralelismo; solo se deja la estructura lista para que la siguiente etapa agregue OpenMP o cualquier estrategia paralela sin mezclar logica de diferentes versiones.

La regla general es:

- src/secuencial/: codigo funcional de la version secuencial
- src/paralela/: espacio reservado para la version paralela
- include/: cabeceras compartidas
- docs/: documentacion de desarrollo y decisiones de arquitectura
- raiz del proyecto: solo carpetas, Makefile, README y binario generado

---

## 3. Estructura del proyecto

```text
PRY1-CPD/
├── Makefile
├── README.md
├── Screensaver_seq
├── docs/
│   └── Dev.md
├── include/
│   ├── campos.h
│   ├── comun.h
│   ├── configuracion.h
│   ├── fuentes.h
│   ├── render.h
│   ├── solver.h
│   └── utilidades.h
└── src/
    ├── secuencial/
    │   ├── campos.c
    │   ├── configuracion.c
    │   ├── fuentes.c
    │   ├── main.c
    │   ├── render.c
    │   ├── Screensaver_seq.c
    │   ├── solver.c
    │   └── utilidades.c
    └── paralela/
        └── (reservado para la siguiente etapa)
```

### Convencion de organizacion

- No se permiten archivos .c en la raiz.
- La logica principal de cada version vive bajo src/secuencial/ o src/paralela/.
- Las cabeceras comunes y reutilizables viven en include/.
- La documentacion viva en docs/.
- El binario generado queda en la raiz, como se indica en el Makefile.

---

## 4. Descripcion funcional del proyecto

El programa simula un flujo de fluidos incompresible en 2D sobre una malla regular. El sistema se resuelve con un enfoque de pasos temporales:

1. Inyeccion de tinta y fuerza desde fuentes aleatorias.
2. Difusion de la tinta y la velocidad.
3. Adveccion del campo usando transporte semi-Lagrangiano.
4. Proyeccion de Hodge para mantener la incompresibilidad del fluido.

El render visualiza el campo de tinta como una textura que se muestra en una ventana SDL2.

Las variables clave son:

- malla_n: resolucion de la malla
- fuentes: numero de puntos de inyeccion de tinta
- dt: paso de tiempo
- viscosidad: control del amortiguamiento del fluido
- difusion: control de la dispersion del color

---

## 5. Responsabilidad de cada modulo

### include/comun.h
Contiene definiciones globales del proyecto:

- constantes fisicas y de validacion
- macro IX() para indexar la malla
- codigos de frontera
- valores por defecto del sistema

### include/configuracion.h
Define la estructura de configuracion del programa, asi como la interfaz para parsear argumentos desde la linea de comandos:

- -n
- -f
- -W
- -H
- -s
- -t
- -v
- -d
- -h

### include/campos.h
Define la estructura principal del estado del fluido:

- velocidades x/y
- campos de tinta RGB
- buffers temporales para los pasos numericos
- metadatos de la malla

### include/fuentes.h
Define el modelo de las fuentes de tinta y la forma en que se inyectan:

- posicion
- fuerza
- direccion angular
- caudal de tinta

### include/solver.h
Expone el nucleo numerico del algoritmo:

- difusion
- adveccion
- proyeccion
- intercambio de buffers
- solucion lineal de Gauss-Seidel

### include/render.h
Define la API del render para asignar la textura SDL y convertir la informacion del campo de tinta a imagen visible.

### include/utilidades.h
Contiene utilidades generales como:

- acotado de valores
- entrega de numeros aleatorios
- validacion de argumentos
- helpers de apoyo para la simulacion

---

## 6. Archivos fuente actuales

### src/secuencial/main.c
Es el punto de entrada del programa. Aqui se:

- parsean argumentos
- inicializa la memoria de los campos
- crea la ventana SDL
- ejecuta el ciclo principal
- actualiza la simulacion por frame
- renderiza la textura
- mide FPS y muestra el titulo de la ventana

### src/secuencial/configuracion.c
Implementa la validacion y lectura de las opciones del programa. Mantiene la logica de parseo separada y aislada del resto del flujo principal.

### src/secuencial/campos.c
Gestiona la reserva, liberacion y limpieza de los campos del fluido. Esta es la capa de memoria del sistema.

### src/secuencial/fuentes.c
Genera las fuentes de tinta y las inyecta durante la simulacion. Aqui se controla la variacion temporal y el giro de la direccion del chorro.

### src/secuencial/solver.c
Implementa la numerica principal: la difusion, adveccion, frontera, proyeccion y relajacion lineal. Es el bloque mas importante para la futura version paralela.

### src/secuencial/render.c
Conecta el estado numerico con SDL para dibujar la simulacion en pantalla.

### src/secuencial/utilidades.c
Centraliza funciones complementarias que no pertenecen al nucleo fisico ni a la capa visual.

### src/secuencial/Screensaver_seq.c
Es una version del programa en un solo archivo de referencia para comparacion historica. Actualmente se ha dejado dentro de src/secuencial para evitar la presencia de .c en la raiz y mantener una estructura coherente.

---

## 7. Regla para la version paralela

La carpeta src/paralela/ esta reservada para la siguiente fase. La idea es mantener una separacion clara:

- secuencial: referencia y base funcional
- paralela: implementacion optimizada con OpenMP o con otra tecnica

Se recomienda no mezclar codigo de la version secuencial con la paralela. La mejor estrategia es:

1. mantener la API de modulos sin cambios para la version secuencial
2. implementar una version paralela basada en los mismos archivos conceptuales
3. hacer comparaciones en rendimiento bajo la misma configuracion
4. documentar diferencias de comportamiento y resultados

---

## 8. Compilacion y ejecucion

Desde la raiz del proyecto:

```bash
make
./Screensaver_seq -n 128 -f 6
```

Para limpiar:

```bash
make clean
```

La compilacion usa SDL2 y gcc con soporte C11.

---

## 9. Recomendaciones para continuar con la fase paralela

Cuando se empiece con OpenMP, lo ideal es seguir este orden:

1. Analizar los kernels hot-spot en solver.c.
2. Identificar elementos directamente paralelizables: adveccion, inyeccion, disipacion, bucles de proyecto fuera del solver lineal.
3. Tratar con cuidado la relajacion de Gauss-Seidel, ya que introduce dependencias de datos entre celdas vecinas.
4. Mantener la misma API de entrada/salida para realizar comparacion fair entre secuencial y paralela.
5. Medir FPS y speedup con una bateria de pruebas reproducibles.

---

## 10. Resumen

La estructura actual cumple con los criterios de orden y mantenibilidad esperados para una fase de paralelizacion:

- codigo organizado por capas
- separacion entre version secuencial y futura paralela
- includes centralizadas
- documentacion tecnica en docs/
- raiz limpia y sin archivos .c
- Makefile funcional y reproducible

Esto deja el proyecto listo para evolucionar hacia la version paralela sin perder trazabilidad ni claridad en el desarrollo.
