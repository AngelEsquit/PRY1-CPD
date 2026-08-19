# Makefile - Screensaver de fluidos (Navier-Stokes)
# Version secuencial. La carpeta src/paralela queda preparada para la
# siguiente fase de implementacion de la version paralela con OpenMP.

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -Iinclude
SDLFLAGS = $(shell pkg-config --cflags sdl2)
SDLLIBS  = $(shell pkg-config --libs sdl2)
LIBS     = $(SDLLIBS) -lm

SEQ_SRC  = $(filter-out src/secuencial/Screensaver_seq.c,$(wildcard src/secuencial/*.c))
SEQ_BIN  = Screensaver_seq

.PHONY: all clean run

all: $(SEQ_BIN)

$(SEQ_BIN): $(SEQ_SRC)
	$(CC) $(CFLAGS) $(SDLFLAGS) -o $@ $(SEQ_SRC) $(LIBS)

# Ejecucion rapida con parametros de prueba
run: $(SEQ_BIN)
	./$(SEQ_BIN) -n 128 -f 6

clean:
	rm -f $(SEQ_BIN)
