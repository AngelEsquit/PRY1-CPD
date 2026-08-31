# Makefile - Screensaver de fluidos (Navier-Stokes)
# Version secuencial organizada por modulos funcionales.

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -Iinclude
SDLFLAGS = $(shell pkg-config --cflags sdl2)
SDLLIBS  = $(shell pkg-config --libs sdl2)
LIBS     = $(SDLLIBS) -lm

SEQ_SRC  = $(shell find src/secuencial -path 'src/secuencial/legacy' -prune -o -type f -name '*.c' -print | sort)
SEQ_BIN  = Screensaver

.PHONY: all clean run

all: $(SEQ_BIN)

$(SEQ_BIN): $(SEQ_SRC)
	$(CC) $(CFLAGS) $(SDLFLAGS) -o $@ $(SEQ_SRC) $(LIBS)

# Ejecucion rapida con parametros de prueba
run: $(SEQ_BIN)
	./$(SEQ_BIN) -n 128 -f 6

clean:
	rm -f $(SEQ_BIN)
