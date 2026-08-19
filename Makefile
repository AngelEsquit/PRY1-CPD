# Makefile - Screensaver de fluidos (Navier-Stokes)
# Version secuencial. La regla paralela queda preparada para la siguiente fase.

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2
SDLFLAGS = $(shell pkg-config --cflags sdl2)
SDLLIBS  = $(shell pkg-config --libs sdl2)
LIBS     = $(SDLLIBS) -lm

SEQ_SRC  = Screensaver_seq.c
SEQ_BIN  = Screensaver_seq

.PHONY: all clean run

all: $(SEQ_BIN)

$(SEQ_BIN): $(SEQ_SRC)
	$(CC) $(CFLAGS) $(SDLFLAGS) -o $@ $< $(LIBS)

# Ejecucion rapida con parametros de prueba
run: $(SEQ_BIN)
	./$(SEQ_BIN) -n 128 -f 6

clean:
	rm -f $(SEQ_BIN)
