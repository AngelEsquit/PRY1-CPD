# Makefile - Fluid screensaver (Navier-Stokes)
# Sequential version organized into functional modules.

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -Iinclude -fopenmp
SDLFLAGS = $(shell pkg-config --cflags sdl2)
SDLLIBS  = $(shell pkg-config --libs sdl2)
LIBS     = $(SDLLIBS) -lm -fopenmp

SRC      = $(shell find src -type f -name '*.c' | sort)
BIN      = ss

.PHONY: all clean run

all: $(BIN)

$(BIN): $(SRC)
	$(CC) $(CFLAGS) $(SDLFLAGS) -o $@ $(SRC) $(LIBS)

# Quick run with test parameters
run: $(BIN)
	./$(BIN) -n 128 -f 6

clean:
	rm -f $(BIN)
