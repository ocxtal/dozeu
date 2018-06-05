
CC=gcc
CFLAGS=-std=c99 -O3 -march=native

all: example

example: example.c dozeu.h
	$(CC) $(CFLAGS) -o example example.c

clean:
	rm example

