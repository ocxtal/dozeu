
CC=gcc
CFLAGS=-O3 -march=native

all: example

example:
	$(CC) $(CFLAGS) -o example example.c

