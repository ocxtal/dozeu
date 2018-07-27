
CC=gcc
OFLAGS=-O3
CFLAGS=$(OFLAGS) -std=c99 -Wall -Wno-unused-variable -Wno-unused-function -march=native

all: example example.2bit example.protein

example: example.c dozeu.h
	$(CC) $(CFLAGS) -o example example.c
example.2bit: example.2bit.c dozeu.h
	$(CC) $(CFLAGS) -o example.2bit example.2bit.c
example.protein: example.protein.c dozeu.h
	$(CC) $(CFLAGS) -o example.protein example.protein.c

test: all
	(./example && ./example.2bit && ./example.protein && echo "succeeded") || echo "failed"

clean:
	rm -f example example.2bit example.protein

