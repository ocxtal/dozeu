
CC = gcc
OFLAGS = -g
CFLAGS = $(OFLAGS) -std=c99 -Wall -Wno-unused-variable -Wno-unused-function -Wno-constant-conversion -march=native


SRCS = example.ascii.c example.2bit.c example.4bit.c example.protein.c
TGTS = $(SRCS:.c=)

all: $(TGTS)

$(TGTS): $(SRCS) dozeu.h
	$(CC) $(CFLAGS) -o $@ $@.c

test: all
	(for t in $(TGTS); do ./$$t; done && echo "succeeded") || echo "failed"

clean:
	rm -f $(TGTS)
