
CC = gcc
OFLAGS = -g
CFLAGS = $(OFLAGS) -std=c99 -Wall -Wno-unused-variable -Wno-unused-function -Wno-constant-conversion -march=native


SRCS = example.ascii.c example.2bit.c example.protein.c
TGTS = $(SRCS:.c=)

all: $(TGTS)

.c.o:
	$(CC) $(CFLAGS) -o $(<:c=o) $<

test: all
	(for t in $(TGTS); do ./$$t; done && echo "succeeded") || echo "failed"

clean:
	rm -f $(TGTS)

$(TGTS): dozeu.h

