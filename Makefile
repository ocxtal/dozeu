
CC = gcc
OFLAGS = -O3
WFLAGS = -Wall -Wextra -Wshadow
NWFLAGS = $(shell bash -c "if [[ $(CC) = icc* ]]; then echo '-Wno-unused-function'; else echo '-Wno-unused-function -Wno-unused-label -Wno-constant-conversion -Wno-implicit-fallthrough -Wno-missing-field-initializers'; fi")
CFLAGS = -std=c99 -march=native $(WFLAGS) $(NWFLAGS) -fsanitize=address
GFLAGS = -g -DDEBUG # -fsanitize=address -fsanitize=leak


SRCS = example.ascii.c example.2bit.c example.4bit.c example.protein.c
TGTS = $(SRCS:.c=)

all: $(TGTS)

$(TGTS): $(SRCS) dozeu.h
	$(CC) $(CFLAGS) $(OFLAGS) -o $@ $@.c

test: all
	(for t in $(TGTS); do ./$$t; done && echo "succeeded") || echo "failed"

clean:
	rm -f $(TGTS)
