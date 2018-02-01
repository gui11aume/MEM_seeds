P= do

CC= gcc
CFLAGS= -std=gnu99 -Wall
LDLIBS= -lm

all: CFLAGS += -O3
all: $(P)

debug: CFLAGS += -DDEBUG -g -O0
debug: $(P)

profile: CFLAGS += -g -pg -O3 -fprofile-arcs -ftest-coverage
profile: $(P)


$(P): memprob.c
	gcc $(CFLAGS) memprob.c -lm -o $@
	gcc $(CFLAGS) kpoly.c -lm -o kpoly

clean:
	rm do
