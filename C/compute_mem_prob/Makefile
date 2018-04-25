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


$(P): compute_mem_prob.c
	gcc $(CFLAGS) compute_mem_prob.c -lm -o $@

clean:
	rm do
