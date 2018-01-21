SHELL= bash
CC= gcc
CFLAGS= -std=gnu99 -O3

EXEC= sim

SIMULATIONS=                       \
	out.txt

FIGS=


all: $(FIGS) $(SIMULATIONS)

sim: sim.c mt.o
	$(CC) $(CFLAGS) sim.c mt.o -o sim

out.txt: sim
	./sim 18 17 0.01 0.05 0 1000000
	./sim 18 17 0.01 0.05 100 1000000

clean:
	rm -f $(SIMULATIONS) $(EXEC) $(FIGS)
