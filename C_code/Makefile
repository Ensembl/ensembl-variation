CC = gcc

all: calc_genotypes

install: calc_genotypes
	mkdir -p ../bin
	cp calc_genotypes ../bin/

calc_genotypes: calc_genotypes.o
	$(CC) -o calc_genotypes calc_genotypes.c
calc_genotypes.o:
	$(CC) -c calc_genotypes.c
clean:
	\rm -f *.o calc_genotypes
