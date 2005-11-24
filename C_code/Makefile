CC = gcc

all: calc_genotypes

calc_genotypes: calc_genotypes.o
	$(CC) -o $@ $<	
calc_genotypes.o:
	$(CC) -c $<			
clean:
	\rm -f *.o calc_genotypes
