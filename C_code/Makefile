CC = gcc

ifndef HTSLIB_DIR
  $(error HTSLIB_DIR is undefined, see README.txt for details)
endif

all: calc_genotypes ld_vcf

install: calc_genotypes
	mkdir -p ../bin
	cp calc_genotypes ../bin/
	cp ld_vcf ../bin/

calc_genotypes: calc_genotypes.o
	$(CC) -o calc_genotypes calc_genotypes.c
calc_genotypes.o: calc_genotypes.c
	$(CC) -c calc_genotypes.c
ld_vcf: ld_vcf.c
	$(CC) -Wall -O3 ld_vcf.c -I $(HTSLIB_DIR)/htslib -o ld_vcf -L$(HTSLIB_DIR) -Wl,-rpath,$(HTSLIB_DIR) -lhts

clean:
	\rm -f *.o calc_genotypes ld_vcf
