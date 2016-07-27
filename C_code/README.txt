-----------------------------------
INSTRUCTIONS TO RUN LD CALCULATIONS
-----------------------------------

In order to run LD calculations, either using the API or the webcode, you will
need:

1 - Download and compile htslib. Note if you have already installed the
    ensembl-io module you may already have htslib:
  
  cd [your_local_sw_dir]
  git clone --branch master --depth 1 https://github.com/samtools/htslib.git
  cd htslib
  make
  cd -
  export HTSLIB_DIR=[your_local_sw_dir]/htslib

2 - Compile the C code present in this directory, called calc_genotypes.c and ld_vcf.c.

You can compile the code doing:

  make
   
if the gcc compiler is installed in your system. Otherwise, you will need to
find out which C compiler you have installed and modify the CC variable in the
Makefile pointing to your compiler (icc,...)

You need to copy the created binary files (calc_genotypes and ld_vcf) to one of the
directories in your path in order for the API code to find it or add this directory
(ensembl-variation/C_code) in your path doing something like

  set path = ($path (your_path_to_ensembl_code)/ensembl-variation/C_code)

