-----------------------------------
INSTRUCTIONS TO RUN LD CALCULATIONS
-----------------------------------

In order to run LD calculations, either using the API or the webcode, you will
need:

1 - Install the IPC::Run module. You can get it from CPAN at: 

http://search.cpan.org/~rsod/IPC-Run-0.79/lib/IPC/Run.pm

2 - Compile the C code present in this directory, called calc_genotypes.c. You
need to copy the binary file to one of the directories in your path in order for
the API code to find it or add this directory (ensembl-variation/C_code in your
path doing something like set path = ($path (your_path_to_ensembl_code)/ensembl-variation/C_code). You can compile the code doing:

   make calc_genotypes
   
if the gcc compiler is installed in your system. Otherwise, you will need to
find out which C compiler you have installed and modify the CC variable in the
Makefile pointing to your compiler (icc,...)

