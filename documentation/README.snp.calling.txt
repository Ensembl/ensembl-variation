This document refers to creating a variation database from the following
resources:

1) A set of traces from the Ensembl Trace Archive (query sequences)
2a) An Ensembl core database of the relevant organism OR
2b) A FASTA file of the organism's reference sequence by chromosome (target)

NB most if not all of these commands require large amounts of memory and/or
runtime, and hence should be run on the farm wherever possible.


1.0 Preparing trace files
=========================

Trace files can be obtained from the Ensembl trace archive:

http://trace.ensembl.org/

For example, the files for Taeniopygia guttata (Zebra finch) can be found at:

ftp://ftp.ensembl.org/pub/traces/taeniopygia_guttata/

Two sets of files are required - the trace files themselves (in FASTA) format in
the fasta/ directory, and the corresponding quality files in the qual/
directory. Each trace file has a corresponding quality file; the two files need
to be combined into a fastq file.

This is done using the compFastq_all script in
ensembl-variation/scripts/ssahaSNP/

a) Put all fasta and fastq files in the same directory

b) Unzip both sets

c) Edit compFastq_all such that the directory used by the while() loop in sub
compFastq is the directory containing your files

d) run:

perl compFastq_all



2.0 Preparing the target files
==============================

If available, a set of target FASTA files should be generated from the relevant
Ensembl core database. This is done using the get_target_dna.pl script in
ensembl-variation/scripts/import/

a) Edit the ensembl.registry file such that the entry for the core database for
your species is pointing to the correct database. For example, here is the entry
used for Taeniopygia guttata:

Bio::EnsEMBL::DBSQL::DBAdaptor->new
  ( '-species' => "Taeniopygia_guttata",
    '-group'   => "core",
    '-port'    => 3306,
    '-host'    => 'ens-staging',
    '-user'    => 'ensro',
    '-pass'    => '',
    '-dbname'  => 'taeniopygia_guttata_core_54_1',
  );

b) run e.g.:

perl get_target_dna.pl --target_dir [target_dir]/ --species Taeniopygia_guttata

changing the --target_dir flag to the directory where you want the target
sequences, and the --species flag to the relevant species.

c) Concatenate the files together:

cat [target_dir]/*.fa > [target_dir]/zfinch.fa



3.0 Running pileup
==================

The process of running pileup is multi-staged and is well documented in the
"ssaha_pileup-readme" file in the pileup archive which is currently available at

http://www.sanger.ac.uk/Software/analysis/SSAHA2/

Choose the appropriate section according to the sequencing technology used to
generate the traces and follow the steps detailed therein.

The final two steps (get_seqreads and ssaha_pileup) require very large amounts
of memory (~30gb for Taeniopygia guttata), and hence need to be run on Turing.
Files need to be copied to an area on Turing's own filesystem using scp before
running.



4.0 Create a variation database
===============================

The next few stages of the process involve parsing the output from pileup into a
variation database. The first step is to create a database on e.g.
ens-genomics2:

mysql -uensadmin -pensembl -hens-genomics2 -e"create database
will_zebrafinch_var_54"
mysql -uensadmin -pensembl -hens-genomics2 will_zebrafinch_var_54 < table.sql

table.sql is found in ensembl-variation/sql/

Your ensembl.registry file should then be modified such that the variation
database for your species points to this newly created database:

Bio::EnsEMBL::DBSQL::DBAdaptor->new
  ( '-species' => "Taeniopygia_guttata",
    '-group'   => "variation",
    '-port'    => 3306,
    '-host'    => 'ens-genomics2',
    '-user'    => '******',
    '-pass'    => '******',
    '-dbname'  => 'will_zebrafinch_var_54',
  );

The next few stages involve using the script run_pileup.pl, found in
ensemb-variation/scripts/ssahaSNP/. This script contains several
methods/subroutines; they are run by uncommenting the appropriate line out of
the following:

# read_cigar_file();
# read_match_file();
# make_pileup_reads_file();
# parse_pileup_snp_file();
# merge_pileup_snp();
# create_vdb();
# PAR_regions();
# read_coverage();

e.g. to run the sub-routine create_vdb(), uncomment the line with create_vdb()
on it. In theory the subroutines could be run sequentially by uncommenting
several at once, though in practice it is simpler to run one at a time to aid
debugging of any errors that occur.

The script needs to be run with several parameters:

perl run_pileup.pl --species Taeniopygia_guttata --input_file genome_SNP.out
--tmpdir [tmp_dir] --tmpfile temp --strain_name tg1
--cigar_file genome_cigar.dat

where:

--tmpdir is a temporary directory (preferably on lustre since the generated
files are very large)
--tmpfile is the name for a temporary file used by the database import library
(Import.pm)
--strain_name is the strain name of the sequenced individual
--cigar_file is the processed cigar file produced when running the pileup
pipeline

a) Run using parse_pileup_snp_file() - this parses the output file from pileup
(in the running example this is named genome_SNP.out) into various tables in the
variation database.

b) Run using merge_pileup_snp() - the line

my $variation_name = "TEMP";

in this subroutine should be changed to something appropriate for the species;
for our running example use "temptgu".

c) Run using create_vdb() - some of the code in this sub-routine needs to be
changed before running:

my %rec_strain("tg1", => 1);

should be modified such that each sequenced strain has a name and entry in the
hash.

my $variation_name = "ENSTNISNP";

should be changed to something appropriate for the species.

d) Run using read_coverage() - this makes read coverage calculations based on
the genome_cigar.dat file produced by pileup and writes them to a file in the
directory specified by --tmpdir; these are used later in post processing so
should not be deleted/misplaced.



5.0 Post processing
===================

Various post-processing stages need to be applied to the database; these are run
by various sub-scripts controlled by a central dispatching script,
parallel_post_process.pl, found in ensembl-variation/scripts/import/

The following command should be run:

perl parallel_post_process.pl -species Taeniopygia_guttata -tmpdir
[tmp_dir] -tmpfile zfinch.txt -num_processes 10 -variation_feature
-flanking_sequence -transcript_variation -top_level 1

This will dispatch several processes to the farm, and populate the
variation_feature, flanking_sequence and transcript_variation tables.

A separate script needs to be run for read coverage, parallel_read_coverage.pl:

perl parallel_read_coverage.pl --species Taeniopygia_guttata --tmpdir
[tmp_dir] --readdir [read_dir] --maxlevel 2

where --readdir is the directory containing the previously generated read
coverage files