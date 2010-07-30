########################
#                      #
# SNP Effect Predictor #
#                      #
########################

Copyright (c) 1999-2010 The European Bioinformatics Institute and
Genome Research Limited.  All rights reserved.

This software is distributed under a modified Apache license.
For license details, please see

http://www.ensembl.org/info/about/code_licence.html

Please email comments or questions to the public Ensembl
developers list at <ensembl-dev@ebi.ac.uk>.

Questions may also be sent to the Ensembl help desk at
<helpdesk@ensembl.org>
  


1.0 Requirements
================

This script requires the Ensembl Core and Variation APIs and their relevant
dependencies to be installed. See
http://www.ensembl.org/info/docs/api/index.html for details. No explicit
installation is required.



2.0 Introduction
================

Given a set of variant positions and alleles, it is often useful to know how
these variants may affect transcripts and other genetic factors. Ensembl
provides such annotation for the variants stored in its Ensembl Variation
database; using this script it is possible to do the same for novel variants.



3.0 Running the script
======================

The script is run on the command line as follows:

perl snp_effect_predictor.pl [options]

where [options] represent a set of flags and options to the script. These can be
listed using the flag --help:

perl snp_effect_predictor.pl --help

By default the script connects to the public Ensembl database server at
ensembldb.ensembl.org; other connection options are available.

3.1 Command line options
------------------------

Long form shown in parentheses.

--help : display help message and quit

-s (--species) : species for your data. This can be the latin name
   e.g. "homo_sapiens" or any Ensembl alias e.g. "mouse". Default = "human"
   
-i (--input_file) : input file name. If not specified, the script will attempt
   to read from STDIN.

-f (--format) : input file format. By default, the script auto-detects the input
   file format. Using this option you can force the script to read the input
   file as VCF or pileup format. Not used by default.
   
-o (--output_file) : output file name. Default = "snp_effect_output.txt"

-b (--buffer_size) : sets the internal buffer size, corresponding to the number
   of variations that are used to query the database simultaneously. Set this
   lower to use less memory at the expense of longer run time, and higher to use
   more memory with a faster run time. Default = 500
   
--check_ref : force the script to check the supplied reference allele against
   the sequence stored in the Ensembl Core database. Lines that do not match are
   skipped. Not used by default.
   
--check_exisiting [0|1] : by default the script checks for the existence of
   variants that are co-located with your input. Disabling this by specifying
   --check_existing=0 will bypass this and provide a speed boost. Default : 1
   
--hgnc : adds the HGNC gene identifer (where available) to the output e.g.
    ENSG00000001626;CFTR. Not used by default
   
-d (--db_host) : manually define the database host to connect to.
   Default = "ensembldb.ensembl.org"
   
-u (--user) : manually define the database username used when connecting to the
   database host. Default = "anonymous"
   
-p (--password) : manually define the database password used when connecting to
   the database host. Not used by default
   
-r (--registry_file) : defining a registry file (see 4.3 below) overwrites other
   connection settings and uses those found in the specified registry file to
   connect.
   
3.2 Examples
------------

perl snp_effect_predictor.pl -i snps.txt -o snps_consequences.txt

perl snp_effect_predictor.pl -i snps.txt -o snps_consequences.txt -b 1000 -hgnc

perl snp_effect_predictor.pl -i snps.txt -r ensembl.registry

perl snp_effect_predictor.pl -i snps.txt -d mydbserver -u user -p password



4.0 File formats
================

The SNP Effect Predictor script uses plain text files both as input and output.

4.1 Input file
--------------

The default input file format consists of five columns; these can be comma, tab
(or any whitespace) separated. The columns are:

- Chromosome name : the name of the chromosome to which the variant maps e.g. 1
  for chromosome 1. This can also be the name of a contig or other seq_region
  for species with unassembled genomes.
- Start position
- End position : along with start position, defines the location of the variant.
  These are 1-indexed (i.e. the first base of the chromosome is base 1). For a
  SNP, start should be the same as end. For a multi-basepair substitution, start
  should be the first and end should be the last base affected (e.g. for a 3
  base substitution, end = start + 2). For an insertion relative to the
  reference, start = end + 1, regardless of its length. This means that in this
  instance the start coordinate will be greater than the end coordinate; end
  represents the base immediately 5' of the insertion site, and start the base
  immediately 3'. For a deletion relative to the reference, follow the same
  rules for a multi-basepair substitution.
- Allele string : a "/"-separated string of alleles. The first of these is
  assumed to be the reference allele (although this is not interpreted or
  checked by the script, so if unknown can be set to anything).
- Strand : the strand to which the variant maps. Possible values are 1 and -1;
  + and - can also be used.
- Variant name : optionally, a sixth column can be specified containing an
  identifier for the variant. If not specified, the name is derived from the
  coordinates as described below.

Examples:

- Substitution of C for T at base 100 on forward strand of chromosome 1:
1 100 100 T/C +

- Substitution of CAT for GGC at base 100 on reverse strand of chromosome Y:
Y 100 102 GGC/CAT -

- Deletion of bases 100, 101 and 102 (CTG) on forward strand of chromosome 3:
3 100 102 CTG/- +

- Insertion of ACCG between bases 100 and 101 on reverse strand of chromosome X:
X 101 100 -/ACCG -

Other input formats are also supported; VCF (see
http://1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcfv3.2) and
pileup. These should be auto-detected by the script; if they are not, you can
force the script to use either with the -f|--format command line flags.


4.2 Output file
---------------

The output file contains one line for each predicted consequence; if a variant
falls in more than one transcript it may produce more than one line. The output
file is a tab-delimited file with the following columns:

- Uploaded Variation : a temporary name assigned to the variant, based on the
  position and alleles, or the optional name if specified
- Location : string corresponding to the variant's location, represented as
  Chr:Start-End
- Gene : the Ensembl stable ID of the gene where the variant is located (where
  applicable), with the HGNC name if specified using the --hgnc flag
- Transcript : the Ensembl stabled ID of the transcript where the variant is
  located (where applicable)
- Consequence : the consequence type predicted. See
  http://www.ensembl.org/info/docs/variation/index.html for a description of
  these
- Position in cDNA : the position of the variant in the transcript's cDNA
  sequence (if applicable)
- Position in protein : the position of the variant in the resultant protein
- Amino acid change : a "/"-separated string of possible amino acids generated
  by the variant, with the reference amino acid first
- Corresponding Variation : if an existing variant is found in the Ensembl
  Variation database whose position matches exactly that of the entered
  variation, its identifier appears here

4.3 Registry file
-----------------

It is possible to configure the databases that the script connects to using an
Ensembl Registry file. This should be considered an option for advanced users -
most users should use the default connection options. The file is a perl file
containing an object for each database connection. An example is shown below:

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;

Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  '-species' => "Homo_sapiens",
  '-group'   => "core",
  '-port'    => 5306,
  '-host'    => 'ensembldb.ensembl.org',
  '-user'    => 'anonymous',
  '-pass'    => '',
  '-dbname'  => 'homo_sapiens_core_59_37d'
);


Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(
  '-species' => "Homo_sapiens",
  '-group'   => "variation",
  '-port'    => 5306,
  '-host'    => 'ensembldb.ensembl.org',
  '-user'    => 'anonymous',
  '-pass'    => '',
  '-dbname'  => 'homo_sapiens_variation_59_37d'
);

The minimum requirement is a connection to a core database for your species,
such as homo_sapiens_core_59_37d in the example above. A variation database is
required to find existing co-located variations, and a function genomics
database is required to assess if variations fall in regulatory regions.



5.0 Notes
=========

Run time of the script is generally proportional to the number of variants given
in the input file. Other dependent factors include:

- route of connection to the database : if using the public Ensembl database
  server at ensembldb.ensembl.org, the connection may be slow since data needs
  to be sent through the internet. If using large numbers of variants, it may be
  worth installing a copy of the Ensembl databases locally. They can be
  downloaded from the Ensembl FTP site - see
  http://www.ensembl.org/info/data/ftp/index.html for details.

- size of core database : the speed of retrieving results will also scale
  approximately linearly with the number of transcripts annotated in your
  species of interest.

- buffer size : if large amounts of memory are avaiable, increasing the buffer
  size using the "-b" flag may provide significant speed benefits for large
  numbers of variations. Using a larger buffer size increases the number of
  simultaneous database transactions and reduces the number of transactions
  overall.
  
- checking for existing variations. If not interested in existing co-located
  variations, this process can be bypassed using the --check_existing=0 command
  line flag.
