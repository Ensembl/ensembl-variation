############################
#                          #
# Variant Effect Predictor #
#                          #
############################

Copyright (c) 1999-2011 The European Bioinformatics Institute and
Genome Research Limited.  All rights reserved.

This software is distributed under a modified Apache license.
For license details, please see

http://www.ensembl.org/info/about/code_licence.html

Please email comments or questions to the public Ensembl
developers list at <ensembl-dev@ebi.ac.uk>.

Questions may also be sent to the Ensembl help desk at
<helpdesk@ensembl.org>
  
  

0.0 What's new
==============

Version 2.0 of the Variant Effect Predictor script (VEP) constitutes a complete
overhaul of both the script and the API behind it. It requires at least version
62 of the Ensembl API to function. Here follows a summary of the changes:

- support for SIFT, PolyPhen and Condel non-synonymous predictions in human

- per-allele and compound consequence types

- support for Sequence Ontology (SO) and NCBI consequence terms

- modified output format
  - support for new output fields in Extra column
  - header section containing information on database and software versions
  - codon change shown in output
  - CDS position shown in output
  - option to output Ensembl protein identifiers
  - option to output HGVS nomenclature for variants
  
- support for gzipped input files
  
- enhanced configuration options, including the ability to read configuration
  from a file

- verbose output now much more useful

- whole-genome mode now more stable

- finding existing co-located variations now ~5x faster



1.0 Requirements
================

Version 2.0 of the script requires at least version 62 of the Ensembl Core and
Variation APIs and their relevant dependencies to be installed. See
http://www.ensembl.org/info/docs/api/index.html for details. No explicit
installation of the script is required.



2.0 Introduction
================

Given a set of variant positions and alleles, it is often useful to know how
these variants may affect transcripts and other genetic factors. Ensembl
provides such annotation for the variants stored in its Ensembl Variation
database; using this script it is possible to do the same for novel variants.



3.0 Running the script
======================

The script is run on the command line as follows:

perl variant_effect_predictor.pl [options]

where [options] represent a set of flags and options to the script. These can be
listed using the flag --help:

perl variant_effect_predictor.pl --help

By default the script connects to the public Ensembl database server at
ensembldb.ensembl.org; other connection options are available.


3.1 Basic options
-----------------

--help : display help message and quit

-v (--verbose) : Output status messages as the script runs. Not used by default

-q (--quiet) : Suppress status and warning messages. Not used by default

-c (--config) : Load configuration options from a config file. The config file
should consist of whitespace-separated pairs of option names and settings e.g.:

output_file   my_output.txt
species       mus_musculus
format        vcf
whole_genome  1

This is useful if you find yourself using the same configuration options each
time. You can create a quick version file of this by setting the flags as normal
and running the script in verbose (-v) mode. This will output lines that can be
copied to a config file that can be loaded in on the next run using -c. Note
that any options specified in the normal way overwrite those in the config file.
Not used by default


3.2 Input options
-----------------

-s (--species) : species for your data. This can be the latin name e.g.
"homo_sapiens" or any Ensembl alias e.g. "mouse". Specifying the latin name can
speed up initial database connection as the registry does not have to load all
available database aliases on the server. Default = "human"

-i (--input_file) : input file name. If not specified, the script will attempt
to read from STDIN.

--format : input file format. By default, the script auto-detects the input
file format. Using this option you can force the script to read the input
file as VCF or pileup format. Not used by default.
 
-o (--output_file) : output file name. Default = "variant_effect_output.txt"


3.3 Database options
--------------------

-d (--db_host) : manually define the database host to connect to.
Default = "ensembldb.ensembl.org"

-u (--user) : manually define the database username used when connecting to the
database host. Default = "anonymous"

-p (--password) : manually define the database password used when connecting to
the database host. Not used by default

-r (--registry_file) : defining a registry file (see 4.3 below) overwrites other
connection settings and uses those found in the specified registry file to
connect.
  
-g (--genomes) : override the default connection settings with those for the
Ensembl Genomes public MySQL server. Required when using any of the Ensembl
Genomes species. Not used by default

--db_version [version]: force the script to connect to a specific version of the
Ensembl databases. Not recommended as there will usually be conflicts between
software and database versions. Not used by default


3.4 Output options
------------------

-t (--terms) [ensembl|so|ncbi] : the type of consequence terms to output. The
Ensembl terms are described at
http://www.ensembl.org/info/docs/variation/index.html. The sequence ontology
(http://www.sequenceontology.org/) is a joint effort by genome annotation
centres to standardise descriptions of biological sequences. The NCBI terms are
those used by dbSNP (http://www.ncbi.nlm.nih.gov/projects/SNP/), and are the
least complete set - where no NCBI term is available, the script will output the
Ensembl term. Default : ensembl

--sift [p|s|b] : SIFT (http://sift.jcvi.org/) predicts whether an amino acid
substitution affects protein function based on sequence homology and the
physical properties of amino acids. The VEP can output the [p]rediction term,
[s]core or [b]oth. Using this may increase the run time of the script. Not used
by default. Human only

--polyphen [p|s|b] : PolyPhen (http://genetics.bwh.harvard.edu/pph/) is a tool
which predicts possible impact of an amino acid substitution on the structure
and function of a human protein using straightforward physical and comparative
considerations. The VEP can output the [p]rediction term, [s]core or [b]oth.
Using this may increase the run time of the script. Not used by default. Human
only

--condel [p|s|b] : Condel (https://bg.upf.edu/forge/wiki/condel/Condel) computes
a weighed average of the scores (WAS) of several computational tools aimed at
classifying missense mutations as likely deleterious or likely neutral. The VEP
currently presents a Condel WAS from SIFT and PolyPhen. The VEP can output the
[p]rediction term, [s]core or [b]oth. Using this may increase the run time of
the script. Not used by default. Human only

--hgvs : add HGVS (http://www.hgvs.org/mutnomen/) nomenclature based on Ensembl
stable identifiers to the output. Both coding and protein sequence names are
added where appropriate. Not used by default

--protein : add the Ensembl protein identifier to the output where appropriate.
Not used by default

--hgnc : adds the HGNC gene identifer (where available) to the output. Not used
by default


3.5 Filtering and QC options
----------------------------

--check_ref : force the script to check the supplied reference allele against
the sequence stored in the Ensembl Core database. Lines that do not match are
skipped. Not used by default.

--coding_only : only return consequences that fall in the coding regions of
transcripts. Not used by default

--check_existing [0|1] : by default the script checks for the existence of
variants that are co-located with your input. Disabling this by specifying
--check_existing=0 will bypass this and provide a large speed boost. Default : 1

--failed [0|1] : when checking for co-located variants, by default the script
will include variants that have been flagged as failed. Set this flag to 0 to
exclude such variants. Default : 1


3.6 Whole-genome mode options
-----------------------------

-w (--whole_genome) : Setting this flag forces the script to run in whole-genome
mode. This should only be used for appropriate datasets - see section on
whole-genome mode for more details. Not used by default

-b (--buffer_size) : sets the internal buffer size, corresponding to the number
of variations that are used to query the database simultaneously. Set this
lower to use less memory at the expense of longer run time, and higher to use
more memory with a faster run time. Default = 5000
  
--chunk_size : Sets the chunk size of the internal data structure used in
whole-genome mode. Only change this if you know what you are doing!
Default : 50kb
  

3.7 Examples
------------

- Read input from STDIN, run in verbose mode

perl variant_effect_predictor.pl -v


- Input file variants.vcf.txt, input file format VCF, add HGNC gene identifiers,
  output SO consequence terms

perl variant_effect_predictor.pl -i variants.vcf.txt -f vcf -hgnc -t so


- Output file variants_output.txt, don't check for existing co-located variants,
  output only coding sequence consequences, output HGVS names
  
perl variant_effect_predictor.pl -i variants.txt -o variants_output.txt \
-check_existing 0 -coding_only -hgvs


- Specify DB connection parameters in registry file ensembl.registry, add SIFT
  score and prediction, PolyPhen prediction, Condel score

perl variant_effect_predictor.pl -i variants.txt -r ensembl.registry -sift b \
-polyphen p -condel s


- Connect to Ensembl Genomes db server for A.thaliana, run in whole-genome mode
  with buffer size of 10000

perl variant_effect_predictor.pl -i variants.txt -genomes -species \
arabidopsis_thaliana -w -b 10000


- Load config from ini file, run in quiet mode

perl variant_effect_predictor.pl -c vep.ini -i variants.vcf.txt -q



4.0 File formats
================

The Variant Effect Predictor script uses plain text files both as input and
output. Input files can be gzip compressed - the zcat utility must be in your
path to use gzipped files.

4.1 Input file
--------------

The script now supports VCF version 4 as input - see the URL below for details:

http://www.1000genomes.org/wiki/Analysis/variant-call-format

The default input file format consists of five columns; these can be comma, tab
(or any whitespace) separated. The columns are:

- Chromosome name : the name of the chromosome to which the variant maps e.g. 1
  for chromosome 1. This can also be the name of a contig or other seq_region
  for species with unassembled genomes.
  
- Start position

- End position : along with start position, defines the location of the variant.
  These are 1-indexed (i.e. the first base of the chromosome is base 1), and
  should be defined according to the following rules:
  
  - For a SNP, start should be the same as end
  - For a multi-basepair substitution, start should be the first and end should
    be the last base affected (e.g. for a 3 base substitution, end = start + 2).
    This also applies to any unbalanced substitution - for example replacing 3
    bases with 5, the coordinates should reflect the 3 bases spanned by the
    reference allele.
  - For an insertion relative to the reference, start = end + 1, regardless of
    its length. This means that in this instance the start coordinate will be
    greater than the end coordinate; end represents the base immediately 5' of
    the insertion site, and start the base immediately 3'.
  - For a deletion relative to the reference, follow the same rules for a
    multi-basepair substitution.
  
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

- Substitution of CATTCC for GGC at base 100 on reverse strand of chromosome Y:
Y 100 102 GGC/CATTCC -

- Deletion of bases 100, 101 and 102 (CTG) on forward strand of chromosome 3:
3 100 102 CTG/- +

- Insertion of ACCG between bases 100 and 101 on reverse strand of chromosome X:
X 101 100 -/ACCG -

Other input formats are also supported; VCF (see
http://1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcfv3.2) and
pileup. These should be auto-detected by the script; if they are not, you can
force the script to use either with the -f|--format command line flag.


4.2 Output file
---------------

4.2.1 Header
''''''''''''

The output file as of version 2.0 contains a header section denoted by lines
starting with "## ". This header contains the following information (where
appropriate):

- the script version
- the time and date that the output file was started (not finished!)
- the name of the database and server connected to
- the API and database versions
- the key descriptions for the key/value pairs found in the Extra column

4.2.2 Main section
''''''''''''''''''

The output file contains one line for each predicted consequence; if a variant
falls in more than one transcript it may produce more than one line. The column
headers are given in line starting "#". The output file is a tab-delimited file
with the following columns:

- Uploaded_variation : a temporary name assigned to the variant, based on the
  position and alleles, or the optional name if specified

- Location : string corresponding to the variant's location, represented as
  Chr:Start-End
  
- Allele : the variant allele used to calculate the consequence. Generally only
  relevant when the input variant had multiple alternate alleles.

- Gene : the Ensembl stable ID of the gene where the variant is located (where
  applicable)

- Transcript : the Ensembl stabled ID of the transcript where the variant is
  located (where applicable)

- Consequence : the consequence type predicted. By default this will be the
  Ensembl (http://www.ensembl.org/info/docs/variation/index.html) term. It is
  possible to configure the script to output SO
  (http://www.sequenceontology.org/) or NCBI
  (http://www.ncbi.nlm.nih.gov/projects/SNP/) terms in place using the
  -t|--terms command line flag. Compound consequences separated by commas are
  given where appropriate. For example, a variant falls in the intron of a
  non-coding gene, it will be represented (in Ensembl terms) as
  INTRONIC,WITHIN_NON_CODING_GENE

- cDNA_position : the position of the variant in the transcript's unspliced cDNA
  sequence (if applicable)

- CDS_position : the position of the variant in the transcript's coding
  sequence (if applicable)

- Protein_position : the position of the variant in the resultant protein

- Amino_acids : a "/"-separated string of possible amino acids generated
  by the variant, with the reference amino acid first
  
- Codons : a "/"-separated string of possible codons generated by the variant,
  with the reference codon first. The variant allele itself is in upper case
  with the other bases of the codon in lower case text

- Existing_variation : if an existing variant is found in the Ensembl
  Variation database whose position matches exactly that of the entered
  variation, its identifier appears here. Where variants from multiple sources
  are found, dbSNP identifiers are given preference.
  
- Extra : this column contains a ";"-separated list of key-value pairs of the
  form "key=value;", containing any of the following pieces of information:
  
  HGNC : the HGNC gene identifier
  ENSP : the Ensembl protein identifier of the affected transcript
  HGVSc : the HGVS coding sequence name
  HGVSp : the HGVS protein sequence name
  SIFT : the SIFT prediction and/or score, with both given as prediction(score)
  PolyPhen : the PolyPhen prediction and/or score
  Condel : the Condel consensus prediction and/or score


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
such as homo_sapiens_core_62_37g in the example above. A variation database is
required to find existing co-located variations.



5.0 Whole-genome mode
=====================

Whole-genome mode allows the analysis of dense datasets that cover an entire
genome, chromosome, gene or set of genes, such as those produced from
resequencing projects.

In "standard" mode the script, for each variant, finds any transcripts that
overlap and works out the consequences of the variant on each of these
transcripts. Users with data more sparsely scattered across the genome should
generally stick to using the default run mode of the VEP.

In whole-genome mode, the script takes the opposite approach - it retrieves all
transcripts for a chromosome (or contig as appropriate), then finds overlapping
variants and works out the consequences. This approach is much faster for dense
datasets since much of the data calculated on the fly per transcript can be
cached and reused for each variant that overlaps it. Hence, users with data
covering just one gene, for example, will find using whole-genome mode many many
times faster than standard mode.

The restriction of this approach is that it depends on having a database-like
approach to retrieving the variants that overlap each transcript. To this end,
the script implements a chunked, or stratified, internal data structure to store
the variants. It is also strongly recommended that the input file(s) be sorted
by chromosome and position before running the script.

By default, the script stores variants in chunks of 50kb (roughly the average
length of a transcript in human). This means that generally only between 1 and 3
"chunks" of variants need to be checked for overlap per transcript. This chunk
size can be modified via the --chunk_size option.

The buffer size setting (--buffer_size) controls how many variants are stored in
memory before scanning the genome for transcripts that overlap them. Setting the
buffer size higher than the default of 5000 will yield faster performance but
will use more memory.

Using a locally installed core database will also dramatically improve
performance versus connecting to the public Ensembl MySQL server.

By default, checking for existing variants is disabled, as is output of the
Ensembl Gene identifier.

Summary of recommendations for using whole-genome mode:

- sort your data by chromosome and position

- install a local copy of the relevant core database to connect to

- ensure your system has a large amount of free memory (1GB+ recommended)

- if possible, run in parallel - divide the data by chromosome and run one
  parallel process per chromosome (considering memory requirements as
  appropriate!)



6.0 Notes
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
  
- checking for existing variations. If not interested in existing co-located
  variations, this process can be bypassed using the --check_existing=0 command
  line flag. Disabling this option saves considerable time.

- outputting HGVS names, as well as SIFT, PolyPhen and Condel predictions can
  significantly increase run time for datasets with large numbers of
  non-synonymous variants.
