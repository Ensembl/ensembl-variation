#!/usr/bin/env perl

# tidy_gvf_dumps.pl - Tidy up the GVF dump directory before a release 

use strict;
use warnings;

use Getopt::Long;

my $toplevel_dir;

GetOptions(
    "output_dir|o=s" => \$toplevel_dir,
);

die "Usage: $0 --output_dir DIR\n" unless $toplevel_dir;

my @species_dirs = `ls $toplevel_dir`;

# slurp in the README text

my $readme;

{
    local $/ = undef;
    $readme = <DATA>;
}

# remove any .out or .err files and print the README in
# each species directory

for my $species (@species_dirs) {
    chomp $species;
    next unless $species;
    my $dir = "$toplevel_dir/$species";
    die "$dir doesn't look like a GVF dump directory\n" 
        unless grep {/\.gvf\.gz/} `ls $dir`;
    `rm -f $dir/*.out`;
    `rm -f $dir/*.err`;
    open README, ">$dir/README" or die "Failed to create README file\n";
    print README $readme, "\n";
}

__DATA__
This directory contains a gzip compressed GVF (Genome Variation Format) file
containing all germline variations from the current Ensembl release for this 
species, named <species>.gvf.gz. If this species has any variations failed 
by the Ensembl QC checks there will also be a file called 
<species>_failed.gvf.gz identifying these. The consequences of the variants, 
as called by the Ensembl variation consequence pipeline, will be included in 
the GVF files as of release 62, scheduled for April 2011.

For human we also provide a file containing all somatic mutations in the
database, files with germline variations observed in the Watson and Venter 
genomes along with their genotypes, and files containing allele frequencies
from several of the HapMap and 1000 genomes pilot study populations.

Please note that depending on the amount of variation data available for 
this species the uncompressed file may be very large (e.g. the entire 
germline file for human is ~3GB).

The data contained in these files is presented in GVF format, this is a
simple tab-delimited format which shows the location of each variant along
with the reference and variant sequences, an identifier for the source
of the data (typically a dbSNP rsID), and other relevant information (e.g. 
genotypes and allele freqeuncies), a short example is presented below. For 
more details about GVF please refer to:

Reese, M.G. et al. A standard variation file format for human genome sequences.
Genome Biology. 2010;11(8):R88 PMID: 20796305

and:

http://www.sequenceontology.org/gvf.html

Questions about these files can be addressed to the Ensembl helpdesk: 
helpdesk@ensembl.org, or to the developer's mailing list: dev@ensembl.org.

-----

Example content from the human germline GVF dump is shown below:

##gff-version 3
##file-date 2011-01-31
##genome-build ensembl GRCh37
##gvf-version 1.05
##feature-ontology http://song.cvs.sourceforge.net/viewvc/song/ontology/so.obo?revision=1.283
##data-source Source=ensembl;version=61;url=http://e61.ensembl.org/Homo_sapiens
##file-version 61
##sequence-region 11 1 135006516
11	dbSNP	SNV	61554	61554	.	+	.	ID=1;Variant_seq=C;Dbxref=dbSNP:rs77355429;Reference_seq=A
11	dbSNP	SNV	61645	61645	.	+	.	ID=2;Variant_seq=C;Dbxref=dbSNP:rs61869610;Reference_seq=G
11	dbSNP	SNV	61868	61868	.	+	.	ID=3;Variant_seq=A;Dbxref=dbSNP:rs365553;Reference_seq=G
11	dbSNP	SNV	67521	67521	.	+	.	ID=4;Variant_seq=G;Dbxref=dbSNP:rs76077193;Reference_seq=A
11	dbSNP	SNV	70073	70073	.	+	.	ID=5;Variant_seq=C;Dbxref=dbSNP:rs1099703;Reference_seq=T
11	dbSNP	SNV	70113	70113	.	+	.	ID=6;Variant_seq=T;Dbxref=dbSNP:rs76956518;Reference_seq=C
11	dbSNP	SNV	70133	70133	.	+	.	ID=7;Variant_seq=T;Dbxref=dbSNP:rs112541484;Reference_seq=C
11	dbSNP	SNV	70135	70135	.	-	.	ID=8;Variant_seq=A;Dbxref=dbSNP:rs4120101;Reference_seq=C
11	dbSNP	SNV	70146	70146	.	-	.	ID=9;Variant_seq=C;Dbxref=dbSNP:rs4120100;Reference_seq=A

