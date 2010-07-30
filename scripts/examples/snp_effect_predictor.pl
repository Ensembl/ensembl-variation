#!/usr/bin/perl

=head1 LICENSE

  Copyright (c) 1999-2010 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

SNP Effect Predictor - a script to predict the consequences of variations

by Will McLaren (wm2@ebi.ac.uk)
=cut

use strict;
use Getopt::Long;
use FileHandle;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;

# get command-line options
my ($in_file, $out_file, $buffer_size, $species, $registry_file, $help, $host, $user, $password);

our ($most_severe, $check_ref, $check_existing, $hgnc, $input_format);

my $args = scalar @ARGV;

GetOptions(
	'input_file=s'     => \$in_file,
	'output_file=s'    => \$out_file,
	'species=s'		   => \$species,
	'buffer_size=s'	   => \$buffer_size,
	'registry=s'	   => \$registry_file,
	'db_host=s'		   => \$host,
	'user=s'		   => \$user,
	'password=s'	   => \$password,
	'most_severe'	   => \$most_severe,
	'check_ref'        => \$check_ref,
	'check_existing=i' => \$check_existing,
	'hgnc'             => \$hgnc,
	'help'			   => \$help,
	'format=s'         => \$input_format,
);

# set defaults
$out_file ||= "snp_effect_output.txt";
$species ||= "human";
$buffer_size ||= 500;
$host ||= 'ensembldb.ensembl.org';
$user ||= 'anonymous';
$check_existing = 1 unless defined $check_existing;

# check file format
if(defined $input_format) {
	die "ERROR: Unrecognised input format specified $input_format\n" unless $input_format =~ /pileup|vcf/i;
}

$input_format ||= 'guess';

# print usage message if requested or no args supplied
if(defined($help) || !$args) {
	&usage;
	exit(0);
}


# get registry
my $reg = 'Bio::EnsEMBL::Registry';

# load DB options from registry file if given
if(defined($registry_file)) {
	$reg->load_all($registry_file);
}

# otherwise manually connect to DB server
else {
	if(defined($password)) {
		$reg->load_registry_from_db(-host => $host,-user => $user, -pass => $password);
	}
	
	else {
		$reg->load_registry_from_db(-host => $host,-user => $user);
	}
}

# get a meta container adaptors to check version
#my $core_mca = $reg->get_adaptor($species, 'core', 'metacontainer');
#my $var_mca = $reg->get_adaptor($species, 'variation', 'metacontainer');
#
#print
#	"Connected to core version ", $core_mca->get_schema_version, " database ",
#	"and variation version ", $var_mca->get_schema_version, " database\n";


# define the filehandle to read input from
my $in_file_handle = new FileHandle;

if(defined($in_file)) {
	
	# check defined input file exists
	die("ERROR: Could not find input file ", $in_file, "\n") unless -e $in_file;
	
	$in_file_handle->open($in_file) or die("ERROR: Could not read from input file ", $in_file, "\n");
}

# no file specified - try to read data off command line
else {
	$in_file_handle->open(<STDIN>);
}

# define filehandle to write to
my $out_file_handle = new FileHandle;
$out_file_handle->open(">$out_file") or die("ERROR: Could not write to output file ", $out_file, "\n");

# write the header line
print $out_file_handle
	"Uploaded Variation\tLocation\tGene\tTranscript\t",
	"Consequence\tPosition in cDNA\tPosition in protein\t",
	"Amino acid change\tCorresponding Variation\n";



# get adaptors
my $vfa = $reg->get_adaptor($species, 'variation', 'variationfeature');
my $tva = $reg->get_adaptor($species, 'variation', 'transcriptvariation');

# get fake ones for species with no var DB
if(!defined($vfa)) {
	$vfa = Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor->new_fake($species);
}

if(!defined($tva)) {
	$tva = Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor->new_fake($species);
}

my $sa = $reg->get_adaptor($species, 'core', 'slice');
our $ga = $reg->get_adaptor($species, 'core', 'gene');

# check we got slice adaptor - can't continue without a core DB
die("ERROR: Could not connect to core database\n") unless defined $sa and defined $ga;



# create a hash to hold slices so we don't get the same one twice
my %slice_hash = ();
my @new_vfs;

my $line_number = 0;

# read the file
while(<$in_file_handle>) {
  chomp;
  
  $line_number++;
  
  # header line?
  next if /^\#/;
  
  my ($chr, $start, $end, $allele_string, $strand, $var_name) = parse_line($_);
  
  # non-variant line from VCF
  next if $chr eq 'non-variant';
  
  #print "$line_number ", (join " ", ($chr, $start, $end, $allele_string, $strand, $var_name));
  #print "\n";
  
  # fix inputs
  $chr =~ s/chr//ig;
  $strand = ($strand =~ /\-/ ? "-1" : "1");
  
  # sanity checks
  unless($start =~ /^\d+$/ && $end =~ /^\d+$/) {
	warn("WARNING: Start $start or end $end coordinate invalid on line $line_number\n");
	next;
  }
  
  unless($allele_string =~ /([ACGT-]+\/*)+/) {
	warn("WARNING: Invalid allele string $allele_string on line $line_number\n");
	next;
  }
  
  # now get the slice
  my $slice;
 
  # check if we have fetched this slice already
  if(defined $slice_hash{$chr}) {
    $slice = $slice_hash{$chr};
  }
 
  # if not create a new one
  else {
	
	# first try to get a chromosome
    eval { $slice = $sa->fetch_by_region('chromosome', $chr); };
	
	# if failed, try to get any seq region
	if(!defined($slice)) {
		$slice = $sa->fetch_by_region(undef, $chr);
	}
	
	# if failed, warn and skip this line
	if(!defined($slice)) {
		warn("WARNING: Could not fetch slice named $chr on line $line_number\n");
		next;
	}	
	
	# store the hash
	$slice_hash{$chr} = $slice;
  }
  
  # check reference allele if requested
  if($check_ref) {
	my $ref_allele = (split /\//, $allele_string)[0];
	
	my $ok = 0;
	my $slice_ref_allele;
	
	# insertion, therefore no ref allele to check
	if($ref_allele eq '-') {
		$ok = 1;
	}
	else {
		my $slice_ref = $slice->sub_Slice($start, $end, $strand);
		
		if(!defined($slice_ref)) {
			warn "WARNING: Could not fetch sub-slice from $start\-$end\($strand\) on line $line_number";
		}
		
		else {
			$slice_ref_allele = $slice_ref->seq;
			$ok = ($slice_ref_allele eq $ref_allele ? 1 : 0);
		}
	}
	
	if(!$ok) {
		warn
			"WARNING: Specified reference allele $ref_allele ",
			"could not be matched",
			($slice_ref_allele ? " to $slice_ref_allele" : ""),
			" on line $line_number";
		next;
	}
  }
 
  # create a new VariationFeature object
  my $new_vf = Bio::EnsEMBL::Variation::VariationFeature->new(
    -start => $start,
    -end => $end,
    -slice => $slice,           # the variation must be attached to a slice
    -allele_string => $allele_string,
    -strand => $strand,
    -map_weight => 1,
    -adaptor => $vfa,           # we must attach a variation feature adaptor
    -variation_name => (defined $var_name ? $var_name : $chr.'_'.$start.'_'.$allele_string),
  );
  
  push @new_vfs, $new_vf;
  
  # if the array is "full" or there are no more items in @list
  if(scalar @new_vfs == $buffer_size) {
	
	# get consequences
	# results are stored attached to reference VF objects
	# so no need to capture return value here
	$tva->fetch_all_by_VariationFeatures(\@new_vfs);
	
	# now print the results to the file handle
	&print_consequences(\@new_vfs, $out_file_handle);
	
	# clear the array
	@new_vfs = ();
  }
}

# clean up any remaining
if(scalar @new_vfs) {
	$tva->fetch_all_by_VariationFeatures(\@new_vfs);
	&print_consequences(\@new_vfs, $out_file_handle);
}

sub print_consequences {
	
	my $vfs = shift;
	my $out_file_handle = shift;
	
	foreach my $new_vf(@$vfs) {
		
		# find any co-located existing VFs
		my $existing_vf = "-";
		
		if(defined($new_vf->adaptor->db) && $check_existing == 1) {
			my $fs = $new_vf->feature_Slice;
			if($fs->start > $fs->end) {
				($fs->{'start'}, $fs->{'end'}) = ($fs->{'end'}, $fs->{'start'});
			}
			foreach my $existing_vf_obj(@{$new_vf->adaptor->fetch_all_by_Slice($fs)}) {
				$existing_vf = $existing_vf_obj->variation_name
					if $existing_vf_obj->seq_region_start == $new_vf->seq_region_start
					and $existing_vf_obj->seq_region_end == $new_vf->seq_region_end;
			}
		}
		
		# the get_all_TranscriptVariations here now just retrieves the
		# objects that were attached above - it doesn't go off and do
		# the calculation again		
		foreach my $con(@{$new_vf->get_all_TranscriptVariations}) {
			foreach my $string(@{$con->consequence_type}) {
			  
			  if($con->cdna_start > $con->cdna_end) {
				($con->{'cdna_start'}, $con->{'cdna_end'}) = ($con->{'cdna_end'}, $con->{'cdna_start'});
			  }
			  
			  if($con->translation_start > $con->translation_end) {
				($con->{'translation_start'}, $con->{'translation_end'}) = ($con->{'translation_end'}, $con->{'translation_start'});
			  }
			  
			  my $gene = $ga->fetch_by_transcript_stable_id($con->transcript->stable_id);
			  
			  my $hgnc_name;
			  if($hgnc) {
				my @entries = grep {$_->database eq 'HGNC'} @{$gene->get_all_DBEntries()};
				$hgnc_name = $entries[0]->display_id || undef;
			  }
			  
			  print $out_file_handle
				$new_vf->variation_name, "\t",
				$new_vf->seq_region_name, ":",
				$new_vf->seq_region_start,
				($new_vf->seq_region_end == $new_vf->seq_region_start ? "" : "-".$new_vf->seq_region_end), "\t".
				($con->transcript ? $gene->stable_id.(defined $hgnc_name ? ";$hgnc_name" : "") : "-"), "\t",
				($con->transcript ? $con->transcript->stable_id : "-"), "\t",
				$string, "\t",
				($con->cdna_start ? $con->cdna_start.($con->cdna_end eq $con->cdna_start ? "" : "-".$con->cdna_end) : "-"), "\t",
				($con->translation_start ? $con->translation_start.($con->translation_end eq $con->translation_start ? "" : "-".$con->translation_end) : "-"), "\t",
				($con->pep_allele_string ? $con->pep_allele_string : "-"), "\t",
				$existing_vf, "\n";
			}
		}
	}
}


sub parse_line {
	my $line = shift;
	
	my @data = (split /\s+/, $_);
	
	# pileup: chr1 60 T A
	if(($input_format =~ /pileup/i) || ($data[0] =~ /(chr)?\w+/ && $data[1] =~ /\d+/ && $data[2] =~ /^[ACGTN-]+$/ && $data[3] =~ /^[ACGTN-]+$/)) {
		return ($data[0], $data[1], $data[1], $data[2]."/".$data[3], 1, (defined $data[4] ? $data[4] : undef));
	}
	
	# VCF: 20      14370   rs6054257 G     A      29    0       NS=58;DP=258;AF=0.786;DB;H2          GT:GQ:DP:HQ
	elsif(($input_format =~ /vcf/i) || ($data[0] =~ /(chr)?\w+/ && $data[1] =~ /\d+/ && $data[3] =~ /^[ACGTN-]+$/ && $data[4] =~ /^([\.ACGTN-]+\,?)+$/)) {
		
		# non-variant line in VCF, return dummy line
		if($data[4] eq '.') {
			return ('non-variant');
		}
		
		my ($start, $end, $alt) = ($data[1], $data[1], $data[4]);
		
		# indel?
		if($alt =~ /D|I/) {
			if($alt =~ /\,/) {
				warn "WARNING: can't deal with multiple different indel types in one variation";
				return ('non-variant');
			}
			
			# deletion
			if($alt =~ /D/) {
				my $num_deleted = $alt;
				$num_deleted =~ s/\D+//g;
				$end += $num_deleted - 1;
				$alt = "-";
			}
			
			# insertion
			else {
				$data[3] = '-';
				$alt =~ s/^I//g;
				$start++;
			}
		}
		
		$alt =~ s/\,/\//g;
		return ($data[0], $start, $end, $data[3]."/".$alt, 1, ($data[2] eq '.' ? undef : $data[2]));
	}
	
	# our format
	else {
		# we allow commas as delimiter so re-split
		@data = (split /\s+|\,/, $_);
		return @data;
	}
}

sub usage {
	my $usage =<<END;
#----------------------#
# SNP EFFECT PREDICTOR #
#----------------------#

By Will McLaren (wm2\@ebi.ac.uk)

Usage:
perl snp_effect_predictor.pl [arguments]

Options
--help                 Display this message and quit

-i | --input_file      Input file - if not specified, reads from STDIN
-o | --output_file     Output file [default: "snp_effect_output.txt"]
-f | --format          Input file format - one of "pileup", "vcf" [default: auto-detect]

-s | --species         Species to use [default: "human"]
-b | --buffer_size     Sets the number of SNPs sent in each batch [default: 500]
                       Increasing buffer size will retrieve results more quickly
                       but requires more memory.
--check_ref            If specified, checks supplied reference allele against stored
                       entry in Ensembl Core database [default: not used]
--check_existing=[0|1] If set to 1, checks for existing co-located variations in the
                       Ensembl Variation database [default: 1]
--hgnc                 If specified, HGNC gene identifiers are output alongside the
                       Ensembl Gene identifier [default: not used]

-d | --db_host         Manually define database host [default: "ensembldb.ensembl.org"]
-u | --user            Database username [default: "anonymous"]
-p | --password        Database password [default: not used]
-r | --registry_file   Registry file to use defines DB connections [default: not used]
                       Defining a registry file overrides above connection settings.
END

	print $usage;
}
