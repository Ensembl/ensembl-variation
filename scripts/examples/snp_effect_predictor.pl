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

GetOptions(
	'input_file=s'  => \$in_file,
	'output_file=s' => \$out_file,
	'species=s'		=> \$species,
	'buffer_size=s'	=> \$buffer_size,
	'registry=s'	=> \$registry_file,
	'db_host=s'		=> \$host,
	'user=s'		=> \$user,
	'password=s'	=> \$password,
	'help'			=> \$help,
);

# set defaults
$out_file ||= "snp_effect_output.txt";
$species ||= "human";
$buffer_size ||= 500;
$host ||= 'ensembldb.ensembl.org';
$user ||= 'anonymous';

if(defined($help)) {
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



# read the file
while(<$in_file_handle>) {
  chomp;
  
  my ($chr, $start, $end, $allele_string, $strand) = split /\s+|\,/, $_;
  
  # fix inputs
  $chr =~ s/chr//g;
  $strand = ($strand =~ /\-/ ? "-1" : "1");
  
  # sanity checks
  die("ERROR: Start $start or end $end coordinate invalid\n") unless $start =~ /^\d+$/ && $end =~ /^\d+$/;
  die("ERROR: Invalid allele string $allele_string\n") unless $allele_string =~ /([ACGT-]+\/*)+/;  
  
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
	
	# if failed, die
	if(!defined($slice)) {
		die("ERROR: Could not fetch slice named $chr\n");
	}	
	
	# store the hash
	$slice_hash{$chr} = $slice;
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
    -variation_name => $chr.'_'.$start.'_'.$allele_string,
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

&print_consequences(\@new_vfs, $out_file_handle);

sub print_consequences {
	
	my $vfs = shift;
	my $out_file_handle = shift;
	
	foreach my $new_vf(@$vfs) {
		
		# find any co-located existing VFs
		my $existing_vf = "N/A";
		
		if(defined($new_vf->adaptor->db)) {
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
			  
			  print $out_file_handle
				$new_vf->variation_name, "\t",
				$new_vf->seq_region_name, ":",
				$new_vf->seq_region_start,
				($new_vf->seq_region_end == $new_vf->seq_region_start ? "" : "-".$new_vf->seq_region_end), "\t".
				($con->transcript ? $ga->fetch_by_transcript_stable_id($con->transcript->stable_id)->stable_id : "N/A"), "\t",
				($con->transcript ? $con->transcript->stable_id : "N/A"), "\t",
				$string, "\t",
				($con->cdna_start ? $con->cdna_start.($con->cdna_end eq $con->cdna_start ? "" : "-".$con->cdna_end) : "N/A"), "\t",
				($con->translation_start ? $con->translation_start.($con->translation_end eq $con->translation_start ? "" : "-".$con->translation_end) : "N/A"), "\t",
				($con->pep_allele_string ? $con->pep_allele_string : "N/A"), "\t",
				$existing_vf, "\n";
			}
		}
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
-h | --help           Display this message and quit

-s | --species        Species to use [default: "human"]
-i | --input_file     Input file - if not specified, reads from STDIN
-o | --output_file    Output file [default: "snp_effect_output.txt"]
-b | --buffer_size    Sets the number of SNPs sent in each batch [default: 500]
                      Increasing buffer size will retrieve results more quickly
					  but requires more memory.

-d | --db_host        Manually define database host [default: "ensembldb.ensembl.org"]
-u | --user           Database username [default: "anonymous"]
-p | --password       Database password [default: not used]
-r | --registry_file  Registry file to use defines DB connections [default: not used]
                      Defining a registry file overrides above connection settings.
END

	print $usage;
}