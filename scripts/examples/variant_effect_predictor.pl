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

Variant Effect Predictor - a script to predict the consequences of genomic variants

by Will McLaren (wm2@ebi.ac.uk)
=cut

use strict;
use Getopt::Long;
use FileHandle;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);

# debug
#use Devel::Size qw(total_size);
use Time::HiRes qw(gettimeofday tv_interval);
our %times;



# get command-line options
my ($in_file, $out_file, $buffer_size, $species, $registry_file, $help, $host, $user, $password, $tmpdir, $db_version, $regulation, $include_failed);

our ($most_severe, $check_ref, $check_existing, $hgnc, $input_format, $whole_genome, $chunk_size);

my $args = scalar @ARGV;

GetOptions(
	'input_file=s'     => \$in_file,
	'output_file=s'    => \$out_file,
	'species=s'		   => \$species,
	'buffer_size=i'	   => \$buffer_size,
	'registry=s'	   => \$registry_file,
	'db_host=s'		   => \$host,
	'user=s'		   => \$user,
	'password=s'	   => \$password,
	'most_severe'	   => \$most_severe,
	'check_ref'        => \$check_ref,
	'check_existing=i' => \$check_existing,
	'failed=i'         => \$include_failed,
	'hgnc'             => \$hgnc,
	'help'			   => \$help,
	'format=s'         => \$input_format,
	'whole_genome'     => \$whole_genome,
	'tmp_dir=s'        => \$tmpdir,
	'version=i'        => \$db_version,
	'chunk_size=s'     => \$chunk_size,
);

# set defaults
$out_file    ||= "variant_effect_output.txt";
$species     ||= "human";
$buffer_size ||= 500;
$chunk_size  ||= '50kb';
$host        ||= 'ensembldb.ensembl.org';
$user        ||= 'anonymous';
$tmpdir      ||= '/tmp';

$include_failed = 1 unless defined $include_failed;
$check_existing = 1 unless (defined $check_existing || defined $whole_genome);
$chunk_size =~ s/mb?/000000/i;
$chunk_size =~ s/kb?/000/i;

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
	$reg->load_registry_from_db(
		-host       => $host,
		-user       => $user,
		-pass       => $password,
		-db_version => $db_version,
	);
}

## get a meta container adaptors to check version
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
	$in_file_handle = 'STDIN';
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
our $vfa = $reg->get_adaptor($species, 'variation', 'variationfeature');
our $tva = $reg->get_adaptor($species, 'variation', 'transcriptvariation');

# get fake ones for species with no var DB
if(!defined($vfa)) {
	$vfa = Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor->new_fake($species);
}

if(!defined($tva)) {
	$tva = Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor->new_fake($species);
}

$vfa->db->include_failed_variations($include_failed) if $vfa->db->can('include_failed_variations');

our $sa = $reg->get_adaptor($species, 'core', 'slice');
our $ga = $reg->get_adaptor($species, 'core', 'gene');

# check we got slice adaptor - can't continue without a core DB
die("ERROR: Could not connect to core database\n") unless defined $sa and defined $ga;

# regulatory stuff
our ($fsa, $mfa, $fss);
if($regulation) {
	$fsa = $reg->get_adaptor($species, 'funcgen', 'featureset');
	$mfa = $reg->get_adaptor($species, 'funcgen', 'motiffeature');
	
	if(defined $fsa and defined $mfa) {
		$fss = $fsa->fetch_all();
	}
	
	else {
		warn("WARNING: Could not get functgenomics adaptors");
		$regulation = undef;
	}
}


# create a hash to hold slices so we don't get the same one twice
my %slice_hash = ();
my @new_vfs;
my %vf_hash;

our $transcript_cache;

my $line_number = 0;
my $vf_count;

# read the file
while(<$in_file_handle>) {
  chomp;
  
  $line_number++;
  
  # header line?
  next if /^\#/;
  
  # some lines (pileup) may actually parse out into more than one variant)
  foreach my $sub_line(@{parse_line($_)}) {
	
	# get the sub-line into named variables
	my ($chr, $start, $end, $allele_string, $strand, $var_name) = @{$sub_line};
	
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
	
	$vf_count++;
	
	if($whole_genome) {
		push @{$vf_hash{$chr}{int($start / $chunk_size)}{$start}}, $new_vf;
	}
	else {
		push @new_vfs, $new_vf;
	}
	
	# process if buffer "full"
	if($vf_count == $buffer_size) {
	  
	  if($whole_genome) {
		&print_consequences(&whole_genome_fetch(\%vf_hash), $out_file_handle);
		%vf_hash = ();
	  }
	  
	  else {
		# get consequences
		# results are stored attached to reference VF objects
		# so no need to capture return value here
		$tva->fetch_all_by_VariationFeatures(\@new_vfs);
		
		# now print the results to the file handle
		&print_consequences(\@new_vfs, $out_file_handle);
		
		# do regulation stuff
		if($regulation) {
			foreach my $vf(@new_vfs) {
				foreach my $mf(@{$mfa->fetch_all_by_Slice_FeatureSets($vf->feature_Slice, $fss)}) {
					print $vf->variation_name, " ", $mf->display_label, "\n";
				}
			}
		}
		
		# clear the array
		@new_vfs = ();
	  }
	  
	  $vf_count = 0;
	}
  }
}

# clean up any remaining
if(scalar @new_vfs) {
	$tva->fetch_all_by_VariationFeatures(\@new_vfs);
	&print_consequences(\@new_vfs, $out_file_handle);
}

# if in whole-genome mode
if($whole_genome && %vf_hash) {
	&print_consequences(&whole_genome_fetch(\%vf_hash), $out_file_handle);
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
			  
			  my $gene = ($con->transcript ? $ga->fetch_by_transcript_stable_id($con->transcript->stable_id) : undef) unless $whole_genome;
			  
			  my $hgnc_name;
			  if($hgnc && $gene) {
				my @entries = grep {$_->database eq 'HGNC'} @{$gene->get_all_DBEntries()};
				$hgnc_name = (scalar @entries ? $entries[0]->display_id : undef);
			  }
			  
			  print $out_file_handle
				$new_vf->variation_name, "\t",
				$new_vf->seq_region_name, ":",
				$new_vf->seq_region_start,
				($new_vf->seq_region_end == $new_vf->seq_region_start ? "" : "-".$new_vf->seq_region_end), "\t".
				($gene ? $gene->stable_id.(defined $hgnc_name ? ";$hgnc_name" : "") : "-"), "\t",
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
	if(($input_format =~ /pileup/i) || ($data[0] =~ /(chr)?\w+/ && $data[1] =~ /\d+/ && $data[2] =~ /^[ACGTN-]+$/ && $data[3] =~ /^[ACGTNRYSWKM*+\/-]+$/)) {
		my @return = ();
		
		if($data[2] ne "*"){
			(my $var = unambiguity_code($data[3])) =~ s/$data[2]//ig;
			if(length($var)==1){
				push @return, [$data[0], $data[1], $data[1], $data[2]."/".$var, 1, undef];
			}
			else{
				for my $nt(split //,$var){
					push @return, [$data[0], $data[1], $data[1], $data[2]."/".$nt, 1, undef];
				}
			}
		}
		else{ #indel
			my @genotype=split /\//,$data[3];
			foreach my $allele(@genotype){
				if(substr($allele,0,1) eq "+") { #ins
					push @return, [$data[0], $data[1]+1, $data[1], "-/".substr($allele,1), 1, undef];
				}
				elsif(substr($allele,0,1) eq "-"){ #del
					push @return, [$data[0], $data[1], $data[1]+length($data[3])-4, substr($allele,1)."/-", 1, undef];			
				}
				elsif($allele ne "*"){
					warn("WARNING: invalid pileup indel genotype: $line\n");
					push @return, ['non-variant'];
				}
			}
		}
		return \@return;
	}
	
	# VCF: 20      14370   rs6054257 G     A      29    0       NS=58;DP=258;AF=0.786;DB;H2          GT:GQ:DP:HQ
	elsif(($input_format =~ /vcf/i) || ($data[0] =~ /(chr)?\w+/ && $data[1] =~ /\d+/ && $data[3] =~ /^[ACGTN-]+$/ && $data[4] =~ /^([\.ACGTN-]+\,?)+$/)) {
		
		# non-variant line in VCF, return dummy line
		if($data[4] eq '.') {
			return [['non-variant']];
		}
		
		# get relevant data
		my ($start, $end, $ref, $alt) = ($data[1], $data[1], $data[3], $data[4]);
		
		# adjust end coord
		$end += (length($ref) - 1);
		
		# find out if any of the alt alleles make this an insertion or a deletion
		my ($is_indel, $is_sub, $ins_count, $total_count);
		foreach my $alt_allele(split /\,/, $alt) {
			$is_indel = 1 if $alt_allele =~ /D|I/;
			$is_indel = 1 if length($alt_allele) != length($ref);
			$is_sub = 1 if length($alt_allele) == length($ref);
			$ins_count++ if length($alt_allele) > length($ref);
			$total_count++;
		}
		
		# multiple alt alleles?
		if($alt =~ /\,/) {
			if($is_indel) {
				
				my @alts;
				
				if($alt =~ /D|I/) {
					foreach my $alt_allele(split /\,/, $alt) {
						# deletion (VCF <4)
						if($alt_allele =~ /D/) {
							push @alts, '-';
						}
						
						elsif($alt_allele =~ /I/) {
							$alt_allele =~ s/^I//g;
							push @alts, $alt_allele;
						}
					}
				}
				
				else {
					$ref = substr($ref, 1);
					$ref = '-' if $ref eq '';
					$start++;
					
					foreach my $alt_allele(split /\,/, $alt) {
						$alt_allele = substr($alt_allele, 1);
						$alt_allele = '-' if $alt_allele eq '';
						push @alts, $alt_allele;
					}
				}
				
				$alt = join "/", @alts;
			}
			
			else {
				# for substitutions we just need to replace ',' with '/' in $alt
				$ref =~ s/\,/\//;
			}
		}
		
		else {
			if($is_indel) {
				# deletion (VCF <4)
				if($alt =~ /D/) {
					my $num_deleted = $alt;
					$num_deleted =~ s/\D+//g;
					$end += $num_deleted - 1;
					$alt = "-";
					$ref .= ("N" x ($num_deleted - 1)) unless length($ref) > 1;
				}
				
				# insertion (VCF <4)
				elsif($alt =~ /I/) {
					$ref = '-';
					$alt =~ s/^I//g;
					$start++;
				}
				
				# insertion or deletion (VCF 4+)
				else {
					# chop off first base
					$ref = substr($ref, 1);
					$alt = substr($alt, 1);
					
					$start++;
					
					if($ref eq '') {
						# make ref '-' if no ref allele left
						$ref = '-';
						
						# extra adjustment required for Ensembl
						$start++;
					}
					
					# make alt '-' if no alt allele left
					$alt = '-' if $alt eq '';
				}
			}
		}
		
		return [[$data[0], $start, $end, $ref."/".$alt, 1, ($data[2] eq '.' ? undef : $data[2])]];
		
	}
	
	# our format
	else {
		# we allow commas as delimiter so re-split
		@data = (split /\s+|\,/, $_);
		return [\@data];
	}
}


sub whole_genome_fetch {
	my $vf_hash = shift;
	
	my $up_down_size = $Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor::UP_DOWN_SIZE;
	
	my %vf_done;
	my @return;
	
	&start("transcripts");
	
	foreach my $chr(keys %$vf_hash) {
		my $slice = $sa->fetch_by_region('chromosome', $chr);
		
		warn "Analyzing chromosome $chr";
		
		if(defined($transcript_cache->{$chr})) {
			warn " - USING CACHE";
		}
		
		$transcript_cache->{$chr} = $slice->get_all_Transcripts() unless defined($transcript_cache->{$chr});
		
		my $tr_count = scalar @{$transcript_cache->{$chr}};
		
		warn " - fetched $tr_count transcripts";
		
		my $tr_counter;
		
		while($tr_counter < $tr_count) {
			
			my $tr = $transcript_cache->{$chr}->[$tr_counter++];
			
			if($tr_counter =~ /(5|0)00$/) {
				warn " - analysed $tr_counter\/$tr_count transcripts on chr $chr";
				&end("transcripts");
				&start("transcripts");
			}
			
			# do each overlapping VF
			my $chr = $tr->seq_region_name;
			my $s = $tr->seq_region_start - $up_down_size;
			my $e = $tr->seq_region_end + $up_down_size;
			
			# get the chunks this transcript overlaps
			my %chunks;
			$chunks{$_} = 1 for (int($s/$chunk_size)..int($e/$chunk_size));
			map {delete $chunks{$_} unless defined($vf_hash{$chr}{$_})} keys %chunks;
			
			foreach my $chunk(keys %chunks) {				
				foreach my $pos(grep {$_ >= $s && $_ <= $e} keys %{$vf_hash{$chr}{$chunk}}) {
					foreach my $vf(@{$vf_hash{$chr}{$chunk}{$pos}}) {
						$tva->_calc_consequences($slice, $tr, $vf);
						
						if(!defined($vf_done{$vf->{'variation_name'}.'_'.$vf->{'start'}})) {
							push @return, $vf;
							$vf_done{$vf->{'variation_name'}.'_'.$vf->{'start'}} = 1;
						}
					}
				}
			}
		}
		
		# clean hash
		delete $vf_hash{$chr};
	}
	
	&end("transcripts");
	
	return \@return;
}

sub usage {
	my $usage =<<END;
#--------------------------#
# Variant EFFECT PREDICTOR #
#--------------------------#

By Will McLaren (wm2\@ebi.ac.uk)

Usage:
perl variant_effect_predictor.pl [arguments]

Options
--help                 Display this message and quit

-i | --input_file      Input file - if not specified, reads from STDIN
-o | --output_file     Output file [default: "variant_effect_output.txt"]
-f | --format          Input file format - one of "pileup", "vcf" [default: auto-detect]

-s | --species         Species to use [default: "human"]
-b | --buffer_size     Sets the number of variants sent in each batch [default: 500]
                       Increasing buffer size will retrieve results more quickly
                       but requires more memory.
--check_ref            If specified, checks supplied reference allele against stored
                       entry in Ensembl Core database [default: not used]
--check_existing=[0|1] If set to 1, checks for existing co-located variations in the
                       Ensembl Variation database [default: 1]
--failed=[0|1]         If set to 1, includes variations flagged as failed when checking
                       for co-located variations. Only applies in Ensembl 61 or later
					   [default: 1]
					   
--hgnc                 If specified, HGNC gene identifiers are output alongside the
                       Ensembl Gene identifier [default: not used]

-d | --db_host         Manually define database host [default: "ensembldb.ensembl.org"]
-u | --user            Database username [default: "anonymous"]
-p | --password        Database password [default: not used]
-r | --registry_file   Registry file to use defines DB connections [default: not used]
                       Defining a registry file overrides above connection settings.
					   
-w | --whole_genome    EXPERIMENTAL! Run in whole genome mode [default: not used]
                       Should only be used with data covering a whole genome or
					   chromosome e.g. from resequencing. For better performance,
					   set --buffer_size higher (>10000 if memory allows).
					   Disables --check_existing option and gene column by default.
--chunk_size           Sets the chunk size of internal data structure [default: 50kb]
                       Setting this lower may improve speed for variant-dense
					   datasets. Only applies to whole genome mode.
END

	print $usage;
}


sub start {
	my $id = shift;
	$id ||= '-';
	$times{$id} = [gettimeofday];
}

sub end {
	my $id = shift;
	$id ||= '-';
	warn "Time for $id : ", tv_interval($times{$id}), "\n";
}
