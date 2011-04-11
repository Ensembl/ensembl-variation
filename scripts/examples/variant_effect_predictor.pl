#!/usr/bin/perl

=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

Version 2.0

by Will McLaren (wm2@ebi.ac.uk)
=cut

use strict;
use Getopt::Long;
use FileHandle;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);
	
# configure from command line opts
my $config = &configure(scalar @ARGV);

# run the main sub routine
&main($config);

# this is the main sub-routine - it needs the configured $config hash
sub main {
	my $config = shift;
	
	debug("Starting...") if defined $config->{verbose};
	
	my $species = $config->{species};
	
	# get adaptors
	my $vfa = $config->{reg}->get_adaptor($species, 'variation', 'variationfeature');
	$config->{tva} = $config->{reg}->get_adaptor($species, 'variation', 'transcriptvariation');
	
	# get fake ones for species with no var DB
	if(!defined($vfa)) {
		$vfa = Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor->new_fake($species);
	}
	else {
		$vfa->db->include_failed_variations($config->{include_failed}) if defined($vfa->db) && $vfa->db->can('include_failed_variations');
	}
	
	
	$config->{sa} = $config->{reg}->get_adaptor($species, 'core', 'slice');
	$config->{ga} = $config->{reg}->get_adaptor($species, 'core', 'gene');
	
	# check we got slice adaptor - can't continue without a core DB
	die("ERROR: Could not connect to core database\n") unless defined $config->{sa} and defined $config->{ga};
	
	
	# create a hash to hold slices so we don't get the same one twice
	my %slice_hash = ();
	my @new_vfs;
	my %vf_hash;
	
	my $transcript_cache;
	
	my $line_number = 0;
	my $vf_count;
	my $in_file_handle = $config->{in_file_handle};
	
	# read the file
	while(<$in_file_handle>) {
	  chomp;
	  
	  $line_number++;
	  
	  # header line?
	  next if /^\#/;
	  
	  # some lines (pileup) may actually parse out into more than one variant)
	  foreach my $sub_line(@{&parse_line($config, $_)}) {
		
		# get the sub-line into named variables
		my ($chr, $start, $end, $allele_string, $strand, $var_name) = @{$sub_line};
		
		# non-variant line from VCF
		next if $chr eq 'non-variant';
		
		# fix inputs
		$chr =~ s/chr//ig;
		$strand = ($strand =~ /\-/ ? "-1" : "1");
		$allele_string =~ tr/acgt/ACGT/;
		
		# sanity checks
		unless($start =~ /^\d+$/ && $end =~ /^\d+$/) {
		  warn("WARNING: Start $start or end $end coordinate invalid on line $line_number\n") unless defined $config->{quiet};
		  next;
		}
		
		unless($allele_string =~ /([ACGT-]+\/*)+/) {
		  warn("WARNING: Invalid allele string $allele_string on line $line_number\n") unless defined $config->{quiet};
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
		  eval { $slice = $config->{sa}->fetch_by_region('chromosome', $chr); };
		  
		  # if failed, try to get any seq region
		  if(!defined($slice)) {
			  $slice = $config->{sa}->fetch_by_region(undef, $chr);
		  }
		  
		  # if failed, warn and skip this line
		  if(!defined($slice)) {
			  warn("WARNING: Could not fetch slice named $chr on line $line_number\n") unless defined $config->{quiet};
			  next;
		  }	
		  
		  # store the hash
		  $slice_hash{$chr} = $slice;
		}
		
		# check reference allele if requested
		if(defined $config->{check_ref}) {
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
				  warn "WARNING: Could not fetch sub-slice from $start\-$end\($strand\) on line $line_number" unless defined $config->{quiet};
			  }
			  
			  else {
				  $slice_ref_allele = $slice_ref->seq;
				  $ok = ($slice_ref_allele eq $ref_allele ? 1 : 0);
			  }
		  }
		  
		  if(!$ok) {
			  warn
				  "WARNING: Specified reference allele $ref_allele ",
				  "does not match Ensembl reference allele",
				  ($slice_ref_allele ? " $slice_ref_allele" : ""),
				  " on line $line_number" unless defined $config->{quiet};
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
		
		if(defined $config->{whole_genome}) {
			push @{$vf_hash{$chr}{int($start / $config->{chunk_size})}{$start}}, $new_vf;
			$vf_count++;
			
			if($vf_count == $config->{buffer_size}) {
				&print_consequences($config, &whole_genome_fetch($config, \%vf_hash, $transcript_cache));
				%vf_hash = ();
				$vf_count = 0;
			}
		}
		else {
			&print_consequences($config, [$new_vf]);
			$vf_count++;
			debug("Processed $vf_count variants") if $vf_count =~ /0$/ && defined($config->{verbose});
		}
	  }
	}
	
	# if in whole-genome mode, finish off the rest of the buffer
	if(defined $config->{whole_genome} && %vf_hash) {
		&print_consequences($config, &whole_genome_fetch($config, \%vf_hash, $transcript_cache));
	}
	
	debug("Finished!") if defined $config->{verbose};
}

sub print_consequences {
	my $config = shift;
	my $vfs = shift;
	
	my $out_file_handle = $config->{out_file_handle};
	
	foreach my $new_vf(@$vfs) {
		
		# find any co-located existing VFs
		my $existing_vf;
		$existing_vf = &find_existing($new_vf) if defined $config->{check_existing} && $config->{check_existing} == 1;
		$existing_vf ||= '-';
		
		foreach my $tv(@{$new_vf->get_all_TranscriptVariations}) {
			
			next if(defined $config->{coding_only} && !($tv->affects_transcript));
			
			foreach my $tva(@{$tv->get_all_alternate_TranscriptVariationAlleles}) {
				
				my $method_name = $config->{terms}.'_term';
				my $type = join ",", map {$_->$method_name || $_->display_term} @{$tva->get_all_OverlapConsequences};
				
				my $gene = ($tv->transcript ? $config->{ga}->fetch_by_transcript_stable_id($tv->transcript->stable_id) : undef) unless defined $config->{whole_genome};
				
				# extra
				my $extra;
				
				# HGNC
				if(defined $config->{hgnc} && $gene) {
				  my @entries = grep {$_->database eq 'HGNC'} @{$gene->get_all_DBEntries()};
				  if(scalar @entries) {
					$extra .= 'HGNC='.$entries[0]->display_id.';';
				  }
				}
				
				# protein ID
				if(defined $config->{protein} && $tv->transcript->translation) {
					$extra .= 'ENSP='.$tv->transcript->translation->stable_id.';';
				}
				
				# HGVS
				if(defined $config->{hgvs}) {
					$extra .= 'HGVSc='.$tva->hgvs_coding.';' if defined($tva->hgvs_coding);
					$extra .= 'HGVSp='.$tva->hgvs_protein.';' if defined($tva->hgvs_protein);
				}
				
				foreach my $tool (qw(SIFT PolyPhen Condel)) {
					my $lc_tool = lc($tool);
					
					if (my $opt = $config->{$lc_tool}) {
						my $want_pred   = $opt =~ /^p/i;
						my $want_score  = $opt =~ /^s/i;
						my $want_both   = $opt =~ /^b/i;
						
						if ($want_both) {
							$want_pred  = 1;
							$want_score = 1;
						}
						
						next unless $want_pred || $want_score;
						
						my $pred_meth   = $lc_tool.'_prediction';
						my $score_meth  = $lc_tool.'_score';
						
						my $pred = $tva->$pred_meth;
						
						if($pred) {
							$extra .= "$tool=";
							
							if ($want_pred) {
								$pred =~ s/\s+/\_/;
								$extra .= $pred;
							}
								
							if ($want_score) {
								my $score = $tva->$score_meth;
								
								if(defined $score) {
									if($want_pred) {
										$extra .= "($score)";
									}
									else {
										$extra .= $score;
									}
								}
							}
							
							$extra .= ';';
						}
					}
				}
				
				$extra =~ s/\;$//g;
				
				print $out_file_handle
				  $new_vf->variation_name, "\t",
				  $new_vf->seq_region_name, ":",
				  &format_coords($new_vf->start, $new_vf->end), "\t",
				  $tva->variation_feature_seq, "\t",
				  ($gene ? $gene->stable_id : '-'), "\t",
				  ($tv->transcript ? $tv->transcript->stable_id : "-"), "\t",
				  $type, "\t",
				  &format_coords($tv->cdna_start, $tv->cdna_end), "\t",
				  &format_coords($tv->cds_start, $tv->cds_end), "\t",
				  &format_coords($tv->translation_start, $tv->translation_end), "\t",
				  ($tva->pep_allele_string || "-"), "\t",
				  ($tva->display_codon_allele_string || "-"), "\t",
				  $existing_vf, "\t",
				  ($extra || '-'), "\n";
			}
		}
	}
}

sub configure {
	my $args = shift;
	
	my $config = {};
	
	GetOptions(
		$config,
		'help',
		
		# input options,
		'config=s',
		'input_file=s',
		'format=s',
		
		# DB options
		'species=s',
		'registry=s',
		'host=s',
		'user=s',
		'port=s',
		'password=s',
		'db_version=i',
		'genomes',
		
		# runtime options
		'most_severe',
		'buffer_size=i',
		'chunk_size=s',
		'check_ref',
		'check_existing=i',
		'failed=i',
		'whole_genome',
		'tmp_dir=s',
		'gp',
		
		# output options
		'output_file=s',
		'terms=s',
		'verbose',
		'quiet',
		'coding_only',
		'protein',
		'hgnc',
		'hgvs',
		'sift=s',
		'polyphen=s',
		'condel=s',
	);
	
	# print usage message if requested or no args supplied
	if(defined($config->{help}) || !$args) {
		&usage;
		exit(0);
	}
	
	# config file?
	if(defined $config->{config}) {
		
		open CONFIG, $config->{config} or die "ERROR: Could not open config file \"".$config->{config}."\"\n";
		
		while(<CONFIG>) {
			next if /^\#/;
			my ($key, $value) = split /\s+|\=/;
			$key =~ s/^\-//g;
			$config->{$key} = $value unless defined $config->{$key};
		}
		
		close CONFIG;
	}
	
	# check file format
	if(defined $config->{input_format}) {
		die "ERROR: Unrecognised input format specified \"".$config->{input_format}."\"\n" unless $config->{input_format} =~ /pileup|vcf|guess/i;
	}
	
	# output term
	if(defined $config->{terms}) {
		die "ERROR: Unrecognised consequence term type specified \"".$config->{terms}."\" - must be one of ensembl, so, ncbi\n" unless $config->{terms} =~ /ensembl|display|so|ncbi/i;
		if($config->{terms} =~ /ensembl|display/i) {
			$config->{terms} = 'display';
		}
		else {
			$config->{terms} = uc($config->{terms});
		}
	}
	
	# check nsSNP tools
	foreach my $tool(grep {defined $config->{lc($_)}} qw(SIFT PolyPhen Condel)) {
		die "ERROR: Unrecognised option for $tool \"", $config->{lc($tool)}, "\" - must be one of p (prediction), s (score) or b (both)\n" unless $config->{lc($tool)} =~ /^(s|p|b)/;
	}
	
	# summarise options if verbose
	if(defined $config->{verbose}) {
		my $header =<<INTRO;
#----------------------------------#
# ENSEMBL VARIANT EFFECT PREDICTOR #
#----------------------------------#

version 2.0

By Will McLaren (wm2\@ebi.ac.uk)

Configuration options:

INTRO
		print $header;
		
		my $max_length = (sort {$a <=> $b} map {length($_)} keys %$config)[-1];
		
		foreach my $key(sort keys %$config) {
			print $key.(' ' x (($max_length - length($key)) + 4)).$config->{$key}."\n";
		}
		
		print "\n".("-" x 20)."\n\n";
	}
	
	# connection settings for Ensembl Genomes
	if($config->{genomes}) {
		$config->{host} ||= 'mysql.ebi.ac.uk';
		$config->{port} ||= 4157;
	}
	
	# connection settings for main Ensembl
	else {
		$config->{species} ||= "homo_sapiens";
		$config->{host}    ||= 'ensembldb.ensembl.org';
		$config->{port}    ||= 5306;
	}
	
	# set defaults
	$config->{user}         ||= 'anonymous';
	$config->{buffer_size}  ||= 5000;
	$config->{chunk_size}   ||= '50kb';
	$config->{output_file}  ||= "variant_effect_output.txt";
	$config->{tmpdir}       ||= '/tmp';
	$config->{format}       ||= 'guess';
	$config->{terms}        ||= 'display';
	
	$config->{include_failed} = 1 unless defined $config->{include_failed};
	$config->{check_existing} = 1 unless (defined $config->{check_existing} || defined $config->{whole_genome});
	$config->{chunk_size} =~ s/mb?/000000/i;
	$config->{chunk_size} =~ s/kb?/000/i;
	
	# connect to databases
	$config->{reg} = &connect_to_dbs($config);
	
	# get input file handle
	$config->{in_file_handle} = &get_in_file_handle($config);
	
	# configure output file
	$config->{out_file_handle} = &get_out_file_handle($config);
	
	return $config;
}

sub usage {
	my $usage =<<END;
#----------------------------------#
# ENSEMBL VARIANT EFFECT PREDICTOR #
#----------------------------------#

version 2.0

By Will McLaren (wm2\@ebi.ac.uk)

Usage:
perl variant_effect_predictor.pl [arguments]

Options
--help                 Display this message and quit
--verbose              Display verbose output as the script runs [default: off]
--quiet                Suppress status and warning messages [default: off]

--config               Load configuration from file. Any command line options
                       specified overwrite those in the file [default: off]

-i | --input_file      Input file - if not specified, reads from STDIN. Files
                       may be gzip compressed.
--format               Alternative input file format - one of "pileup", "vcf"
-o | --output_file     Output file [default: "variant_effect_output.txt"]

-t | --terms           Type of consequence terms to output - one of "ensembl", "SO",
                       "NCBI" [default: ensembl]
					   
--sift=[p|s|b]         Add SIFT [p]rediction, [s]core or [b]oth [default: off]
--polyphen=[p|s|b]     Add PolyPhen [p]rediction, [s]core or [b]oth [default: off]
--condel=[p|s|b]       Add Condel SIFT/PolyPhen consensus [p]rediction, [s]core or
                       [b]oth [default: off]

NB: SIFT, PolyPhen and Condel predictions are currently available for human only

--hgnc                 If specified, HGNC gene identifiers are output alongside the
                       Ensembl Gene identifier [default: off]
--hgvs                 Output HGVS identifiers (coding and protein) [default: off]
--protein              Output Ensembl protein identifer [default: off]

--coding_only          Only return consequences that fall in the coding region of
                       transcripts [default: off]

--check_ref            If specified, checks supplied reference allele against stored
                       entry in Ensembl Core database [default: off]
--check_existing=[0|1] If set to 1, checks for existing co-located variations in the
                       Ensembl Variation database [default: 1]
--failed=[0|1]         If set to 1, includes variations flagged as failed when checking
                       for co-located variations. Only applies in Ensembl 61 or later
                       [default: 1]
--gp                   If specified, tries to read GRCh37 position from GP field in the
                       INFO column of a VCF file. Only applies when VCF is the input
                       format and human is the species [default: off]

-s | --species         Species to use [default: "human"]
--host                 Manually define database host [default: "ensembldb.ensembl.org"]
-u | --user            Database username [default: "anonymous"]
--port                 Database port [default: 5306]
--password             Database password [default: no password]
--genomes              Sets DB connection params for Ensembl Genomes [default: off]
-r | --registry_file   Registry file to use defines DB connections [default: off]
                       Defining a registry file overrides above connection settings.
--db_version=[number]  Force script to load DBs from a specific Ensembl version. Not
                       advised due to likely incompatibilities between API and DB
					   
-w | --whole_genome    Run in whole genome mode [default: off]
                       Recommended for use with data covering a whole genome,
                       chromosome or gene/set of genes e.g. from resequencing. For
                       better performance, set --buffer_size higher (>10000 if memory
                       allows).
-b | --buffer_size     Sets the number of variants sent in each batch [default: 5000]
                       Increasing buffer size can retrieve results more quickly
                       but requires more memory. Only applies to whole genome mode.
--chunk_size           Sets the chunk size of internal data structure [default: 50kb]
                       Setting this lower may improve speed for variant-dense
                       datasets. Only applies to whole genome mode.
					   
NB: Whole genome mode disables --check_existing option and gene column by default.
END

	print $usage;
}


sub connect_to_dbs {
	my $config = shift;
	
	# get registry
	my $reg = 'Bio::EnsEMBL::Registry';
	
	# load DB options from registry file if given
	if(defined($config->{registry})) {
		debug("Loading DB config from registry file ", $config->{registry});
		$reg->load_all($config->{registry});
	}
	
	# otherwise manually connect to DB server
	else {
		$reg->load_registry_from_db(
			-host       => $config->{host},
			-user       => $config->{user},
			-pass       => $config->{password},
			-port       => $config->{port},
			-db_version => $config->{db_version},
			-species    => $config->{species} =~ /^[a-z]+\_[a-z]+/i ? $config->{species} : undef,
			-verbose    => $config->{verbose},
		);
	}
	
	$reg->set_disconnect_when_inactive();
	
	if($config->{verbose}) {
		# get a meta container adaptors to check version
		my $core_mca = $reg->get_adaptor($config->{species}, 'core', 'metacontainer');
		my $var_mca = $reg->get_adaptor($config->{species}, 'variation', 'metacontainer');
		
		if($core_mca && $var_mca) {
			debug(
				"Connected to core version ", $core_mca->get_schema_version, " database ",
				"and variation version ", $var_mca->get_schema_version, " database"
			);
		}
	}
	
	return $reg;
}


sub get_in_file_handle {
	my $config = shift;

	# define the filehandle to read input from
	my $in_file_handle = new FileHandle;
	
	if(defined($config->{input_file})) {
		
		# check defined input file exists
		die("ERROR: Could not find input file ", $config->{input_file}, "\n") unless -e $config->{input_file};
		
		if($config->{input_file} =~ /\.gz$/){
			$in_file_handle->open("zcat ". $config->{input_file} . " | " ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
		}
		else {
			$in_file_handle->open( $config->{input_file} ) or die("ERROR: Could not read from input file ", $config->{in_file}, "\n");
		}
	}
	
	# no file specified - try to read data off command line
	else {
		$in_file_handle = 'STDIN';
		debug("Reading input from STDIN (or maybe you forgot to specify an input file?)...") unless defined $config->{quiet};
	}
	
	return $in_file_handle;
}


sub get_out_file_handle {
	my $config = shift;
	
	# define filehandle to write to
	my $out_file_handle = new FileHandle;
	$out_file_handle->open(">".$config->{output_file}) or die("ERROR: Could not write to output file ", $config->{output_file}, "\n");
	
	# make header
	my $time = &get_time;
	my $core_mca = $config->{reg}->get_adaptor($config->{species}, 'core', 'metacontainer');
	my $db_string = $core_mca->dbc->dbname." on ".$core_mca->dbc->host if defined $core_mca;
	my $version_string =
		"Using API version ".$config->{reg}->software_version.
		", DB version ".(defined $core_mca && $core_mca->get_schema_version ? $core_mca->get_schema_version : '?');
	
	my $header =<<HEAD;
## ENSEMBL VARIANT EFFECT PREDICTOR v2.0
## Output produced at $time
## Connected to $db_string
## $version_string
## Extra column keys:
## HGNC     : HGNC gene identifier
## ENSP     : Ensembl protein identifer
## HGVSc    : HGVS coding sequence name
## HGVSp    : HGVS protein sequence name
## SIFT     : SIFT prediction
## PolyPhen : PolyPhen prediction
## Condel   : Condel SIFT/PolyPhen consensus prediction
HEAD
	
	# add headers
	print $out_file_handle $header;
	
	# add column headers
	print $out_file_handle join "\t", qw(
		#Uploaded_variation
		Location
		Allele
		Gene
		Transcript
		Consequence
		cDNA_position
		CDS_position
		Protein_position
		Amino_acids
		Codons
		Existing_variation
		Extra
	);
	
	print $out_file_handle "\n";
	
	return $out_file_handle;
}



sub parse_line {
	my $config = shift;
	my $line   = shift;
	
	my @data = (split /\s+/, $_);
	
	# pileup: chr1 60 T A
	if(
	   ($config->{input_format} =~ /pileup/i) ||
	   (
			$data[0] =~ /(chr)?\w+/ &&
			$data[1] =~ /\d+/ &&
			$data[2] =~ /^[ACGTN-]+$/ &&
			$data[3] =~ /^[ACGTNRYSWKM*+\/-]+$/
		)
	) {
		my @return = ();
		
		if($data[2] ne "*"){
			my $var;
			
			if($data[2] =~ /^[A|C|G|T]$/) {
				$var = $data[2];
			}
			else {
				($var = unambiguity_code($data[3])) =~ s/$data[2]//ig;
			}
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
					warn("WARNING: invalid pileup indel genotype: $line\n") unless defined $config->{quiet};
					push @return, ['non-variant'];
				}
			}
		}
		return \@return;
	}
	
	# VCF: 20      14370   rs6054257 G     A      29    0       NS=58;DP=258;AF=0.786;DB;H2          GT:GQ:DP:HQ
	elsif(
		($config->{input_format} =~ /vcf/i) ||
		(
			$data[0] =~ /(chr)?\w+/ &&
			$data[1] =~ /\d+/ &&
			$data[3] =~ /^[ACGTN-]+$/ &&
			$data[4] =~ /^([\.ACGTN-]+\,?)+$/
		)
	) {
		
		# non-variant line in VCF, return dummy line
		if($data[4] eq '.') {
			return [['non-variant']];
		}
		
		# get relevant data
		my ($chr, $start, $end, $ref, $alt) = ($data[0], $data[1], $data[1], $data[3], $data[4]);
		
		if(defined $config->{use_gp}) {
			$chr = undef;
			$start = undef;
			
			foreach my $pair(split /\;/, $data[7]) {
				my ($key, $value) = split /\=/, $pair;
				if($key eq 'GP') {
					($chr,$start) = split /\:/, $value;
					$end = $start;
				}
			}
			
			unless(defined($chr) and defined($start)) {
				warn "No GP flag found in INFO column" unless defined $config->{quiet};
				return [['non-variant']];
			}
		}
		
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
				$alt =~ s/\,/\//;
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
					}
					
					# make alt '-' if no alt allele left
					$alt = '-' if $alt eq '';
				}
			}
		}
		
		return [[$chr, $start, $end, $ref."/".$alt, 1, ($data[2] eq '.' ? undef : $data[2])]];
		
	}
	
	# our format
	else {
		# we allow commas as delimiter so re-split
		@data = (split /\s+|\,/, $_);
		return [\@data];
	}
}


sub whole_genome_fetch {
	my $config = shift;
	my $vf_hash = shift;
	my $transcript_cache = shift;
	
	my $up_down_size = MAX_DISTANCE_FROM_TRANSCRIPT;
	
	my %vf_done;
	my @return;
	
	foreach my $chr(keys %$vf_hash) {
		my $slice;
		
 		# first try to get a chromosome
 		eval { $slice = $config->{sa}->fetch_by_region('chromosome', $chr); };
 		
 		# if failed, try to get any seq region
 		if(!defined($slice)) {
 			$slice = $config->{sa}->fetch_by_region(undef, $chr);
 		}
		
		debug("Analyzing chromosome $chr") unless defined($config->{quiet});
		
		$transcript_cache->{$chr} = $slice->get_all_Transcripts() unless defined($transcript_cache->{$chr});
		
		my $tr_count = scalar @{$transcript_cache->{$chr}};
		
		debug("Fetched $tr_count transcripts") unless defined($config->{quiet});
		
		my $tr_counter;
		
		while($tr_counter < $tr_count) {
			
			my $tr = $transcript_cache->{$chr}->[$tr_counter++];
			
			if($tr_counter =~ /(5|0)00$/) {
				debug("Analysed $tr_counter\/$tr_count transcripts") unless defined($config->{quiet});
			}
			
			# do each overlapping VF
			my $chr = $tr->seq_region_name;
			my $s = $tr->start - $up_down_size;
			my $e = $tr->end + $up_down_size;
			
			# get the chunks this transcript overlaps
			my %chunks;
			$chunks{$_} = 1 for (int($s/$config->{chunk_size})..int($e/$config->{chunk_size}));
			map {delete $chunks{$_} unless defined($vf_hash->{$chr}{$_})} keys %chunks;
			
			foreach my $chunk(keys %chunks) {				
				foreach my $pos(grep {$_ >= $s && $_ <= $e} keys %{$vf_hash->{$chr}{$chunk}}) {
					foreach my $vf(@{$vf_hash->{$chr}{$chunk}{$pos}}) {
						$vf->add_TranscriptVariation(Bio::EnsEMBL::Variation::TranscriptVariation->new(
							-transcript => $tr,
							-variation_feature => $vf,
							-adaptor => $config->{tva},
						));
						
						if(!defined($vf_done{$vf->{'variation_name'}.'_'.$vf->{'start'}})) {
							push @return, $vf;
							$vf_done{$vf->{'variation_name'}.'_'.$vf->{'start'}} = 1;
						}
					}
				}
			}
		}
		
		# clean hash
		delete $vf_hash->{$chr};
	}
	
	debug("Calculating and writing output") unless defined($config->{quiet});
	
	return \@return;
}

# gets time
sub get_time() {
	my @time = localtime(time());

	# increment the month (Jan = 0)
	$time[4]++;

	# add leading zeroes as required
	for my $i(0..4) {
		$time[$i] = "0".$time[$i] if $time[$i] < 10;
	}

	# put the components together in a string
	my $time =
 		($time[5] + 1900)."-".
 		$time[4]."-".
 		$time[3]." ".
		$time[2].":".
		$time[1].":".
		$time[0];

	return $time;
}



# prints debug output with time
sub debug {
	my $text = (@_ ? (join "", @_) : "No message");
	my $time = get_time;
	
	print $time." - ".$text.($text =~ /\n$/ ? "" : "\n");
}

sub format_coords {
	my ($start, $end) = @_;
	
	if(!defined($start)) {
		return '-';
	}
	elsif(!defined($end)) {
		return $start;
	}
	elsif($start == $end) {
		return $start;
	}
	elsif($start > $end) {
		return $end.'-'.$start;
	}
	else {
		return $start.'-'.$end;
	}
}

sub find_existing {
	my $new_vf = shift;
	
	if(defined($new_vf->adaptor->db)) {
		
		my $sth = $new_vf->adaptor->db->dbc->prepare(qq{
		  SELECT variation_name, source_id
		  FROM variation_feature
		  WHERE seq_region_id = ?
		  AND seq_region_start = ?
		  AND seq_region_end = ?
		});
		
		$sth->execute($new_vf->slice->get_seq_region_id, $new_vf->start, $new_vf->end);
		
		my ($name, $source);
		$sth->bind_columns(\$name, \$source);
		
		my %by_source;
		
		push @{$by_source{$source}}, $name while $sth->fetch;
		$sth->finish;
		
		if(scalar keys %by_source) {
			foreach my $s(sort {$a <=> $b} keys %by_source) {
				return shift @{$by_source{$s}};
			}
		}
	}
	
	return undef;
}
