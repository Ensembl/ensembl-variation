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
import_vcf.pl - imports variations from a VCF file into an Ensembl variation DB

by Will McLaren (wm2@ebi.ac.uk)
=cut

use strict;
use Bio::EnsEMBL::Registry;
use Getopt::Long;
use FileHandle;
use Data::Dumper;
use Time::HiRes qw(gettimeofday tv_interval);

# get command-line options
my ($in_file, $species, $registry_file, $help, $host, $user, $password, $source, $population, $flank_size);

my $args = scalar @ARGV;

GetOptions(
	'input_file=s'  => \$in_file,
	'species=s'		=> \$species,
	'registry=s'	=> \$registry_file,
	'db_host=s'		=> \$host,
	'user=s'		=> \$user,
	'password=s'	=> \$password,
	'help'			=> \$help,
	'source=s'      => \$source,
	'population=s'  => \$population,
	'flank=s'       => \$flank_size,
);

# set defaults
$species ||= "human";
$host ||= 'ensembldb.ensembl.org';
$user ||= 'anonymous';
$flank_size ||= 200;

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

# connect to DB
my $vdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation')
    || usage( "Cannot find variation db for $species in $registry_file" );
my $dbVar = $vdba->dbc->db_handle;

# get seq_region_id hash
my %seq_region_ids;
my ($seq_region_id, $chr_name);
my $sth = $dbVar->prepare(qq{SELECT seq_region_id, name FROM seq_region});
$sth->execute;
$sth->bind_columns(\$seq_region_id, \$chr_name);
$seq_region_ids{$chr_name} = $seq_region_id while $sth->fetch;

die("ERROR: seq_region not populated\n") unless scalar keys %seq_region_ids;

# get/set source_id
die("ERROR: no source specified\n") unless defined $source;
our $source_id;

# check existing
$sth = $dbVar->prepare(qq{select source_id from source where name = ?});
$sth->execute($source);
$sth->bind_columns(\$source_id);
$sth->fetch;

if(!defined($source_id)) {
	$sth = $dbVar->prepare(qq{insert into source(name) values(?)});
	$sth->execute($source);
	$source_id = $dbVar->last_insert_id(undef, undef, qw(source source_id));
}



# now do population
die("ERROR: no population specified\n") unless defined $population;
our $pop_id;

# check existing
$sth = $dbVar->prepare(qq{select sample_id from sample where name = ?});
$sth->execute($population);
$sth->bind_columns(\$pop_id);
$sth->fetch;

if(!defined($pop_id)) {
	# insert into sample
	$sth = $dbVar->prepare(qq{insert into sample(name) values(?)});
	$sth->execute($population);
	$pop_id = $dbVar->last_insert_id(undef, undef, qw(sample sample_id));
	
	# insert into population
	$sth = $dbVar->prepare(qq{insert ignore into population(sample_id) values(?)});
	$sth->execute($pop_id);
}




my (%headers, $first_sample_col, @sample_ids);

my $start_time = time();

my $var_counter = 0;
our %times;
&start(5000);

# read the file
while(<$in_file_handle>) {
  chomp;
  
  # header lines
  next if /^##/;
  
  my @split = split /\t/;
  
  # column definition line
  if(/^#/) {
	$headers{$split[$_]} = $_ for(0..$#split);
	$first_sample_col = $headers{FORMAT} + 1;
	
	# populate sample-type tables
	$sth = $dbVar->prepare(qq{INSERT INTO sample(name) VALUES(?)});
	my $sth2 = $dbVar->prepare(qq{INSERT INTO individual_population(population_sample_id, individual_sample_id) VALUES(?,?)});
	my $sth3 = $dbVar->prepare(qq{INSERT INTO individual(sample_id, individual_type_id) VALUES(?,?)});
	my $sth4 = $dbVar->prepare(qq{select sample_id from sample where name = ?});
	
	for my $i($first_sample_col..$#split) {
		
		my $sample_id;
		$sth4->execute($split[$i]);
		$sth4->bind_columns(\$sample_id);
		$sth4->fetch;
		
		if(!$sample_id) {
			$sth->execute($split[$i]);
			$sample_id = $dbVar->last_insert_id(undef, undef, qw(sample sample_id));
			$sth2->execute($pop_id, $sample_id);
			$sth3->execute($sample_id, 3);
		}
		
		push @sample_ids, $sample_id;
	}
  }
  
  # data
  else {
	
	$var_counter++;
	if($var_counter =~ /000$/) {
		warn "COUNTER $var_counter";
		&end(5000);
		&start(5000);
	}
	
	# parse into a hash
	my %data = ();
	$data{$_} = $split[$headers{$_}] for keys %headers;
	
	# skip non-variant lines
	next if $data{ALT} eq '.';
	
	# make a var name if none exists
	if($data{ID} eq '.') {
		$data{ID} = 'test_'.$data{'#CHROM'}.'_'.$data{POS};
	}
	
	# populate variation
	$sth = $dbVar->prepare(qq{INSERT INTO variation(source_id, name) VALUES (?,?);});
	$sth->execute($source_id, $data{ID});
	my $var_id = $dbVar->last_insert_id(undef, undef, qw(variation variation_id));
	
	
	# positions
	my ($start, $end) = ($data{POS}, $data{POS});
	
	# work out if this is an indel or not
	my $is_indel = 0;
	
	if($data{ALT} =~ /D|I/) {
		$is_indel = 1;
		
		if($data{ALT} =~ /\,/) {
			warn "WARNING: can't deal with multiple different indel types in one variation";
			next;
		}
		
		# deletion
		if($data{ALT} =~ /D/) {
			my $num_deleted = $data{ALT};
			$num_deleted =~ s/\D+//g;
			
			$end += $num_deleted - 1;
			
			$data{ALT} = "-";
		}
		
		# insertion
		else {
			$data{REF} = '-';
			
			my $insert = $data{ALT};
			$insert =~ s/^I//g;
			
			$data{ALT} = $insert;
			$start++;
		}
	}
	
	# get alleles
	my @alleles = ($data{REF}, split /\,/, $data{ALT});
	
	
	# populate flanking_sequence
	$sth = $dbVar->prepare(qq{INSERT INTO flanking_sequence(variation_id, up_seq_region_start, up_seq_region_end, down_seq_region_start, down_seq_region_end, seq_region_id, seq_region_strand) VALUES (?,?,?,?,?,?,?);});
	$sth->execute(
		$var_id,
		$start - $flank_size,
		$start - 1,
		$end + 1,
		$end + $flank_size,
		$seq_region_ids{$data{'#CHROM'}},
		1
	);
	
	# populate variation_feature
	$sth = $dbVar->prepare(qq{INSERT INTO variation_feature(variation_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, variation_name, allele_string, map_weight, source_id) VALUES (?,?,?,?,?,?,?,?,?);});
	$sth->execute(
		$var_id,
		$seq_region_ids{$data{'#CHROM'}},
		$data{POS},
		$data{POS},
		1,
		$data{ID},
		(join "/", @alleles),
		1,
		$source_id #needs changing
	);
	
	# parse info column
	my %info;	
	foreach my $chunk(split /\;/, $data{INFO}) {
		my ($key, $val) = split /\=/, $chunk;
		$info{$key} = $val;
	}
	
	# deal with frequencies
	my (@freqs, %gt_freqs, @genotypes);
	
	# may be defined in INFO column
	if(defined($info{AF})) {
		@freqs = split /\,/, $info{AF};
		my $total_alt_freq = 0;
		$total_alt_freq += $_ for @freqs;
		unshift @freqs, 1 - $total_alt_freq;
	}
	
	# otherwise calculate from genotypes
	else {
		my $total_count;
		my %counts;
		
		# get counts
		for my $i($first_sample_col..$#split) {
			my $gt;
			my @bits;
			push @bits, $alleles[$_] for split /\||\/|\\/, (split /\:/, $split[$i])[0];
			$gt = join '|', sort @bits;
			
			# store genotypes for later
			push @genotypes, $gt;
			
			$counts{$gt}++;
			$total_count++;
		}
		
		# now calculate frequencies
		for my $i(0..$#alleles) {
			my ($c_aa, $c_ab);
			my $a = $alleles[$i];
			
			$c_aa = $counts{$a.'|'.$a};
			
			foreach my $gt(keys %counts) {
				$c_ab += $counts{$gt} if $gt =~ /$a/ and $gt ne $a.'|'.$a;
				$gt_freqs{$gt} = $counts{$gt}/$total_count;
			}
			
			push @freqs, (((2*$c_aa) + $c_ab) / (2*$total_count));
		}
	}
	
	# populate allele table
	$sth = $dbVar->prepare(qq{INSERT INTO allele(variation_id, allele, frequency, sample_id) VALUES (?,?,?,?)});
	$sth->execute($var_id, $alleles[$_], $freqs[$_], $pop_id) for (0..$#alleles);
	
	# populate population_genotype table
	if(scalar keys %gt_freqs) {
		$sth = $dbVar->prepare(qq{INSERT INTO population_genotype(variation_id, allele_1, allele_2, frequency, sample_id) VALUES (?,?,?,?,?)});
		
		foreach my $gt(keys %gt_freqs) {
			my @bits = split /\|/, $gt;
			$sth->execute($var_id, $bits[0], $bits[1], $gt_freqs{$gt}, $pop_id);
		}
	}
	
	# now do genotypes
	my @rows;
	
	for my $i($first_sample_col..$#split) {
		my $sample_id = $sample_ids[$i-$first_sample_col];
		
		my @bits;
		
		if(scalar @genotypes) {
			@bits = split /\|/, $genotypes[$i-$first_sample_col];
		}
		
		else {
			push @bits, $alleles[$_] for split /\||\/|\\/, (split /\:/, $split[$i])[0];
		}
		
		push @rows, "(".(join ",", ($var_id, '"'.$bits[0].'"', '"'.$bits[1].'"', $sample_id)).")";
	}
	
	my $vals = join ",", @rows;
	
	my $table = ($is_indel ? 'individual_genotype_multiple_bp' : 'tmp_individual_genotype_single_bp');
	
	$sth = $dbVar->prepare(qq{INSERT INTO $table(variation_id, allele_1, allele_2, sample_id) values$vals});
	$sth->execute;
  }
}

&end(5000);
print "Took ", time() - $start_time, " to run\n";

sub usage {
	my $usage =<<END;
Usage:
perl import_vcf.pl [arguments]

Options
-h | --help           Display this message and quit


-i | --input_file     Input file - if not specified, reads from STDIN
--species             Species to use [default: "human"]
--source              Name of source [required]
--population          Name of population [required]
-f | --flank          Size of flanking sequence [default: 200]

-d | --db_host        Manually define database host [default: "ensembldb.ensembl.org"]
-u | --user           Database username [default: "anonymous"]
--password            Database password [default: not used]
-r | --registry_file  Registry file to use defines DB connections [default: not used]
                      Defining a registry file overrides above connection settings.
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
