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
use ImportUtils qw(debug load);

use constant DISTANCE => 100_000;
use constant MAX_SHORT => 2**16 -1;

my %Printable = ( "\\"=>'\\', "\r"=>'r', "\n"=>'n', "\t"=>'t', "\""=>'"' );

# get command-line options
my ($in_file, $species, $registry_file, $help, $host, $user, $password, $source, $population, $flank_size, $TMP_DIR, $TMP_FILE, $skip_multi, $use_gp, $sample_prefix, $variation_prefix, $disable_keys);

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
	'tmpdir=s'      => \$TMP_DIR,
	'tmpfile=s'     => \$TMP_FILE,
	'skip_multi'    => \$skip_multi,
	'gp'            => \$use_gp,
	'ind_prefix=s'  => \$sample_prefix,
	'var_prefix=s'  => \$variation_prefix,
	'disable_keys'  => \$disable_keys,
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

die "ERROR: tmpdir not specified" unless defined $TMP_DIR;

$TMP_FILE ||= 'compress.txt';


$ImportUtils::TMP_DIR = $TMP_DIR;
$ImportUtils::TMP_FILE = $TMP_FILE;

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
	$in_file_handle = 'STDIN';
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


# disable keys if requested
if($disable_keys) {
	foreach my $table(qw/allele population_genotype individual_genotype_multiple_bp variation_feature/) {
		$sth = $dbVar->do(qq{ALTER TABLE $table DISABLE KEYS;})
	}
}

my (%headers, $first_sample_col, @sample_ids, $region_start, $region_end, $prev_region_end, $prev_seq_region, $genotypes);

my $start_time = time();

my $var_counter = 0;
our %times;
&start(1000);

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
		
		my $sample_name = $sample_prefix.$split[$i];
		my $sample_id;
		$sth4->execute($sample_name);
		$sth4->bind_columns(\$sample_id);
		$sth4->fetch;
		
		if(!$sample_id) {
			$sth->execute($sample_name);
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
		&end(1000);
		&start(1000);
	}
	
	# parse into a hash
	my %data = ();
	$data{$_} = $split[$headers{$_}] for keys %headers;
	
	# skip non-variant lines
	next if $data{ALT} eq '.';
	
	# make a var name if none exists
	if($data{ID} eq '.') {
		$data{ID} =
			($variation_prefix ? $variation_prefix : 'tmp').
			'_'.$data{'#CHROM'}.'_'.$data{POS};
	}
	
	# check if variation exists
	my ($variation_already_exists, $var_id);
	
	$sth = $dbVar->prepare(qq{SELECT variation_id FROM variation WHERE name = ?;});
	$sth->execute($data{ID});
	$sth->bind_columns(\$var_id);
	$sth->fetch;
	
	if(defined($var_id)) {
		$variation_already_exists = 1;
	}
	
	else {
		# populate variation
		$sth = $dbVar->prepare(qq{INSERT INTO variation(source_id, name) VALUES (?,?);});
		$sth->execute($source_id, $data{ID});
		$var_id = $dbVar->last_insert_id(undef, undef, qw(variation variation_id));
	}
	
	
	# positions
	my ($start, $end, $chr, $seq_region);
	
	if($use_gp) {
		foreach my $pair(split /\;/, $data{INFO}) {
			my ($key, $value) = split /\=/, $pair;
			if($key eq 'GP') {
				($chr,$start) = split /\:/, $value;
				$seq_region = $seq_region_ids{$chr};
				$end = $start;
			}
		}
		
		next unless defined($seq_region) and defined($start);
	}
	
	else {
		($start, $end) = ($data{POS}, $data{POS});
		$seq_region = $seq_region_ids{$data{'#CHROM'}};
	}
	
	# work out if this is an indel or multi-bp
	my $is_indel = 0;
	my $is_multi = 0;
	
	$is_multi = 1 if length($data{REF}) > 1;
	foreach my $alt(split /\,/, $data{ALT}) {
		$is_multi = 1 if length($alt) > 1;
	}
	
	# adjust end coord
	$end += (length($data{REF}) - 1);
	
	# find out if any of the alt alleles make this an insertion or a deletion
	foreach my $alt_allele(split /\,/, $data{ALT}) {
		$is_indel = 1 if $alt_allele =~ /D|I/;
		$is_indel = 1 if length($alt_allele) != length($data{REF});
	}
	
	# multiple alt alleles?
	if($data{ALT} =~ /\,/) {
		if($is_indel) {
			
			my @alts;
			
			if($data{ALT} =~ /D|I/) {
				foreach my $alt_allele(split /\,/, $data{ALT}) {
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
				$data{REF} = substr($data{REF}, 1);
				$data{REF} = '-' if $data{REF} eq '';
				$start++;
				
				foreach my $alt_allele(split /\,/, $data{ALT}) {
					$alt_allele = substr($alt_allele, 1);
					$alt_allele = '-' if $alt_allele eq '';
					push @alts, $alt_allele;
				}
			}
			
			$data{ALT} = join "/", @alts;
		}
		
		else {
			# for substitutions we just need to replace ',' with '/' in $alt
			$data{ALT} =~ s/\,/\//;
		}
	}
	
	else {
		if($is_indel) {
			# deletion (VCF <4)
			if($data{ALT} =~ /D/) {
				my $num_deleted = $data{ALT};
				$num_deleted =~ s/\D+//g;
				$end += $num_deleted - 1;
				$data{ALT} = "-";
				$data{REF} .= ("N" x ($num_deleted - 1)) unless length($data{REF}) > 1;
			}
			
			# insertion (VCF <4)
			elsif($data{ALT} =~ /I/) {
				$data{REF} = '-';
				$data{ALT} =~ s/^I//g;
				$start++;
			}
			
			# insertion or deletion (VCF 4+)
			else {
				# chop off first base
				$data{REF} = substr($data{REF}, 1);
				$data{ALT} = substr($data{ALT}, 1);
				
				$start++;
				
				if($data{REF} eq '') {
					# make ref '-' if no ref allele left
					$data{REF} = '-';
					
					# extra adjustment required for Ensembl
					$start++;
				}
				
				# make alt '-' if no alt allele left
				$data{ALT} = '-' if $data{ALT} eq '';
			}
		}
	}
	
	
	# get alleles
	my @alleles = ($data{REF}, split /\,|\//, $data{ALT});
	
	
	unless($variation_already_exists) {
		
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
	}
	
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
			my $gt = (split /\:/, $split[$i])[0];
			foreach my $bit(split /\||\/|\\/, $gt) {
				push @bits, $alleles[$bit] unless $bit eq '.';
			}
			
			if(scalar @bits) {
				$gt = join '|', sort @bits;
				
				# store genotypes for later
				push @genotypes, $gt;
				
				$counts{$gt}++;
				$total_count++;
			}
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
	$sth = $dbVar->prepare(qq{INSERT IGNORE INTO allele(variation_id, allele, frequency, sample_id) VALUES (?,?,?,?)});
	$sth->execute($var_id, $alleles[$_], $freqs[$_], $pop_id) for (0..$#alleles);
	
	# now do genotypes
	my (@multi_rows, $total_count, %counts);
	
	for my $i($first_sample_col..$#split) {
		my $sample_id = $sample_ids[$i-$first_sample_col];
		
		my @bits;
		my $gt = (split /\:/, $split[$i])[0];
		foreach my $bit(split /\||\/|\\/, $gt) {
			push @bits, $alleles[$bit] unless $bit eq '.';
		}
		
		if(scalar @bits && !(scalar keys %gt_freqs)) {
			$counts{$bits[0].$bits[1]}++;
			$total_count++;
		}
		
		# for indels, write directly to temp table
		if($is_indel || $is_multi) {
			push @multi_rows,
				"(".
					(join ",",
						(
							$var_id,
							'"'.($bits[0] || '.').'"',
							'"'.($bits[1] || '.').'"',
							$sample_id
						)
					).
				")";
		}
		
		# otherwise add to compress hash for writing later
		elsif(scalar @bits) {
			
			if (!defined $genotypes->{$sample_id}->{region_start}){
				$genotypes->{$sample_id}->{region_start} = $start;
				$genotypes->{$sample_id}->{region_end} = $end;
			}
			
			# write previous data?
			#compare with the beginning of the region if it is within the DISTANCE of compression
			if (
				(abs($genotypes->{$sample_id}->{region_start} - $start) > DISTANCE()) ||
				(abs($start - $genotypes->{$sample_id}->{region_end}) > MAX_SHORT) ||
				(defined($prev_seq_region) && $seq_region != $prev_seq_region)
			) {
				#snp outside the region, print the region for the sample we have already visited and start a new one
				print_file("$TMP_DIR/compressed_genotype_".$$.".txt",$genotypes, $prev_seq_region, $sample_id);
				delete $genotypes->{$sample_id}; #and remove the printed entry
				$genotypes->{$sample_id}->{region_start} = $start;
			}
			
			# not first genotype
			if ($start != $genotypes->{$sample_id}->{region_start}){
				#compress information
				my $blob = pack ("n",$start - $genotypes->{$sample_id}->{region_end} - 1);
				$genotypes->{$sample_id}->{genotypes} .= &escape($blob).($bits[0] || '-').($bits[1] || '0');
			}
			# first genotype
			else{
				$genotypes->{$sample_id}->{genotypes} = ($bits[0] || '-').($bits[1] || '0');
			}
			
			$genotypes->{$sample_id}->{region_end} = $start;
		}
	}	
	
	# write indel/multibp genotypes
	if(scalar @multi_rows) {
		my $vals = join ",", @multi_rows;	
		$sth = $dbVar->prepare(qq{INSERT IGNORE INTO individual_genotype_multiple_bp(variation_id, allele_1, allele_2, sample_id) values$vals});
		$sth->execute;
	}
	
	if(!(scalar keys %gt_freqs)) {
		foreach my $gt(keys %counts) {
			$gt_freqs{$gt} = $counts{$gt}/$total_count;
		}
	}
	
	# populate population_genotype table
	if(scalar keys %gt_freqs) {
		$sth = $dbVar->prepare(qq{INSERT IGNORE INTO population_genotype(variation_id, subsnp_id, allele_1, allele_2, frequency, sample_id) VALUES (?,NULL,?,?,?,?)});
		
		foreach my $gt(keys %gt_freqs) {
			my @bits = split /\||\/|\\/, $gt;
			$sth->execute($var_id, $bits[0], $bits[1], $gt_freqs{$gt}, $pop_id);
		}
	}
	
	$prev_seq_region = $seq_region;
  }
}

print_file("$TMP_DIR/compressed_genotype_".$$.".txt",$genotypes, $prev_seq_region);
&import_genotypes($dbVar);

$dbVar->do(qq{DELETE FROM individual_genotype_multiple_bp WHERE allele_1 = '.' AND allele_2 = '.';});

# re-enable keys if requested
if($disable_keys) {
	foreach my $table(qw/allele population_genotype individual_genotype_multiple_bp variation_feature/) {
		$sth = $dbVar->do(qq{ALTER TABLE $table ENABLE KEYS;})
	}
}

&end(1000);
print "Took ", time() - $start_time, " to run\n";

sub usage {
	my $usage =<<END;
Usage:
perl import_vcf.pl [arguments]

Options
-h | --help           Display this message and quit


-i | --input_file     Input file - if not specified, reads from STDIN
--tmpdir              Temporary directory to write genotype dumpfile. Required.
--tmpfile             Name for temporary file [default: compress.txt]

--species             Species to use [default: "human"]
--source              Name of source [required]
--population          Name of population [required]
--ind_prefix          Prefix added to sample names [default: not used]
--var_prefix          Prefix added to constructed variation names [default: not used]
-f | --flank          Size of flanking sequence [default: 200]
--gp                  Use GP tag from INFO column to get locations [default: not used]

--disable_keys        Disable MySQL keys during inserts [default: not used]

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


sub print_file{
    my $file = shift;
    my $genotypes = shift;
    my $seq_region_id = shift;
    my $sample_id = shift;

    open( FH, ">>$file") or die "Could not add compressed information: $!\n";
    if (!defined $sample_id){
		#new chromosome, print all the genotypes and flush the hash
		foreach my $sample_id (keys %{$genotypes}){
			print FH join("\t",
				$sample_id,
				$seq_region_id,
				$genotypes->{$sample_id}->{region_start},
				$genotypes->{$sample_id}->{region_end},
				1,
				$genotypes->{$sample_id}->{genotypes}) . "\n";
		}
    }
    else{
		#only print the region corresponding to sample_id
		print FH join("\t",
			$sample_id,
			$seq_region_id,
			$genotypes->{$sample_id}->{region_start},
			$genotypes->{$sample_id}->{region_end},
			1,
			$genotypes->{$sample_id}->{genotypes}) . "\n";
    }
    close FH;
}

# $special_characters_escaped = printable( $source_string );
sub escape ($) {
  local $_ = ( defined $_[0] ? $_[0] : '' );
  s/([\r\n\t\\\"])/\\$Printable{$1}/sg;
  return $_;
}

sub import_genotypes{
    my $dbVar = shift;

    warn("Importing compressed genotype data");
    my $call = "mv $TMP_DIR/compressed_genotype_".$$.".txt $TMP_DIR/$TMP_FILE";
    system($call);
    load($dbVar,qw(compressed_genotype_single_bp sample_id seq_region_id seq_region_start seq_region_end seq_region_strand genotypes));
}