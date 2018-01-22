#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

=head1 NAME
import_vcf_compress.pl - imports genotypes directly into
compressed_genotype_single_bp from a VCF file

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
my ($in_file, $species, $registry_file, $help, $host, $user, $password, $source, $population, $flank_size, $TMP_DIR, $TMP_FILE, $skip_multi, $use_gp, $sample_prefix);

my $args = scalar @ARGV;

GetOptions(
	'input_file=s'  => \$in_file,
	'species=s'		=> \$species,
	'registry=s'	=> \$registry_file,
	'db_host=s'		=> \$host,
	'user=s'		=> \$user,
	'password=s'	=> \$password,
	'help'			=> \$help,
	'population=s'  => \$population,
	'flank=s'       => \$flank_size,
	'tmpdir=s'      => \$TMP_DIR,
	'tmpfile=s'     => \$TMP_FILE,
	'skip_multi'    => \$skip_multi,
	'gp'            => \$use_gp,
	'prefix=s'      => \$sample_prefix,
);

# set defaults
$species ||= "human";
$flank_size ||= 200;
$sample_prefix ||= "";

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
	die "ERROR: Host or user not defined" unless defined $host and defined $user;
	
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


# create tmp table for indel/multibp genotypes
$sth = $dbVar->do(qq{
	CREATE TABLE IF NOT EXISTS `individual_genotype_multiple_bp_named` (
		`name` varchar(255) NOT NULL,
		`allele_1` varchar(255) DEFAULT NULL,
		`allele_2` varchar(255) DEFAULT NULL,
		`sample_id` int(10) unsigned DEFAULT NULL,
		KEY `variation_idx` (`name`),
		KEY `sample_idx` (`sample_id`)
	);
});



my (%headers, $first_sample_col, @sample_ids, $region_start, $region_end, $prev_region_end, $prev_seq_region, $genotypes);

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
	if($var_counter =~ /0000$/) {
		warn "COUNTER $var_counter";
		&end(10000);
		&start(10000);
	}
	
	# parse into a hash
	my %data = ();
	$data{$_} = $split[$headers{$_}] for keys %headers;
	
	# skip non-variant lines
	#next if $data{ALT} eq '.';	
	
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
	foreach my $alt(split /\,/, $data{REF}) {
		$is_multi = 1 if length($alt) > 1;
	}
	
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
	
	next if $skip_multi && ($is_indel || $is_multi);
	
	# get alleles
	my @alleles = ($data{REF}, split /\,/, $data{ALT});
	
	# get genotypes
	my @genotypes;
	
	# now do genotypes
	my @rows;
	
	for my $i($first_sample_col..$#split) {
		my $sample_id = $sample_ids[$i-$first_sample_col];
		
		my @bits;
		my $gt = (split /\:/, $split[$i])[0];
		foreach my $bit(split /\||\/|\\/, $gt) {
			push @bits, $alleles[$bit] unless $bit eq '.';
		}
		
		# for indels, write directly to temp table
		if($is_indel || $is_multi) {
			push @rows,
				"(".
					(join ",",
						(
							'"'.$data{ID}.'"',
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
				print_file("$TMP_DIR/compressed_genotype.txt",$genotypes, $prev_seq_region, $sample_id);
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
	if(scalar @rows) {
		my $vals = join ",", @rows;	
		$sth = $dbVar->prepare(qq{INSERT INTO individual_genotype_multiple_bp_named(name, allele_1, allele_2, sample_id) values$vals});
		$sth->execute;
	}
	
	$prev_seq_region = $seq_region;
  }
}

print_file("$TMP_DIR/compressed_genotype.txt",$genotypes, $prev_seq_region);
&import_genotypes($dbVar);

$dbVar->do(qq{DELETE FROM individual_genotype_multiple_bp_named WHERE allele_1 = '.' AND allele_2 = '.';});

&end(10000);
print "Took ", time() - $start_time, " to run\n";

sub usage {
	my $usage =<<END;
Usage:
perl import_vcf_compress.pl [arguments]

Options
-h | --help           Display this message and quit


-i | --input_file     Input file - if not specified, reads from STDIN
--species             Species to use [default: "human"]
--population          Name of population [required]
--skip_multi          Do not import multi bp or indel genotypes
--gp                  Use the GP info tag to extract location information
--prefix              String to add as a prefix to sample names

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
    my $call = "mv $TMP_DIR/compressed_genotype.txt $TMP_DIR/$TMP_FILE";
    system($call);
    load($dbVar,qw(compressed_genotype_single_bp sample_id seq_region_id seq_region_start seq_region_end seq_region_strand genotypes));
}
