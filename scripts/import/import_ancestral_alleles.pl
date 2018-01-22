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

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;

use Bio::DB::Fasta;
use DBI;
use FileHandle;
use Getopt::Long;
use ImportUtils qw(load);

usage() if (!scalar(@ARGV));

my $config = {};

GetOptions(
  $config,
  'tmp_dir=s',
  'fasta_files_dir=s',
  'mode=s',
  'resume_update=s',
  'version=i',
  'registry=s',
  'species=s',
  'host=s',
  'dbname=s',
  'user=s',
  'pass=s',
  'port=i',
  'help!',
) or die "Error: Failed to parse command line arguments\n";

usage() if ($config->{help});

die ("fasta files dir (--fasta_file_dir) is required") unless (defined($config->{fasta_files_dir}));
die ("tmp dir (--tmp_dir) is required") unless (defined($config->{tmp_dir}));
die ("Mode is required (--mode update or --mode load (try --help))") unless (defined($config->{mode}));
die ("Not a valid value for mode. (--mode update or --mode reload (try --help))") unless ($config->{mode} eq 'update' || $config->{mode} eq 'load');
die ("Release version is required (--version)") unless (defined($config->{version}));
die ("Species is required (--species)") unless (defined($config->{species}));
my $variation_credentials = defined($config->{host}) && defined($config->{port}) && defined($config->{dbname}) && defined($config->{user}) && defined($config->{pass});
my $registry = defined($config->{registry});
die ("Database credentials or a registry file are required (try --help)") unless ($variation_credentials || $registry);
die ("Variation database credentials (--host, --dbname, --user, --pass, --port) are required") if (!$variation_credentials && !$registry);

set_up_db_connections($config);
if ($config->{resume_update} eq 'update_AA') {
	update_ancestral_alleles($config);
} elsif ($config->{resume_update} eq 'assign_AA') {
	assign_ancestral_alleles($config);
	preprocess_ancestral_alleles($config);
	update_ancestral_alleles($config);
} else {
	fetch_variation_ids_without_AA($config);
	fetch_variation_features($config);
	assign_ancestral_alleles($config);
	preprocess_ancestral_alleles($config);
	update_ancestral_alleles($config);
}

sub set_up_db_connections {
    my $config = shift;
    my $registry = 'Bio::EnsEMBL::Registry';
	my $species = $config->{species};
    my ($dbh, $dbc, $vdba);
    if ($config->{registry}) {
        $registry->load_all($config->{registry});
        $vdba = $registry->get_DBAdaptor($species, 'variation');
    } else { 
        $vdba = new Bio::EnsEMBL::Variation::DBSQL::DBAdaptor(
            -host => $config->{host},
            -user => $config->{user},
            -pass => $config->{pass},
            -port => $config->{port},
            -dbname => $config->{dbname},
        ) or die("Could not get a database adaptor for $config->{dbname} on $config->{host}:$config->{port}");
    }
    $dbh = $vdba->dbc->db_handle;
	$dbc = $vdba->dbc;
    $config->{dbc} = $dbc;
    $config->{dbh} = $dbh;
}

sub fetch_variation_ids_without_AA {
	my $config = shift;
	my $version = $config->{version};
	my $TMP_DIR = $config->{tmp_dir};
	my $dbh     = $config->{dbh};

	my $fh = FileHandle->new("$TMP_DIR/unassigned_AA_$version.txt", 'w');

	my $sth = $dbh->prepare(qq{
		SELECT variation_id FROM variation
		WHERE ancestral_allele IS NULL
		OR ancestral_allele = '';}, {mysql_use_result => 1}
	) or die $dbh->errstr;
	$sth->execute() or die $dbh->errstr;
	my $variation_id;
	$sth->bind_columns(\$variation_id);
	while ($sth->fetch) {
		print $fh $variation_id, "\n";
	}
	$sth->finish();
	$fh->close();	
}

sub fetch_variation_features {
	my $config = shift;
	my $version = $config->{version};
	my $TMP_DIR = $config->{tmp_dir};
	my $dbh = $config->{dbh};
	
	my $fh = FileHandle->new("$TMP_DIR/variation_features_$version.txt", 'w');
	
	my $sth = $dbh->prepare(qq{
		SELECT vf.variation_id, sr.name, vf.seq_region_start, vf.seq_region_end, vf.seq_region_strand
		FROM variation_feature vf, seq_region sr
		WHERE vf.seq_region_id = sr.seq_region_id
		ORDER BY vf.variation_id;}, {mysql_use_result => 1}
	) or die $dbh->errstr;
	$sth->execute() or die $dbh->errstr;
	my ($variation_id, $seq_region, $start, $end, $strand);
	$sth->bind_columns(\($variation_id, $seq_region, $start, $end, $strand));
	while ($sth->fetch) {
		print $fh join("\t", ($variation_id, $seq_region, $start, $end, $strand)), "\n";
	}

	$sth->finish();
	$fh->close();
}

sub assign_ancestral_alleles {
	my $config = shift;
	my $fasta_files_dir = $config->{fasta_files_dir};
	my $version = $config->{version};
	my $TMP_DIR = $config->{tmp_dir};

	my $db = Bio::DB::Fasta->new($fasta_files_dir, -reindex => 1); 
	my @sequence_ids = $db->get_all_ids;
	my %sequence_id_2_chr_number;	

	foreach my $sequence_id (@sequence_ids) {
		my @split = split(/:/, $sequence_id);
		$sequence_id_2_chr_number{$split[2]} = $sequence_id;
	}

	my $mode = $config->{mode};
	my $reload = 0;
	my $unassigned_AA = {};
	if ($mode eq 'update') {
		my $fh = FileHandle->new("$TMP_DIR/unassigned_AA_$version.txt", 'r');
		while (<$fh>) {
			chomp;
			my $variation_id = $_;
			$unassigned_AA->{$variation_id} = 1; 
		}
		$fh->close();
	} else {
		$reload = 1;
	}
	
	my $input = FileHandle->new("$TMP_DIR/variation_features_$version.txt", 'r');
	my $fh    = FileHandle->new("$TMP_DIR/variation_features_$version\_AA.txt", 'w');	

	while (<$input>) {
		chomp;
		my ($variation_id, $chrom, $start, $end, $strand) = split /\t/;
		if ($unassigned_AA->{$variation_id} || $reload) {
			if ($start > $end) { # Insertion, we don't assign AA to inserions
				print $fh "$variation_id\t-\n"; 
			} elsif (($end - $start) > 50) { # We don't assign AA to variations of length > 50
				print $fh "$variation_id\t-\n"; 
			} else {
				my $chrom_name = $sequence_id_2_chr_number{$chrom};
				if ($chrom_name && $start && $end) {
					my $AA = $db->seq("$chrom_name:$start,$end");
					if ($AA) {
						print $fh "$variation_id\t$AA\n";
					} else {
						print STDERR "No AA for $variation_id\n";
					}
				} else {
					print STDERR "Incomplete coords for variation_id $variation_id\n";
				}
			}
		}
	}

	$input->close();
	$fh->close();
}

sub preprocess_ancestral_alleles {	
	my $config = shift;
	my $version = $config->{version};
	my $TMP_DIR = $config->{tmp_dir};
	my $IN  = FileHandle->new("$TMP_DIR/variation_features_$version\_AA.txt", 'r');	
	my $OUT = FileHandle->new("$TMP_DIR/preprocessed_AA_$version.txt", 'w');

	my %ancestral_alleles;
	my $prev_id = 0;
	push @{$ancestral_alleles{$prev_id}}, '-';

	while (<$IN>) {
		chomp;
		my ($id, $AA) = split /\t/;
		if ($id ne $prev_id) {
			if ($prev_id > 0) {
				print $OUT $prev_id, "\t", join(',', @{$ancestral_alleles{$prev_id}}), "\n";
			}
			%ancestral_alleles = ();
		}
		push @{$ancestral_alleles{$id}}, $AA;
		$prev_id = $id;
	}	
	print $OUT $prev_id, "\t", join(',', @{$ancestral_alleles{$prev_id}}), "\n";
	$IN->close();
	$OUT->close();
}

sub update_ancestral_alleles {
	my $config = shift;
	my $version   =	$config->{version};
	my $TMP_DIR   = $config->{tmp_dir};
	my $tmp_file  = "update_AA_$version.txt";
	# consolidate information from preprocess step and write file: variation_id to AA

	my $IN = FileHandle->new("$TMP_DIR/preprocessed_AA_$version.txt", 'r');
	my $fh = FileHandle->new("$TMP_DIR/$tmp_file", 'w');
    my %assigned_alleles;
    while (<$IN>) {
        chomp;
        my ($id, $aas) = split/\t/;
        my @values = split(',', $aas);
        my %map;
        foreach my $value (@values) {
            if ($value =~ m/^[ACGTacgt]$/) {
                my $uc_value = uc $value;
                $map{$uc_value} = 1;
            }
        }
        my @final = keys %map;
        if (scalar @final == 1) {
            if ($final[0] =~ m/^[ACGT]$/) {
                my $AA = $final[0];
                $assigned_alleles{$AA} = 1;
				print $fh "$id\t$AA\n";
            }
        }
    }
	print STDERR "Assigned ancestral_alleles:\n";
    for my $aa (keys %assigned_alleles) {
        print STDERR $aa, "\n";
    }
	$IN->close();
	$fh->close();

	# load data into variation database

	$ImportUtils::TMP_DIR  = $TMP_DIR;
	$ImportUtils::TMP_FILE = $tmp_file;

	my $dbh = $config->{dbh};
	my $dbc = $config->{dbc};
	$dbh->do(qq{DROP TABLE IF EXISTS variation_id_AA}) or die $dbh->errstr;
	$dbh->do(qq{
    	CREATE TABLE `variation_id_AA` (
        	`variation_id` int(10) unsigned NOT NULL DEFAULT '0',
            `ancestral_allele` varchar(255) DEFAULT NULL,
            UNIQUE KEY `variation_idx` (`variation_id`));
   		}) or die $dbh->errstr;

	load($dbc, qw(variation_id_AA variation_id ancestral_allele));	

	# update variation table

	$dbh->do(qq{
		UPDATE variation_id_AA vaa JOIN variation v ON (vaa.variation_id = v.variation_id)
		SET v.ancestral_allele = vaa.ancestral_allele; 
	}) or die $dbh->errstr;

}

sub usage {

  print qq{
  Usage: perl import_ancestral_alleles.pl [OPTION]

  Import ancestral alleles into a Variation database

  Options:

      -help               Print this message

    Either provide all the necessary database credentials or a registry file:
    Database credentials are specified on the command line

      -host      Variation database host name (Required)
      -dbname    Variation database name (Required)
      -user      Variation database user (Required)
      -pass      Variation database password (Required)
      -port      Variation database port (Default: 3306)

      -registry  Location of registry file with database connections to variation database

	Species. Define name of the species for which to import ancestral alleles. Use latin name
	or alias. But make sure the alias is stored in the registry file.

	  -species    species name (Required)

	Version. Release number is used for creating tmp files. Could be necessary to compare results
	between releases.
	  -version    release number (Required)

    Tmp directory. The directory is needed for temporary files created during the import.

      -tmp_dir          Location of a directory to save temporary results (Required)

	  -fasta_files_dir  Location of fasta files storing the ancestral genome (Required)

	Mode. The script can be run in two modes: update or load. Load is run after a fresh dbSNP
	import. Update is run if only a few new variants have been imported.

	  -mode           update, load (Required)

	Resume update. If loading ancestral alleles into the database failed, the script can be
	resumed at the update step.

	  -resume_update  (Otional)


  } . "\n";
  exit(0);
}


