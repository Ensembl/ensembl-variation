# Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

use warnings;

use DBI;
use FindBin qw($Bin);
use Getopt::Long;
use Pod::Usage;

my ($tmp_host, $tmp_port, $tmp_user, $tmp_password, $tmp_db_name, $host, $port, $user, $password, $db_name, $TMP_DIR, $processes, $help);

GetOptions(
	'tmp_host=s'    => \$tmp_host,
	'tmp_port=s'    => \$tmp_port,
	'tmp_user=s'    => \$tmp_user,
	'tmp_password'  => \$tmp_password,
	'tmp_db_name=s' => \$tmp_db_name,
	'host=s'        => \$host,
	'port=i'        => \$port,
	'user=s'        => \$user,
	'password=s'    => \$password,
	'db_name=s'     => \$db_name,
	'TMP_DIR=s'     => \$TMP_DIR,
	'processes=i'   => \$processes,
	'help|h'        => \$help,
); or die pod2usage(1);
pod2usage(1) if $help;

die "Argument list is not complete, try --help for usage" unless ($tmp_host && $tmp_port && $tmp_user && $tmp_password && $tmp_db_name && $host && $port && $user && $password && $db_name && $TMP_DIR && $processes);

my $count = 1;

while ($count <= $processes) {
	my $input  = "$TMPDIR/ancestral_allele_output.$count.txt";
	my $output = "$TMPDIR/ancestral_allele_db_import.$count.txt";
	open INPUT, "<$input" or die "Can't open $input";
	open OUTPUT, ">$output" or die "Can't open $output";
	while (<INPUT>) {
		chomp;
		my ($variation_id, $seq_region_name, $seq_region_start, $seq_region_end, $seq_region_strand, $ancestral_allele) = split /\t/;
		print OUTPUT "$variation_id\t$ancestral_allele\n";
	}
	close INPUT;
	close OUTPUT;
	$count++;
}
# commands
my $tmp_db_connect_2_mysql = "--host=$tmp_host --port=$tmp_port --user=$tmp_user --password=$tmp_password";
my $tmp_db_connect_command = "mysql $tmp_db_connect_2_mysql $tmp_db_name";

# create database drop old database if exists before
my $drop_database = "mysql $tmp_db_connect_2_mysql -e \"DROP DATABASE IF EXISTS $tmp_db_name;\";";
print $drop_database, "\n";
system($drop_database);

my $create_database = "mysql $tmp_db_connect_2_mysql -e \"CREATE DATABASE $tmp_db_name;\";";
print $create_database, "\n";
system($create_database);

# create table
my $create_table = "$tmp_db_connect_command < $Bin/ancestral_allele_table.sql;";
print $create_table, "\n";
system($create_table);

# insert data into table
$count = 1;
while ($count <= $processes) {
	my $insert_ancestral_allele = "$tmp_db_connect_command -e \"LOAD DATA LOCAL INFILE '$TMPDIR/ancestral_allele_db_import.$count.txt' INTO TABLE ancestral_alleles (variation_id, ancestral_allele)\";";
	print $insert_ancestral_allele, "\n";
	system($insert_ancestral_allele);
	$count++;
}

# sort by variation id
my $dbh = DBI->connect("DBI:mysql:database=$tmp_db_name;host=$tmp_host;port=$tmp_port;user=$tmp_user;password=$tmp_password", {RaiseError => 1});
my $sth;

$sth = $dbh->prepare(qq{select count(*) from ancestral_alleles});
$sth->execute();
my $ancestral_alleles = $sth->fetch();
$sth->finish();

$sth = $dbh->prepare(qq{select variation_id, ancestral_allele 
						from ancestral_alleles
                        order by variation_id
                        }, {mysql_use_result => 1});

$sth->execute();
$count = 0;
my $process = 1;
my $num_processes = 10;
my $sub_ancestral_alleles = int($ancestral_alleles->[0] / $num_processes);
my $previous_variation_id = 0;
my $current_variation_id;

my $file = $TMPDIR . "/sorted_variation_id.$process.txt";
open FILE, ">$file" or die ("Could not open $file");
  
while (my $row = $sth->fetch()){
	$count++;
	$current_variation_id = $row->[0];

	if ($current_variation_id ne $previous_variation_id) {
		if ((($sub_ancestral_alleles * $process) <= $count) && ($process < $num_processes)) {
			close(FILE);
			$process++;
			$file = $TMPDIR . "/sorted_variation_id.$process.txt";
			open FILE, ">$file" or die ("Could not open $file");
		}
	}
	my @a = map {defined($_) ? $_ : '-'} @$row;
	print FILE join("\t", @a) . "\n";
	$previous_variation_id = $current_variation_id;
}
close(FILE);

# assign ancestral allele to variation
# consolidate ancestral alleles in case of several variation features
my $i;
my %alleles;
my $previous = -1;
my ($variation_id, $ancestral_allele);

for ($i = 1; $i <= 10; $i++) {
	open INPUT, "<$TMPDIR/sorted_variation_id.$i.txt" or die "Could not open $TMPDIR/sorted_variation_id.$i.txt";
	open OUTPUT, ">$TMPDIR/update_ancestral_allele_in_variation_table.$i.txt" or die "Could not open $TMPDIR/final_ancestral_alleles.$i.txt";
	while (<INPUT>) {
		chomp;
		($variation_id, $ancestral_allele) = split /\t/;
		$alleles{$variation_id}{$ancestral_allele} = 1;
	
		if (scalar keys %alleles > 1) {
			foreach my $key (keys %alleles) {
				if ($key != $variation_id) {
					if (scalar (keys %{$alleles{$key}}) > 1) { # different ancestral alleles for variation
						print OUTPUT "update ignore variation set ancestral_allele=NULL where variation_id=$key;\n";
					} else {
						foreach my $allele (keys %{$alleles{$key}}) {
							if ($allele =~ /(A-Z)*/) {
								print OUTPUT "update ignore variation set ancestral_allele=\"$allele\" where variation_id=$key;\n";
							} else {
								print OUTPUT "update ignore variation set ancestral_allele=NULL where variation_id=$key;\n";
							}
						} 
					}
				}
			}
			%alleles = ();
			$alleles{$variation_id}{$ancestral_allele} = 1;
		}
	} # end while
	# take care of last variation id in file (could be more than one line)
	if (scalar (keys %{$alleles{$variation_id}}) > 1) {
		print OUTPUT "update ignore variation set ancestral_allele=NULL where variation_id=$variation_id;\n";
	} else {
		foreach my $allele (keys %{$alleles{$variation_id}}) {
			if ($allele =~ /[A-Z]*/) {
				print OUTPUT "update ignore variation set ancestral_allele=\"$allele\" where variation_id=$variation_id;\n";
			} else {
				print OUTPUT "update ignore variation set ancestral_allele=NULL where variation_id=$variation_id;\n";
			}
		} 
	}
	close INPUT;
	close OUTPUT;
}

# update ancestral alleles in working variation database 
$db_connect_2_mysql = "--host=$host --port=$port --user=$user --password=$password";
$db_connect_command = "mysql $db_connect_2_mysql $db_name";

my $update_table;
for (my $i = 1; $i <= 10; $i++) {
	$update_table = "$db_connect_command < $TMPDIR/update_ancestral_allele_in_varition_table.$i.txt";
	print $update_table, "\n";
	system($update_table);
}

__END__

=head1 NAME

ancestral_allele_post_processing.pl

=head1 DESCRIPTION

# steps:
# redistribute results
# sort by variation id and distribute on several files
# consolidate ancestral alleles for variations with several features
# import results into variation database
# transform ancestral_allele_output to database update file

=head1 SYNOPSIS

ancestral_allele_post_processing.pl [arguments]

=head1 ARGUMENTS

=over 4

=item B<--tmp_host NAME>

Connection parameter to connect to server. Used for building a tmp db for sorting.

=item B<--tmp_port>

Connection parameter to connect to server. Used for building a tmp db for sorting.

=item B<--tmp_user>

Connection parameter to connect to server. Used for building a tmp db for sorting.

=item B<--tmp_password>

Connection parameter to connect to server. Used for building a tmp db for sorting.

=item B<--tmp_db_name>

Connection parameter to connect to server. Used for building a tmp db for sorting.

=item B<--host>

Connection parameter to connect to server. Used to update ancestral alleles.

=item B<--port>

Connection parameter to connect to server. Used to update ancestral alleles.

=item B<--user>

Connection parameter to connect to server. Used to update ancestral alleles.

=item B<--password>

Connection parameter to connect to server. Used to update ancestral alleles.

=item B<--db_name>

Connection parameter to connect to server. Used to update ancestral alleles.

=item B<--TMP_DIR>

Directory where intermediate results for computation of ancestral alleles are saved.

=item B<--processes NUMBER>

Number of files containg computed ancestral alleles.

=item B<--help>

Show help

=head1

For help with this script address questions to dev@ensembl.org
