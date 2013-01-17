use strict;

=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

use warnings;

use DBI;
use File::Copy;
use File::Path qw(make_path remove_tree);
use Getopt::Long;
use Pod::Usage;

my ($host, $port, $user, $db_name, $TMP_DIR, $fasta_dir, $prefix, $help);
my $postfix = ".fa";
my (%seq_region_id_2_name, %seq_region_name_2_id, $seq_region_id, $seq_region_name, %load_map, $count, $sum_variation_feature);

GetOptions(
	'TMP_DIR=s'   => \$TMP_DIR,
	'db_name=s'   => \$db_name,
	'host=s'      => \$host,
	'user=s'      => \$user,
	'port=s'      => \$port,
	'fasta_dir=s' => \$fasta_dir,
	'prefix=s'    => \$prefix,
	'help|h'      => \$help,
) or die pod2usage(1);
pod2usage(1) if $help;

die "Argument list is not complete, try --help for usage" unless ($TMP_DIR && $db_name && $host && $user && $port && $fasta_dir && $prefix);

my $logfile = "logfile_ancestral_allele_pre_processing.txt";
open LOG_FILE, ">$TMP_DIR" . $logfile or die("Could not open log file: " . $TMP_DIR . $logfile); 
my @parameters = ($host, $port, $user, $db_name, $TMP_DIR, $fasta_dir, $prefix);
print LOG_FILE join("\n", @parameters);

my $dbh = DBI->connect("DBI:mysql:database=$db_name;host=$host;port=$port;user=$user", {RaiseError => 1});
my $sth;

# fetch seq_region id and name
print LOG_FILE "Fetch and build mapping for seq_region id to name and name to id\n";
$sth = $dbh->prepare(qq{select seq_region_id, name from seq_region;});
$sth->execute() or die "Could not execute statement: " . $sth->errstr;
$sth->bind_columns(\$seq_region_id, \$seq_region_name);
while ($sth->fetch) {
	$seq_region_id_2_name{$seq_region_id} = $seq_region_name;
	$seq_region_name_2_id{$seq_region_name} = $seq_region_id;
}
$sth->finish();

# fetch variation_feature count for each seq_region
print LOG_FILE "Fetch variation_feature count for each seq_region id\n";
$sth = $dbh->prepare(qq{select count(*) from variation_feature where seq_region_id=?;});

foreach my $seq_region_id (keys %seq_region_id_2_name) {
	$sth->execute($seq_region_id);
	$sth->bind_columns(\$count);
	$sth->fetch;
	if ($count) {
		$load_map{$seq_region_id_2_name{$seq_region_id}} = $count;
		$sum_variation_feature += $count;
	}
}

$sth->finish;

# load balance distribute over files
print LOG_FILE "Compute a load balance scheme: Distribute variation_feature over a given number of files\n";
my %load;
my $processes = 10;
my $avg_load = ($sum_variation_feature / $processes);
my $tmp_load = 0;
my $file_number = 1;

foreach my $seq_region_name (keys %load_map) {
	if ($tmp_load < $avg_load) {
		push @{$load{$file_number}} , $seq_region_name;
		$tmp_load += $load_map{$seq_region_name};
	} else {
		$file_number++;
		push @{$load{$file_number}}, $seq_region_name;
		$tmp_load = $load_map{$seq_region_name};	
	}
}

my $load_file = $TMP_DIR . "load_schema_". $db_name . ".txt";
print LOG_FILE "Save load balance schema in $load_file\n";
open OUTPUT, ">" . $load_file or die $!;
foreach my $file_number (sort keys %load) {
	print OUTPUT $file_number, "\t", join(";", sort @{$load{$file_number}}), "\n";
}
close OUTPUT;

#fetch variation feature by seq_region 
print LOG_FILE "Fetch variation_feature by seq_region and save in file containing given seq_regions\n";
$sth = $dbh->prepare(qq{select seq_region_id, seq_region_start, seq_region_end, seq_region_strand, variation_id
				 from variation_feature 
                 where seq_region_id=?;});

foreach my $file_number (sort keys %load) {
	open OUTPUT, ">" . $TMP_DIR . "variation_feature.$file_number.txt" or die $!;
	foreach my $seq_region_name (@{$load{$file_number}}) {
		$sth->execute($seq_region_name_2_id{$seq_region_name});
		while (my $row = $sth->fetch) {
			my @a = map {defined($_) ? $_ : 'NULL'} @$row;
			print OUTPUT $seq_region_name, "\t", join("\t", @a), "\n";
		}
	}
	close OUTPUT;
}
$sth->finish;

# redistribute fasta files
print LOG_FILE "Redistribute fasta files asccording to load balance schema.\n";
opendir(DIR, $fasta_dir);
my @fasta_files = readdir(DIR);
closedir(DIR);

my %fasta_file_2_folder_number;
my @sub_processes;
open INPUT, "<" . $load_file or die $!;
while (<INPUT>) {
	chomp;
	my ($folder_number, $ids) = split /\t/;
	push @sub_processes, $folder_number;
	if (-d $fasta_dir . $folder_number) {
		remove_tree($fasta_dir . $folder_number);
		print LOG_FILE "Remove recursively old folder " . $fasta_dir . $folder_number . "\n";
	}
	make_path($fasta_dir . $folder_number);
	print LOG_FILE "Create directory " . $fasta_dir . $folder_number . "\n";
	my @seq_regions = split(/;/, $ids);
	foreach my $seq_region (@seq_regions) {
		# map which fasta file belongs in which folder according to the load file
		$fasta_file_2_folder_number{$prefix . $seq_region . $postfix} = $folder_number;
	}
}
close INPUT; 

print LOG_FILE "Copy fasta files in sub directories as given by load balance schema.\n";
foreach my $fasta_file (@fasta_files) {
	if (defined $fasta_file_2_folder_number{$fasta_file}) {
		copy($fasta_dir . $fasta_file, $fasta_dir . $fasta_file_2_folder_number{$fasta_file} . "/" . $fasta_file);
	}
}

print LOG_FILE "Pre processing is finished.\n";
close LOG_FILE;

__END__

=head1 NAME

ancestral_allele_pre_processing.pl

=head1 DESCRIPTION

# steps:
fetch seq_region id and name
# fetch for each seq_region number of variation_feature on that region
# compute a load balance: store wich seq_region is contained in which file
# finally fetch variation_feature by seq_region and store in files according to load balance

=head1 SYNOPSIS

ancestral_allele_pre_processing.pl [arguments]

=head1 ARGUMENTS

=over 4

=item B<--TMP_DIR DIR>

Save all intermediate results here.

=item B<--db_name NAME>

Database name of variation database to which ancestral allele data should be added.

=item B<--host NAME>

Connection parameter to connect to server.

=item B<--user NAME>

Connection parameter to connect to server.

=item B<--port NUMBER>

Connection parameter to connect to server.

=item B<--fasta_dir DIR>

Directory with fasta files for the ancestral genome.

=item B<--prefix NAME>

Prefix of fasta file. The part before the actual seq_region name.

=item B<--help>

Show help

=head1

For help with this script address questions to dev@ensembl.org

