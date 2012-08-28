#!/usr/bin/perl

use DBI;
use Getopt::Long;
use ImportUtils qw(load);

my $config = {};

my @exclude = qw(transcript_variation_id hgvs_genomic hgvs_protein hgvs_transcript somatic codon_allele_string);

GetOptions(
    $config,
	'tmpdir=s',
	'tmpfile=s',
	'host|h=s',
	'user|u=s',
	'port|P=i',
	'password|p=s',
	'db|d=s',
	'table|t=s'
) or die "ERROR: Could not parse command line options\n";


foreach my $opt(qw(tmpdir host user port password db)) {
	die "ERROR: --$opt not specified\n";
}


$config->{table} ||= 'transcript_variation';
my $source_table = $config->{table};
my $table = 'MTMP_'.$source_table.'_consequence_types';

my $TMP_DIR = $config->{tmpdir};
my $TMP_FILE = $config->{tmpfile};

$ImportUtils::TMP_DIR = $TMP_DIR;
$ImportUtils::TMP_FILE = $TMP_FILE;

# check options
for(qw(host user port password db tmpdir tmpfile)) {
	die("ERROR: $_ not defined, use --$_\n") unless defined $config->{$_};
}

my $dbc = DBI->connect(
    sprintf(
        "DBI:mysql(RaiseError=>1):host=%s;port=%s;db=%s",
        $config->{host},
        $config->{port},
        $config->{db},
    ), $config->{user}, $config->{password}
);

print "Getting column definition\n";

# get column definition from transcript_variation
my $sth = $dbc->prepare(qq{
	SHOW CREATE TABLE $source_table
});
$sth->execute();

my $create_sth = $sth->fetchall_arrayref->[0]->[1];
$sth->finish;

# convert set to enum
$create_sth =~ s/^set/enum/;

# rename table
$create_sth =~ s/TABLE \`$source_table\`/TABLE \`$table\`/;

# filter out some columns
$create_sth =~ s/\`?$_.+?,// for @exclude;

# filter out some indices
$create_sth =~ s/AUTO_INCREMENT=\d+//;
$create_sth =~ s/$_.+,// for ('PRIMARY KEY', 'KEY `somatic', 'KEY `cons');
$create_sth =~ s/\`somatic\`\)//;

# remove final comma
$create_sth =~ s/,(\s+\))/$1/;

#print "$create_sth\n";

print "Creating table $table\n";

# create a new table
$dbc->do($create_sth);

# get columns of new table
$sth = $dbc->prepare(qq{
	DESCRIBE $table
});
$sth->execute();
my @cols = map {$_->[0]} @{$sth->fetchall_arrayref};
$sth->finish;


print "Populating table $table\n";


$sth = $dbc->prepare(qq{
	SELECT count(*)
	FROM $source_table
});
$sth->execute();
my $row_count = $sth->fetchall_arrayref->[0]->[0];
$sth->finish;

# populate it
$sth = $dbc->prepare(qq{
	SELECT *
	FROM $source_table
	#LIMIT ?, ?
}, {mysql_use_result => 1});

open OUT, ">$TMP_DIR/$TMP_FILE";

$sth->execute();
	
while(my $row = $sth->fetchrow_hashref()) {
	my $cons = $row->{consequence_types};
	
	foreach my $con(split /\,/, $cons) {
		my %vals = %{$row};
		$vals{consequence_types} = $con;
		delete $vals{$_} for @exclude;
		my @values = map {defined $vals{$_} ? $vals{$_} : '\N'} @cols;
		
		print OUT join("\t", @values);
		print OUT "\n";
		
		#$sth2->execute(@values);
	}
}

$sth->finish;

close OUT;

load($dbc, $table);
