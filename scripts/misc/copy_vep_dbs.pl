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


use Getopt::Long;
use DBI;

my $config = {};

GetOptions(
  $config,
  
  'help|h',
  'version|v=i',
  
	'target_host|h=s',
	'target_port|P=i',
	'target_password|p=s',
  
	'source_hosts|sh=s',
	'source_user|su=s',
	'source_port|sP=i',
	'source_password|sp=s',
  
  'copy_script|c=s',
  'force|f',
  'tables|t=s',
  
  'pattern|p=s',
) or die "ERROR: Could not parse command line options\n";

if(defined($config->{help}) || !@ARGV) {
  print qq{
This script creates variation DBs with a minimal set of tables required for the
VEP web interface to work. By default it copies all variation DBs found on the
staging servers.

Usage:
perl copy_vep_dbs.pl [options]

# Required options

--version | -v          Database version to copy

--target_host | -h     Target DB host
--target_port | -P     Target DB port
--target_pass | -p     Target DB password

NB target DB user is hard-coded to ensadmin in CopyDBoverServer.pl


# Options with defaults that shouldn't need to be overwritten

--source_host | -sh      Source DB host(s) (default: ens-staging1,ens-staging2)
--source_user | -su      Source DB user (default: ensro)
--source_port | -sP      Source DB port (default: 3306)
--source_pass | -sp      Source DB password

--copy_script | -c      Path to CopyDBoverServer.pl
--tables | -t           Comma-separated list of tables to copy


# Other

--pattern | -p          Select only databases that match '%pattern%'
};
  exit(0);
}

# set defaults
$config->{source_hosts}    ||= 'ens-staging,ens-staging2';
$config->{source_user}     ||= 'ensro';
$config->{source_port}     ||= 3306;
$config->{target_port}     ||= 3306;
$config->{copy_script}     ||= '../../../ensembl/misc-scripts/CopyDBoverServer.pl';
$config->{tables}          ||= join(',', qw(
  associate_study
  attrib
  attrib_type
  failed_variation
  failed_structural_variation
  meta
  meta_coord
  seq_region
  source
  structural_variation
  structural_variation_feature
  study
  variation
  variation_feature
  variation_synonym
));

die("ERROR: Could not find DB copy script ".$config->{copy_script}."\n") unless defined($config->{copy_script}) && -e $config->{copy_script};

for(qw(host port password)) {
  die("ERROR: No target $_ given - use --target_$_ [$_]") unless defined($config->{'target_'.$_});
}

# check version defined
die "ERROR: No Ensembl DB version defined - use --version [version]\n" unless defined($config->{version});

# parse hosts
$config->{source_hosts} = [split /\,/, $config->{source_hosts}];

foreach my $host(@{$config->{source_hosts}}) {
	debug("Getting database list for host $host");
	my $db_list = get_database_list($config, $host);
	debug("Found ".(scalar @$db_list)." valid databases\n");
	
	foreach my $db(@$db_list) {
		debug("Copying data from $db");
    
    copy_database($config, $host, $db);
	}
}

sub get_database_list {
	my $config = shift;
	my $host   = shift;
  
  my $dbc = connect_to_host($host, $config->{source_port}, $config->{source_user}, $config->{source_password});
	
	my $version = $config->{version};

	my $sth = $dbc->prepare(qq{
		SHOW DATABASES LIKE '%\_variation\_$version%'
	});
	$sth->execute();
	
	my $db;
	$sth->bind_columns(\$db);
	
	my @dbs;
	push @dbs, $db while $sth->fetch;
	$sth->finish;
  
  # don't include master
  @dbs = grep {$_ !~ /master/} @dbs;
  
  my $pattern = $config->{pattern};
  @dbs = grep {$_ =~ /$pattern/} @dbs if $pattern;

	return \@dbs;
}

sub connect_to_host {
  my ($host, $port, $user, $pass) = @_;
  
	my $connection_string = sprintf(
    "DBI:mysql(RaiseError=>1):host=%s;port=%s;db=mysql",
    $host,
    $port
  );
	
	# connect to DB
	my $dbc = DBI->connect(
	  $connection_string, $user, $pass
	);
  
  return $dbc;
}

sub copy_database {
  my $config = shift;
  my $host = shift;
  my $db = shift;
  
  my $cmd = sprintf(
    'perl %s --source %s@%s:%i --target %s@%s:%i --pass %s --only_tables %s%s',
    $config->{copy_script},
    $db, $host, $config->{source_port},
    $db, $config->{target_host}, $config->{target_port},
    $config->{target_password},
    $config->{tables},
    $config->{force} ? ' --force' : ''
  );
  
  print "$cmd\n";
  die("ERROR: Copy script failed") if system($cmd);
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
