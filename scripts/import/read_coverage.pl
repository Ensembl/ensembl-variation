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


use strict;
use warnings;

use Getopt::Long;
use DBH;

use Bio::EnsEMBL::Mapper::RangeRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use ImportUtils qw(debug);

my ($TMP_DIR, $TMP_FILE); #global variables for the tmp files and folder
my $range_registry = []; #reference to an array containing an rr for each of the possible levels (from 1-6)
#hash that, for a given individual_name, returns the individual_id in the database, if present. The last 4 samples are defined individuals in dbSNP
#with a predefined id that do not change
my %individual_id = ();
my ($chost, $cport, $cdbname, $cuser, $cpass,
    $vhost, $vport, $vdbname, $vuser, $vpass,
    $read_file, $max_level,
    $ind_file);

  GetOptions('chost=s'     => \$chost,
             'cuser=s'     => \$cuser,
             'cpass=s'     => \$cpass,
             'cport=i'     => \$cport,
             'cdbname=s'   => \$cdbname,
	     'vhost=s'     => \$vhost,
             'vuser=s'     => \$vuser,
             'vpass=s'     => \$vpass,
             'vport=i'     => \$vport,
             'vdbname=s'   => \$vdbname,
             'tmpdir=s'    => \$ImportUtils::TMP_DIR,
             'tmpfile=s'   => \$ImportUtils::TMP_FILE,
	     'readfile=s'  => \$read_file,
	     'maxlevel=i'  => \$max_level,
	     'indfile=s'   => \$ind_file);

#added default options
$chost    ||= 'ecs2';
$cuser    ||= 'ensro';
$cport    ||= 3364;

usage('-vdbname argument is required') if(!$vdbname);
usage('-cdbname argument is required') if(!$cdbname);
usage('-readfile argument is required (name of file with read information)') if (!$read_file);
usage('-maxlevel argument is required (max level coverage calculated)') if (!$max_level);

my $dbCore = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (-host   => $chost,
     -user   => $cuser,
     -pass   => $cpass,
     -port   => $cport,
     -dbname => $cdbname);

my $dbVar = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new
    (-host   => $vhost,
     -user   => $vuser,
     -pass   => $vpass,
     -port   => $vport,
     -dbname => $vdbname
     );


$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

if ($ind_file)
{
    open IND, $ind_file;

    my ($ind_name, $ind_id);

    foreach my $ind_line (<IND>)
    {
        chomp $ind_line;

        ($ind_name, $ind_id) = split /\t/, $ind_line;

        $individual_id{$ind_name} = $ind_id;
    }

    close IND;
}

&load_individuals($dbVar,\%individual_id); #first of all, load the hash with the individual_id
&initialize_range_registry($range_registry, $max_level);
my $pair; #reference to an array with id,start,end format
$read_file =~ /^.*\/(.*)\.mapped$/;  #extract the chromosome from the file name
my $region = $1;print "the region is $region\n";
my $individuals = {}; #reference to a hash containing all individuals present in the chromosome
open IN, "$read_file" or die "Could not open file $read_file with read information in the format seq_region_id\tstart\tend:$!\n";
while (<IN>){
    chomp;  #remove the last \n    
    #line format
    #LIB     ID_read start   end
    
    #NA17109 7495470 1833    2533
    #TSC     8093917 617     957
    ($pair) = [split /\t/];
#    splice @{$pair},1,1;    #we need to get rid of the second column, ID_read
    if ($pair->[0] eq ''){$pair->[0] = 'Unknown'}
    &register_range_level($range_registry,$pair,1,$individuals,$max_level);
}
close IN;
    
    
#get the dbID for the region where getting the reads
my $slice_adaptor = $dbCore->get_SliceAdaptor;
my $slice = $slice_adaptor->fetch_by_region('toplevel',$region);
my $seq_region_id = $slice->get_seq_region_id();
    
&import_data(1,$range_registry,"read_coverage",$individuals,$seq_region_id,\%individual_id);
&import_data($max_level,$range_registry,"read_coverage",$individuals,$seq_region_id,\%individual_id);
&update_meta_table($dbVar,$max_level);
    
debug("File $read_file finished");

sub initialize_range_registry{
    my $range_registry = shift;
    my $max_level = shift;
    foreach my $level (1..$max_level){
	$range_registry->[$level] = Bio::EnsEMBL::Mapper::RangeRegistry->new();
    }

    return;
}

#for a given range and the one covered, returns the inverted
sub invert_pair{
    my $range = shift; #initial range of the array
    my $pairs = shift; #listref with the pairs that have been added to the range

    my @inverted_pairs;
    my $inverted;

    my $rr = Bio::EnsEMBL::Mapper::RangeRegistry->new();

    foreach my $pair (@{$pairs}){
	$rr->check_and_register(1,$pair->[0],$pair->[1]);
    }
    return $rr->check_and_register(1,$range->[0],$range->[1]); #register again the range
}

sub register_range_level{
    my $range_registry = shift;
    my $range = shift;
    my $level = shift;
    my $individuals = shift;
    my $max_level = shift;

    my $individual_name = shift @{$range}; #get the individual from the array
    return if ($level > $max_level);
    $individuals->{$individual_name}++;
    my $rr = $range_registry->[$level];
    my $pair = $rr->check_and_register($individual_name,$range->[0],$range->[1]);
    my $pair_inverted = &invert_pair($range,$pair);
    return if (!defined $pair_inverted);
    foreach my $inverted_range (@{$pair_inverted}){
	unshift @{$inverted_range},$individual_name; #add again the name of the individual
	&register_range_level($range_registry,$inverted_range,$level+1,$individuals,$max_level);
    }
}

#loads a table in the database with the read coverage data
sub import_data{
    my $coverage_level = shift;
    my $range_registry = shift;
    my $table = shift;
    my $individuals = shift; #reference to a hash containing the individuals in the chromosome
    my $seq_region_id = shift;
    my $individual_id = shift;
    
    open FH,">>$TMP_DIR/$TMP_FILE";
    foreach my $individual_name (keys %{$individuals}){
	if ($individual_name ne '' && defined $range_registry->[$coverage_level]->get_ranges($individual_name)){
	    foreach my $range (@{$range_registry->[$coverage_level]->get_ranges($individual_name)}){
		print FH join("\t",$seq_region_id,$range->[0],$range->[1],$coverage_level,$individual_id->{$individual_name}),"\n";
	    }
	}
    }
    close FH;
}

#will be necessary to load the hash with the relation individual_name => individual_id at the beginning of the script. Some of the individuals in
#Sanger data are already defined in the database, and the ID is permanent. Others, will be necessary to be added (if not already) in the database
sub load_individuals{
    my $dbVar = shift;
    my $individual_id = shift;
    
    #first, get the individuals that need to be loaded from the database
    foreach my $individual_name (keys %{$individual_id}){
	#only load ones without individual_id
	if ($individual_id->{$individual_name} eq ''){
	    my $individualAdaptor = $dbVar->get_IndividualAdaptor();
	    my $individual = ($individualAdaptor->fetch_all_by_name($individual_name))->[0];
	    my ($sample_id);
	    #the individual is not in the database, will need to create it
	    if (!defined $individual){
		#insert the sample
		$dbVar->dbc()->do(qq{INSERT INTO sample (name,description) VALUES ("$individual_name",'Individuals used to get read information')});
		$sample_id = $dbVar->dbc()->db_handle->{'mysql_insertid'}; #get the id inserted
		#insert the individual
		$dbVar->dbc()->do(qq{INSERT INTO individual (sample_id,gender) VALUES ("$sample_id",'Unknown')});
		#get the population. We will used Unknown, since the origin it is not clear
		my $population_adaptor = $dbVar->get_PopulationAdaptor();
		my $population = $population_adaptor->fetch_by_name('UNKNOWN');
		my $population_id = $population->dbID();
		#insert the relation individual-population
		$dbVar->dbc()->do(qq{INSERT INTO individual_population (individual_sample_id,population_sample_id) VALUES ("$sample_id","$population_id")});
	    }
	    else{
		$sample_id = $individual->dbID();
	    }
	    $individual_id->{$individual_name} = $sample_id;
	}
    }
}

#
# updates the meta table
#
sub update_meta_table {
  my $dbVar  = shift;
  my $max_level = shift;

  #find out if the meta table has been updated with the coverage levels calculated
  my $rca = $dbVar->get_ReadCoverageAdaptor();
  my $levels = $rca->get_coverage_levels();
  
  
  my $sth = $dbVar->dbc->prepare
    ('INSERT INTO meta (meta_key,meta_value) VALUES (?,?)');

  foreach (1,$max_level){
      my $cur_level = $_;
      if ((grep {$cur_level eq $_} @{$levels}) == 0){
	  $sth->execute('read_coverage.coverage_level',$_);
      }
  }
  
  $sth->finish();

  return;
}

sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: perl read_coverage.pl <options>

options:
    -chost <hostname>    hostname of core Ensembl MySQL database (default = ecs2)
    -cuser <user>        username of core Ensembl MySQL database (default = ensro)
    -cpass <pass>        password of core Ensembl MySQL database
    -cport <port>        TCP port of core Ensembl MySQL database (default = 3364)
    -cdbname <dbname>    dbname of core Ensembl MySQL database
    -vhost <hostname>    hostname of variation MySQL database to write to
    -vuser <user>        username of variation MySQL database to write to (default = ensadmin)
    -vpass <pass>        password of variation MySQL database to write to
    -vport <port>        TCP port of variation MySQL database to write to (default = 3306)
    -vdbname <dbname>    dbname of variation MySQL database to write to
    -tmpdir <dir>        temp directory to use (with lots of space!)
    -tmpfile <filename>  name of temp file to use
    -readfile <filename> name of file with read information
    -maxlevel <number>   max level of coverage calculated
EOF

  die("\n$msg\n\n");
}
