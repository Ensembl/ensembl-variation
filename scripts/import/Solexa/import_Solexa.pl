#!/usr/bin/env perl
# Copyright 2013 Ensembl
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
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


use strict;
use warnings;

use Getopt::Long;
use ImportUtils qw(debug load dumpSQL);
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use DBH;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use FindBin qw( $Bin );
use DBI qw(:sql_types);
use Data::Dumper;

my $species = 'human';
my $file;

GetOptions('tmpdir=s'  => \$ImportUtils::TMP_DIR,
	   'tmpfile=s' => \$ImportUtils::TMP_FILE,
	   'species=s' => \$species,
	   'solexa=s'  => \$file
	   );

warn("Make sure you have a updated ensembl.registry file!\n");
die("Need to enter a valid file with the Consensus data") if(! -e $file);
my $registry_file ||= $Bin . "/ensembl.registry";

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $dbVar = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');

my $TMP_DIR  = $ImportUtils::TMP_DIR;
my $TMP_FILE = $ImportUtils::TMP_FILE;

my $last_sample_id = get_last_table_id($dbVar,"sample"); #get last sample in the database
my $last_source_id = get_last_table_id($dbVar,"source");
my $last_variation_id = get_last_table_id($dbVar,"variation");
my $last_allele_id = get_last_table_id($dbVar,"allele");
my $last_variation_feature_id = get_last_table_id($dbVar,"variation_feature");
my $chr = 6;
my $start_BAC = 30537344;

import_Sample_table($dbVar, $last_sample_id);
import_Source_table($dbVar,$last_source_id);
import_Population_table($dbVar,$last_sample_id);
parse_Consensus_file($dbVar,$file,$last_variation_id,$last_allele_id, $last_source_id,$last_sample_id,$last_variation_feature_id,$chr,$start_BAC);

sub import_Sample_table{
    my $dbVar = shift;
    my $last_sample_id = shift;
    #write a new entry in the table, for the BAC sequence
    write_file("sample.txt",$last_sample_id+1,"Human BAC","First human sequence from Solexa");
    rename("$TMP_DIR/sample.txt", "$TMP_DIR/$TMP_FILE");
    load($dbVar->dbc,"sample","sample_id","name","description");
 }

sub import_Source_table{
    my $dbVar = shift;
    my $last_source_id = shift;
    #write a new entry in the table, for the BAC sequence
    write_file("source.txt",$last_source_id+1,"Solexa");
    rename("$TMP_DIR/source.txt", "$TMP_DIR/$TMP_FILE");
    load($dbVar->dbc,"source","source_id","name");
 }

sub import_Population_table{
    my $dbVar = shift;
    my $last_sample_id = shift;
    #write a new entry in the table, for the BAC sequence
    write_file("population.txt",$last_sample_id+1,1);
    rename("$TMP_DIR/population.txt", "$TMP_DIR/$TMP_FILE");
    load($dbVar->dbc,"population","sample_id","is_strain");
 }

sub parse_Consensus_file{
    my $dbVar = shift;
    my $file = shift;   
    my $last_variation_id = shift;
    my $last_allele_id = shift;
    my $last_source_id = shift;
    my $last_sample_id = shift;
    my $last_variation_feature_id = shift;
    my $chr = shift;
    my $start_BAC = shift;

    my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor('human','core','Slice');
    my $slice = $slice_adaptor->fetch_by_region('chromosome',$chr);
    my @line;
    my $position;
    my $ref_base;
    my $call_base;
    my $score;
    my $new_variation_id;
    my $new_allele_id;
    my $new_variation_feature_id;
    open IN, "<$file" || die "Could not open $file: $!\n";
    while(<IN>){
	next if(/^\#/); #skip comments
	chomp;
	@line = split/\t/;
	$position = $line[0];
	$ref_base = $line[1];
	$call_base = $line[6];
	$score = $line[7];
	next if ($ref_base eq $call_base); 
	next if ($score < 0.9);
	#we now have a variation, calculate the position and add to the hash    
	$new_variation_id = $last_variation_id + 1;
	$new_allele_id = $last_allele_id +1;
	$new_variation_feature_id = $last_variation_feature_id +1;

	write_file("variation.txt",$new_variation_id,$last_source_id+1,"SOL$new_variation_id"); #write the Variation file
	write_file("allele.txt",$new_allele_id,$new_variation_id,$last_sample_id+1,$call_base,1); #write the Allele file
	write_file("flanking_sequence.txt",$new_variation_id,$slice->get_seq_region_id,1,$start_BAC + $position - 1 - 26, $start_BAC + $position -2,$start_BAC+$position,$start_BAC+$position+25); #write the flanking sequence file
	write_file("variation_feature.txt",$new_variation_feature_id,$new_variation_id,$slice->get_seq_region_id,$start_BAC + $position -1,$start_BAC + $position -1,1,"SOL$new_variation_id","$ref_base/$call_base",1,$last_source_id +1); #write the variation_feature file

	$last_variation_id = $new_variation_id;
	$last_allele_id = $new_allele_id;
	$last_variation_feature_id = $new_variation_feature_id;
	
    }
    rename("$TMP_DIR/variation.txt", "$TMP_DIR/$TMP_FILE");
    load($dbVar->dbc,"variation","variation_id","source_id","name");
    rename("$TMP_DIR/allele.txt", "$TMP_DIR/$TMP_FILE");
    load($dbVar->dbc,"allele","allele_id","variation_id","sample_id","allele","frequency");
    rename("$TMP_DIR/flanking_sequence.txt", "$TMP_DIR/$TMP_FILE");
    load($dbVar->dbc,"flanking_sequence","variation_id","seq_region_id","seq_region_strand","up_seq_region_start","up_seq_region_end","down_seq_region_start","down_seq_region_end");
    rename("$TMP_DIR/variation_feature.txt", "$TMP_DIR/$TMP_FILE");
    load($dbVar->dbc,"variation_feature","variation_feature_id","variation_id","seq_region_id","seq_region_start","seq_region_end","seq_region_strand","variation_name","allele_string","map_weight","source_id");
    close IN;
}


#function to return the last id used in the table tablename (the id must be called "tablename_id")
sub get_last_table_id{
    my $dbSanger = shift;
    my $tablename = shift;

    my $max_id;
    my $sth = $dbSanger->dbc->prepare(qq{SELECT MAX($tablename\_id) from $tablename});
    $sth->execute();
    $sth->bind_columns(\$max_id);
    $sth->fetch;
    $sth->finish();

    return $max_id if (defined $max_id);
    return 0 if (!defined $max_id);
}

sub write_file{
    my $filename = shift;
    my @values = @_;

    open FH, ">>$TMP_DIR/$filename" || die "Could not open file with information: $!\n";
    my @a = map {defined($_) ? $_ : '\N'} @values; #to replace undefined values by \N in the file
    print FH join("\t", @a), "\n";
    close FH || die "Could not close file with information: $!\n";
    
}
