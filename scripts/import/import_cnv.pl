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


#to run it type perl import_cnv.pl data_file.txt
use strict;
use warnings;

use Getopt::Long;
use FindBin qw( $Bin );
use ImportUtils qw(debug load);
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use Bio::EnsEMBL::Registry;
use Data::Dumper;

warn("Make sure you have a updated ensembl.registry file!\n");

my $cnv_file; #file containing CNV data
GetOptions('tmpdir=s'  => \$ImportUtils::TMP_DIR,
	   'tmpfile=s' => \$ImportUtils::TMP_FILE,
	   'cnv_file=s' => \$cnv_file
	   );

my $registry_file ||= $Bin . "/ensembl.registry";
my $species = 'human';
Bio::EnsEMBL::Registry->load_all( $registry_file );

my $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $dbVar = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');

my $TMP_DIR  = $ImportUtils::TMP_DIR;
my $TMP_FILE = $ImportUtils::TMP_FILE;

#parse the file to get the cnv position
#the file is given as an argument
my @line;
my $source_id_ensembl = 1;
my $source_id_cnv = 2;
my $variation_id = 1;
my $slice_adaptor = $dbCore->get_SliceAdaptor();
my ($chr,$start,$end,$seq_region_start,$name,$seq_region_id);
open IN, "<$cnv_file" || die "Could not open file with CNV data: $!\n";
write_file("$TMP_DIR/source.txt",$source_id_ensembl,'ENSEMBL');
write_file("$TMP_DIR/source.txt",$source_id_cnv,'Database of genomics Variants'); #create the source file
<IN>;# don't read the first line, the header
while (<IN>){
    chomp;
    @line = split/\s+/; 
    $name = $line[0]; #get variation name from line
    $chr = substr($line[2],3); #get chromsome name
    $start = $line[3];
    $end = $line[4];
    $seq_region_id = get_seq_region_id($slice_adaptor,$chr);
    write_file("$TMP_DIR/variation.txt",$variation_id,$source_id_ensembl,'ENSCNV' . $variation_id);
    write_file("$TMP_DIR/variation_synonym.txt",$variation_id,$variation_id,$source_id_cnv,$name);
    write_file("$TMP_DIR/variation_feature.txt",$variation_id,$seq_region_id,$start,$end,1,$variation_id,'CNV','ENSCNV' . $variation_id,1,$source_id_ensembl);
    write_file("$TMP_DIR/flanking_sequence.txt",$variation_id,$start - 100,$start - 1,$end + 1, $end + 100, $seq_region_id,1);
    $variation_id++
}

close IN;

#time to import the data in the database
rename "$TMP_DIR/source.txt","$TMP_DIR/$TMP_FILE";
load($dbVar->dbc,"source","source_id","name");
unlink "$TMP_DIR/source.txt";
rename "$TMP_DIR/variation.txt","$TMP_DIR/$TMP_FILE";
load($dbVar->dbc,"variation","variation_id","source_id","name");
unlink "$TMP_DIR/variation.txt";
rename "$TMP_DIR/variation_synonym.txt","$TMP_DIR/$TMP_FILE";
load($dbVar->dbc,"variation_synonym","variation_synonym_id","variation_id","source_id","name");
unlink "$TMP_DIR/variation_synonym.txt";
rename "$TMP_DIR/variation_feature.txt","$TMP_DIR/$TMP_FILE";
load($dbVar->dbc,"variation_feature","variation_feature_id","seq_region_id","seq_region_start","seq_region_end","seq_region_strand","variation_id","allele_string","variation_name","map_weight","source_id");
unlink "$TMP_DIR/variation_feature.txt";
rename "$TMP_DIR/flanking_sequence.txt","$TMP_DIR/$TMP_FILE";
load($dbVar->dbc,"flanking_sequence","variation_id","up_seq_region_start","up_seq_region_end","down_seq_region_start","down_seq_region_end","seq_region_id","seq_region_strand");
unlink "$TMP_DIR/flanking_sequence.txt";

#get the internal seq_region_id from the core database
sub get_seq_region_id{
    my $slice_adaptor = shift;
    my $chr = shift;

    my $slice = $slice_adaptor->fetch_by_region('chromosome',$chr);

    return $slice->get_seq_region_id();
}

#method to write in a file a list of values
sub write_file{
    my $table_name = shift;
    my @values = @_;

    open FH, ">>$table_name" || die "Could not open file with information: $!\n";
    my @a = map {defined($_) ? $_ : '\N'} @values; #to replace undefined values by \N in the file
    print FH join("\t", @a), "\n";
    close FH || die "Could not close file with information: $!\n";
    
}
