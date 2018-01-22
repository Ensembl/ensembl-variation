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

# this script read mapping data ALL_DIFF into variation_feature table

use strict;
use warnings;
use DBI;
use DBH;
use Getopt::Long;
use Benchmark;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBConnection;
use ImportUtils qw(dumpSQL debug create_and_load load );
use FindBin qw( $Bin );

my ($TAX_ID, $LIMIT_SQL, $CONTIG_SQL, $TMP_DIR, $TMP_FILE);

my $dbSNP;
my $dbVar;
my $dbCore;

{
  my ($species,$mapping_file,$limit);
  
  GetOptions('species=s' => \$species,
	     'mapping_file=s' => \$mapping_file,
             'tmpdir=s'  => \$ImportUtils::TMP_DIR,
             'tmpfile=s' => \$ImportUtils::TMP_FILE,
             'limit=i'   => \$limit);
  
  usage('-mapping_file argument is required.') if(!$mapping_file);

  $TMP_DIR  = $ImportUtils::TMP_DIR;
  $TMP_FILE = $ImportUtils::TMP_FILE;


  $LIMIT_SQL = ($limit) ? " LIMIT $limit " : '';

  my $registry_file;
  $registry_file ||= $Bin . "/ensembl.registry";
  Bio::EnsEMBL::Registry->load_all( $registry_file );

  my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
  my $vdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
  $dbCore = $cdb->dbc;
  $dbVar = $vdb->dbc;


  variation_feature($mapping_file);
}

sub variation_feature {
  
  my $mapping_file = shift;

  my (%rec, %rec_pos, %rec_line, %rec_seq_region);


  debug("Dumping Variation");
  
  my $sth = $dbVar->prepare (qq{SELECT v.variation_id, v.name, v.source_id, f.allele_string, v.validation_status 
	                              FROM variation v LEFT JOIN variation_feature f ON v.variation_id=f.variation_id;
                             });									 
  $sth->execute();

  while(my ($variation_id, $name, $source_id, $allele_string, $validation_status) = $sth->fetchrow_array()) {
		next if ($rec{$name});
		
		my $h={};
    $h->{name} = $name;
    $h->{variation_id} = $variation_id;
    $h->{source_id} = $source_id;
    $h->{allele_string} = defined $allele_string ?  $allele_string : '\N';
    $h->{validation} = defined $validation_status ? $validation_status : '\N';
    $rec{$name} = $h;
  }

  $sth->finish();

  debug("Reading Mapping file...");

  open (IN, "$mapping_file");
  open ( FH, ">$TMP_DIR/$TMP_FILE" );

  while (<IN>) {
    chomp;
    next if /^more |^PARSING/;
    s/^MORE_HITS\s+//;
    my ($ref_id, $slice_name, $start, $end, $strand, $ratio) =split;
    push @{$rec_line{$ref_id}}, $_;
  }
  
  foreach my $key (keys %rec_line) {
    #next if @{$rec_line{$key}} >3;
    foreach my $line (@{$rec_line{$key}}) {
      my ($ref_id, $slice_name, $start, $end, $strand, $ratio) =split /\s+/,$line;
      next if $ratio <0.7;

      $strand = ($strand eq "+") ? 1 : -1;
      my ($coord_sys,$assembly,$seq_region_name,$seq_region_start,$seq_region_end,$version);

      ($seq_region_name,$seq_region_start,$seq_region_end) = split /\-/, $slice_name
	if $slice_name =~ /\-/;

      ($coord_sys,$assembly,$seq_region_name,$seq_region_start,$seq_region_end,$version) = split /\:/, $slice_name
	if $slice_name =~ /\:/;

      my $sth;

      if (!$rec_seq_region{$seq_region_name}) {
	$sth = $dbCore->prepare (qq{SELECT seq_region_id from seq_region where name = ?});
	$sth->execute("$seq_region_name");
     
	my $seq_region_id = $sth->fetchrow_array();
	$rec_seq_region{$seq_region_name}=$seq_region_id;
      
	if (!$seq_region_id) {
	  warn "There is no seq_region_id for $ref_id\n";
	  next;
	}
      }
      my $seq_region_id = $rec_seq_region{$seq_region_name};
      my $new_seq_region_start = $seq_region_start + $start -1 if ($seq_region_start);
      my $new_seq_region_end = $seq_region_start + $end -1 if ($seq_region_start);
    
      if (!$rec_pos{$ref_id}{$seq_region_id}{$new_seq_region_start}{$new_seq_region_end}) {
				print FH join ("\t", $seq_region_id,$new_seq_region_start,$new_seq_region_end,$strand,$rec{$ref_id}->{variation_id},$ref_id,$rec{$ref_id}->{allele_string},$rec{$ref_id}->{source_id},$rec{$ref_id}->{validation})."\n";
				$rec_pos{$ref_id}{$seq_region_id}{$new_seq_region_start}{$new_seq_region_end}=1;
      }
    }
  }
  
  $sth->finish();
  
  close IN;
  close FH;


  $dbVar->do(qq{CREATE TABLE IF NOT EXISTS variation_feature_mapping_MT LIKE variation_feature});
  load($dbVar, "variation_feature_mapping_MT","seq_region_id","seq_region_start","seq_region_end",
               "seq_region_strand","variation_id","variation_name","allele_string","source_id", "validation_status");
}

###
###example of usage###
###./load_mapping2vf.pl -cdbname danio_rerio_core_24_4 -vdbname yuan_zfish_snp -alldiff /ecs2/scratch4/yuan/zfish/ds_chNotOndir/ALL_DIFF -vpass xxxx
