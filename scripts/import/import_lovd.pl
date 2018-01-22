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

# Import LOVD data into an Ensembl variation schema database

use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(SO_variation_class);

my $USE_DB = 1;
my $ADD_SETS = 1;
my $VERBOSE = 1;

my $import_file;
my $registry_file;
my $version;
my $help;

GetOptions(
    "import|i=s"    => \$import_file,
    "registry|r=s"  => \$registry_file,
    "verbose|v"     => \$VERBOSE,
    "version=s"     => \$version,
    "test|t"        => \$USE_DB,
    "help|h"        => \$help,
);

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_all($registry_file);

my $sa = $registry->get_adaptor(
    'human', 'core', 'slice'
);

my $dbh = $registry->get_adaptor(
    'human', 'variation_private', 'variation'
)->dbc->db_handle;

open my $INPUT, "<$import_file" or die "Can't open '$import_file'";    
    
# check to see if we already have the LOVD source

my $src_sth = $dbh->prepare(qq{
    SELECT  source_id
    FROM    source
    WHERE   name LIKE "LOVD%"
});

$src_sth->execute;

my $existing_src = $src_sth->fetchrow_arrayref;

my $source_id;
my $not_flipped_col = 'not_flipped_allele';

if ($existing_src) {
    $source_id = $existing_src->[0];
    
    print "Found existing source_id: $source_id\n";
    my $sth = $dbh->prepare(qq{  UPDATE source SET version=? WHERE source_id=? });
    $sth->execute($version,$source_id);
}
else {
    # if not, add it
    my $sth = $dbh->prepare(qq{
        INSERT INTO source (name, description, url, somatic_status) 
        VALUES (
            'LOVD', 
            'Leiden Open (source) Variation Database', 
            'http://www.lovd.nl',
            'germline'
        );
    });

    $sth->execute;

    $source_id = $dbh->last_insert_id(undef, undef, undef, undef);

    print "New source_id: $source_id\n";
}



my $attr_type_sth = $dbh->prepare(qq{ SELECT attrib_type_id FROM attrib_type WHERE code='SO_term'});
$attr_type_sth->execute;
my $attr_type_id = ($attr_type_sth->fetchrow_array)[0];
$attr_type_sth->finish;
die "No attribute type found for the entry 'SO_term'" if (!defined($attr_type_id));

my $SO_terms = get_SO_term_list($attr_type_id);

my $find_existing_var_sth = $dbh->prepare(qq{
    SELECT variation_id FROM variation WHERE name = ?
});

my $add_var_sth = $dbh->prepare(qq{
    INSERT IGNORE INTO variation (source_id, name, flipped, class_attrib_id) VALUES (?,?,?,?)
});

my $add_vf_sth = $dbh->prepare(qq{
    INSERT IGNORE INTO variation_feature (variation_id, seq_region_id, seq_region_start,
        seq_region_end, seq_region_strand, variation_name, allele_string, map_weight, 
        source_id, class_attrib_id, $not_flipped_col)
    VALUES (?,?,?,?,?,?,?,?,?,?,?)
});


my $find_existing_svar_sth = $dbh->prepare(qq{
    SELECT structural_variation_id FROM structural_variation WHERE variation_name = ?
});

my $add_svar_sth = $dbh->prepare(qq{
    INSERT INTO structural_variation (source_id, variation_name, class_attrib_id) VALUES (?,?,?)
});

my $add_svf_sth = $dbh->prepare(qq{
    INSERT IGNORE INTO structural_variation_feature (structural_variation_id, seq_region_id, seq_region_start,
        seq_region_end, seq_region_strand, variation_name, allele_string, source_id, class_attrib_id)
    VALUES (?,?,?,?,?,?,?,?,?)
});


pre_process_variation_feature();

# loop over the input file
my $sv_prefix = 'structural_';

MAIN_LOOP : while(<$INPUT>) {

    next unless /^chr/;
  
    my $line = $_;
    
    my (
        $chr,
        $start, 
        $stop,
        $accession,
        $url_param,
        $description
    ) = split "\t", $_;

    $chr = (split('r',$chr))[1];
    $accession =~ s/ //g;
      
    $start ++; # Conversion from BED file
  
    my ($type, $al_ref, $al_alt);
    my $strand = 1;
    my $failed_mapping = 0;
    
    # Choose between Variation and Structural variation
    my $prefix = '';
    
    ## COMPLEX
    if ($accession =~ /;/i) {
      $al_alt = 'SEQUENCE ALTERATION';
      $type = 'sequence_alteration';
    }
    ## SUBSTITUTION
    elsif ($accession =~ /([ATGC]+)\>([ATGC]++)\]?\)?\??$/i) {
      $al_ref = uc($1);
      $al_alt = uc($2);
      $type = 'SNV';
    }
    ## INDEL
    elsif ($accession =~ /delins([ATGC]+)\]?\)?\??$/) {
      $al_ref = 'DELETION';
      $al_alt = $1;
      $type = 'indel';
    }
    elsif ($accession =~ /delins(\d+)\]?\)?\??$/) {
      $al_alt = "$1 BP INSERTION";
      $type = 'indel';
    }
    elsif ($accession =~ /del([ATGC]+)ins([ATGC]+)\]?\)?\??$/) {
      $al_ref = $1;
      $al_alt = $2;
      $type = 'indel';
    }
    elsif ($accession =~ /del(\d+)ins(\d+)\]?\)?\??$/) {
      $al_ref = "$1 BP DELETION";
      $al_alt = "$2 BP INSERTION";
      $type = 'indel';
    }
    elsif ($accession =~ /del([ATGC]+)ins(\d+)\]?\)?\??$/) {
      $al_ref = "$1";
      $al_alt = "$2 BP INSERTION";
      $type = 'indel';
    }
    ## DELETION
    elsif ($accession =~ /del([ATGC]+)\]?\)?\??$/) {
      $al_ref = $1;
      $al_alt = '-';
      $type = 'deletion';
    }
    elsif ($accession =~ /del(\d+)\]?\)?\??$/) {
      $al_alt = "$1 BP DELETION";
      $type = 'deletion';
    }
    elsif ($accession =~ /del\]?\)?\??$/) {
      $al_alt = "DELETION";
      $type = 'deletion';
    }
    ## DUPLICATION
    elsif ($accession =~ /dup([ATGC]+)\]?\)?\??$/) {
      $al_ref = $1;
      $al_alt = "$al_ref$al_ref";
      $type = 'duplication';
    }
    elsif ($accession =~ /dup(\d+)\]?\)?\??$/) {
      $al_alt = "$1 BP DUPLICATION";
      $type = 'duplication';
    }
    elsif ($accession =~ /dup\]?\)?\??$/) {
      $al_alt = 'DUPLICATION';
      $type = 'duplication';
    }
    elsif ($accession =~ /\[(\d+)\]?\)?\??$/) {
      my $count = $1 - 1;
      $al_alt = "DUPLICATION $count times";
      $type = 'duplication';
    }
    ## INVERSION
    elsif ($accession =~ /inv[ATGC]*\]?\)?\??$/) {
      $al_alt = "INVERSION";
      $type = 'inversion';
    }
    ## INSERTION
    elsif ($accession =~ /ins([ATGC]+)\]?\)?\??$/) {
      $al_ref = '-';
      $al_alt = $1;
      $type = 'insertion';
      if ($stop != $start+1) {
        $failed_mapping = 2;
      }
      else {
        my $tmp_start = $start;
        $start = $stop;
        $stop  = $tmp_start;
      }
    }
    elsif ($accession =~ /ins(\d+)\]?\)?\??$/) {
      $al_ref = '-';
      $al_alt = "$1 BP INSERTION";
      $type = 'insertion';
      if ($stop != $start+1) {
        $failed_mapping = 2;
      }
      else {
        my $tmp_start = $start;
        $start = $stop;
        $stop  = $tmp_start;
      }
    }
    else {
      print STDERR "Type for the variant '$accession' not found\n";
      next;
    }
    
    # Choose between Variation and Structural variation
    
    ####################################
    # Add tandem_duplication > 50 ???? #
    ####################################
    if ($type !~ /insertion/i && $stop-$start > 50) {
      $prefix = $sv_prefix;
    }

    my $slice = $sa->fetch_by_region('chromosome', $chr, $start, $stop);

    if ($slice) {
      my $seq = $slice->seq;
      
      if (defined($al_ref) && $al_ref =~ /^[ATGC]+$/) {
        if ($seq ne $al_ref) {
          reverse_comp(\$seq);
          if ($seq ne $al_ref) {
            $failed_mapping = 1;
          } 
          else {
            $strand = -1;
          }
        }
        else {
          $strand = 1;
        }
      }
      elsif (defined($al_alt) && $al_alt =~ /^[ATGC]+$/ && $type ne 'insertion' && $type ne 'indel' && $type ne 'sequence_alteration') {
        if ($seq ne $al_alt) {
          reverse_comp(\$seq);
          if ($seq ne $al_alt) {
            $failed_mapping = 1;
          } 
          else {
            $strand = -1;
          }
        }
        else {
          $strand = 1;
        }  
      }
    }
    else {
      print STDERR "Coordinates of the variant '$accession' can't be mapped to the genome ($chr:$start-$stop)\n";
      next;
    }
    
    
    my $allele_string = (defined($al_ref)) ? "$al_ref/$al_alt" : $al_alt;
    
    if ($failed_mapping != 0) {
      if ($failed_mapping == 1) {
        print STDERR "Alleles of the variant '$accession' can't be mapped to the genome ($allele_string) => ($chr:$start-$stop)\n";
      }
      elsif ($failed_mapping == 2) {
        print STDERR "Coordinates of the variant '$accession' ($chr:$start-$stop) don't correspond to the type of variant '$type'\n";
      }
      next;
    }
    
    my $flipped = 0;
    my $not_flipped_allele;
    if ($strand == -1 and $prefix eq '') {
      if ($al_ref =~ /^[ATGC]+$/) {
        reverse_comp(\$al_ref);
        $flipped = 1;
      }
      if ($al_alt =~ /^[ATGC]+$/) {
        reverse_comp(\$al_alt);
        $flipped = 1;
      }
      if ($flipped == 1) {
        $strand = 1;
        $not_flipped_allele = $allele_string;
      }
      $allele_string = (defined($al_ref)) ? "$al_ref/$al_alt" : $al_alt;
    }
    
    
    my $class_attrib_id = $SO_terms->{$type};
    my $seq_region_id = $sa->get_seq_region_id($slice);
    
    if ($prefix eq '') {
      # Import the Variation
        $find_existing_var_sth->execute($accession);

        my $existing_var = $find_existing_var_sth->fetchrow_arrayref;

        my $variation_id;
        
        if ($existing_var) {
          $variation_id = $existing_var->[0];
        }
        else {
          $add_var_sth->execute($source_id, $accession, $flipped, $class_attrib_id);

          $variation_id = $dbh->last_insert_id(undef, undef, undef, undef);
        }
      # Import the Variation Feature
      $add_vf_sth->execute(
         $variation_id,
         $seq_region_id,
         $start,
         $stop,
         $strand,
         $accession,
         $allele_string,
         1,
         $source_id,
         $class_attrib_id,
         $not_flipped_allele
      );
    }
    else {
      # Import the Structural Variation
      $find_existing_svar_sth->execute($accession);

      my $existing_svar = $find_existing_svar_sth->fetchrow_arrayref;

      my $structural_variation_id;
        
      if ($existing_svar) {
        $structural_variation_id = $existing_svar->[0];
      }
      else {
        $add_svar_sth->execute($source_id, $accession, $class_attrib_id);
       
        $structural_variation_id = $dbh->last_insert_id(undef, undef, undef, undef);
      }
      
      # Import the Structural Variation Feature
      $add_svf_sth->execute(
         $structural_variation_id,
         $seq_region_id,
         $start,
         $stop,
         $strand,
         $accession,
         $allele_string,
         $source_id,
         $class_attrib_id
      );
    }
#    print "$chr\t$start\t$stop\t$strand\t$accession\t$allele_string\n";
}

# Post processing
post_process_variation_feature();


sub get_SO_term_list {
  my $type_id = shift;
  my %SO_list;
  my $attr_list_sth = $dbh->prepare(qq{ SELECT attrib_id,value FROM attrib WHERE attrib_type_id=$type_id});
  $attr_list_sth->execute;
  while(my @row = $attr_list_sth->fetchrow_array) {
   $SO_list{$row[1]} = $row[0];
  }
  return \%SO_list;
}

sub pre_process_variation_feature {
  # Add tmp column for allele_string before flipping (variation_feature)
  if ($dbh->do(qq{show columns from variation_feature like '$not_flipped_col';}) != 1){
    $dbh->do(qq{ALTER TABLE variation_feature ADD COLUMN $not_flipped_col varchar(255) DEFAULT NULL;});
  }  
  
  # Prepare the variation_feature table
  if ($dbh->do(qq{SHOW KEYS FROM variation_feature WHERE Key_name='name_key';}) < 1){
    $dbh->do(qq{ALTER TABLE variation_feature ADD CONSTRAINT UNIQUE KEY `name_key` (`variation_name`,`seq_region_id`,`seq_region_start`)});
  }
  # Prepare the structural_variation_feature table
  if ($dbh->do(qq{SHOW KEYS FROM structural_variation_feature WHERE Key_name='name_key';}) < 1){
    $dbh->do(qq{ALTER TABLE structural_variation_feature ADD CONSTRAINT UNIQUE KEY `name_key` (`variation_name`,`seq_region_id`,`seq_region_start`)});
  }
}  
  

sub post_process_variation_feature {

  $dbh->do(qq{ALTER TABLE variation_feature DROP INDEX `name_key`});
  $dbh->do(qq{ALTER TABLE structural_variation_feature DROP INDEX `name_key`});

  # Populate the map_weight column in variation_feature;
  $dbh->do(qq{DROP TABLE IF EXISTS tmp_vf_map_weight});
  $dbh->do(qq{
    CREATE table tmp_vf_map_weight (
      variation_id int(10), vf_count int(10),
      PRIMARY KEY (variation_id)
    );
  });
  $dbh->do(qq{
    INSERT INTO tmp_vf_map_weight (variation_id, vf_count) 
    SELECT variation_id, count(variation_feature_id) FROM variation_feature GROUP BY variation_id
  });
  $dbh->do(qq{
    UPDATE tmp_vf_map_weight w, variation_feature f SET f.map_weight=w.vf_count WHERE w.variation_id=f.variation_id AND w.vf_count>1;
  });
  $dbh->do(qq{ DROP table tmp_vf_map_weight});
  
  $dbh->do(qq{
    UPDATE variation_feature SET allele_string=$not_flipped_col,seq_region_strand=-1 WHERE $not_flipped_col is not null and map_weight>1;
  });
  
  #$dbh->do(qq{ALTER TABLE variation_feature DROP COLUMN $not_flipped_col});
}


