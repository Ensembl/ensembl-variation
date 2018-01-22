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

use Getopt::Long;

use ImportUtils qw(debug);
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(SO_variation_class);
use JSON;

our ($species, $input_file, $source_name, $version, $registry_file, $debug);

GetOptions('species=s'         => \$species,
           'source_name=s'     => \$source_name,
           'input_file=s'      => \$input_file,
           'version=i'         => \$version,
           'registry=s'        => \$registry_file,
           'debug!'            => \$debug,
          );
$registry_file ||= "./ensembl.registry";
$source_name   ||= 'DECIPHER';

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_all($registry_file);

my $sa = $registry->get_adaptor('human', 'core', 'slice');

my $dbh = $registry->get_adaptor('human', 'variation_private', 'variation')->dbc->db_handle;


my $source_id = source();

my $type = 'Variation';
my $short_variant_type = 'SNV';

# SO terms
my $attr_type_sth = $dbh->prepare(qq{ SELECT attrib_type_id FROM attrib_type WHERE code='SO_term'});
$attr_type_sth->execute;
my $attr_type_id = ($attr_type_sth->fetchrow_array)[0];
$attr_type_sth->finish;
die "No attribute type found for the entry 'SO_term'" if (!defined($attr_type_id));
my $SO_terms = get_attrib_list($attr_type_id);

# Clinical significance
my $cs_attr_type_sth = $dbh->prepare(qq{ SELECT attrib_type_id FROM attrib_type WHERE code='clinvar_clin_sig'});
$cs_attr_type_sth->execute;
my $cs_attr_type_id = ($cs_attr_type_sth->fetchrow_array)[0];
$cs_attr_type_sth->finish;
die "No attribute type found for the entry 'clinvar_clin_sig'" if (!defined($cs_attr_type_id));
my $clin_sign_types = get_attrib_list($cs_attr_type_id);

# Inheritance type
my $inheritance_attrib_type = 'inheritance_type';
my $stmt = qq{ SELECT attrib_type_id FROM attrib_type WHERE code='$inheritance_attrib_type'};
my $inheritance_attrib_type_id = ($dbh->selectall_arrayref($stmt))->[0][0];

my $find_existing_var_sth = $dbh->prepare(qq{
    SELECT variation_id FROM variation WHERE name = ?
});

my $add_var_sth = $dbh->prepare(qq{
    INSERT INTO variation (source_id, name, flipped, class_attrib_id, clinical_significance) VALUES (?,?,?,?,?)
});

my $find_existing_vf_sth = $dbh->prepare(qq{
    SELECT variation_feature_id FROM variation_feature WHERE variation_name = ? AND seq_region_id = ? 
    AND seq_region_start = ? AND seq_region_end = ?
});

my $add_vf_sth = $dbh->prepare(qq{
    INSERT IGNORE INTO variation_feature (variation_id, seq_region_id, seq_region_start,
        seq_region_end, seq_region_strand, variation_name, allele_string, map_weight, 
        source_id, class_attrib_id, clinical_significance)
    VALUES (?,?,?,?,?,?,?,?,?,?,?)
});

my $select_phe_sth = $dbh->prepare(qq{
    SELECT phenotype_id FROM phenotype WHERE description = ?
});

my $add_phe_sth = $dbh->prepare(qq{
    INSERT IGNORE INTO phenotype (description)
    VALUES (?)
});

my $find_existing_pf_sth = $dbh->prepare(qq{
    SELECT phenotype_feature_id FROM phenotype_feature WHERE object_id = ? AND seq_region_id = ? 
    AND seq_region_start = ? AND seq_region_end = ? AND phenotype_id = ?
});

my $add_pf_sth = $dbh->prepare(qq{
    INSERT IGNORE INTO phenotype_feature (seq_region_id, seq_region_start,
        seq_region_end, seq_region_strand, object_id, type, 
        source_id, phenotype_id)
    VALUES (?,?,?,?,?,?,?,?)
});

my $add_pfa_sth = $dbh->prepare(qq{
    INSERT IGNORE INTO phenotype_feature_attrib (phenotype_feature_id, attrib_type_id, value)
    SELECT phenotype_feature_id, ?, ? FROM phenotype_feature WHERE source_id=? AND object_id=? AND type=? AND phenotype_id=?
});

my $add_failed_sth = $dbh->prepare(qq{
    INSERT IGNORE INTO failed_variation (variation_id, failed_description_id) VALUES (?,?)
});


debug(localtime()." Parse file and insert data into the database");

open my $INPUT, "<$input_file" or die "Can't open '$input_file'";  

MAIN_LOOP : while(<$INPUT>) {

  next if /^track/;
  chomp $_;
    
  my $data = parse_line($_);
    
  next if !$data;
  
  my $slice = $data->{slice};
    
  my $class_attrib_id = $SO_terms->{$data->{SO_term}};
  my $seq_region_id = $sa->get_seq_region_id($slice);
  my $clin_sign = $data->{clin_sign} if ($clin_sign_types->{$data->{clin_sign}});
    
    
  # Import Variation
  $find_existing_var_sth->execute($data->{ID});

  my $existing_var = $find_existing_var_sth->fetchrow_arrayref;

  my $variation_id;
        
  if ($existing_var) {
    $variation_id = $existing_var->[0];
  }
  else {
    $add_var_sth->execute($source_id, $data->{ID}, $data->{flipped}, $class_attrib_id, $clin_sign);

    $variation_id = $dbh->last_insert_id(undef, undef, undef, undef);
  }
    
    
  if (!$data->{failed}) {
    
    # Import Variation Feature
    $find_existing_vf_sth->execute($data->{ID},$seq_region_id, $data->{start}, $data->{end});

    my $existing_vf = $find_existing_vf_sth->fetchrow_arrayref;
      
    if (!$existing_vf) {
      $add_vf_sth->execute(
        $variation_id,
        $seq_region_id,
        $data->{start},
        $data->{end},
        $data->{strand},
        $data->{ID},
        $data->{ref}.'/'.$data->{alt},
        1,
        $source_id,
        $class_attrib_id,
        $clin_sign
      );
    }
      
      
    if ($data->{phenotype} && $data->{phenotype} ne '') {
        
      my @phenotypes = split(/\|/,$data->{phenotype});
      foreach my $phe (keys(%{$data->{phenotype}})) {

        # Import Phenotype
        $select_phe_sth->execute($phe);
        my $phe_id = ($select_phe_sth->fetchrow_array)[0];
        if (!$phe_id) {
          $add_phe_sth->execute($phe);
          $phe_id = $dbh->last_insert_id(undef, undef, undef, undef);
        }

        # Import Phenotype Feature
        $find_existing_pf_sth->execute($data->{ID},$seq_region_id, $data->{start}, $data->{end}, $phe_id);

        my $existing_pf = $find_existing_pf_sth->fetchrow_arrayref;
        
        if (!$existing_pf) {
          $add_pf_sth->execute(
            $seq_region_id,
            $data->{start},
            $data->{end},
            $data->{strand},
            $data->{ID},
            $type,
            $source_id,
            $phe_id
          );

          # Import Phenotype Feature Attrib
          if ($data->{inheritance} && $data->{inheritance} ne '') {

            $add_pfa_sth->execute(
              $inheritance_attrib_type_id,
              $data->{inheritance},
              $source_id,
              $data->{ID},
              $type,
              $phe_id,
            );
          }
        }
      }
    }
  }
  else {
    $add_failed_sth->execute($variation_id,$data->{failed});
  }
}

sub get_attrib_list {
  my $type_id = shift;
  my %attrib_list;
  my $attr_list_sth = $dbh->prepare(qq{ SELECT attrib_id,value FROM attrib WHERE attrib_type_id=$type_id});
  $attr_list_sth->execute;
  while(my @row = $attr_list_sth->fetchrow_array) {
   $attrib_list{$row[1]} = $row[0];
  }
  return \%attrib_list;
}


sub source{
  debug(localtime()." Inserting into source table");
  
  my $name = $source_name;
  my $url  = 'http://decipher.sanger.ac.uk';
  my $desc = 'Database of Chromosomal Imbalance and Phenotype in Humans Using Ensembl Resources';
  # Check if the DGVa source already exists, else it create the entry
  if ($dbh->selectrow_arrayref(qq{SELECT source_id FROM source WHERE name='$name';})) {
    $dbh->do(qq{UPDATE IGNORE source SET description='$desc',url='$url',version=$version where name='$name';});
  }
  else {
    $dbh->do(qq{INSERT INTO source (name,description,url,version) VALUES ('$name','$desc','$url',$version);});
  }
  my @source_id = @{$dbh->selectrow_arrayref(qq{SELECT source_id FROM source WHERE name='$name';})};
  return $source_id[0];
}

sub parse_line {
  my $line = shift;
  my @data = split(/\t/, $line);
  my $info;
  
  # Extract data
  $info->{chr}         = $data[0];
  $info->{start}       = $data[1]+1;
  $info->{end}         = $data[2];
  $info->{subject}     = $data[3];
  $info->{strand}      = ($data[5] eq '+') ? 1 : -1;
  
  $info->{chr} =~ s/chr//i;
  $info->{chr} = 'MT' if ($info->{chr} eq 'M');
  
  my $json_text = decode_json($data[13]);
  
  return if ($json_text->{'variant_type'} ne $short_variant_type);
  
  $info->{ref}         = $json_text->{ref_allele};
  $info->{alt}         = $json_text->{alt_allele};
  $info->{inheritance} = $json_text->{inheritance};
  $info->{clin_sign}   = lc($json_text->{pathogenicity}) if ($json_text->{pathogenicity});
  foreach my $phe (@{$json_text->{phenotypes}}) {
    my $p_name = $phe->{name};
    $info->{phenotype}{$p_name} = 1;
  }
  $info->{ID}          = "DEC".$info->{subject}."_".$info->{chr}."_".$info->{start}."_".$info->{end};
  $info->{flipped}     = 0;

  
  # Strand / alleles
  my $slice = $sa->fetch_by_region('chromosome', $info->{chr},$info->{start},$info->{end},1);
  my $seq   = $slice->seq;
  my $ref   = $info->{ref};
  my $alt   = $info->{alt};
  
  if ($seq ne $ref) {
    if ($info->{strand} == 1) {
      $info->{failed} = 15; # Mapped position is not compatible with reported alleles
      print STDERR "Variant ".$info->{ID}.": Mapped position is not compatible with reported alleles\n";
    }
    else {
      reverse_comp(\$ref);
      if ($ref eq $seq) {
        reverse_comp(\$alt);
        $info->{ref} = $ref;
        $info->{alt} = $alt;
        $info->{flipped} = 1;
      }
      else {
        $info->{failed} = 2; # None of the variant alleles match the reference allele
        print STDERR "Variant ".$info->{ID}.": None of the variant alleles match the reference allele\n";
      }
    }
  }
  $info->{slice} = $slice;

  # SO term
  if ($info->{start} == $info->{end}) {
    if (length($info->{ref}) == length($info->{alt})) {
      $info->{SO_term} = 'SNV';
    }
    else {
      $info->{SO_term} = 'insertion';
      $info->{start} ++;
    }
  }
  else {
    $info->{SO_term} = (length($info->{ref}) > length($info->{alt})) ? 'deletion' : 'indel';
  }

  return $info;
}


