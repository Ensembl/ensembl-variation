=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation;

use strict;
use warnings;

use DBI qw(:sql_types);
use String::Approx qw(amatch adist);
use Algorithm::Diff qw(diff);

#TODO: check that all libraries from import_script are imported. Any missing? Were they actually used?
use Data::Dumper;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');
#use base qw(Exporter);

#our @EXPORT_OK = qw($core_dba $variation_dba $core_dba);

my %special_characters = (
  'Å' => 'A',
  'Ä' => 'A',
  'Ö' => 'O',
  'Ü' => 'U',
  'ö' => 'o',
  'ü' => 'u',
  'ä' => 'a',
  'í' => 'i',
  'é' => 'e',
  'è' => 'e',
  'ë' => 'e',
  'ç' => 'c',
  '<' => ' less than ',
  '>' => ' more than ',
  'Đ' => 'D',
);

my $pubmed_prefix = 'PMID:';

my $prev_prog;

our $debug; #TODO: difference between our and my in this context?
#our $core_dba;
#our $variation_dba; #TODO: info: this was $db_adaptor in previous script
#our $phenotype_dba;
our $skip_synonyms = 0; #TODO: Q does this have to be a parameter to the entrire pipeline? script has it on for impc, nhgri
our $skip_phenotypes = 0; #TODO; as above, script has it on for uniprot
our $skip_sets = 0; #TODO: this is never set in the script BUT is input param to script

our %phenotype_cache;

sub set_skip_synonyms {
  $skip_synonyms = shift;
}

sub get_special_characters {
  return \%special_characters;
}

sub add_phenotypes {
  my $data = shift;
  my $source_info = shift;
  my $db_adaptor = shift;

  my $phenotypes = $data->{phenotypes};
  my $coords = $data->{coords};
  my @attrib_types = @{$data->{attrib_types}};
  my $variation_ids = $data->{variation_ids};

  my $object_type = $source_info->{object_type};
  my $source_id = $source_info->{source_id};
  my $threshold = $source_info->{threshold};

  # Prepared statements
  my $st_ins_stmt = qq{
    INSERT INTO
      study ( source_id, external_reference, study_type, description
    )
    VALUES (
      $source_id, ?, ?, ?
    )
  };

  my $extra_cond = '';
  my $left_join = '';

  # add special joins for checking certain sources
  if ($source_info->{source} =~ m/uniprot/i) {
    $left_join = qq{
      LEFT JOIN
      (
        phenotype_feature_attrib pfa
        JOIN attrib_type at
        ON pfa.attrib_type_id = at.attrib_type_id
      )
      ON pf.phenotype_feature_id = pfa.phenotype_feature_id
    };

    $extra_cond = 'AND at.code = "variation_names" and pfa.value = ? ';
  }
  elsif ($source_info->{source} =~ m/omim/i) {
    $left_join = qq{
      LEFT JOIN
      (
        phenotype_feature_attrib pfa1
        JOIN attrib_type at1
        ON pfa1.attrib_type_id = at1.attrib_type_id
      )
      ON pf.phenotype_feature_id = pfa1.phenotype_feature_id

      LEFT JOIN
      (
        phenotype_feature_attrib pfa2
        JOIN attrib_type at2
        ON pfa2.attrib_type_id = at2.attrib_type_id
      )
      ON pf.phenotype_feature_id = pfa2.phenotype_feature_id
    };

    $extra_cond = qq{
      AND at1.code = "risk_allele"
      AND pfa1.value = ?
      AND at2.code = "associated_gene"
      AND pfa2.value = ?
    };
  }
  elsif ($source_info->{source} =~ m/nhgri/i) {
    $left_join = qq{
      LEFT JOIN
      (
        phenotype_feature_attrib pfa
        JOIN attrib_type at
        ON pfa.attrib_type_id = at.attrib_type_id
      )
      ON pf.phenotype_feature_id = pfa.phenotype_feature_id
    };

    $extra_cond = qq{ AND at.code = "p_value" AND pfa.value = ? };
  }

  my $pf_check_stmt = qq{
    SELECT
      pf.phenotype_feature_id
    FROM
      phenotype_feature pf
    $left_join
    WHERE
      pf.object_id = ? AND
      pf.type = ? AND
      pf.phenotype_id = ? AND
      pf.source_id = ? AND
      (pf.study_id = ? OR pf.study_id IS NULL)
      $extra_cond
    LIMIT 1
  };

  my $pf_ins_stmt = qq{
    INSERT INTO phenotype_feature (
      phenotype_id, source_id, study_id,
      type, object_id, is_significant,
      seq_region_id, seq_region_start, 
      seq_region_end, seq_region_strand
    )
    VALUES (
      ?, ?, ?,
      ?, ?, ?,
      ?, ?, 
      ?, ?
    )
  };

  my $attrib_ins_stmt = qq{
    INSERT INTO phenotype_feature_attrib (
      phenotype_feature_id,
      attrib_type_id,
      value
    )
    SELECT ?, at.attrib_type_id, ?
    FROM attrib_type at
    WHERE at.code = ?
  };

  my $attrib_ins_cast_stmt = qq{
    INSERT INTO phenotype_feature_attrib (
      phenotype_feature_id,
      attrib_type_id,
      value
    )
    SELECT ?, at.attrib_type_id, CAST(? AS CHAR)
    FROM attrib_type at
    WHERE at.code = ?
  };

  my $attrib_id_ext_stmt = qq{
    SELECT attrib_id, value from attrib, attrib_type
    WHERE attrib_type.code = 'ontology_mapping'
    AND attrib.attrib_type_id = attrib_type.attrib_type_id
  };

  my $ontology_accession_ins_stmt = qq{
    INSERT IGNORE INTO phenotype_ontology_accession
    (phenotype_id, accession, mapped_by_attrib, mapping_type)
    values (?,?,?,?)
   };

  my $st_ins_sth   = $db_adaptor->dbc->prepare($st_ins_stmt);
  my $pf_check_sth = $db_adaptor->dbc->prepare($pf_check_stmt);
  my $pf_ins_sth   = $db_adaptor->dbc->prepare($pf_ins_stmt);
  my $attrib_ins_sth = $db_adaptor->dbc->prepare($attrib_ins_stmt);
  my $attrib_ins_cast_sth = $db_adaptor->dbc->prepare($attrib_ins_cast_stmt);
  my $ontology_accession_ins_sth = $db_adaptor->dbc->prepare($ontology_accession_ins_stmt);

  ## get the attrib id for the type of description to ontology term linking
  my $attrib_id_ext_sth = $db_adaptor->dbc->prepare($attrib_id_ext_stmt);
  $attrib_id_ext_sth->execute();
  my $ont_attrib_type = $attrib_id_ext_sth->fetchall_hashref("value");

  if ($source_info->{source_name} eq 'RGD'){ #required as db source:name is RGD and attrib_type_id 509 for ontology_mapping is Rat Genome Database
    $source_info->{source_attrib_type} = 'Rat Genome Database';
  }
  my $mapped_by = 'Data source';
  if (exists $source_info->{source_mapped_attrib_type} && exists $ont_attrib_type->{$source_info->{source_mapped_attrib_type}}){
    $mapped_by = $source_info->{source_mapped_attrib_type};
  }

  # First, sort the array according to the phenotype description
  my @sorted = sort {($a->{description} || $a->{name}) cmp ($b->{description} || $b->{name})} @{$phenotypes};
  my $study_count = 0;
  my $phenotype_feature_count = 0;

  my $total = scalar @sorted;
  my $i = 0;

  while (my $phenotype = shift(@sorted)) {
    progress($i++, $total);

    $object_type = $phenotype->{type} if defined($phenotype->{type});

    my $object_id = $phenotype->{"id"};
    # If the rs could not be mapped to a variation id, skip it
    next if $object_type =~ /Variation/ && (!defined($variation_ids->{$object_id}));

    $object_id = $variation_ids->{$object_id}[1] if ($object_type eq 'Variation');

    # if we have no coords, skip it
    my $study_id;

    if(defined($phenotype->{study}) || defined($phenotype->{study_description}) || defined($phenotype->{study_type})) {

      my $sql_study = '= ?';
      my $sql_type  = '= ?';

      # To avoid duplication of study entries
      if (!defined $phenotype->{"study"}) {$sql_study = 'IS NULL'; }
      if (!defined $phenotype->{"study_type"}) {$sql_type = 'IS NULL'; }

      my $st_check_stmt = qq{
        SELECT study_id
        FROM study
        WHERE
        source_id = $source_id AND
        external_reference $sql_study AND
        study_type $sql_type
        LIMIT 1
      };

      my $st_check_sth = $db_adaptor->dbc->prepare($st_check_stmt);
      my $param_num = 1;

      if (defined $phenotype->{"study"}) {
        $st_check_sth->bind_param($param_num,$phenotype->{"study"},SQL_VARCHAR);
        $param_num++;
      }
      if (defined $phenotype->{"study_type"}) {
        $st_check_sth->bind_param($param_num,$phenotype->{"study_type"},SQL_VARCHAR);
        $param_num++;
      }
      $st_check_sth->execute();
      $st_check_sth->bind_columns(\$study_id);
      $st_check_sth->fetch();

      if (!defined($study_id)) {
        if (length($phenotype->{"study"}) > 255) {
          warn "WARNING study external_references truncated FROM:>",$phenotype->{"study"}, "<\n";
          $phenotype->{"study"} = substr($phenotype->{"study"}, 0, 254);
          $phenotype->{"study"} = substr($phenotype->{"study"}, 0,rindex($phenotype->{"study"}, ",PMID"));
          warn "WARNING study external_references truncated TO  :>",$phenotype->{"study"}, "<\n";
        }
        $st_ins_sth->bind_param(1,$phenotype->{"study"},SQL_VARCHAR);
        $st_ins_sth->bind_param(2,$phenotype->{"study_type"},SQL_VARCHAR);
        $st_ins_sth->bind_param(3,$phenotype->{"study_description"},SQL_VARCHAR);
        $st_ins_sth->execute();
        $study_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
        $study_count++;
      }
    }

    # Remove special characters from the phenotype description
    foreach my $char (keys(%special_characters)) {
      my $new_char = $special_characters{$char};
      $phenotype->{description} =~ s/$char/$new_char/g;
    }

    # get phenotype ID
    my $phenotype_id = get_phenotype_id($phenotype, $db_adaptor);

    foreach my $acc (@{$phenotype->{accessions}}){
      $acc =~ s/\s+//g;
      $acc = iri2acc($acc) if $acc =~ /^http/;
      $ontology_accession_ins_sth->execute( $phenotype_id, $acc,  $ont_attrib_type->{$mapped_by}->{'attrib_id'}, $phenotype->{ontology_mapping_type} ) ||die "Failed to import phenotype accession\n";
    }

    if ($phenotype->{"associated_gene"}) {
      $phenotype->{"associated_gene"} =~ s/\s//g;
      $phenotype->{"associated_gene"} =~ s/;/,/g;
    }

    $phenotype->{"p_value"} = convert_p_value($phenotype->{"p_value"}) if (defined($phenotype->{"p_value"}));

    # Check if this phenotype_feature already exists for this variation and source, in that case we probably want to skip it
    my $pf_id;

    $pf_check_sth->bind_param(1,$object_id,SQL_VARCHAR);
    $pf_check_sth->bind_param(2,$object_type,SQL_VARCHAR);
    $pf_check_sth->bind_param(3,$phenotype_id,SQL_INTEGER);
    $pf_check_sth->bind_param(4,$source_id,SQL_INTEGER);
    $pf_check_sth->bind_param(5,$study_id,SQL_INTEGER);
    # For uniprot data
    if ($source_info->{source} =~ m/uniprot/i) {
      $pf_check_sth->bind_param(6,$phenotype->{"variation_names"},SQL_VARCHAR);
    }
    # For nhgri-ebi gwas data
    elsif ($source_info->{source} =~ m/nhgri/i) {
      $pf_check_sth->bind_param(6,$phenotype->{"p_value"},SQL_VARCHAR);
    }
    # For omim data
    elsif ($source_info->{source} =~ m/omim/i) {
      $pf_check_sth->bind_param(6,$phenotype->{"risk_allele"},SQL_VARCHAR);
      $pf_check_sth->bind_param(7,$phenotype->{"associated_gene"},SQL_VARCHAR);
    }

    $pf_check_sth->execute();
    $pf_check_sth->bind_columns(\$pf_id);
    $pf_check_sth->fetch();
    next if (defined($pf_id));

    my $is_significant = defined($threshold) ? ($phenotype->{"p_value"} && $phenotype->{"p_value"} < $threshold ? 1 : 0) : 1;

    # Else, insert this variation annotation
    foreach my $coord(@{$coords->{$object_id}}) {
      $pf_ins_sth->bind_param(1,$phenotype_id,SQL_INTEGER);
      $pf_ins_sth->bind_param(2,$source_id,SQL_INTEGER);
      $pf_ins_sth->bind_param(3,$study_id,SQL_INTEGER);
      $pf_ins_sth->bind_param(4,$object_type,SQL_VARCHAR);
      $pf_ins_sth->bind_param(5,$object_id,SQL_VARCHAR);
      $pf_ins_sth->bind_param(6,$is_significant,SQL_INTEGER);
      $pf_ins_sth->bind_param(7,$coord->{seq_region_id},SQL_INTEGER);
      $pf_ins_sth->bind_param(8,$coord->{seq_region_start},SQL_INTEGER);
      $pf_ins_sth->bind_param(9,$coord->{seq_region_end},SQL_INTEGER);
      $pf_ins_sth->bind_param(10,$coord->{seq_region_strand},SQL_INTEGER);
      $pf_ins_sth->execute();
      $phenotype_feature_count++;

      # get inserted ID
      $pf_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};

      # add attribs
      foreach my $attrib_type(grep {defined($phenotype->{$_}) && $phenotype->{$_} ne ''} @attrib_types) {
        my $value = $phenotype->{$attrib_type};
        my $sth = $value =~ m/^\d+(\.\d+)?$/ ? $attrib_ins_cast_sth : $attrib_ins_sth;
        $sth->bind_param(1,$pf_id,SQL_INTEGER);
        $sth->bind_param(2,$value,SQL_VARCHAR);
        $sth->bind_param(3,$attrib_type,SQL_VARCHAR);
        $sth->execute();
      }
    }
  }
  end_progress();
  print STDOUT "$study_count new studies added\n" if ($debug);
  print STDOUT "$phenotype_feature_count new phenotype_features added\n" if ($debug);
}

sub add_synonyms { #TODO: test this method when I have synonyms
  my $synonyms = shift;
  my $variation_ids = shift;
  my $source_id = shift;
  my $db_adaptor = shift;
  
  # If we actually didn't get any synonyms, just return
  return if (!defined($synonyms) || !scalar(keys(%{$synonyms})));
  
  # Some prepeared statements needed for inserting the synonyms into database
  my $ins_stmt = qq{
    INSERT IGNORE INTO
      variation_synonym (
      variation_id,
      source_id,
      name
      )
    VALUES (
      ?,
      $source_id,
      ?
    )
  };
  my $ins_sth = $db_adaptor->dbc->prepare($ins_stmt);
  
  my $alt_count = 0;
  my $variation_count = 0;
  
  foreach my $rs_id (keys %{$variation_ids}) {
    
    my $var_id = $variation_ids->{$rs_id}[0];
    
    # If we have a variation id, we can proceed
    if (defined($var_id)) {
      
      $variation_count++;
      
      $ins_sth->bind_param(1,$var_id,SQL_INTEGER);
      
      # Handle all synonym ids for this rs_id
      while (my $alt_id = shift(@{$synonyms->{$rs_id}})) {
      
        # Add the id as synonym, if it is already present, it will just be ignored
        $ins_sth->bind_param(2,$alt_id,SQL_VARCHAR);
        $ins_sth->execute();
        $alt_count++;
      }
    }
  }
  
  print STDOUT "Added $alt_count synonyms for $variation_count rs-ids\n" if ($debug);
}

sub add_set {
  my $set = shift;
  my $source_id = shift;
  my $db_adaptor = shift;
  my $set_from = shift;
  
  return if (!defined($set));
  return if (!defined($set_from) && ($set_from ne 'phenotype' || $set_from ne 'synonym')); 
   
  my $variation_set_id;
  
  # Get variation_set_id
  my $select_set_stmt = qq{
  SELECT v.variation_set_id
  FROM variation_set v, attrib a
  WHERE v.short_name_attrib_id=a.attrib_id 
  AND a.value = ?
  };
  my $sth1 = $db_adaptor->dbc->prepare($select_set_stmt);
  $sth1->bind_param(1,$set,SQL_VARCHAR);
  $sth1->execute();
  $sth1->bind_columns(\$variation_set_id);
  $sth1->fetch();
  return if (!defined($variation_set_id));
  
  # Insert into variation_set_variation
  my $insert_set_stmt = qq{
    INSERT IGNORE INTO variation_set_variation (variation_id,variation_set_id)
    SELECT distinct v.variation_id, ? 
  };
  if ($set_from eq 'phenotype') {
    $insert_set_stmt .= qq{ 
      FROM phenotype_feature pf, variation v WHERE 
        v.name=pf.object_id AND
        pf.type='Variation' AND
        pf.source_id=?
    };
  }
  elsif ($set_from eq 'synonym') {
    $insert_set_stmt .= qq{ 
      FROM variation_synonym vs, variation v WHERE 
        v.variation_id=vs.variation_id AND
        vs.source_id=?
    };
  }
  
  my $sth2 = $db_adaptor->dbc->prepare($insert_set_stmt);
  $sth2->bind_param(1,$variation_set_id,SQL_INTEGER);
  $sth2->bind_param(2,$source_id,SQL_INTEGER);
  $sth2->execute();
}

sub convert_p_value {
  my $pval = shift;
  
  my $sci_pval = '';
  # If a scientific format is not found, then ...
  if ($pval !~ /^\d+.*e.+$/i) {  
    # If a range format is found (e.g. 10^-2 > p > 10^-3)
    if ($pval =~ /^\d+\^(-\d+)/) {
      if (length("$1")==1) { $1 = "0$1"; } 
      $sci_pval = "1.00e$1"; # e.g 10^-2 > p > 10^-3 => 1.00e-2
    }
    # If a decimal format is found (e.g. 0.0023)
    elsif ($pval =~ /(\d+.*)/){
      $sci_pval = $1;
    #$sci_pval = sprintf("%.2e",$pval); # e.g. 0.002 => 2,30e-3
    }
    elsif ($pval =~ /^\w+/) {
      $sci_pval = "NULL";
    }
  }
  else {
    $pval =~ tr/E/e/;
    if ($pval =~ /^(\d+)(e-?\d+)/) {
      $pval="$1.00$2";  
    }
    if ($pval =~ /^(\d+\.\d{1})(e-?\d+)/) {
      $pval="$1"."0$2";  
    }
    $sci_pval = $pval;
  }
  return $sci_pval;
}

## format ontology accessions
sub iri2acc{

  my $iri = shift;
  my @a = split/\//, $iri;
  my $acc = pop @a;
  $acc =~ s/\_/\:/;

  return $acc;
}

sub get_seq_region_ids {
  my ($self, $db_adaptor) = @_;

  my $sth = $db_adaptor->dbc->prepare(qq{
    SELECT seq_region_id, name
    FROM seq_region
  });
  $sth->execute;

  my (%seq_region_ids, $id, $name);
  $sth->bind_columns(\$id, \$name);
  $seq_region_ids{$name} = $id while $sth->fetch();
  $sth->finish;

  return \%seq_region_ids;
}

sub get_coords {
  my $ids = shift;
  my $variation_ids = shift;
  my $type = shift;
  my $db_adaptor = shift;
  
  my $tables;
  my $where_clause;
  
  my @object_ids = ($type eq 'Variation') ? map { $variation_ids->{$_}[1] } keys(%$variation_ids): @$ids;

  if($type eq 'Variation') {
    $tables = 'variation_feature f, variation v';
    $where_clause = 'v.variation_id = f.variation_id AND v.name = ?';
  }
  elsif($type =~ /Structural/) {
    $tables = 'structural_variation_feature f, structural_variation v';
    $where_clause = 'v.structural_variation_id = f.structural_variation_id AND v.variation_name = ?';
  }
  elsif($type eq 'Gene') {
    $tables = 'gene f';
    $where_clause = 'f.stable_id = ?';
  }
  
  my $sth = $db_adaptor->dbc->prepare(qq{
    SELECT
      f.seq_region_id, f.seq_region_start, f.seq_region_end, f.seq_region_strand
    FROM
      $tables
    WHERE
      $where_clause
  });
  
  my $coords = {};
  my ($sr_id, $start, $end, $strand);
  
  foreach my $id (@object_ids) {
    $sth->bind_param(1,$id,SQL_VARCHAR);
    $sth->execute();
    $sth->bind_columns(\$sr_id, \$start, \$end, \$strand);
    
    push @{$coords->{$id}}, {
      seq_region_id => $sr_id,
      seq_region_start => $start,
      seq_region_end => $end,
      seq_region_strand => $strand
    } while $sth->fetch();
  }
  
  $sth->finish();
  
  return $coords;
}

sub get_dbIDs {
  my $rs_ids = shift;
  my $db_adaptor = shift;

  my $id_stmt = qq{
    SELECT DISTINCT
      v.variation_id, v.name
    FROM
      variation v,
      variation_feature vf
    WHERE
      v.name = ? AND
      v.variation_id=vf.variation_id
    LIMIT 1
  };
  my $syn_stmt = qq{
    SELECT DISTINCT
      v.variation_id, v.name
    FROM
      variation_feature vf,
      variation_synonym vs JOIN
      variation v ON vs.variation_id = v.variation_id
    WHERE
      vs.name = ? AND
      v.variation_id=vf.variation_id
    LIMIT 1
  };
  my $id_sth = $db_adaptor->dbc->prepare($id_stmt);
  my $syn_sth = $db_adaptor->dbc->prepare($syn_stmt);

  my %mapping;

  foreach my $rs_id (@{$rs_ids}) {
    $id_sth->bind_param(1,$rs_id,SQL_VARCHAR);
    $id_sth->execute();
    my ($var_id,$var_name);
    $id_sth->bind_columns(\$var_id,\$var_name);
    $id_sth->fetch();

    # If we couldn't find the rs_id, look in the synonym table
    if (!defined($var_id)) {
      $syn_sth->bind_param(1,$rs_id,SQL_VARCHAR);
      $syn_sth->execute();
      $syn_sth->bind_columns(\$var_id,\$var_name);
      $syn_sth->fetch();
    }
    warn "$rs_id - no mapping found in variation db\n" unless $var_id ;
    $mapping{$rs_id} = [$var_id,$var_name] if $var_id && $var_name;
  }

  return \%mapping;
}

sub get_attrib_types {
  my $db_adaptor = shift;

  my $sth = $db_adaptor->dbc->prepare(qq{
    SELECT code
    FROM attrib_type
  });
  $sth->execute();

  my ($attrib_type, @tmp_types);
  $sth->bind_columns(\$attrib_type);
  while ($sth->fetch()) {
    push @tmp_types, $attrib_type if ($attrib_type ne 'description');
  }
  $sth->finish;

  return \@tmp_types;
}

sub get_or_add_source {
  my $self = shift;

  my $source_info = shift;
  my $db_adaptor  = shift;
  
#  my $source_name = shift;
#  my $source_description = shift;
#  my $source_url    = shift;
#  my $source_status = shift;
#  my $source_version = shift;
#  my $db_adaptor    = shift;
#  my $debug     = shift;

  my $stmt = qq{
    SELECT
      source_id
    FROM
      source
    WHERE
      name = '$source_info->{source_name}'
    LIMIT 1
  };
  my $sth = $db_adaptor->dbc->prepare($stmt);
  $sth->execute();
  my $source_id;
  $sth->bind_columns(\$source_id);
  $sth->fetch();

  if (!defined($source_id)) {
    $stmt = qq{
      INSERT INTO
        source (name, description, url, somatic_status, version )
      VALUES (
        '$source_info->{source_name}',
        '$source_info->{source_description}',
        '$source_info->{source_url}',
        '$source_info->{source_status}',
        $source_info->{source_version}
      )
    };
    $db_adaptor->dbc->do($stmt);
    $source_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};

    print STDOUT "Added source for $source_info->{source_name} (source_id = $source_id)\n" if ($debug);
  }
  else {
    $stmt = qq{
      UPDATE
        source
      SET name=?,
        description=?,
        url=?,
        version=?
      WHERE
        source_id=?
    };
    my $update_source_sth =$db_adaptor->dbc->prepare($stmt);
    $update_source_sth->execute($source_info->{source_name},
                                $source_info->{source_description},
                                $source_info->{source_url},
                                $source_info->{source_version},
                                $source_id);
  }

  return $source_id;
}

sub get_phenotype_id {
  my $phenotype = shift;
  my $db_adaptor = shift;

  my ($name, $description);
  $name = $phenotype->{name};
  $description = $phenotype->{description};

  # Clean up
  $description =~ s/^\s+|\s+$//g; # Remove spaces at the beginning and the end of the description
  $description =~ s/\n//g; # Remove 'new line' characters

  # Check phenotype description in the format "description; name"
  if (!defined($name) || $name eq '') {
    my ($p_desc,$p_name) = split(";",$description);
    if ($p_name) {
      $p_name =~ s/ //g;
      if ($p_name =~ /^\w+$/) {
        $description = $p_desc;
        $name = $p_name;
      }
    }
  }

  if(scalar keys %phenotype_cache) {
    
    # check cache first
    return $phenotype_cache{$description} if defined $phenotype_cache{$description};
    
    my @tmp = keys %phenotype_cache;
    
    # lc everything
    my $description_bak = $description;
    $description = lc($description);
    
    # store mapped
    my %mapped;
    $mapped{lc($_)} = $_ for @tmp;
    @tmp = keys %mapped;
    
    # check if it matches lc
    if(defined($mapped{$description})) {
      $phenotype_cache{$description_bak} = $phenotype_cache{$mapped{$description}};
      return $phenotype_cache{$description_bak};
    }
    
    # try a fuzzy match using String::Approx
    my @matches = amatch($description, [10], @tmp);
  
    if(@matches) {
      
      # we only want the best match
      my $best = scalar @matches == 1 ? $matches[0] : (sort {abs(adist($description, $a)) <=> abs(adist($description, $b))} @matches)[0];
      print STDERR "\nPHENOTYPE:\n\tINPUT: $description\n\tBEST:  $best\n\tDIST: ".adist($description, $best)."\n";
      
      my $skip = 0;

      # Assuming a perfect match check has been done in the previous lines, e.g. if ($phenotype_cache{$description})
      $skip = 1 if (adist($description, $best) == 0);

      # find characters that differ
      my $diff = diff([split(//,$description)], [split(//,$best)]);
      
      # skip if mismatch is anything word-like
      my $diff_string = '';
      my $previous_diff_pos;
      foreach(map {@$_} @$diff) {
        my $diff_pos  = $_->[1];
        my $diff_char = $_->[2];

        $skip = 1 if $diff_char =~ /\w/i;

        if (!$previous_diff_pos) {
          $diff_string .= "[$diff_char";
        }
        elsif ($diff_pos == $previous_diff_pos) {
          $diff_string .= "/$diff_char]";
        }
        elsif ($diff_pos == ($previous_diff_pos+1)) {
          $diff_string .= $diff_char;
        }
        else {
          $diff_string .= "] [$diff_char";
        }
        $previous_diff_pos = $diff_pos;
      }
      if ( $diff_string ne '') {
        $diff_string .= ']' if ($diff_string !~ /\]$/);
        print STDERR "\tDIFFERENCE(S): $diff_string\n" if ( $diff_string ne '');
      }

      # cache this match so we don't have to fuzz again
      $phenotype_cache{$description_bak} = $phenotype_cache{$mapped{$best}};
      unless ($skip) {
        print STDERR "\tUSED (with diff): $best\n";
        return $phenotype_cache{$mapped{$best}};
      }
    }
    
    # restore from backup before inserting
    $description = $description_bak;
  }


  # finally if no match, do an insert
  my $sth = $db_adaptor->dbc->prepare(qq{
    INSERT IGNORE INTO phenotype ( name, description ) VALUES ( ?,? )
  });
  
  $sth->bind_param(1,$name,SQL_VARCHAR);
  $sth->bind_param(2,$description,SQL_VARCHAR);
  $sth->execute();
  my $phenotype_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
  
  # update cache
  $phenotype_cache{$description} = $phenotype_id;
  
  return $phenotype_id;
}

#TODO: Q: do I still need a progress bar?
# update or initiate progress bar
sub progress {
  my ($i, $total) = @_;
  
  my $width = 60;
  my $percent = int(($i/$total) * 100);
  my $numblobs = int((($i/$total) * $width) - 2);
  
  # this ensures we're not writing to the terminal too much
  return if defined($prev_prog) && $numblobs.'-'.$percent eq $prev_prog;
  $prev_prog = $numblobs.'-'.$percent;
  
  printf("\r% -${width}s% 1s% 10s", '['.('=' x $numblobs).($numblobs == $width - 2 ? '=' : '>'), ']', "[ " . $percent . "% ]");
}

# end progress bar
sub end_progress {
  progress(1,1);
  print "\n";
}

sub run {

    my $self = shift;

    #TODO: figure out if I need this
}

=head2 save_phenotypes

    Arg [1]              : hashref $source_info
    Arg [2]              : hashref $input_data
    Arg [3]              : Bio::EnsEMBL::DBSQL::DBAdaptor $core_dba
    Arg [4]              : Bio::EnsEMBL::Variation::DBSQL::DBAdaptor $variation_dba
    Example              : $self->save_phenotypes(\%source_info, $input_data, $core_dba, $variation_dba);
    Description          : Saves the phenotype data ($input_data) for the given source ($source_info) using the given core and variation adaptors.
    Returntype           : none
    Exceptions           : none
    Caller               : general
    Status               : Stable

=cut
sub save_phenotypes {
  my $self        = shift;

  my $source_info = shift;
  my $input_data      = shift;
  my $core_dba    = shift;
  my $variation_dba = shift;

  my $phenotype_dba =  $variation_dba->get_PhenotypeAdaptor; #TODO: maybe work a way to have these global, OR move lines 577-579 to general one
#  $variation_dba  = $self->get_species_adaptor('variation'); #TODO: why did the init -> Base class -> this not work?
  my $debug = $self->param('debug_mode');

  my $set = defined($source_info->{set}) ? $source_info->{set} : undef;
  $source_info->{source_status} = 'germline' unless defined($source_info->{source_status});
  
  my $initial_phenotype_count = $phenotype_dba->generic_count;
  my @attrib_types = @{get_attrib_types($variation_dba)};
  %phenotype_cache = map {$_->description() => $_->dbID()} @{$phenotype_dba->fetch_all};
  
  my @ids;
  my %synonym;
  my @phenotypes;
  if (exists($input_data->{'synonyms'})) {
    %synonym = %{$input_data->{'synonyms'}};
    # To get all the ids of the source (Uniprot)
    @ids = keys(%synonym);
  }
  if (exists($input_data->{'phenotypes'})) {
    @phenotypes = @{$input_data->{'phenotypes'}};
  }
  print STDOUT "Got ".(scalar @phenotypes)." objects\n" if ($debug);

  # Get internal variation ids for the rsIds
  print STDOUT "Retrieving internal variation IDs\n" if ($debug);
  if (scalar @ids == 0) {
    @ids = map {$_->{'id'}} @phenotypes;
  }
  my $variation_ids = $source_info->{object_type} =~ /Variation/ ? get_dbIDs(\@ids,$variation_dba) : {};

  # Get coordinates of objects
  my $coords;

  print STDOUT "Retrieving object coordinates\n" if ($debug);

  # might be able to copy them from data (QTLs and Genes)
  if(defined($phenotypes[0]->{seq_region_id})) {
    foreach my $p(@phenotypes) {
      my $coord = {};
      
      foreach my $key(qw(id start end strand)) {
        $coord->{'seq_region_'.$key} = $p->{'seq_region_'.$key};
      }
      
      push @{$coords->{$p->{id}}}, $coord;
    }
  }

  $coords ||= get_coords(\@ids,$variation_ids, $source_info->{object_type}, $source_info->{object_type} eq 'Gene' ? $core_dba : $variation_dba);

  # uniquify coords
  foreach my $id (keys %$coords) {
    my %tmp = map {
      $_->{seq_region_id}."_".
      $_->{seq_region_start}."_".
      $_->{seq_region_end}."_".
      $_->{seq_region_strand}."_" => $_
    } @{$coords->{$id}};
    
    $coords->{$id} = [values %tmp];
  }

  # Get or add a source
  my $source_id = $self->get_or_add_source($source_info,$variation_dba);
  print STDOUT "$source_info->{source} source_id is $source_id\n" if ($debug);

  # Add the synonyms if required
  unless ($skip_synonyms) {
    print STDOUT "Adding synonyms\n" if ($debug);
    add_synonyms(\%synonym,$variation_ids,$source_id,$variation_dba);
  }

  # Now, insert phenotypes
  unless ($skip_phenotypes) {
    die("ERROR: No phenotypes or objects retrieved from input\n") unless scalar @phenotypes;

    print STDOUT "Adding phenotypes\n" if ($debug);
    $source_info->{source_id}=$source_id;
    my %data = (phenotypes => \@phenotypes, coords => $coords, 
                variation_ids => $variation_ids, attrib_types=> \@attrib_types);
    add_phenotypes(\%data,$source_info,$variation_dba);

    print STDOUT "$initial_phenotype_count initial phenotypes\n" if ($debug);
    my $added_phenotypes = $phenotype_dba->generic_count - $initial_phenotype_count;
    print STDOUT "$added_phenotypes new phenotypes added\n" if ($debug);
  }

  # Add the variation sets if required
  unless ($skip_sets) {
    if (%synonym && !$skip_synonyms && $skip_phenotypes) {
      print STDOUT "Adding variation sets for synonyms\n" if ($debug);
      add_set($set,$source_id,$variation_dba,'synonym');
    }
    elsif (@phenotypes && !$skip_phenotypes) {
      print STDOUT "Adding variation sets for phenotypes\n" if ($debug);
      add_set($set,$source_id,$variation_dba,'phenotype');
    }
  }
}

1;
