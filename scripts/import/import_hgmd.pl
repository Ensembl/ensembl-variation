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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use DBI qw(:sql_types);
use Getopt::Long;
our ($registry_file);

GetOptions('registry=s' => \$registry_file);

$registry_file ||= "./ensembl.registry";
my $species = 'human';

my $var_table = 'HGMD_PUBLIC_variation';
my $vf_table  = 'HGMD_PUBLIC_variation_feature';
my $va_table  = 'HGMD_PUBLIC_variation_annotation';
my $short_set = 'ph_hgmd_pub';

Bio::EnsEMBL::Registry->load_all( $registry_file );
my $vdb2 = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbh = $vdb2->dbc->db_handle;

my $select_source_sth = $dbh->prepare(qq{
  SELECT source_id FROM source WHERE name='HGMD-PUBLIC';
});
$select_source_sth->execute();
my $source_id = ($select_source_sth->fetchrow_array)[0];

die ("Source not found") if (!defined($source_id));
die ("HGMD tables not found in the database!\nBe sure you ran the script 'map_hgmd_coord.pl' before running the current script.\n") if ($dbh->do(qq{show tables like "$var_table";}) != 1);


# Main
add_variation();
add_variation_feature();
add_annotation(); # phenotype_feature & phenotype_feature_attrib
update_features();
add_attrib();
add_set();



# Methods
sub add_variation {

  my $select_vh_sth = $dbh->prepare(qq{
    SELECT distinct name FROM $var_table;
  });

  my $insert_v_sth = $dbh->prepare(qq{
    INSERT IGNORE INTO variation (name,source_id) VALUES (?,?);
  });

  my $select_v_sth = $dbh->prepare(qq{
    SELECT variation_id FROM variation WHERE name=? and source_id=?;
  });

  my $update_vh_sth = $dbh->prepare(qq{
    UPDATE $var_table SET new_var_id=? WHERE name=?;
  });
  
  my $select_not_existing_sth = $dbh->prepare(qq{
    SELECT name FROM variation WHERE source_id=? AND name NOT IN (SELECT name FROM $var_table)
  });


  if ($dbh->do(qq{show columns from $var_table like 'new_var_id';}) != 1){
    $dbh->do(qq{ALTER TABLE $var_table ADD COLUMN new_var_id int(10);});
    $dbh->do(qq{ALTER TABLE $var_table ADD INDEX new_var_idx (new_var_id);});
  }


  $select_vh_sth->execute();
  while (my @res = $select_vh_sth->fetchrow_array) {
    $insert_v_sth->execute($res[0],$source_id);
  
    $select_v_sth->execute($res[0],$source_id);
    my $new_id = ($select_v_sth->fetchrow_array)[0];
    if (defined($new_id)) {
      $update_vh_sth->execute($new_id,$res[0]);
    }
  }
  
  # Check if some HGMD entries have beeen remove from the previous version
  my @not_existing;
  $select_not_existing_sth->execute($source_id);
  while (my @res = $select_not_existing_sth->fetchrow_array) {
    push(@not_existing,$res[0]);
  }
  if (scalar @not_existing) {
    print "WARNING: ".scalar(@not_existing)." entries are not anymore in the HGMD database!\n";
    print "Please, remove the following HGMD mutations from the variation database:\n";
    print join("\n",@not_existing)."\n";
  }
}


sub add_variation_feature {
  my $insert_vf_sth = $dbh->prepare(qq{
    INSERT IGNORE INTO 
      variation_feature (
        seq_region_id,
        seq_region_start,
        seq_region_end,
        seq_region_strand,
        allele_string,
        map_weight,
        variation_id,
        variation_name,
        source_id
      ) 
      SELECT 
        vf.seq_region_id,
        vf.seq_region_start,
        vf.seq_region_end,
        1,
        vf.allele_string,
        vf.map_weight,
        v.new_var_id,
        v.name,
        ?
        FROM $vf_table vf, $var_table v 
        WHERE v.variation_id=vf.variation_id
        AND NOT EXISTS ( SELECT * 
                         FROM variation_feature vf2 
                         WHERE vf2.variation_id=v.new_var_id
                       );
  });
  $insert_vf_sth->execute($source_id);
}

sub add_annotation {
  
  # Phenotype
  my $select_phe_sth = $dbh->prepare(qq{
    SELECT distinct phenotype_id FROM phenotype where name='HGMD_MUTATION' limit 1;
  });
  $select_phe_sth->execute();
  my $phenotype_id = ($select_phe_sth->fetchrow_array)[0];
  if (!defined($phenotype_id)) {
    $dbh->do(qq{INSERT INTO phenotype (name,description) 
                VALUES ('HGMD_MUTATION','Annotated by HGMD but no phenotype description is publicly available')
               });
    $phenotype_id = $dbh->{'mysql_insertid'};
  }
  if (!defined($phenotype_id)) { print "Phenotype not found\n"; exit(0); }

  # Phenotype feature
  my $insert_pf_sth = $dbh->prepare(qq{
    INSERT IGNORE INTO 
      phenotype_feature (
        phenotype_id,
        object_id,
        source_id,
        type,
        seq_region_id,
        seq_region_start,
        seq_region_end,
        seq_region_strand
      ) 
      SELECT 
        ?,
        v.name,
        ?,
        'Variation',
        vf.seq_region_id,
        vf.seq_region_start,
        vf.seq_region_end,
        vf.seq_region_strand 
      FROM $va_table va, $var_table v , $vf_table vf
      WHERE v.variation_id=va.variation_id AND v.variation_id=vf.variation_id
      AND NOT EXISTS ( SELECT * FROM phenotype_feature pf2
                       WHERE pf2.object_id=v.name AND pf2.type='Variation' 
                     );
  });
  $insert_pf_sth->execute($phenotype_id,$source_id);
  
  # Attribute type "associated_gene"
  my $select_attr_sth = $dbh->prepare(qq{
    SELECT distinct attrib_type_id FROM attrib_type where code='associated_gene' limit 1;
  });
  $select_attr_sth->execute();
  my $attrib_type_id = ($select_attr_sth->fetchrow_array)[0];
  die("Variation annotation - Can't find an attribute type ID for 'associated_gene in the table attrib_type'") if (!defined($attrib_type_id));
  
  # Phenotype feature attrib
  my $insert_pfa_sth = $dbh->prepare(qq{
    INSERT IGNORE INTO 
      phenotype_feature_attrib (
        phenotype_feature_id,
        attrib_type_id,
        value
      ) 
      SELECT DISTINCT
        pf.phenotype_feature_id,
        ?,
        va.associated_gene
      FROM $var_table v, $va_table va, phenotype_feature pf
      WHERE v.variation_id=va.variation_id 
        AND v.name=pf.object_id
        AND pf.type='Variation'
  });
  $insert_pfa_sth->execute($attrib_type_id);
}


sub update_features {
	# Variation feature
	my $update_vf_sth = $dbh->prepare(qq{
	  UPDATE $vf_table hvf, $var_table hv, variation_feature vf SET
		  vf.seq_region_id=hvf.seq_region_id,
			vf.seq_region_start=hvf.seq_region_start,
			vf.seq_region_end=hvf.seq_region_end
		WHERE
		  vf.variation_id=hv.new_var_id AND
			hv.variation_id=hvf.variation_id AND 
			vf.source_id = ? AND
			(vf.seq_region_id!=hvf.seq_region_id OR
			 vf.seq_region_start!=hvf.seq_region_start OR
			 vf.seq_region_end!=hvf.seq_region_end) 
	});
	$update_vf_sth->execute($source_id);
	
	# Phenotype feature
	my $update_pf_sth = $dbh->prepare(qq{
	  UPDATE phenotype_feature pf, $vf_table vf SET
		  pf.seq_region_id=vf.seq_region_id,
			pf.seq_region_start=vf.seq_region_start,
			pf.seq_region_end=vf.seq_region_end
		WHERE
		  pf.object_id=vf.variation_name AND
			pf.source_id = ? AND
			(vf.seq_region_id!=pf.seq_region_id OR
			 vf.seq_region_start!=pf.seq_region_start OR
			 vf.seq_region_end!=pf.seq_region_end) 
	});
	$update_pf_sth->execute($source_id);
}


sub add_attrib {
  my %attrib = ('M' => 'SNV',
                'D' => 'deletion',
                'I' => 'insertion',
                'X' => 'indel',
                'P' => 'indel',
                'R' => 'sequence_alteration',
                'S' => 'sequence_alteration'
               );

  my %class = ();

  my $select_a_sth = $dbh->prepare(qq{
    SELECT a.attrib_id FROM attrib a, attrib_type att 
    WHERE a.attrib_type_id = att.attrib_type_id AND att.code = 'SO_term' AND a.value = ?;
  });


  my $select_v_sth = $dbh->prepare(qq{
    SELECT DISTINCT new_var_id,type FROM $var_table;
  });

  my $update_v_sth = $dbh->prepare(qq{
    UPDATE variation SET class_attrib_id = ? WHERE variation_id = ?;
  });

  my $update_vf_sth = $dbh->prepare(qq{
    UPDATE variation_feature vf, variation v SET vf.class_attrib_id = v.class_attrib_id 
    WHERE v.variation_id = vf.variation_id AND v.source_id=?;
  });

  while (my ($k,$v) = each (%attrib)) {
    $select_a_sth->execute($v);
    $class{$k} = ($select_a_sth->fetchrow_array)[0];
    print "$k: ".$class{$k}."\n";
  }

  $select_v_sth->execute();
  while (my @res = $select_v_sth->fetchrow_array) {
    my $att = $class{$res[1]};
    if (defined $att) {
      $update_v_sth->execute($att,$res[0]) or die $!;
    }
  }

  $update_vf_sth->execute($source_id) or die $!;
}


sub add_set {
  
  my $variation_set_id;
  
  # Get variation_set_id
  my $select_set_stmt = qq{
        SELECT v.variation_set_id
        FROM variation_set v, attrib a
        WHERE v.short_name_attrib_id=a.attrib_id 
          AND a.value = ?
  };
  my $sth1 = $dbh->prepare($select_set_stmt);
  $sth1->bind_param(1,$short_set,SQL_VARCHAR);
  $sth1->execute();
  $sth1->bind_columns(\$variation_set_id);
  $sth1->fetch();
  return if (!defined($variation_set_id));
  
  # Insert into variation_set_variation
  my $insert_set_stmt = qq{ 
    INSERT IGNORE INTO variation_set_variation (variation_id,variation_set_id)
      SELECT distinct variation_id, ? 
      FROM variation WHERE source_id=?
  };
  my $sth2 = $dbh->prepare($insert_set_stmt);
  $sth2->bind_param(1,$variation_set_id,SQL_INTEGER);
  $sth2->bind_param(2,$source_id,SQL_INTEGER);
  $sth2->execute();
}
