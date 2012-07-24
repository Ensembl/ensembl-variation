use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use DBI qw(:sql_types);

my $registry_file = "./ensembl.registry";
my $species = 'human';

my $var_table = 'HGMD_PUBLIC_variation';
my $vf_table  = 'HGMD_PUBLIC_variation_feature';
my $fs_table  = 'HGMD_PUBLIC_flanking_sequence';
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
add_feature();
add_flanking();
add_allele();
add_annotation();
add_attrib();
add_set();



# Methods
sub add_variation {

	my $select_vh_sth = $dbh->prepare(qq{
		SELECT distinct name FROM $var_table;
	});

	my $insert_v_sth = $dbh->prepare(qq{
		INSERT INTO variation (name,source_id) VALUES (?,?);
	});

	my $select_v_sth = $dbh->prepare(qq{
		SELECT variation_id FROM variation WHERE name=? and source_id=?;
	});

	my $update_vh_sth = $dbh->prepare(qq{
		UPDATE $var_table SET new_var_id=? WHERE name=?;
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
}


sub add_feature {
	my $insert_vf_sth = $dbh->prepare(qq{
		INSERT INTO 
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
	});
	$insert_vf_sth->execute($source_id);
}


sub add_flanking {
	my $insert_fs_sth = $dbh->prepare(qq{
		INSERT INTO 
			flanking_sequence (
				variation_id,
				up_seq_region_start, 
				up_seq_region_end,
				down_seq_region_start,
				down_seq_region_end,
				seq_region_id,
				seq_region_strand
			) 
			SELECT 
				v.new_var_id,
				fs.up_seq_region_start, 
				fs.up_seq_region_end,
				fs.down_seq_region_start,
				fs.down_seq_region_end,
				fs.seq_region_id,
				1
			FROM $fs_table fs, $var_table v 
			WHERE v.variation_id=fs.variation_id
	});
	$insert_fs_sth->execute();
}


sub add_allele {
	my $insert_al_sth;
	my $allele = 'HGMD_MUTATION';
	
	if ($dbh->do(qq{show columns from allele like 'allele';}) != 1){
		$insert_al_sth = $dbh->prepare(qq{
			INSERT INTO 
				allele (
					variation_id,
					allele_code_id
				) 
				SELECT DISTINCT
					v.new_var_id,
					a.allele_code_id
				FROM $var_table v, allele_code a 
				WHERE a.allele="$allele";
		});
	} else {

		$insert_al_sth = $dbh->prepare(qq{
			INSERT INTO 
				allele (
					variation_id,
					allele
				) 
				SELECT DISTINCT
					new_var_id,
					"$allele"
				FROM $var_table
		});
	}
	$insert_al_sth->execute();
}


sub add_annotation {
	
	# Study
	my $select_study_sth = $dbh->prepare(qq{
		SELECT study_id FROM study WHERE source_id=?;
	});
	$select_study_sth->execute($source_id);
	my $study_id = ($select_study_sth->fetchrow_array)[0];
	if (!defined($study_id)) { 
		$dbh->do(qq{INSERT INTO study (name,source_id) VALUES ('HGMD-PUBLIC',$source_id)});
		$study_id = $dbh->{'mysql_insertid'};
	}
	if (!defined($study_id)) { print "Study not found\n"; exit(0); }
	
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

	my $insert_va_sth = $dbh->prepare(qq{
		INSERT INTO 
			variation_annotation (
				associated_gene,
				variation_id,
				variation_names,
				study_id,
				phenotype_id
			) 
			SELECT 
				va.associated_gene,
				v.new_var_id,
				v.name,
				?,
				?
			FROM $va_table va, $var_table v 
			WHERE v.variation_id=va.variation_id
	});
	$insert_va_sth->execute($study_id, $phenotype_id);
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
			FROM variation_annotation WHERE study_id IN
				(SELECT study_id FROM study WHERE source_id=?)
	};
	my $sth2 = $dbh->prepare($insert_set_stmt);
	$sth2->bind_param(1,$variation_set_id,SQL_INTEGER);
  $sth2->bind_param(2,$source_id,SQL_INTEGER);
  $sth2->execute();
}
