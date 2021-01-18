=head1 LICENSE
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
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

use strict;
use warnings;

package Bio::EnsEMBL::Variation::Pipeline::Remapping::FinishVariationFeatureQC;

use base ('Bio::EnsEMBL::Variation::Pipeline::Remapping::BaseRemapping');

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Utils::RemappingUtils qw(qc_mapped_vf);
use FileHandle;

sub fetch_input {
  my $self = shift;
}

sub run {
  my $self = shift;
  my $qc_update_features_dir = $self->param('qc_update_features_dir');
  my $qc_failure_reasons_dir = $self->param('qc_failure_reasons_dir');
  my $feature_table = $self->param('feature_table') . '_mapping_results';

  $self->add_flip_column($feature_table);
  $self->run_qc_updates($qc_update_features_dir);
#  $self->backup_tables(['allele', 'population_genotype', 'tmp_sample_genotype_single_bp', 'failed_variation']);

  my $failed_variations_after_remapping = $self->failed_variations_after_remapping($qc_failure_reasons_dir);
  my $previous_unmapped_variants = $self->previous_unmapped_variants(); # did not map in previous assembly
  my $unmapped_variants = $self->unmapped_variants($feature_table); # couldn't be mapped to new assembly

  foreach my $variant_id (keys %$unmapped_variants) {
    if ($previous_unmapped_variants->{$variant_id}) {
      $failed_variations_after_remapping->{$variant_id}->{5} = 1; # Variant does not map to the genome
    } else {
      $failed_variations_after_remapping->{$variant_id}->{17} = 1; # Variant can not be re-mapped to the current assembly
    }
  }
  $self->cleanup_failed_variants($feature_table);
  $self->load_failed_variants($failed_variations_after_remapping);

  $self->init_flip_features($feature_table);
  $self->flip_features();

  $self->update_variation_set_variation;
  $self->update_display($feature_table);
  $self->cleanup_mapped_feature_table;
}

sub add_flip_column {
  my $self = shift;
  my $feature_table = shift;
  my $vdba = $self->get_newasm_variation_database_connection;
  my $dbh = $vdba->dbc->db_handle;
  if ($self->column_exists($vdba, $feature_table, 'flip')) {
    $dbh->do("ALTER TABLE $feature_table DROP COLUMN flip;")
  }
  $dbh->do("ALTER TABLE $feature_table ADD COLUMN flip INT DEFAULT 0;")
}

sub run_qc_updates {
  my $self = shift;
  my $qc_update_features = shift;
  my $vdba = $self->get_newasm_variation_database_connection;
  my $dbh = $vdba->dbc->db_handle;
  opendir(DIR, $qc_update_features) or die $!;
  while (my $file = readdir(DIR)) {
    next unless (-f "$qc_update_features/$file");
    next unless ($file =~ m/\.txt$/);
    my $fh = FileHandle->new("$qc_update_features/$file", 'r');
    while(<$fh>) {
      chomp;
      $dbh->do("$_") or die $!;
    }
    $fh->close;
  }
  closedir(DIR);
}

sub init_flip_features {
  my $self = shift;
  my $feature_table = shift;
  my $pipeline_dir = $self->param('pipeline_dir');

  my $failed_variations = $self->failed_variations($feature_table);
  my $failed_alleles = $self->failed_alleles($feature_table); 
 
  # empty files 
  open my $fh, ">$pipeline_dir/update_flip_features.txt" or die $!;
  open my $fh_err, ">$pipeline_dir/errors_update_flip_features.txt" or die $!;

  $self->dump_features_for_flipping(qq{
  SELECT vf.seq_region_id, vf.seq_region_start, vf.seq_region_end, vf.seq_region_strand, vf.variation_id, vf.variation_name, vf.allele_string, vf.map_weight, a.allele_id, a.allele_code_id
  FROM $feature_table vf
  INNER JOIN allele a ON vf.variation_id = a.variation_id
  WHERE vf.flip = 1
  AND vf.map_weight = 1;
  }, "$pipeline_dir/alleles_for_flipping.txt");
  $self->flip_alleles($failed_alleles, $failed_variations, "$pipeline_dir/alleles_for_flipping.txt", $fh, $fh_err);

  $self->dump_features_for_flipping(qq{
  SELECT vf.seq_region_id, vf.seq_region_start, vf.seq_region_end, vf.seq_region_strand, vf.variation_id, vf.variation_name, vf.allele_string, vf.map_weight, pg.population_genotype_id, pg.genotype_code_id
  FROM $feature_table vf
  INNER JOIN population_genotype pg ON vf.variation_id = pg.variation_id
  WHERE vf.flip = 1
  AND vf.map_weight = 1;
  }, "$pipeline_dir/population_genotypes_for_flipping.txt");
  $self->flip_population_genotypes("$pipeline_dir/population_genotypes_for_flipping.txt", $fh, $fh_err);

  $self->dump_features_for_flipping(qq{
  SELECT vf.seq_region_id, vf.seq_region_start, vf.seq_region_end, vf.seq_region_strand, vf.variation_id, vf.variation_name, vf.allele_string, vf.map_weight, sg.subsnp_id, sg.allele_1, sg.allele_2, sg.sample_id
  FROM $feature_table vf
  INNER JOIN tmp_sample_genotype_single_bp sg ON vf.variation_id = sg.variation_id
  WHERE vf.flip = 1
  AND vf.map_weight = 1;
  }, "$pipeline_dir/sample_genotypes_for_flipping.txt");
  $self->flip_sample_genotypes("$pipeline_dir/sample_genotypes_for_flipping.txt", $fh, $fh_err);

  close $fh;
  close $fh_err;
}

sub flip_features {
  my $self = shift;
  my $pipeline_dir = $self->param('pipeline_dir');
  my $vdba = $self->get_newasm_variation_database_connection;
  my $dbh = $vdba->dbc->db_handle;
  my $fh = FileHandle->new("$pipeline_dir/update_flip_features.txt", 'r');
  while(<$fh>) {
    chomp;
    $dbh->do("$_") or die $!;
  }
  $fh->close;
}

sub backup_tables {
  my $self = shift;
  my $tables = shift;
  my $vdba = $self->get_newasm_variation_database_connection;
  my $dbh = $vdba->dbc->db_handle;

  foreach my $table (@$tables) {
    my $sth = $dbh->prepare(qq{
      SHOW TABLES LIKE 'before_remapping_qc_$table';
    });
    $sth->execute() or die 'Could not execute statement ' . $sth->errstr;
    my @row = $sth->fetchrow_array;
    $sth->finish();
    if (!@row) {
      $dbh->do(qq{CREATE TABLE before_remapping_qc_$table LIKE $table;}) or die $dbh->errstr; 
      $dbh->do(qq{INSERT INTO before_remapping_qc_$table SELECT * FROM $table;}) or die $dbh->errstr; 
    }
  }
}

sub failed_alleles {
  my $self = shift;
  my $feature_table = shift;
  my $failed_alleles = {};
  my $vdba = $self->get_newasm_variation_database_connection;
  my $dbh = $vdba->dbc->db_handle;
  my $sth = $dbh->prepare(qq{
  SELECT fa.allele_id, fa.failed_description_id, a.variation_id
  FROM $feature_table vf
  INNER JOIN allele a ON vf.variation_id = a.variation_id
  INNER JOIN failed_allele fa ON a.allele_id = fa.allele_id
  WHERE vf.flip = 1
  AND vf.map_weight = 1;
  }, {mysql_use_result => 1});

  $sth->execute();
  while (my $row = $sth->fetchrow_arrayref) {
    my ($allele_id, $failed_description_id, $variation_id) = @$row;
    $failed_alleles->{$variation_id}->{$allele_id}->{$failed_description_id} = 1;
  }
  $sth->finish();
  return $failed_alleles;
}

sub failed_variations {
  my $self = shift;
  my $feature_table = shift;
  my $failed_variations = {};
  my $vdba = $self->get_newasm_variation_database_connection;
  my $dbh = $vdba->dbc->db_handle;
  my $sth = $dbh->prepare(qq{
  SELECT fv.variation_id, fv.failed_description_id
  FROM $feature_table vf
  INNER JOIN allele a ON vf.variation_id = a.variation_id
  INNER JOIN failed_variation fv ON a.variation_id = fv.variation_id
  WHERE vf.flip = 1
  AND vf.map_weight = 1;
  }, {mysql_use_result => 1});

  $sth->execute();
  while (my $row = $sth->fetchrow_arrayref) {
    my ($variation_id, $failed_description_id) = @$row;
    $failed_variations->{$variation_id}->{$failed_description_id} = 1;
  }
  $sth->finish();
  return $failed_variations;
}

sub flip_alleles {
  my $self = shift;
  my $failed_alleles = shift;
  my $failed_variations = shift;
  my $dumped_features_file = shift;
  my $fh_out = shift;
  my $fh_err = shift;
  my $vdba = $self->get_newasm_variation_database_connection;
  my $allele_adaptor = $vdba->get_adaptor('Allele');
  
  my %allele_id_2_string = reverse %{$allele_adaptor->_cache_allele_codes};
  my $fh = FileHandle->new($dumped_features_file, 'r');

  while (<$fh>) {
    chomp;
    my ($seq_region_id, $start, $end, $strand, $variation_id, $variation_name, $allele_string, $map_weight, $allele_id, $allele_code_id) = split/\t/;
    # next if failed_allele
    next if ($failed_alleles->{$variation_id}->{$allele_id}); # don't bother with previously failed alleles
    next if ($failed_variations->{$variation_id}); # don't bother with previousely failed variations
    my $allele = $allele_id_2_string{$allele_code_id};
    next if (contained_in_allele_string($allele_string, $allele));
    my $rev_comp_allele = reverse_comp_allele_string($allele);
    if (contained_in_allele_string($allele_string, $rev_comp_allele)) {
      my $allele_code_id = $allele_adaptor->_allele_code($rev_comp_allele);
      print $fh_out "UPDATE allele set allele_code_id = $allele_code_id WHERE allele_id = $allele_id;\n";
    }
  }
  $fh->close;
}

sub flip_population_genotypes {
  my $self = shift;
  my $dumped_features_file = shift;
  my $fh_out = shift;
  my $fh_err = shift;
  my $vdba = $self->get_newasm_variation_database_connection;
  my $sgta = $vdba->get_SampleGenotypeAdaptor;
  my $gtca = $vdba->get_GenotypeCodeAdaptor;
  my $gtcs = $gtca->fetch_all;
  my $gtc_id_2_string = {};
  foreach my $gtc (@$gtcs) {
    my $gtc_dbID = $gtc->dbID;
    my $alleles = join('/', @{$gtc->genotype});
    $gtc_id_2_string->{$gtc_dbID} = $alleles;
  }
  my $fh = FileHandle->new("$dumped_features_file", 'r');
  while (<$fh>) {
    chomp;
    my ($seq_region_id, $start, $end, $strand, $variation_id, $variation_name, $allele_string, $map_weight, $population_genotype_id, $genotype_code) = split/\t/;
    my $genotype_string = $gtc_id_2_string->{$genotype_code};
    next if (contained_in_allele_string($allele_string, $genotype_string));
    my $rev_comp_genotype_string = reverse_comp_allele_string($genotype_string);
    if (contained_in_allele_string($allele_string, $rev_comp_genotype_string)) {
      my $genotype_code_id = $sgta->_genotype_code([split('/', $rev_comp_genotype_string)]);
      print $fh_out "UPDATE population_genotype SET genotype_code_id = $genotype_code_id WHERE population_genotype_id = $population_genotype_id;\n ";
    }
  }
  $fh->close();
}

sub flip_sample_genotypes {
  my $self = shift;
  my $dumped_features_file = shift;
  my $fh_out = shift;
  my $fh_err = shift;
  my $fh = FileHandle->new($dumped_features_file, 'r');
  while (<$fh>) {
    chomp;
    my ($seq_region_id, $start, $end, $strand, $variation_id, $variation_name, $allele_string, $map_weight, $subsnp_id, $allele_1, $allele_2, $sample_id) = split/\t/;
    # deal with N alleles
    my @alleles = ();
    foreach my $allele ($allele_1, $allele_2) {
      if ($allele eq 'N') {
        push @alleles, $allele;
      } else {
        if (contained_in_allele_string($allele_string, $allele)) {
          push @alleles, $allele;
        } else {
          my $rev_comp_allele = reverse_comp_allele_string($allele);
          if (contained_in_allele_string($allele_string, $rev_comp_allele)) {
            push @alleles, $rev_comp_allele;
          }
        }
      }
    }
    if (scalar @alleles == 2 ) {
      @alleles = sort @alleles;
      my $new_allele_1 = $alleles[0];
      my $new_allele_2 = $alleles[1];
      print $fh_out "UPDATE tmp_sample_genotype_single_bp SET allele_1 = '$new_allele_1' WHERE variation_id = $variation_id AND subsnp_id = $subsnp_id AND sample_id = $sample_id;\n";
      print $fh_out "UPDATE tmp_sample_genotype_single_bp SET allele_2 = '$new_allele_2' WHERE variation_id = $variation_id AND subsnp_id = $subsnp_id AND sample_id = $sample_id;\n";
    } else {
      print $fh_err "not enough alleles for $variation_name ", join("/", @alleles), " old alleles $allele_1 $allele_2\n";
    }
  }
  $fh->close();
}

sub reverse_comp_allele_string {
  my $allele_string = shift;
  my @allele_string_rev_comp = split('/', $allele_string);
  foreach my $allele (@allele_string_rev_comp) {
    reverse_comp(\$allele);
  }
  return join('/', sort @allele_string_rev_comp);
}

sub contained_in_allele_string {
  my $allele_string = shift;
  my $alleles = shift;
  my %lookup = map {$_ => 1} split/\//, $allele_string;
  foreach my $allele (split/\//, $alleles) {
    if (!$lookup{$allele}) {
      return 0;
    }
  }
  return 1;
}

sub dump_features_for_flipping {
  my $self = shift;
  my $query = shift;
  my $file_name = shift;
  my $vdba = $self->get_newasm_variation_database_connection;
  my $dbh = $vdba->dbc->db_handle;
  my $fh = FileHandle->new($file_name, 'w');
  my $sth = $dbh->prepare($query, {mysql_use_result => 1});
  $sth->execute();
  while (my $row = $sth->fetchrow_arrayref) {
    my @values = map { defined $_ ? $_ : '\N' } @$row;
    print $fh join("\t", @values), "\n";
  }
  $sth->finish();
  $fh->close();
}

sub column_exists {
  my $self = shift;
  my $vdba = shift;
  my $feature_table = shift;
  my $column = shift;
  my $dbh = $vdba->dbc->db_handle;
  my $dbname = $vdba->dbc->dbname();
  my $sth = $dbh->prepare(qq{
      SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS
      WHERE TABLE_SCHEMA = '$dbname'
      AND TABLE_NAME = '$feature_table';
      });
  $sth->execute();

  my @column_names = ();
  while (my @name = $sth->fetchrow_array) {
    if ($name[0] eq $column) {
      return 1;
    }
  }
  return 0;
}

sub failed_variations_after_remapping {
  my $self = shift;
  my $qc_failure_reasons_dir = shift;
  my $failed_variants = {};
  opendir(DIR, $qc_failure_reasons_dir) or die $!;
  while (my $file = readdir(DIR)) {
    next unless (-f "$qc_failure_reasons_dir/$file");
    next unless ($file =~ m/\.txt$/);
    my $fh = FileHandle->new("$qc_failure_reasons_dir/$file", 'r');
    while(<$fh>) {
      chomp;
      my ($variation_id, $failed_description_id) = split/\t/, $_;
      $failed_variants->{$variation_id}->{$failed_description_id} = 1;
    }
    $fh->close;
  }
  closedir(DIR);
  return $failed_variants;
}

sub previous_unmapped_variants {
  my $self = shift;
  my $vdba = $self->get_oldasm_variation_database_connection;
  my $dbh = $vdba->dbc->db_handle;
  my $previous_unmapped_variants = {};
  my $sth = $dbh->prepare(qq{
    SELECT variation_id FROM failed_variation WHERE failed_description_id=5;
  }, {mysql_use_result => 1});
  $sth->execute();
  while (my $row = $sth->fetchrow_arrayref) {
    $previous_unmapped_variants->{$row->[0]} = 1;
  }
  $sth->finish();
  return $previous_unmapped_variants;
}

sub unmapped_variants {
  my $self = shift;
  my $feature_table_mapping_results = shift;
  my $vdba = $self->get_newasm_variation_database_connection;
  my $dbh = $vdba->dbc->db_handle;

  my $unmapped_variants = {};
  my $sth = $dbh->prepare(qq{
  SELECT v.variation_id
  FROM variation v
  LEFT JOIN $feature_table_mapping_results ftmr ON v.variation_id = ftmr.variation_id
  WHERE ftmr.variation_id IS NULL;
  }, {mysql_use_result => 1});
  $sth->execute();
  while (my $row = $sth->fetchrow_arrayref) {
    $unmapped_variants->{$row->[0]} = 1;
  }
  $sth->finish();

  return $unmapped_variants; 
}

sub cleanup_failed_variants {
  my $self = shift;
  my $feature_table = shift;
    
  my $vdba = $self->get_newasm_variation_database_connection;
  my $dbh = $vdba->dbc->db_handle;

  $dbh->do("DELETE fv FROM failed_variation as fv LEFT JOIN $feature_table ft ON fv.variation_id = ft.variation_id AND ft.variation_id IS NOT NULL ;") or die $!;
}

sub load_failed_variants {
  my $self = shift;
  my $failed_variants = shift;

  my $vdba = $self->get_newasm_variation_database_connection;
  my $dbh = $vdba->dbc->db_handle;

  foreach my $variation_id (keys %$failed_variants) {
    foreach my $failed_description_id (keys %{$failed_variants->{$variation_id}}) {
      $dbh->do("INSERT INTO failed_variation(variation_id, failed_description_id) VALUES($variation_id, $failed_description_id);") or die $!;
    }
  }
}

sub update_variation_set_variation {
  my $self = shift;
  my $vdba = $self->get_newasm_variation_database_connection;
  my $dbh = $vdba->dbc->db_handle;

  my $sth = $dbh->prepare(qq{
    SELECT variation_set_id FROM variation_set WHERE name='All failed variations';
  }, {mysql_use_result => 1});
  $sth->execute();
  my $row = $sth->fetchrow_arrayref;
  $sth->finish();

  my $variation_set_id = $row->[0];
  die "Could not fetch variation_set_id from All failed variations" if (!$variation_set_id);

  $dbh->do("DELETE FROM variation_set_variation WHERE variation_set_id=$variation_set_id;") or die $!;
  $dbh->do("INSERT INTO variation_set_variation(variation_id, variation_set_id) SELECT DISTINCT variation_id, $variation_set_id FROM failed_variation;") or die $!;
}

sub update_display {
  my $self = shift;
  my $feature_table = shift;
  my $vdba = $self->get_newasm_variation_database_connection;
  my $dbh = $vdba->dbc->db_handle;

  $dbh->do("UPDATE variation SET display=1;") or die $!;
  $dbh->do("UPDATE variation v, failed_variation fv SET v.display = 0 WHERE fv.variation_id = v.variation_id;") or die $!;
  $dbh->do("UPDATE variation v, variation_citation vc SET v.display = 1 WHERE vc.variation_id = v.variation_id;") or die $!;
  $dbh->do("UPDATE variation v, $feature_table vf SET vf.display = v.display WHERE v.variation_id = vf.variation_id;") or die $!;
}

sub cleanup_mapped_feature_table {
  my $self = shift;
  my $vdba = $self->get_newasm_variation_database_connection;
  my $dbh = $vdba->dbc->db_handle;

  my $feature_table = $self->param('feature_table');
  my $feature_table_mapping_results = $self->param('feature_table') . '_mapping_results';

  $dbh->do("RENAME TABLE $feature_table to before_remapping_$feature_table;") or die $!;
  $dbh->do("RENAME TABLE $feature_table_mapping_results to $feature_table;") or die $!;
  foreach my $column (qw/flip variation_feature_id_old/) {
    if ($self->column_exists($vdba, $feature_table, $column)) {
      $dbh->do("ALTER TABLE $feature_table DROP COLUMN $column;")
    }
  }
}

1;
