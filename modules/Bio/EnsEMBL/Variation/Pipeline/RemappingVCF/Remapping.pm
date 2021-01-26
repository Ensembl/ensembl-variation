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
 developers list at <https://lists.ensembl.org/mailman/listinfo/dev>.
 Questions may also be sent to the Ensembl help desk at
 <https://www.ensembl.org/Help/Contact>.
=cut
package Bio::EnsEMBL::Variation::Pipeline::RemappingVCF::Remapping;

use strict;
use warnings;

use FileHandle;
use Bio::EnsEMBL::IO::Parser::VCF4Tabix;
use Bio::EnsEMBL::IO::Parser::VCF4;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(trim_sequences);
use Array::Utils qw(:all);
use Bio::EnsEMBL::Variation::Pipeline::RemappingVCF::Mapping;

use base ('Bio::EnsEMBL::Hive::Process');

use constant EQUAL => 1;
use constant EQUAL_AFTER_REV_COMP => 2;
use constant NOT_EQUAL => 0;

sub run {
  my $self = shift;
  my $species = $self->param('species');
  my $chrom = $self->param('chrom');
  my $vcf_chrom = $self->param('vcf_chrom');
  my $seq_region_end = $self->param('seq_region_end');
  my $vcf_file = $self->param('vcf_file');
  my $pipeline_dir = $self->param('pipeline_dir');
  my $registry_file = $self->param('registry_file_newasm');
  my $population = $self->param('population');
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_all($registry_file);
  my $sa = $registry->get_adaptor($species, 'core', 'slice');
  my $dbh = $registry->get_DBAdaptor($species, 'variation')->dbc->db_handle;

  my $fh_out = FileHandle->new("$pipeline_dir/$population.chrom$chrom.new_assembly", 'w');
  my $fh_no_mapping = FileHandle->new("$pipeline_dir/$population.chrom$chrom.no_mapping", 'w');
  my $fh_err = FileHandle->new("$pipeline_dir/$population.chrom$chrom.err", 'w');

  my $parser = Bio::EnsEMBL::IO::Parser::VCF4Tabix->open($vcf_file);
  $parser->seek($vcf_chrom, 1, $seq_region_end);
  while ($parser->next) {
    my $reference = $parser->get_reference;
    my @alternatives = @{$parser->get_alternatives};
    my $allele_string = join('/', $reference, @alternatives);
    my $seq_name = $parser->get_seqname;
    $seq_name =~ s/chr//;
    my $raw_start = $parser->get_raw_start;
    my $start = $parser->get_start; # or raw start?
    my $raw_IDs = $parser->get_raw_IDs;

    my $mappings = $self->update_mappings($dbh, $seq_name, $raw_start, $allele_string);
    my $number_of_mappings = scalar @$mappings;

    if ( $number_of_mappings == 0) {
      print $fh_no_mapping "No mappings for $seq_name $raw_start $raw_IDs count $number_of_mappings\n";
    } else {
      my $qual = $parser->get_raw_score;
      my $filter = $parser->get_raw_filter_results;
      my $info = 'has been removed';
      my $sample_genotypes = $parser->get_samples_genotypes;
      my $samples = $parser->get_samples;
      my $format = 'GT';
      foreach my $mapping (@$mappings) {
        my $is_mapping = compare_allele_strings($mapping);
        if ($is_mapping) {
          my $mapped_sample_genotypes = $self->map_sample_genotypes($sample_genotypes, $mapping, $is_mapping);
          my $vcf_line = $self->get_vcf_line($mapping, $raw_IDs, $qual, $filter, $info, $format, $samples, $mapped_sample_genotypes, $population, $sa);
          print $fh_out "$vcf_line\n";
        } else {
          my $allele_string_old = $mapping->allele_string_old;
          my $allele_string_new = $mapping->allele_string_new;
          print $fh_err "$seq_name $raw_start $raw_IDs $is_mapping $allele_string_old $allele_string_new\n";
        }
      }
    }
  }
  $fh_out->close;
  $fh_err->close;
  $fh_no_mapping->close;
}

sub compare_allele_strings {
  my $mapping = shift;

  my $allele_string_old = $mapping->allele_string_old;
  my $allele_string_new = $mapping->allele_string_new;
  my $allele_string_padded_old = $mapping->allele_string_padded_old;
  my $allele_string_padded_new = $mapping->allele_string_padded_new;

  my $sorted_allele_string_old = sort_allele_string($allele_string_old);
  my $sorted_allele_string_new = sort_allele_string($allele_string_new);

  if ($sorted_allele_string_old eq $sorted_allele_string_new) {
    $mapping->alleles_are_equal(1);
    return EQUAL;
  }
  my $rev_comp_sorted_allele_string_old = rev_comp_allele_string($sorted_allele_string_old);
  if ($rev_comp_sorted_allele_string_old eq $sorted_allele_string_new) {
    $mapping->alleles_are_equal(1);
    return EQUAL_AFTER_REV_COMP;
  }

  if (first_allele_string_contained_in_second_allele_string($allele_string_new, $allele_string_old)) {
    return EQUAL;
  } 

  if (first_allele_string_contained_in_second_allele_string($allele_string_old, $allele_string_new)) {
    return EQUAL;
  } 

  if (first_allele_string_contained_in_second_allele_string($allele_string_new, $rev_comp_sorted_allele_string_old)) {
    return EQUAL_AFTER_REV_COMP;
  } 

  if (first_allele_string_contained_in_second_allele_string($rev_comp_sorted_allele_string_old, $allele_string_new)) {
    return EQUAL_AFTER_REV_COMP;
  } 

  return NOT_EQUAL;
}

sub first_allele_string_contained_in_second_allele_string {
  my $first_allele_string = shift;
  my $second_allele_string = shift;
  my @first_array = split('/', $first_allele_string);
  my @second_array = split('/', $second_allele_string);

  # get items from array @a that are not in array @b
  # my @minus = array_minus( @a, @b );
  my @not_in_second_array = array_minus( @first_array, @second_array);
  return !@not_in_second_array;
}

sub rev_comp_allele_string {
  my $allele_string = shift;
  my @rev_comp_alleles = split('/', $allele_string);
  foreach my $allele (@rev_comp_alleles) {
    reverse_comp(\$allele);
  }
  return join('/', @rev_comp_alleles);
}

sub sort_allele_string {
  my $allele_string = shift;
  return join('/', sort split('/', $allele_string));
}

sub trimmed_alleles {
  my $allele_string = shift;
  my @alleles = split('/', $allele_string);
  my $ref = shift @alleles; 
  my @trimmed_alleles = ();
  my $trimmed_ref = '';
  foreach my $allele (@alleles) {

    my ($_ref, $_alt, $_start) = @{trim_sequences($ref, $allele, 1, 1, 1)};
    $trimmed_ref = $_ref;
    push @trimmed_alleles, $_alt;
  }
  unshift @trimmed_alleles, $trimmed_ref;
  return join('/', @trimmed_alleles);
}

sub get_vcf_line {
  my ($self, $mapping, $rawIDs, $qual, $filter, $info, $format, $samples, $mapped_sample_genotypes, $population, $sa) = @_;

  my $chrom = $mapping->seq_region_name_new;
  my $location = $mapping->seq_region_start_new;
  my $new_allele_string = $mapping->allele_string_new;
  my @alleles = split('/', $new_allele_string);

  if ($mapping->is_indel) {
    @alleles = split('/', $new_allele_string); 
    $location--;
    my $slice = $sa->fetch_by_toplevel_location("$chrom:$location-$location");
    my $prev_base = $slice->seq() || 'N';
    for my $i (0..$#alleles) {
      $alleles[$i] =~ s/\-//g;
      $alleles[$i] = $prev_base . $alleles[$i];
    }
  }

  my @default_alleles = ();
  foreach my $allele (@alleles) {
    push @default_alleles, '.';
  }
  my $default_gt = join('/', @default_alleles);

  my $ref_allele = shift @alleles;
  my $alt_alleles = join(',', @alleles);

  my @genotypes = ();

  foreach my $sample (@$samples) {
    my $gt = $mapped_sample_genotypes->{$sample};

    # deal with ./. genotypes!
    if (! defined $gt) {
      push @genotypes, $default_gt;
    } else {
      push @genotypes, $mapped_sample_genotypes->{$sample};
    }
  }
  #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT genotypes
  $info = allele_frequency($new_allele_string, $mapped_sample_genotypes, $population); 
  if (!$info) {
    $info = 'No INFO';
  }
  my $line = join("\t", $chrom, $location, $rawIDs, $ref_allele, $alt_alleles, $qual, $filter, $info, $format, @genotypes);
  return $line;
}

sub allele_frequency {
  my $allele_string = shift;
  my $sample_genotypes = shift;
  my $population = shift;
  my $allele_map = get_allele_map($allele_string);
  ##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
  ###INFO=<ID=AF,Number=A,Type=Float,Description="Estimated Allele Frequencies">
  ###INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
  ###INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">

  my $counts = {};
  my $sample_count = 0;
  my $total_number_of_alleles = 0;
  foreach my $sample (keys %$sample_genotypes) {
    $sample_count++;
    foreach my $index (split(/\||\//, $sample_genotypes->{$sample})) {
      $total_number_of_alleles++;
      $counts->{$index}++;  
    }
  } 

  my @alleles = split('/', $allele_string);
  my @ac = ();
  my @af = ();
  my $ref = shift @alleles;
  if ($total_number_of_alleles == 0) {
    return undef;
  }
  foreach my $allele (@alleles) {
    my $allele_count = $counts->{$allele_map->{$allele}} // 0;
    my $frequency = sprintf("%.4f", $allele_count / $total_number_of_alleles);  
    push @af, $frequency;
    push @ac, $allele_count;
  }
  return "AN\_$population=$total_number_of_alleles;AF\_$population=" . join(',', @af) . ";AC\_$population=" . join(',', @ac) . ";NS\_$population=$sample_count";
}

sub map_sample_genotypes {
  my $self = shift;
  my $sample_genotypes_old = shift;
  my $mapping = shift;
  my $is_mapping = shift; # rev_comp?

  my $rev_comp = 0;
  if ($is_mapping == EQUAL_AFTER_REV_COMP) {
    $rev_comp = 1;
  }
  my $sample_genotypes_new = {};

  # if alleles from old allele string (from VCF file) are missing from new allele string (from database) add missing alleles to end of new allele string  
  if (!$mapping->alleles_are_equal) {
    $self->add_missing_alleles($mapping, $rev_comp);
  }

  # This will be new entry in VCF file. After remapping we know correct allele order for the new allele string
  my $new_allele_map = get_allele_map($mapping->allele_string_new);
  my $old_rev_comp_map = {};

  if ($rev_comp) {
    $old_rev_comp_map = get_rev_comp_allele_map($mapping->allele_string_old); 
  } 

  foreach my $sample (keys %$sample_genotypes_old) {
    my $old_gt =  $sample_genotypes_old->{$sample};
    my $is_phased = ($old_gt =~ m/\|/);
    my $seperator = ($is_phased) ? '|' : '/';
    my @alleles = split(/\||\//, $old_gt);
    my @new_gt = ();
    foreach my $allele (@alleles) {
      my $allele_index = undef;
      if ($rev_comp) {
        $allele_index = $new_allele_map->{$old_rev_comp_map->{$allele}};
      } else {
        $allele_index = $new_allele_map->{$allele};
      }
      if (!defined $allele_index) {
        my $new_start = $mapping->seq_region_start_new;
        my $new_allele_string = $mapping->allele_string_new;
        my $old_allele_string = $mapping->allele_string_old;
        my $seq_region_name_new = $mapping->seq_region_name_new;
      }
      push @new_gt, $allele_index;

    }
    $sample_genotypes_new->{$sample} = join($seperator, @new_gt);
  }
  return $sample_genotypes_new;
}

sub add_missing_alleles {
  my $self = shift;
  my $mapping = shift;
  my $rev_comp = shift;

  my $new_allele_string = $mapping->allele_string_new;
  my $old_allele_string = $mapping->allele_string_old;

  my $new_allele_map = get_allele_map($new_allele_string);
  my @new_alleles = split('/', $new_allele_string); 
  my @old_alleles = split('/', $old_allele_string);
  foreach my $allele (@old_alleles) {
    if (! defined $new_allele_map->{$allele}) {
      push @new_alleles, $allele;
    }
  }
  $new_allele_string = join('/', @new_alleles);
  my $variation_id = $mapping->variation_id;
  my $new_allele_string_before_change = $mapping->allele_string_new;
  $mapping->allele_string_new($new_allele_string);
}

# allele to index
sub get_allele_map {
  my $allele_string = shift;
  my $allele_map = {};
  my @alleles = split('/', $allele_string);
  for (my $i = 0 ; $i < scalar(@alleles); $i++) {
    $allele_map->{$alleles[$i]} = $i;
  }
  return $allele_map;
}
# allele to rev_allele
sub get_rev_comp_allele_map {
  my $allele_string = shift;
  my @alleles = split('/', $allele_string);
  my @rev_comp = split('/', $allele_string);
  foreach my $allele (@rev_comp) {
    reverse_comp(\$allele);
  }
  my $allele_map = {};
  for my $i (0 .. $#alleles) {
    $allele_map->{$alleles[$i]} = $rev_comp[$i];
  }
  return $allele_map;
}

sub update_mappings {
  my $self = shift;
  my $dbh = shift;
  my $seq_region_name = shift;
  my $raw_start = shift;
  my $allele_string = shift;
  my @mappings = ();

  my $sth = $dbh->prepare(qq{
    SELECT vcf_variation_id, seq_region_name_old, seq_region_id_old, seq_region_start_old, seq_region_start_padded_old, allele_string_old, allele_string_padded_old, vcf_id, variation_id, seq_region_name_new, seq_region_id_new, seq_region_start_new, allele_string_new 
    FROM vcf_variation
    WHERE seq_region_name_old = ?
    AND seq_region_start_old = ?
    AND allele_string_old = ?
    AND variation_id IS NOT NULL
    AND seq_region_name_new IS NOT NULL;
  }, {mysql_use_result => 1});

  $sth->execute($seq_region_name, $raw_start, $allele_string);
  my ($vcf_variation_id, $seq_region_name_old, $seq_region_id_old, $seq_region_start_old, $seq_region_start_padded_old, $allele_string_old, $allele_string_padded_old, $vcf_id, $variation_id, $seq_region_name_new, $seq_region_id_new, $seq_region_start_new, $allele_string_new); 
  $sth->bind_columns(\($vcf_variation_id, $seq_region_name_old, $seq_region_id_old, $seq_region_start_old, $seq_region_start_padded_old, $allele_string_old, $allele_string_padded_old, $vcf_id, $variation_id, $seq_region_name_new, $seq_region_id_new, $seq_region_start_new, $allele_string_new));
  while ($sth->fetch) {
    
    my $mapping = Bio::EnsEMBL::Variation::Pipeline::RemappingVCF::Mapping->new(
      -vcf_variation_id => $vcf_variation_id,
      -seq_region_name_old => $seq_region_name_old,
      -seq_region_id_old => $seq_region_id_old,
      -seq_region_start_old => $seq_region_start_old,
      -seq_region_start_padded_old => $seq_region_start_padded_old || undef,
      -allele_string_old => $allele_string_old,
      -allele_string_padded_old => $allele_string_padded_old || undef,
      -vcf_id => $vcf_id,
      -variation_id => $variation_id,
      -seq_region_name_new => $seq_region_name_new,
      -seq_region_id_new => $seq_region_id_new,
      -seq_region_start_new => $seq_region_start_new,
      -allele_string_new => $allele_string_new,
    );
    push @mappings, $mapping;
  }
  $sth->finish;
  my $count = scalar @mappings;
  return \@mappings;
}

1;
