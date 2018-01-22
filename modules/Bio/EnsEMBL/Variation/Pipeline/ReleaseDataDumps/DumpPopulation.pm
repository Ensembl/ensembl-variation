=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::DumpPopulation;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');
use FileHandle;

sub run {
  my $self = shift;
  my $file_type = $self->param('file_type');
  if ($file_type eq 'gvf') {
    $self->dump_gvf;
  }
  if ($file_type eq 'vcf') {
    $self->dump_vcf;
  }
}

sub dump_gvf {
  my $self = shift;

  my $human_gvf_file = $self->param('human_gvf_file');
  my $variation_ids_dir = $self->param('variation_ids_dir');
  my $prefetched_frequencies_dir = $self->param('prefetched_frequencies'); 

  my $file_name = $self->param('file_name');
  my $population_gvf_file = $self->param('population_gvf_file');
  my $population_name = $self->param('population_name');
  my $short_name = $self->param('short_name');
  my $group = $self->param('group'); # HAPMAP, 1000G, ESP

  my $add_to_header = "##population $population_name\n";
  my $added_name = 0;

  my $prev_seq_id = -1;
  my $cache = {};
  my $fh_gvf = FileHandle->new($human_gvf_file, 'r');
  my $fh_out = FileHandle->new($population_gvf_file, 'w');
  while (<$fh_gvf>) {
    chomp;
    if (/^#/) {
      if (/^##sequence-region/) {
        unless ($added_name) {
          print $fh_out $add_to_header;
          $added_name = 1;
        }   
        print $fh_out $_, "\n";
      } else {
        print $fh_out $_, "\n";
      }
    } else {
      my $line = $_;
      my $gvf_line = get_gvf_line($line);
      my $seq_id = $gvf_line->{seq_id};
      unless ("$seq_id" eq "$prev_seq_id") {
        $cache = update_cache("$variation_ids_dir/$seq_id.txt", "$prefetched_frequencies_dir/$group/$short_name.txt"  );
        $prev_seq_id = $seq_id;
      }
      my $allele_string = $gvf_line->{attributes}->{allele_string};
      my $variation_id  = $gvf_line->{attributes}->{variation_id};
      next unless ($cache->{$variation_id});

      my @vf_alleles = split(',', $allele_string);
      my @alleles = ();
      my @freqs = ();
      if (my $allele_hash = $cache->{$variation_id}) {
        for my $allele (keys %$allele_hash) {
          my $freq = $allele_hash->{$allele};
          if ($allele eq $vf_alleles[0]) {
            $gvf_line->{attributes}->{Reference_seq} = $allele;
            if ($freq == 1) {
              if ((scalar keys %$allele_hash) == 1) {
                push @alleles, $vf_alleles[1];
                push @freqs, '0.0';
              }
            }
          } else {
            push @alleles, $allele;
            push @freqs, $freq;
          }
        }
      }

      if ((scalar @alleles) == 0) {
        $self->warning("No Variant_seq for $variation_id");
      }
      $gvf_line->{attributes}->{Variant_seq} = join(',', @alleles);
      if (grep {$_ =~ /\d/} @freqs) {
        $gvf_line->{attributes}->{Variant_freq} = join(',', @freqs);
      }
      delete $gvf_line->{attributes}->{variation_id};
      delete $gvf_line->{attributes}->{allele_string};
      delete $gvf_line->{attributes}->{global_minor_allele_frequency};

      $line = join("\t", map {$gvf_line->{$_}} (
        'seq_id',
        'source',
        'type',
        'start',
        'end',
        'score',
        'strand',
        'phase'));
      my $attributes = join(";", map{"$_=$gvf_line->{attributes}->{$_}"} keys %{$gvf_line->{attributes}});
      print $fh_out $line, "\t", $attributes, "\n";
    }
  } # while

  system("gzip $population_gvf_file");
}

sub get_gvf_line {
  my $line = shift;
  my $gvf_line = {};
  my @header_names = qw/seq_id source type start end score strand phase/;
  my @header_values = split(/\t/, $line);
  my $attrib = pop @header_values;

  for my $i (0 .. $#header_names) {
    $gvf_line->{$header_names[$i]} = $header_values[$i];
  }

  my @attributes = split(';', $attrib);
  foreach my $attribute (@attributes) {
    my ($key, $value) = split('=', $attribute);
    if ($value) {
      $gvf_line->{attributes}->{$key} = $value;
    }
  }
  return $gvf_line;
}

sub update_cache {
  my $variation_ids_file = shift;
  my $cache_file = shift;
  my $fh_var_ids = FileHandle->new($variation_ids_file, 'r');
  my $variation_ids = {};
  while (<$fh_var_ids>) {
    chomp;
    $variation_ids->{$_} = 1;
  }
  $fh_var_ids->close();
  my $cache = {};
  my $fh = FileHandle->new($cache_file, 'r');
  while (<$fh>) {
    chomp;
    my ($variation_id, $allele, $frequency) = split(/\t/);
    if ($variation_ids->{$variation_id}) {
      $frequency ||= 0;
      $cache->{$variation_id}->{$allele} = $frequency;
    }
  }
  $fh->close();
  return $cache;
}

1;
