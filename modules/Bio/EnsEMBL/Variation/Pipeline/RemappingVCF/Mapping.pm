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

use strict;
use warnings;

package Bio::EnsEMBL::Variation::Pipeline::RemappingVCF::Mapping;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

our @ISA = ('Bio::EnsEMBL::Storable');

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($vcf_variation_id,
      $seq_region_name_old,
      $seq_region_id_old,
      $seq_region_start_old,
      $seq_region_start_padded_old,
      $allele_string_old,
      $allele_string_padded_old,
      $vcf_id,
      $variation_id,
      $seq_region_name_new,
      $seq_region_id_new,
      $seq_region_start_new,
      $seq_region_start_padded_new,
      $allele_string_new,
      $allele_string_padded_new ) = 
    rearrange([
      'vcf_variation_id',
      'seq_region_name_old',
      'seq_region_id_old',
      'seq_region_start_old',
      'seq_region_start_padded_old',
      'allele_string_old',
      'allele_string_padded_old',
      'vcf_id',
      'variation_id',
      'seq_region_name_new',
      'seq_region_id_new',
      'seq_region_start_new',
      'seq_region_start_padded_new',
      'allele_string_new',
      'allele_string_padded_new'
    ], @_);
  my $self = bless {
      vcf_variation_id => $vcf_variation_id,
      seq_region_name_old => $seq_region_name_old,
      seq_region_id_old => $seq_region_id_old,
      seq_region_start_old => $seq_region_start_old,
      seq_region_start_padded_old => $seq_region_start_padded_old,
      allele_string_old => $allele_string_old,
      allele_string_padded_old => $allele_string_padded_old,
      vcf_id => $vcf_id,
      variation_id => $variation_id,
      seq_region_name_new => $seq_region_name_new,
      seq_region_id_new => $seq_region_id_new,
      seq_region_start_new => $seq_region_start_new,
      seq_region_start_padded_new => $seq_region_start_padded_new,
      allele_string_new => $allele_string_new,
      allele_string_padded_new => $allele_string_padded_new, 
  }, $class;

  return $self;
}

sub vcf_variation_id {
  my $self = shift;
  $self->{vcf_variation_id} = shift if @_;
  return $self->{vcf_variation_id};
}

sub seq_region_name_old {
  my $self = shift;
  $self->{seq_region_name_old} = shift if @_;
  return $self->{seq_region_name_old};
}

sub seq_region_name_new {
  my $self = shift;
  $self->{seq_region_name_new} = shift if @_;
  return $self->{seq_region_name_new};
}

sub seq_region_id_old {
  my $self = shift;
  $self->{seq_region_id_old} = shift if @_;
  return $self->{seq_region_id_old};
}

sub seq_region_id_new {
  my $self = shift;
  $self->{seq_region_id_new} = shift if @_;
  return $self->{seq_region_id_new};
}

sub seq_region_start_old {
  my $self = shift;
  $self->{seq_region_start_old} = shift if @_;
  return $self->{seq_region_start_old};
}

sub seq_region_start_new {
  my $self = shift;
  $self->{seq_region_start_new} = shift if @_;
  return $self->{seq_region_start_new};
}

sub seq_region_start_padded_old {
  my $self = shift;
  $self->{seq_region_start_padded_old} = shift if @_;
  return $self->{seq_region_start_padded_old};
}

sub seq_region_start_padded_new {
  my $self = shift;
  $self->{seq_region_start_padded_new} = shift if @_;
  return $self->{seq_region_start_padded_new};
}

sub allele_string_old {
  my $self = shift;
  $self->{allele_string_old} = shift if @_;
  return $self->{allele_string_old};
}

sub allele_string_new {
  my $self = shift;
  $self->{allele_string_new} = shift if @_;
  return $self->{allele_string_new};
}

sub allele_string_padded_old {
  my $self = shift;
  $self->{allele_string_padded_old} = shift if @_;
  return $self->{allele_string_padded_old};
}

sub allele_string_padded_new {
  my $self = shift;
  $self->{allele_string_padded_new} = shift if @_;
  return $self->{allele_string_padded_new};
}

sub is_indel {
  my $self = shift;
  return defined($self->allele_string_padded_old);
}

sub variation_id {
  my $self = shift;
  $self->{variation_id} = shift if @_;
  return $self->{variation_id};
}

sub alleles_are_equal {
  my $self = shift;
  $self->{alleles_are_equal} = shift if @_;
  return $self->{alleles_are_equal};
}

=begin
CREATE TABLE `vcf_variation` (
  `vcf_variation_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_name_old` varchar(50) DEFAULT NULL,
   seq_region_id_old int(10) unsigned DEFAULT NULL, 
  `seq_region_start_old` int(11) DEFAULT NULL,
   seq_region_start_padded_old int(11) DEFAULT NULL,
  `allele_string_old` varchar(50) DEFAULT NULL,
   allele_string_padded_old varchar(50) DEFAULT NULL,
  `vcf_id` int(10) unsigned DEFAULT NULL,
  `variation_id` int(10) unsigned DEFAULT NULL,
  `seq_region_name_new` varchar(50) DEFAULT NULL,
   seq_region_id_new int(10) unsigned DEFAULT NULL,
  `seq_region_start_new` int(11) DEFAULT NULL,
   seq_region_start_padded_new int(11) DEFAULT NULL,
  `allele_string_new` varchar(50) DEFAULT NULL,
  allele_string_padded_new varchar(50) DEFAULT NULL,
  PRIMARY KEY (`vcf_variation_id`),
  KEY `variation_idx` (`variation_id`),
  KEY `seq_region_name_old_idx` (`seq_region_name_old`),
  KEY `seq_region_id_old_idx` (`seq_region_id_old`),
  KEY `seq_region_start_old_idx` (`seq_region_start_old`)
)
=end
=cut


1;
