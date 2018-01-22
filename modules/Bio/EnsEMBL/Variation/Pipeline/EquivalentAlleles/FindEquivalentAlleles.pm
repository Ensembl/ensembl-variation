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



=head1 NAME Bio::EnsEMBL::Variation::Pipeline::EquivalentAlleles::FindEquivalentAlleles

=head1 DESCRIPTION

Find Variants with equivalent alleles by checking for identical HGVS notation

=cut

package Bio::EnsEMBL::Variation::Pipeline::EquivalentAlleles::FindEquivalentAlleles;


use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub run {

  my $self = shift;

  my $core_dba = $self->get_species_adaptor('core');
  my $var_dba  = $self->get_species_adaptor('variation');

  my $sa       = $core_dba->get_SliceAdaptor();
  my $vfa      = $var_dba->get_VariationFeatureAdaptor();


  my $location = $self->required_param('location');
$self->warning( 'starting location '. $location);
  my $slice    = $sa->fetch_by_location($location, 'toplevel');
  die "No slice found for $location\n" unless defined $slice;  

  my $length = $slice->length();
  my ($seqname, $start, $end)  = split/\:|\-/, $location ;
  my $size   = $self->required_param('bin_size'); 

  my %hgvs; ## store variant name by hgvs genomic description  

  while($start < $end){

    my $bin_end = $start + $size;
    my $constraint = " seq_region_start > $start and seq_region_start < $bin_end "; 
    my $vfs = $vfa->fetch_all_by_Slice_constraint($slice, $constraint);

    foreach my $vf(@{$vfs}){
      ## SNVs don't 3' shift
      next if $vf->class_SO_term() eq 'SNV';

      my $hgvsg = $vf->hgvs_genomic($slice); 

      my $alleles = $vf->allele_string();
      my @al = split/\//, $alleles;

      foreach my $al (@al){
        next unless $hgvsg->{$al};
        push @{$hgvs{$hgvsg->{$al}}}, $vf->variation_name();
      }
    }
    $start = $bin_end -  $self->required_param('overlap');
  }

  $self->store(\%hgvs);
}


sub store{

  my $self = shift;
  my $vars = shift;

  my $var_dba  = $self->get_species_adaptor('variation');
  my $va       = $var_dba->get_VariationAdaptor();

  ## Two multi-allelic variants may have multiple sets of matching HGVS currently
  ## Enter attrib only once
  my %done;

  foreach my $hgvs(keys %{$vars}){
    next if scalar @{$vars->{$hgvs}} == 1;

    my @u = unique(@{$vars->{$hgvs}});
    next if scalar @u ==1; ## shouldn't happen

    foreach my $name(@u){
      my $var = $va->fetch_by_name($name);
      die "$name not found\n" unless defined $var; ## shouldn't happen
  
      foreach my $m(@u){

        ## don't store the self matches
        next if $m eq $name;

        ## don't store same match twice
        next if $done{$name}{$m};
        $done{$name}{$m} = 1;

        $var->update_attributes( {"co-located allele"  => $m }) ;
        $va->store_attributes($var);
      }
    }
  }
}
sub unique{
  my %a;
  map { $a{$_} = 1; } @_;
  return sort keys %a;
}


1;
