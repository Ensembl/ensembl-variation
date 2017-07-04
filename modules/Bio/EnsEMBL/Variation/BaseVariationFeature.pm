=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

# Ensembl module for Bio::EnsEMBL::Variation::VariationFeature
#
#


=head1 NAME

Bio::EnsEMBL::Variation::BaseVariationFeature - Abstract base class for variation features

=head1 SYNOPSIS

None

=head1 DESCRIPTION

Abstract base class representing variation features. Should not be instantiated
directly.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::BaseVariationFeature;

use Bio::EnsEMBL::Feature;

our @ISA = ('Bio::EnsEMBL::Feature');

=head2 consequence_type

  Arg [1]    : (optional) String $term_type
  Description: Get a list of all the unique consequence terms of this 
               VariationFeature. By default returns Ensembl display terms
               (e.g. 'NON_SYNONYMOUS_CODING'). $term_type can also be 'label'
               (e.g. 'Non-synonymous coding'), 'SO' (Sequence Ontology, e.g.
               'non_synonymous_codon') or 'NCBI' (e.g. 'missense')
  Returntype : listref of strings
  Exceptions : none
  Status     : Stable

=cut

sub consequence_type {
    
    my $self = shift;
    my $term_type = shift;
    
    my $method_name;
    
    # delete cached term
    if(defined($term_type)) {
        delete $self->{consequence_types};
        $method_name = $term_type.($term_type eq 'label' ? '' : '_term');
        $method_name = 'SO_term' unless @{$self->get_all_OverlapConsequences} && $self->get_all_OverlapConsequences->[0]->can($method_name);
    }
    
    $method_name ||= 'SO_term';

    if (exists($self->{current_consequence_method}) && $self->{current_consequence_method} ne $method_name) {
        delete $self->{consequence_type};
    }

    unless ($self->{consequence_types}) {

        # work out the terms from the OverlapConsequence objects
        
        $self->{consequence_types} = 
            [ map { $_->$method_name } @{ $self->get_all_OverlapConsequences } ];
    }

    $self->{current_consequence_method} = $method_name;
    
    return $self->{consequence_types};
}

=head2 display_consequence

  Arg [1]    : (optional) String $term_type
  Description: Get the term for the most severe consequence of this 
               VariationFeature. By default returns the Sequence Ontology term
               (e.g. 'missense_variant'). $term_type can also be 'label'
               (e.g. 'Missense variant'), 'Ensembl' ((e.g. 'NON_SYNONYMOUS_CODING') 
               or 'NCBI' (e.g. 'missense')
  Returntype : string
  Exceptions : none
  Status     : Stable

=cut

sub display_consequence {
    my $self = shift;
    my $term_type = shift;
    
    my $method_name;
    
    # delete cached term
    if(defined($term_type)) {
        $method_name = $term_type.($term_type eq 'label' ? '' : '_term');
        $method_name = 'SO_term' unless @{$self->get_all_OverlapConsequences} && $self->get_all_OverlapConsequences->[0]->can($method_name);
    }
    
    $method_name ||= 'SO_term';
    
    return $self->most_severe_OverlapConsequence->$method_name;
}

=head2 most_severe_OverlapConsequence

  Description: Get the OverlapConsequence considered (by Ensembl) to be the most severe 
               consequence of all the alleles of this VariationFeature 
  Returntype : Bio::EnsEMBL::Variation::OverlapConsequence
  Exceptions : none
  Status     : At Risk

=cut

sub most_severe_OverlapConsequence {
    my $self = shift;
    
    unless ($self->{_most_severe_consequence}) {
        
        my $highest;
        
        for my $cons (@{ $self->get_all_OverlapConsequences }) {
            $highest ||= $cons;
            if ($cons->rank < $highest->rank) {
                $highest = $cons;
            }
        }
        
        $self->{_most_severe_consequence} = $highest;
    }
    
    return $self->{_most_severe_consequence};
}

sub _get_transcript_key {
  my $self = shift;
  my $tr = shift;

  my $key = $tr->stable_id;

  if(my $dbID = $tr->dbID) {
    $key .= '_'.$dbID;
  }

  return $key;
}


=head2 _get_prev_base
  
  Arg 1      : (optional) int $strand
  Example    : $base = $bvf->get_prev_base();
  Description: Get the base preceding the given variant's position. Will
               use FASTA or database; returns "N" if sequence retrieval fails.
               Defaults to retrieving sequence for the slice's strand;
               strand may be specified with the first argument
  Returntype : string
  Exceptions : none
  Caller     : to_VCF_record(), 
  Status     : Stable

=cut

sub _get_prev_base {
  my $self = shift;
  my $strand = shift;

  # we need the ref base before the variation
  # default to N in case we cant get it
  my $prev_base = 'N';

  if(my $slice = $self->{slice}) {
    my $sub_slice = $slice->sub_Slice($self->start - 1, $self->start - 1, $strand);
    $prev_base = $sub_slice->seq if defined($sub_slice);
  }

  return $prev_base;
}

1;
