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


=head1 NAME

Bio::EnsEMBL::Variation::AlleleSynonym - Ensembl representation of an allele synonym

=head1 SYNOPSIS

    # AlleleSynonym
    my $allele_synonym = Bio::EnsEMBL::Variation::AlleleSynonym->new
         (-hgvs_genomic => 'NC_000003.12:g.125439671T>A',
          -name         => 'CA83132407',
          -variation    => $v);
    
    print join("\t", $allele_synonym->variation()->name(),
                     $allele_synonym->name(),
                     $allele_synonym->hgvs_genomic()), "\n";


=head1 DESCRIPTION

This is a class representing an allele synonym for a variation

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::AlleleSynonym;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);
use Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor;

our @ISA = ('Bio::EnsEMBL::Storable');

=head2 new

  Arg [-dbID] :
    see superclass constructor

  Arg [-ADAPTOR] :
    see superclass constructor

  Arg [-VARIATION] :
    int - the variation object associated with this allele_synonym

  Arg [-_VARIATION_ID] :
    int - the internal id of the variation object associated with this
    identifier. This may be provided instead of a variation object so that
    the variation may be lazy-loaded from the database on demand.

  Arg [-HGVS_GENOMIC] :
    string - hgvs genomic
  
  Arg [-NAME] :
    string - allele synonym name
 
  Example    :

    $asa = $reg->get_adaptor('homo_sapiens', 'variation', 'allelesynonym');
    $allele_synonym  = Bio::EnsEMBL::Variation::AlleleSynonym->new(
        -_variation_id   => 176456403,
        -adaptor         => $asa,
        -hgvs_genomic    => 'NC_000003.12:g.125439671T>A',
        -name            => 'CA83132407'
    );

  Description: Constructor. Instantiates a Allele Synonym object.
  Returntype : Bio::EnsEMBL::Variation::AlleleSynonym
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut 

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);
  my ($dbID, $adaptor, $variation, $variation_id, $hgvs_genomic, $name) = 
      rearrange([qw(dbID ADAPTOR VARIATION _VARIATION_ID HGVS_GENOMIC NAME)], @_);

  return bless {
    'dBID'          => $dbID,
    'adaptor'       => $adaptor,
    'variation'     => $variation,
    '_variation_id' => $variation_id,
    'hgvs_genomic'  => $hgvs_genomic,
    'name'          => $name}, $class;
}

=head2 name

  Arg [1]    : string $newval (optional)
               The new value to set the name attribute to
  Example    : $name = $allele_synonym->name()
  Description: Getter/Setter for the name attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub name {
  my $self = shift;
  return $self->{'name'} = shift if(@_);
  return $self->{'name'};
}

=head2 hgvs_genomic

  Arg [1]    : string $newval (optional)
               The new value to set the hgvs_genomic attribute to
  Example    : $hgvs_genomic = $allele_synonym->hgvs_genomic()
  Description: Getter/Setter for the hgvs_genomic attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub hgvs_genomic {
  my $self = shift;
  return $self->{'hgvs_genomic'} = shift if(@_);
  return $self->{'hgvs_genomic'};
}

=head2 variation

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Variation $var
  Example    : $v = $allele_synonym->variation();
  Description: Getter/Setter for the variation associated with this allele synonym
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk

=cut

sub variation {
  my $self = shift;
  
  if(@_) {
    my $v = shift;
    if(defined($v) && (!ref($v) || !$v->isa('Bio::EnsEMBL::Variation::Variation'))) {
      throw('Bio::EnsEMBL::Variation::Variation argument expected');
    }
    return $self->{variation} = $v;
  }
  
  if(!defined($self->{variation}) && defined($self->{_variation_id})) {
    my $va = $self->adaptor->db->get_VariationAdaptor;
    if(defined($va)) {
      my $v = $va->fetch_by_dbID($self->{_variation_id});
      if(defined($v)) {
        $self->{variation} = $v;
      }
    }
  }
  return $self->{'variation'};
}

1;
