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

# Ensembl module for Bio::EnsEMBL::Variation::Genotype
#
#


=head1 NAME

Bio::EnsEMBL::Variation::Genotype - Abstract base class representing a genotype

=head1 SYNOPSIS

    print $genotype->variation()->name(), "\n";
    print $genotype->allele1(), '/', $genotype->allele2(), "\n";

=head1 DESCRIPTION

This is an abstract base class representing a genotype.  Specific types of
genotype are represented by subclasses such as IndividualGenotype and
PopulationGenotype.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::Genotype;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(strain_ambiguity_code);

use vars qw(@ISA $AUTOLOAD);

@ISA = qw(Bio::EnsEMBL::Storable);


sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}


=head2 allele
  Args       : int $index
               string $new_allele (optional)
  Examples   : $allele1 = $genotype->allele(1);
               $allele2 = $genotype->allele2();
  Description: Getter/Setter for one of the alleles that compose this genotype.
               Can be called as $genotype->allele(1), or via AUTOLOAD as
               $genotype->allele1()
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub allele {
  my $self = shift;
  my $index = shift;
  my $allele = shift;
  $index = 1 unless defined($index) && $index >= 1;
  
  $index--;
  
  $self->{genotype}->[$index] = $allele if defined($allele);
  
  return defined($self->{genotype}->[$index]) ?  $self->{genotype}->[$index] : undef;
}



=head2 genotype
  Examples   : @alleles = @{$genotype->genotype};
  Description: Getter for the genotype as an arrayref of alleles
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub genotype {
  return $_[0]->{genotype}
}



=head2 genotype_string

  Arg [1]    : (optional) bool $sort
  Examples   : $genotype_string = $genotype->genotype_string;
  Description: Gets the genotype as a '|'-separated string. Pass "1" as first
               argument to alphabetically sort genotype.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub genotype_string {  
  my $self = shift;
  my $sort = shift;

  my @gt = @{$self->genotype || []};
  @gt = sort @gt if defined($sort);

  return join '|', @gt;
}

=head2 variation

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Variation $var
  Example    : $var = $genotype->variation();
  Description: Getter/Setter for the Variation as
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

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


=head2 subsnp

  Arg [1]    : string $newval (optional) 
               The new value to set the subsnp attribute to
  Example    : print $genotype->subsnp();
  Description: Getter/Setter for the subsnp attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub subsnp {
  my $self = shift;
  if (@_) {
    $self->{'subsnp'} = shift;
  }
  
  my $ssid = $self->{'subsnp'};
  if (defined($ssid)) {
	  $ssid = 'ss'. $ssid unless $ssid =~ /^ss/;
  }
  
  return $ssid;
}

=head2 subsnp_handle

  Arg [1]    : string $newval (optional) 
               The new value to set the subsnp_handle attribute to
  Example    : print $genotype->subsnp_handle();
  Description: Getter/Setter for the subsnp_handle attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub subsnp_handle {
  my $self = shift;
  my $handle = shift;

  # if changing handle
  if (defined($handle)) {
    $self->{'subsnp_handle'} = $handle;
  }
  elsif (!defined($self->{'subsnp_handle'})) {
    # Check that this allele has an adaptor attached
    assert_ref($self->adaptor(), 'Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor');
    $self->{'subsnp_handle'} = $self->adaptor->get_subsnp_handle($self);
  }

  return $self->{'subsnp_handle'};
}


=head2 ambiguity_code

  Example    : print $genotype->ambiguity_code();
  Description: Get the ambiguity code for this genotype
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub ambiguity_code {
  return strain_ambiguity_code($_[0]->genotype_string);
}

=head2 phased

  Example    : $p = $genotype->phased()
  Description: Getter for the phased status of this genotype.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub phased {
  return $_[0]->{phased};
}

sub AUTOLOAD {
  my $self = shift;
  
  my $method = $AUTOLOAD;
  $method =~ s/.*:://;
  
  if ($method =~ /(allele)\_?(\d+)/) {
    $method = $1;
    unshift @_, $2;
  }
  else {
    return;
  }
  
  return $self->$method(@_);
}

1;
