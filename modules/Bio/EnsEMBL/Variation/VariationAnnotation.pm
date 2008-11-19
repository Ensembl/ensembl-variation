# Ensembl module for Bio::EnsEMBL::Variation::VariationFeature
#
# Copyright (c) 2004 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::VariationAnnotation - A genotype phenotype annotation for a nucleotide variation.

=head1 SYNOPSIS

    # Variation Annotation is associated with a variation object
    $va = Bio::EnsEMBL::Variation::VariationAnnotation->new
       (_variation_id   => 8,
        -phenotype_name     => 'BD',
        -phenotype_description => 'Bipolar Disorder',
        -source_name  => 'EGA',
        -study_type   => 'GWAS',
        -local_stable_id => 'EGAS00000000001',
        -variation => $v);

    ...

    print $va->phenotype_name(),'-',$va->phenotype_description,"\n";
    print "From source ",$va->source_name,'-',$va->local_stable_id,"\n";
    print " With study_type ", $va->study_type(),"\n";

    # Get the Variation object which this annotation represents
    # If not already retrieved from the DB, this will be
    # transparently lazy-loaded
    my $v = $va->variation();

=head1 DESCRIPTION

This is a class representing the genotype-phenotype annotation of a variation
from the ensembl-variation database.  The actual variation information is
represented by an associated Bio::EnsEMBL::Variation::Variation object. 

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::VariationAnnotation;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);
use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Storable;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);

=head2 new

  Arg [-dbID] :
    int - unique internal identifier for variation_annotation
  Arg [-ADAPTOR] :
    Bio::EnsEMBL::Variation::DBSQL::VariationAnnotationAdaptor
  Arg [-PHENOTYPE_NAME] :
    string - name of the phenotype
  Arg [-PHENOTYPE_DESCRIPTION] :
    string - description of the phenotype
  Arg [-SOURCE_NAME] :
    string - name of the source
  Arg [-VARIATION_NAME] :
    string - name of the variation
  Arg [-VARIATION] :
    int - the variation object associated with this annotation.

  Arg [_VARIATION_ID] :
    int _ the internal id of the variation object associated with this
    identifier. This may be provided instead of a variation object so that
    the variation may be lazy-loaded from the database on demand.

  Example    :
    $va = Bio::EnsEMBL::Variation::VariationAnnotation->new
       (-phenotype_name => 'BD',
        -phenotype_description => 'Bipolar Disorder', 
        -souce_name  => 'EGA',
        -variation_name => 'rs123',
        -local_stable_id => 'EGAS00000000001',
        _variation_id => 10,
        -variation => $v);

  Description: Constructor. Instantiates a new VariationAnnotation object.
  Returntype : Bio::EnsEMBL::Variation::VariationAnnotation
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($dbID,$adaptor,$phenotype_name,$phenotype_description,$source_name,$study_type,$local_stable_id,$variation_id,$variation_name,$variation) =
    rearrange([qw(dbID ADAPTOR PHENOTYPE_NAME PHENOTYPE_DESCRIPTION SOURCE_NAME
                  STUDY_TYPE LOCAL_STABLE_ID _VARIATION_ID VARIATION_NAME VARIATION)],@_); 

  $self->{'dbID'} = $dbID;
  $self->{'adaptor'}    = $adaptor;
  $self->{'phenotype_name'}   = $phenotype_name;
  $self->{'phenotype_description'}  = $phenotype_description;
  $self->{'local_stable_id'} = $local_stable_id;
  $self->{'variation'}        = $variation;
  $self->{'_variation_id'}    = $variation_id;
  $self->{'source_name'}      = $source_name;
  $self->{'variation_name'}   = $variation_name;
  $self->{'study_type'}  = $study_type;
 
  return $self;
}



sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}


=head2 phenotype_name

  Arg [1]    : string phenotype_name (optional)
               The new value to set the phenotype_name attribute to
  Example    : $phenotype_name = $obj->phenotype_name()
  Description: Getter/Setter for the phenotype_name attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub phenotype_name{
  my $self = shift;
  return $self->{'phenotype_name'} = shift if(@_);
  return $self->{'phenotype_name'};
}

=head2 phenotype_description

  Arg [1]    : string phenotype_description (optional)
               The new value to set the phenotype_description attribute to
  Example    : $phenotype_description = $obj->phenotype_description()
  Description: Getter/Setter for the phenotype_description attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub phenotype_description{
  my $self = shift;
  return $self->{'phenotype_description'} = shift if(@_);
  return $self->{'phenotype_description'};
}

=head2 source_name

  Arg [1]    : string source_name (optional)
               The new value to set the source_name attribute to
  Example    : $source_name = $obj->source_name()
  Description: Getter/Setter for the source_name attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub source_name{
  my $self = shift;
  return $self->{'source_name'} = shift if(@_);
  return $self->{'source_name'};
}

=head2 study_type  

  Arg [1]    : string study_type (optional)               
               The new value to set the study_type attribute to  
  Example    : $study_type = $obj->study_type()  
  Description: Getter/Setter for the study_type attribute.  
  Returntype : string  
  Exceptions : none  
  Caller     : general  
  Status     : At Risk

=cut

sub study_type{  
  my $self = shift;  
  return $self->{'study_type'} = shift if(@_);  
  return $self->{'study_type'};
}

=head2 local_stable_id  

  Arg [1]    : string local_stable_id (optional)               
               The new value to set the local_stable_id attribute to  
  Example    : $local_stable_id = $obj->local_stable_id()  
  Description: Getter/Setter for the local_stable_id attribute.  
  Returntype : string  
  Exceptions : none  
  Caller     : general  
  Status     : At Risk
  
=cut

sub local_stable_id{  

  my $self = shift;  
  return $self->{'local_stable_id'} = shift if(@_);  
  return $self->{'local_stable_id'};
  
}

=head2 variation

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Variation $variation
  Example    : $v = $va->variation();
  Description: Getter/Setter for the variation associated with this annotation.
               If not set, and this VariationAnnotation has an associated adaptor
               an attempt will be made to lazy-load the variation from the
               database.
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Stable

=cut

sub variation {
  my $self = shift;

  if(@_) {
    if(!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::Variation')) {
      throw("Bio::EnsEMBL::Variation::Variation argument expected");
    }
    $self->{'variation'} = shift;
  }
  elsif(!defined($self->{'variation'}) && $self->{'adaptor'} &&
        defined($self->{'_variation_id'})) {
    # lazy-load from database on demand
    my $va = $self->{'adaptor'}->db()->get_VariationAdaptor();
    $self->{'variation'} = $va->fetch_by_dbID($self->{'_variation_id'});
  }

  return $self->{'variation'};
}

=head2 variation_name

  Arg [1]    : string $newval (optional)
               The new value to set the variation_name attribute to
  Example    : $variation_name = $obj->variation_name()
  Description: Getter/Setter for the variation_name attribute.  This is the
               name of the variation associated with this feature.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub variation_name{
  my $self = shift;
  return $self->{'variation_name'} = shift if(@_);
  return $self->{'variation_name'};
}


1;
