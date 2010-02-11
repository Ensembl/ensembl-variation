# Ensembl module for Bio::EnsEMBL::Variation::VariationAnnotation
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
        -study        => 'pubmed/17293876',
        -associated_gene => 'HHEX',
        -associated_variant_risk_allele => 'rs13266634-C',
        -variation_names => 'rs13266634',
        -risk_allele_freq_in_controls => '0.3',
        -p_value  => '6.00E-08',
        -variation => $v);

    ...

    print $va->phenotype_name(),'-',$va->phenotype_description,"\n";
    print "From source ",$va->source_name,'-',$va->local_stable_id,"\n";
    print " With study_type ", $va->study_type(),"\n";

    ...
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
  Arg [-VARIATION_NAMES] :
    string - name of the associated variations
  Arg [-VARIATION] :
    int - the variation object associated with this annotation.
  Arg [_VARIATION_ID] :
    int _ the internal id of the variation object associated with this
    identifier. This may be provided instead of a variation object so that
    the variation may be lazy-loaded from the database on demand.
  Arg [-STUDY] :
    string - the pubmed/ids or project/study names
  Arg [-ASSOCIATED_GENE] :
    string - the gene names associated with this annotation/variant.
  Arg [-ASSOCIATED_VARIANT_RISK_ALLELE] :
    string - the variants-risk alleles associated with this annotation.
  Arg [-RISK_ALLELE_FREQ_IN_CONTROLS] :
    string - the risk allele frequency in controls associated with this annotation.
  Arg [-VARIATION] :
    string - the p_value associated with this annotation.

  Example    :
    $va = Bio::EnsEMBL::Variation::VariationAnnotation->new
       (-phenotype_name => 'BD',
        -phenotype_description => 'Bipolar Disorder', 
        -souce_name  => 'EGA',
        -variation_names => 'rs123',
        -local_stable_id => 'EGAS00000000001',
        _variation_id => 10,
        -study        => 'pubmed/17293876',
        -associated_gene => 'HHEX',
        -associated_variant_risk_allele => 'rs13266634-C',
        -risk_allele_freq_in_controls => '0.3',
        -p_value  => '6.00E-08',
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

  my ($dbID,$adaptor,$phenotype_id,$phenotype_name,$phenotype_description,$source_name,$study_type,$local_stable_id,$variation_id,$variation_names,$variation,$study,$associated_gene,$associated_variant_risk_allele,$risk_allele_freq_in_controls,$p_value) =
    rearrange([qw(dbID ADAPTOR _PHENOTYPE_ID PHENOTYPE_NAME PHENOTYPE_DESCRIPTION SOURCE_NAME
                  STUDY_TYPE LOCAL_STABLE_ID _VARIATION_ID VARIATION_NAMES VARIATION
		  STUDY ASSOCIATED_GENE ASSOCIATED_VARIANT_RISK_ALLELE
		 RISK_ALLELE_FREQ_IN_CONTROLS P_VALUE)],@_); 

  $self->{'dbID'} = $dbID;
  $self->{'adaptor'}    = $adaptor;
  $self->{'_phenotype_id'} = $phenotype_id;
  $self->{'phenotype_name'}   = $phenotype_name;
  $self->{'phenotype_description'}  = $phenotype_description;
  $self->{'local_stable_id'} = $local_stable_id;
  $self->{'variation'}        = $variation;
  $self->{'_variation_id'}    = $variation_id;
  $self->{'source_name'}      = $source_name;
  $self->{'variation_names'}   = $variation_names;
  $self->{'study_type'}  = $study_type;
  $self->{'study'} = $study;
  $self->{'associated_gene'} = $associated_gene;
  $self->{'associated_variant_risk_allele'} = $associated_variant_risk_allele;
  $self->{'risk_allele_freq_in_controls'} = $risk_allele_freq_in_controls;
  $self->{'p_value'} = $p_value;
 
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

=head2 variation_names

  Arg [1]    : string $newval (optional)
               The new value to set the variation_names attribute to
  Example    : $variation_names = $obj->variation_names()
  Description: Getter/Setter for the variation_names attribute.  This is the
               names of the variation associated with this study.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub variation_names{
  my $self = shift;
  return $self->{'variation_names'} = shift if(@_);
  return $self->{'variation_names'};
}

=head2 study

  Arg [1]    : string $newval (optional)
               The new value to set the study attribute to
  Example    : $variation_names = $obj->study()
  Description: Getter/Setter for the study attribute.  This is the
               pubmed/id or project name associated with this study.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub study{
  my $self = shift;
  return $self->{'study'} = shift if(@_);
  return $self->{'study'};
}

=head2 associated_gene

  Arg [1]    : string $newval (optional)
               The new value to set the associated_gene attribute to
  Example    : $associated_gene = $obj->associated_gene()
  Description: Getter/Setter for the associated_gene attribute.  This is the
               gene names associated with this study.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub associated_gene{
  my $self = shift;
  return $self->{'associated_gene'} = shift if(@_);
  return $self->{'associated_gene'};
}

=head2 associated_variant_risk_allele

  Arg [1]    : string $newval (optional)
               The new value to set the associated_variant_risk_allele attribute to
  Example    : $associated_variant_risk_allele = $obj->associated_variant_risk_allele()
  Description: Getter/Setter for the associated_variant_risk_allele attribute.  This is the
               associated_variant_risk_allele associated with this study.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub associated_variant_risk_allele{
  my $self = shift;
  return $self->{'associated_variant_risk_allele'} = shift if(@_);
  return $self->{'associated_variant_risk_allele'};
}

=head2 risk_allele_freq_in_controls

  Arg [1]    : string $newval (optional)
               The new value to set the risk_allele_freq_in_controls attribute to
  Example    : $risk_allele_freq_in_controls = $obj->risk_allele_freq_in_controls()
  Description: Getter/Setter for the risk_allele_freq_in_controls attribute.  This is the
               risk_allele_freq_in_controls associated with this study.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub risk_allele_freq_in_controls{
  my $self = shift;
  return $self->{'risk_allele_freq_in_controls'} = shift if(@_);
  return $self->{'risk_allele_freq_in_controls'};
}

=head2 p_value

  Arg [1]    : string $newval (optional)
               The new value to set the p_value attribute to
  Example    : $p_value = $obj->p_value()
  Description: Getter/Setter for the p_value attribute.  This is the
               p_value associated with this study.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub p_value{
  my $self = shift;
  return $self->{'p_value'} = shift if(@_);
  return $self->{'p_value'};
}


1;
