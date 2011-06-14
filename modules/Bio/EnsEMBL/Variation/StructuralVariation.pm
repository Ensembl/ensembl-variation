=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

# Ensembl module for Bio::EnsEMBL::Variation::StructuralVariation
#
# Copyright (c) 2004 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::StructuralVariation - A genomic position for a structural variation.

=head1 SYNOPSIS

    # Structural variation feature representing a CNV
    $svf = Bio::EnsEMBL::Variation::StructuralVariation->new
       (-start   => 100,
        -end     => 200,
        -strand  => 1,
        -slice   => $slice,
        -variation_name => 'esv25480',
				-class_so_term => 'structural_variant',
				-source => 'DGVa',
				-source_description => 'Database of Genomic Variants Archive',
				-study_name => 'estd20',
				-study_description => 'Conrad 2009 "Origins and functional impact of copy number variation in the human genome." PMID:19812545 [remapped from build NCBI36]',
				-study_url => 'ftp://ftp.ebi.ac.uk/pub/databases/dgva/estd20_Conrad_et_al_2009',
				-external_reference => 'pubmed/19812545');

    ...

    # a variation feature is like any other ensembl feature, can be
    # transformed etc.
    $svf = $svf->transform('supercontig');

    print $svf->start(), "-", $svf->end(), '(', $svf->strand(), ')', "\n";

    print $svf->name(), ":", $svf->class();

=head1 DESCRIPTION

This is a class representing the genomic position of a structural variation
from the ensembl-variation database.  A StructuralVariation behaves as any other
Ensembl feature. See B<Bio::EnsEMBL::Feature> and
B<Bio::EnsEMBL::Variation::Variation>.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::StructuralVariation;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);
use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Variation::Utils::Constants qw(%VARIATION_CLASSES); 

our @ISA = ('Bio::EnsEMBL::Feature');

=head2 new

  Arg [-dbID] :
    see superclass constructor

  Arg [-ADAPTOR] :
    see superclass constructor

  Arg [-START] :
    see superclass constructor
  Arg [-END] :
    see superclass constructor

  Arg [-STRAND] :
    see superclass constructor

  Arg [-SLICE] :
    see superclass constructor
	
  Arg [-INNER_START] :
	int - the 5'-less coordinate of the underlying structural variation
	
  Arg [-INNER_END] :
	int - the 3'-less coordinate of the underlying structural variation

  Arg [-VARIATION_NAME] :
    string - the name of the variation this feature is for (denormalisation
    from Variation object).

	Arg [-CLASS_SO_TERM] :
		string - the sequence ontology term defining the class of the structural variation.
		
  Arg [-SOURCE] :
    string - the name of the source where the variation comes from
	
  Arg [-SOURCE_DESCRIPTION] :
	string - description of the source

  Arg [-TYPE] :
     string - the class of structural variation e.g. 'SV'

  Arg [-STUDY_NAME] :
    string - the name of the study where the variation comes from
	
  Arg [-STUDY_DESCRIPTION] :
	string - description of the study
	
  Arg [-STUDY_URL] :
	string - url of the database/file where the data are stored

  Arg [-EXTERNAL_REFERENCE] :
	string - the pubmed/ids or project/study names
	
	Arg [-VALIDATION_STATUS] :
	string - the status of the structural variation (e.g. validated, not validated, ...)
	
	
  Example    :
    $svf = Bio::EnsEMBL::Variation::StructuralVariation->new
       (-start   => 100,
        -end     => 200,
        -strand  => 1,
        -slice   => $slice,
        -variation_name => 'esv25480',
		-class_so_term => 'structural_variant',
		-source => 'DGVa',
		-source_description => 'Database of Genomic Variants Archive',
		-study_name => 'estd20',
		-study_description => 'Conrad 2009 "Origins and functional impact of copy number variation in the human genome." PMID:19812545 [remapped from build NCBI36]',
		-study_url => 'ftp://ftp.ebi.ac.uk/pub/databases/dgva/estd20_Conrad_et_al_2009',
		-external_reference => 'pubmed/19812545');

  Description: Constructor. Instantiates a new StructuralVariation object.
  Returntype : Bio::EnsEMBL::Variation::StructuralVariation
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);
  
  my (
    $var_name, 
    $source, 
    $source_version, 
    $source_description, 
    $class_so_term, 
    $inner_start, 
    $inner_end, 
    $allele_string, 
    $study_name, 
    $study_description, 
    $study_url, 
    $external_reference,
		$validation_status
  ) = rearrange([qw(
          VARIATION_NAME 
          SOURCE 
          SOURCE_VERSION
          SOURCE_DESCRIPTION 
          CLASS_SO_TERM
          INNER_START 
          INNER_END 
          ALLELE_STRING 
          STUDY_NAME 
          STUDY_DESCRIPTION 
          STUDY_URL 
          EXTERNAL_REFERENCE
    )], @_);


  $self->{'variation_name'}     = $var_name;
  $self->{'source'}             = $source;
  $self->{'source_version'}     = $source_version;
  $self->{'source_description'} = $source_description;
  $self->{'class_SO_term'}      = $class_so_term;
  $self->{'inner_start'}        = $inner_start;
  $self->{'inner_end'}          = $inner_end;
  $self->{'allele_string'}      = $allele_string;
  $self->{'study_name'}         = $study_name;
  $self->{'study_description'}  = $study_description;
  $self->{'study_url'}          = $study_url;
  $self->{'external_reference'} = $external_reference;
	$self->{'validation_status'}  = $validation_status;

  return $self;
}



sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}


=head2 display_id

  Arg [1]    : none
  Example    : print $svf->display_id(), "\n";
  Description: Returns the 'display' identifier for this feature. For
               StructuralVariations this is simply the name of the variation
               it is associated with.
  Returntype : string
  Exceptions : none
  Caller     : webcode
  Status     : At Risk

=cut

sub display_id {
  my $self = shift;
  return $self->{'variation_name'} || '';
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
  Status     : Stable

=cut

sub variation_name{
  my $self = shift;
  return $self->{'variation_name'} = shift if(@_);
  return $self->{'variation_name'};
}

=head2 allele_string

  Arg [1]    : string $newval (optional)
               The new value to set the allele_string attribute to
  Example    : $allele_string = $obj->allele_string()
  Description: Getter/Setter for the allele_string attribute. This is the
               genomic sequence represented by this feature.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub allele_string{
  my $self = shift;
  return $self->{'allele_string'} = shift if(@_);
  return $self->{'allele_string'};
}


=head2 get_all_SupportingStructuralVariants

  Example     : $svf->get_all_SupportingStructuralVariants();
  Description : Retrieves all SupportingStructuralVariation associated with this structural variation.
                Return empty list if there are none.
  Returntype  : reference to list of Bio::EnsEMBL::Variation::SupportingStructuralVariation objects
  Exceptions  : None
  Caller      : general
  Status      : At Risk

=cut

sub get_all_SupportingStructuralVariants {
	my $self = shift;
	
	if (defined ($self->{'adaptor'})){
		my $ssv_adaptor = $self->{'adaptor'}->db()->get_SupportingStructuralVariationAdaptor();
		return $ssv_adaptor->fetch_all_by_StructuralVariation($self);
    }
    return [];
}




=head2 get_nearest_Gene

  Example     : $svf->get_nearest_Gene($flanking_size);
  Description : Getter a Gene which is associated to or nearest to the StructuralVariation
  Returntype  : a reference to a list of objects of Bio::EnsEMBL::Gene
  Exceptions  : None
  Caller      : general
  Status      : At Risk

=cut

sub get_nearest_Gene{

    my $self = shift;
    my $flanking_size = shift; #flanking size is optional
    $flanking_size ||= 0;
    my $sa = $self->{'adaptor'}->db()->dnadb->get_SliceAdaptor();
    my $slice = $sa->fetch_by_Feature($self,$flanking_size);
    my @genes = @{$slice->get_all_Genes};
    return \@genes if @genes; #$svf is on the gene

    if (! @genes) { #if $svf is not on the gene, increase flanking size
      warning("flanking_size $flanking_size is not big enough to overlap a gene, increase it by 1,000,000");
      $flanking_size += 1000000;
      $slice = $sa->fetch_by_Feature($self,$flanking_size);
      @genes = @{$slice->get_all_Genes};
    }
    if (@genes) {
      my %distances = ();
      foreach my $g (@genes) {
        if ($g->seq_region_start > $self->start) {
          $distances{$g->seq_region_start-$self->start}=$g;
        }
        else {
          $distances{$self->start-$g->seq_region_end}=$g;
        }
      }
      my @distances = sort {$a<=>$b} keys %distances;
      my $shortest_distance = $distances[0];
      if ($shortest_distance) {
        my $nearest_gene = $distances{$shortest_distance};
        return [$nearest_gene];
      }
    }
    else {
      throw("variation_feature with flanking_size $flanking_size is not overlap with a gene, try a bigger flanking_size");
    }
}


=head2 class

    Args         : None
    Example      : my $sv_class = $svf->class()
    Description  : Getter/setter for the class of structural variation
    ReturnType   : String $sv_class
    Exceptions   : none
    Caller       : General
    Status       : At Risk

=cut


sub class {
	my $self = shift;
    
	unless ($self->{class_display_term}) {
        my $display_term =
            $VARIATION_CLASSES{$self->{class_SO_term}}->{display_term};

        warn "No display term for SO term: ".$self->{class_SO_term} unless $display_term;

        $self->{class_display_term} = $display_term || $self->{class_SO_term};
    }

	return $self->{class_display_term};
}

sub class_SO_term {
	my $self = shift;

	return $self->{class_SO_term};
}



=head2 source

  Arg [1]    : string $source (optional)
               The new value to set the source attribute to
  Example    : $source = $svf->source()
  Description: Getter/Setter for the source attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub source{
  my $self = shift;
  return $self->{'source'} = shift if(@_);
  return $self->{'source'};
}

=head2 source_version

  Arg [1]    : string $source_version (optional)
               The new value to set the source_version attribute to
  Example    : $source_version = $svf->source_version()
  Description: Getter/Setter for the source_version attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub source_version {
  my $self = shift;
  return $self->{'source_version'} = shift if(@_);
  return $self->{'source_version'};
}

=head2 source_description

  Arg [1]    : string $source_description (optional)
               The new value to set the source_description attribute to
  Example    : $source_description = $svf->source_description()
  Description: Getter/Setter for the source_description attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub source_description{
  my $self = shift;
  return $self->{'source_description'} = shift if(@_);
  return $self->{'source_description'};
}


=head2 bound_start



	Arg [1]    : int $bound_start (optional)
				The new value to set the bound_start attribute to
    Example     : my $bound_start = $svf->bound_start();
    Description : Getter/setter for the 5'-most coordinate defined for this StructuralVariation
    ReturnType  : int
    Exceptions  : none
    Caller      : general
    Status      : DEPRECATED - Use the method seq_region_start
=cut

sub bound_start{
  my $self = shift;
	deprecate('Use the seq_region_start method instead.');
  return $self->{'start'} = shift if(@_);
  return $self->{'start'};
}


=head2 bound_end

	Arg [1]    : int $bound_end (optional)
				The new value to set the bound_end attribute to
    Example     : my $bound_end = $svf->bound_end();
    Description : Getter/setter for the 3'-most coordinate defined for this StructuralVariation
    ReturnType  : int
    Exceptions  : none
    Caller      : general
    Status      : DEPRECATED - Use the method seq_region_end
=cut

sub bound_end{
  my $self = shift;
	deprecate('Use the seq_region_end method instead.');
  return $self->{'end'} = shift if(@_);
  return $self->{'end'};
}



=head2 inner_start

	Arg [1]    : int $inner_start (optional)
				The new value to set the inner_start attribute to
    Example     : my $inner_start = $svf->inner_start();
    Description : Getter/setter for the 5'-less coordinate defined for this StructuralVariation
    ReturnType  : int
    Exceptions  : none
    Caller      : general
    Status      : At Risk
=cut

sub inner_start{
  my $self = shift;
  return $self->{'inner_start'} = shift if(@_);
  return $self->{'inner_start'};
}


=head2 inner_end

	Arg [1]    : int $inner_end (optional)
				The new value to set the bound_end attribute to
    Example     : my $inner_end = $svf->inner_end();
    Description : Getter/setter for the 3'-less coordinate defined for this StructuralVariation
    ReturnType  : int
    Exceptions  : none
    Caller      : general
    Status      : At Risk
=cut

sub inner_end{
  my $self = shift;
  return $self->{'inner_end'} = shift if(@_);
  return $self->{'inner_end'};
}


=head2 get_all_validation_states

  Arg [1]    : none
  Example    : my @vstates = @{$v->get_all_validation_states()};
  Description: Retrieves all validation states for this structural variation.  Current
               possible validation statuses are 'validated','not validated',
               'high quality'
  Returntype : reference to list of strings
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub get_all_validation_states {
  my $self = shift;

  return $self->{'validation_status'} || [];
}


=head2 get_reference_sequence

    Args        : none
    Example     : my $seq = $svf->get_reference_sequence
    Description : returns a string containing the reference sequence for the region
				  covered by this StructuralVariation
    ReturnType  : string
    Exceptions  : none
    Caller      : general
    Status      : At Risk
=cut

sub get_reference_sequence{
  my $self = shift;
  
  return $self->feature_Slice->seq();
}


sub transform {
  my $self = shift;
  
  # run the transform method from the parent class
  my $transformed = $self->SUPER::transform(@_);
  
  if(defined $transformed) {
	
	# fit the start and end coords to the new coords
	$transformed->_fix_bounds($self);
  }
  
  return $transformed;
}


sub transfer {
  my $self = shift;
  
  # run the transfer method from the parent class
  my $transferred = $self->SUPER::transfer(@_);
  
  if(defined $transferred) {
	
	# fit the start and end coords to the new coords
	$transferred->_fix_bounds($self);
  }
  
  return $transferred;
}



sub _fix_bounds {
  my $self = shift;
  my $old = shift;
  
  if(defined $old->{'start'}) {
	$self->{'start'} = $self->start - ($old->start - $old->{'start'});
  }
  
  if(defined $old->{'end'}) {
	$self->{'end'} = $self->end + ($old->{'end'} - $old->end);
  }
}


=head2 study_name

  Arg [1]    : string $study (optional)
               The new value to set the study attribute to
  Example    : $study = $svf->study()
  Description: Getter/Setter for the study attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub study_name{
  my $self = shift;
  return $self->{'study_name'} = shift if(@_);
  return $self->{'study_name'};
}



=head2 study_description

  Arg [1]    : string $study_description (optional)
               The new value to set the study_description attribute to
  Example    : $study_description = $svf->study_description()
  Description: Getter/Setter for the study_description attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub study_description{
  my $self = shift;
  return $self->{'study_description'} = shift if(@_);
  return $self->{'study_description'};
}

=head2 study_url

  Arg [1]    : string $newval (optional)
               The new value to set the study_url attribute to
  Example    : $paper = $obj->study_url()
  Description: Getter/Setter for the study_url attribute.This is the link to the website where the data are stored.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub study_url{
  my $self = shift;
  return $self->{'study_url'} = shift if(@_);
  return $self->{'study_url'};
}


=head2 external_reference

  Arg [1]    : string $newval (optional)
               The new value to set the external reference attribute to
  Example    : $paper = $obj->external_reference()
  Description: Getter/Setter for the external reference attribute. This is the
               pubmed/id or project name associated with this study.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub external_reference{
  my $self = shift;
  return $self->{'external_reference'} = shift if(@_);
  return $self->{'external_reference'};
}



=head2 is_supporting_structural_variation
  Example    : $sv = $obj->is_supporting_structural_variation()
  Description: Getter of the structural variation object for which this structural variant 
	             is a supporting evidence. 
  Returntype : A different Bio::EnsEMBL::Variation::StructuralVariation
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub is_supporting_structural_variation{
  my $self = shift;

	my $ssva = $self->{'adaptor'}->db()->get_SupportingStructuralVariationAdaptor();
	my $ssv  = $ssva->fetch_by_name($self->{'variation_name'});
	if (defined($ssv)) {
		return $ssv->get_StructuralVariation;
	}
	else { return undef; }
}
1;
