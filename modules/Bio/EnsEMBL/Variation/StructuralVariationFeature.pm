=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
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

# Ensembl module for Bio::EnsEMBL::Variation::StructuralVariationFeature
#
# Copyright (c) 2011 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::StructuralVariationFeature - A genomic position for a structural variation.

=head1 SYNOPSIS

    # Structural variation feature representing a CNV
    $svf = Bio::EnsEMBL::Variation::StructuralVariationFeature->new
       (-start   => 100,
        -end     => 200,
        -strand  => 1,
        -slice   => $slice,
        -variation_name => 'esv1001',
				-class_so_term => 'copy_number_variation',
				-source => 'DGVa',
				-source_description => 'Database of Genomic Variants Archive',
			 );

    ...

    print $svf->start(), "-", $svf->end(), '(', $svf->strand(), ')', "\n";

    print $svf->variation_name(), ":", $svf->var_class();

=head1 DESCRIPTION

This is a class representing the genomic position of a structural variant
from the ensembl-variation database.  A StructuralVariationFeature behaves as any other
Ensembl feature. See B<Bio::EnsEMBL::Feature> and
B<Bio::EnsEMBL::Variation::Variation>.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::StructuralVariationFeature;

use Bio::EnsEMBL::Variation::BaseVariationFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Variation::Utils::Constants qw(%VARIATION_CLASSES);
use Bio::EnsEMBL::Variation::StructuralVariationOverlap;

our @ISA = ('Bio::EnsEMBL::Variation::BaseVariationFeature');

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
	int - the 5'-greater coordinate of the underlying structural variation
	
  Arg [-INNER_END] :
	int - the 3'-less coordinate of the underlying structural variation

	 Arg [-OUTER_START] :
	int - the 5'-less coordinate of the underlying structural variation
	
  Arg [-OUTER_END] :
	int - the 3'-greater coordinate of the underlying structural variation

  Arg [-VARIATION_NAME] :
    string - the name of the variation this feature is for (denormalisation
    from Variation object).

	Arg [-CLASS_SO_TERM] :
		string - the sequence ontology term defining the class of the structural variation.
		
	Arg [-ALLELE_STRING] :
		string - allele sequence of the structural variation.
		
  Arg [-SOURCE] :
    string - the name of the source where the variation comes from
	
  Arg [-SOURCE_VERSION]:
	string - version number of the source
	
	Arg [-IS_SOMATIC] :
	  int - flag to inform whether the structural variant is a somatic (1) or germline (0).

	Arg [-BREAKPOINT_ORDER] :
	  int - For a structural variant with multiple breakpoints, this gives the predicted order of the breakpoint event.
	
  Example    :
    $svf = Bio::EnsEMBL::Variation::StructuralVariationFeature->new
       (-start   => 100,
        -end     => 200,
        -strand  => 1,
        -slice   => $slice,
        -variation_name => 'esv25480',
		    -class_so_term => 'structural_variant',
		    -source => 'DGVa');

  Description: Constructor. Instantiates a new StructuralVariationFeature object.
  Returntype : Bio::EnsEMBL::Variation::StructuralVariationFeature
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
    $class_so_term, 
    $inner_start, 
    $inner_end,
		$outer_start,
		$outer_end, 
    $allele_string,
		$is_somatic,
		$breakpoint_order
  ) = rearrange([qw(
          VARIATION_NAME 
          SOURCE 
          SOURCE_VERSION
          CLASS_SO_TERM
          INNER_START 
          INNER_END 
					OUTER_START
					INNER_START
          ALLELE_STRING
					IS_SOMATIC
					BREAKPOINT_ORDER
    )], @_);


  $self->{'variation_name'}     = $var_name;
  $self->{'source'}             = $source;
  $self->{'source_version'}     = $source_version;
  $self->{'class_SO_term'}      = $class_so_term;
  $self->{'inner_start'}        = $inner_start;
  $self->{'inner_end'}          = $inner_end;
	$self->{'outer_start'}        = $outer_start;
  $self->{'outer_end'}          = $outer_end;
  $self->{'allele_string'}      = $allele_string;
	$self->{'is_somatic'}         = $is_somatic || 0;
	$self->{'breakpoint_order'}   = $breakpoint_order;

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
               StructuralVariationFeatures this is simply the name of the structural variation
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
               name of the structural variant associated with this feature.
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



=head2 structural_variation

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::StructuralVariation or 
	             Bio::EnsEMBL::Variation::SupportingStructuralVariation $structural_variation
  Example    : $sv = $svf->structural_variation();
  Description: Getter/Setter for the structural variant associated with this feature.
               If not set, and this StructuralVariationFeature has an associated adaptor
               an attempt will be made to lazy-load the structural variation from the
               database.
  Returntype : Bio::EnsEMBL::Variation::StructuralVariation or 
	             Bio::EnsEMBL::Variation::SupportingStructuralVariation
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Stable

=cut

sub structural_variation {
  my $self = shift;

  if(@_) {
    if(!ref($_[0]) || (!$_[0]->isa('Bio::EnsEMBL::Variation::StructuralVariation') &&
		                   !$_[0]->isa('Bio::EnsEMBL::Variation::SupportingStructuralVariation')
		)) {
      throw("Bio::EnsEMBL::Variation::StructuralVariation or Bio::EnsEMBL::Variation::SupportingStructuralVariation argument expected");
    }
    $self->{'structural_variation'} = shift;
  }
  elsif(!defined($self->{'structural_variation'}) && $self->{'adaptor'} &&
        defined($self->{'structural_variation_id'})) {
    # lazy-load from database on demand
		my $sva = $self->{'adaptor'}->db()->get_StructuralVariationAdaptor();
		$self->{'structural_variation'} = $sva->fetch_by_dbID($self->{'structural_variation_id'});
		if (!defined($self->{'structural_variation'})) {
			$sva = $self->{'adaptor'}->db()->get_SupportingStructuralVariationAdaptor();
			$self->{'structural_variation'} = $sva->fetch_by_dbID($self->{'structural_variation_id'});
		}
  }

  return $self->{'structural_variation'};
}





=head2 get_all_VariationSets

    Args        : none
    Example     : my @vs = @{$svf->get_all_VariationSets()};
    Description : returns a reference to a list of all the VariationSets this
                  StructuralVariationFeature is a member of
    ReturnType  : reference to list of Bio::EnsEMBL::Variation::VariationSets
    Exceptions  : if no adaptor is attached to this object
    Caller      : general
    Status      : At Risk
=cut

sub get_all_VariationSets {
    my $self = shift;
    
    if (!$self->adaptor()) {
      throw('An adaptor must be attached in order to get all variation sets');
    }
    my $vs_adaptor = $self->adaptor()->db()->get_VariationSetAdaptor();
    my $variation_sets = $vs_adaptor->fetch_all_by_StructuralVariation($self->structural_variation());
    
    return $variation_sets;
}


=head2 get_nearest_Gene

  Example     : $svf->get_nearest_Gene($flanking_size);
  Description : Getter a Gene which is associated to or nearest to the StructuralVariationFeature
  Returntype  : Listref of objects of Bio::EnsEMBL::Gene
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


=head2 is_somatic

  Arg [1]    : boolean $is_somatic (optional)
               The new value to set the is_somatic flag to
  Example    : $is_somatic = $svf->is_somatic
  Description: Getter/Setter for the is_somatic flag, which identifies this structural variation feature as either somatic or germline
  Returntype : boolean
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub is_somatic {
  my ($self, $is_somatic) = @_;
  $self->{'is_somatic'} = $is_somatic if defined $is_somatic;
  return $self->{'is_somatic'};
}


=head2 breakpoint_order

  Arg [1]    : string $bp_order (optional)
               The new value to set the breakpoint order to
  Example    : $bp_order = $svf->breakpoint_order()
  Description: Getter/Setter for the breakpoint_order attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub breakpoint_order {
  my $self = shift;
  return $self->{'breakpoint_order'} = shift if(@_);
  return $self->{'breakpoint_order'};
}


sub get_all_StructuralVariationOverlaps {
  my $self = shift;
  
  if(!defined($self->{structural_variation_overlaps})) {
	my @svos = ();
	
	foreach my $gene(@{$self->feature_Slice->get_all_Genes}) {
	  my $svo = Bio::EnsEMBL::Variation::StructuralVariationOverlap->new(
		-feature                      => $gene,
		-structural_variation_feature => $self,
	  );
	  
	  push @svos, $svo if defined($svo);
	  
	  foreach my $tr(grep {$_->seq_region_start <= $self->seq_region_end && $_->seq_region_end >= $self->seq_region_start} @{$gene->get_all_Transcripts}) {
		$svo = undef;
		
		$svo = Bio::EnsEMBL::Variation::StructuralVariationOverlap->new(
		  -feature                      => $tr,
		  -structural_variation_feature => $self,
		);
		
		push @svos, $svo if defined $svo;
		
		foreach my $exon(grep {$_->seq_region_start <= $self->seq_region_end && $_->seq_region_end >= $self->seq_region_start} @{$tr->get_all_Exons}) {
		  $svo = undef;
		  
		  $svo = Bio::EnsEMBL::Variation::StructuralVariationOverlap->new(
			-feature                      => $exon,
			-structural_variation_feature => $self,
		  );
		  
		  push @svos, $svo if defined $svo;
		}
	  }
	}
	
	$self->{structural_variation_overlaps} = \@svos;
	
	# sort them
	#$self->_sort_svos;
  }
  
  return $self->{structural_variation_overlaps};
}


=head2 var_class

    Args         : None
    Example      : my $sv_class = $svf->var_class()
    Description  : Getter for the class of structural variation
    ReturnType   : String
    Exceptions   : none
    Caller       : General
    Status       : At Risk

=cut

sub var_class {
	my $self = shift;
    
	unless ($self->{class_display_term}) {
        my $display_term = $VARIATION_CLASSES{$self->{class_SO_term}}->{display_term};

        warn "No display term for SO term: ".$self->{class_SO_term} unless $display_term;

        $self->{class_display_term} = $display_term || $self->{class_SO_term};
    }

	return $self->{class_display_term};
}


=head2 class_SO_term

    Args         : None
    Example      : my $sv_so_term = $svf->class_SO_term()
    Description  : Getter for the class of structural variation, returning the SO term
    ReturnType   : String
    Exceptions   : none
    Caller       : General
    Status       : At Risk

=cut

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

=head2 bound_start

		Args        : None
    Example     : my $bound_start = $svf->bound_start();
    Description : Getter/setter for the 5'-most coordinate defined for this StructuralVariationFeature (outer_start or start)
    ReturnType  : int
    Exceptions  : none
    Caller      : general
    Status      : At risk
=cut

sub bound_start{
  my $self = shift;
	return $self->{'outer_start'} if (defined($self->{'outer_start'}));
	return $self->{'start'};
}


=head2 bound_end

		Args        : None
    Example     : my $bound_end = $svf->bound_end();
    Description : Getter/setter for the 3'-most coordinate defined for this StructuralVariationFeature (outer_end or end)
    ReturnType  : int
    Exceptions  : none
    Caller      : general
    Status      : At risk
=cut

sub bound_end{
  my $self = shift;
	return $self->{'outer_end'} if (defined($self->{'outer_end'}));
	return $self->{'end'};
}


=head2 outer_start

		Arg [1]     : int $outer_start (optional)
								  The new value to set the outer_start attribute to
    Example     : my $outer_start = $svf->outer_start();
    Description : Getter/setter for the 5'-most coordinate defined for this StructuralVariationFeature
    ReturnType  : int
    Exceptions  : none
    Caller      : general
    Status      : At risk
=cut

sub outer_start{
  my $self = shift;
  return $self->{'outer_start'} = shift if(@_);
  return $self->{'outer_start'};
}


=head2 outer_end

	  Arg [1]     : int $outer_end (optional)
				          The new value to set the outer_end attribute to
    Example     : my $outer_end = $svf->outer_end();
    Description : Getter/setter for the 3'-most coordinate defined for this StructuralVariationFeature
    ReturnType  : int
    Exceptions  : none
    Caller      : general
    Status      : At risk
=cut

sub outer_end{
  my $self = shift;
  return $self->{'outer_end'} = shift if(@_);
  return $self->{'outer_end'};
}


=head2 inner_start

	Arg [1]       : int $inner_start (optional)
				          The new value to set the inner_start attribute to
    Example     : my $inner_start = $svf->inner_start();
    Description : Getter/setter for the 5'-less coordinate defined for this StructuralVariationFeature
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

	  Arg [1]     : int $inner_end (optional)
				          The new value to set the inner_end attribute to
    Example     : my $inner_end = $svf->inner_end();
    Description : Getter/setter for the 3'-less coordinate defined for this StructuralVariationFeature
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


=head2 get_reference_sequence

    Args        : none
    Example     : my $seq = $svf->get_reference_sequence
    Description : returns a string containing the reference sequence for the region
				  covered by this StructuralVariationFeature
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
  
  if(defined $old->{'outer_start'}) {
	$self->{'outer_start'} = $self->start - ($old->start - $old->{'outer_start'});
  }
  
  if(defined $old->{'outer_end'}) {
	$self->{'outer_end'} = $self->end + ($old->{'outer_end'} - $old->end);
  }
}

sub _sort_svos {
  my $self = shift;
  
  return unless defined $self->{structural_variation_overlaps};
  
  my @svos = @{$self->{structural_variation_overlaps}};
  
  # define a feature order for sorting
  my %feature_order = (
	'Bio::EnsEMBL::Gene'       => 1,
	'Bio::EnsEMBL::Transcript' => 2,
	'Bio::EnsEMBL::Exon'       => 3,
  );
  
  # sort them nicely by feature type and position
  @svos = sort {
	$feature_order{ref($a->feature)} <=> $feature_order{ref($b->feature)} ||
	$a->feature->start <=> $b->feature->start
  } @svos;
  
  $self->{structural_variation_overlaps} = \@svos;
}

1;
