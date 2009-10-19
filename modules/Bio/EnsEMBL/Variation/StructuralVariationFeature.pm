# Ensembl module for Bio::EnsEMBL::Variation::StructuralVariationFeature
#
# Copyright (c) 2004 Ensembl
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
        -variation_name => 'esv25480',
        -map_weight  => 1,
        -variation => $v,
		-class     => 'CNV');

    ...

    # a variation feature is like any other ensembl feature, can be
    # transformed etc.
    $svf = $svf->transform('supercontig');

    print $svf->start(), "-", $svf->end(), '(', $svf->strand(), ')', "\n";

    print $svf->name(), ":", $svf->class();

    # Get the Variation object which this feature represents the genomic
    # position of. If not already retrieved from the DB, this will be
    # transparently lazy-loaded
    my $v = $svf->variation();

=head1 DESCRIPTION

This is a class representing the genomic position of a structural variation
from the ensembl-variation database.  The actual variation information is
represented by an associated Bio::EnsEMBL::Variation::Variation object. Some
of the information has been denormalized and is available on the feature for
speed purposes.  A StructuralVariationFeature behaves as any other Ensembl feature.
See B<Bio::EnsEMBL::Feature> and B<Bio::EnsEMBL::Variation::Variation>.

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::StructuralVariationFeature;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);
use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Slice;

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
	
  Arg [-BOUND_START] :
	int - the 5'-most coordinate of the underlying structural variation
	
  Arg [-BOUND_END] :
	int - the 3'-most coordinate of the underlying structural variation

  Arg [-VARIATION_NAME] :
    string - the name of the variation this feature is for (denormalisation
    from Variation object).

  Arg [-MAP_WEIGHT] :
    int - the number of times that the variation associated with this feature
    has hit the genome. If this was the only feature associated with this
    variation_feature the map_weight would be 1.

  Arg [-VARIATION] :
    int - the variation object associated with this feature.

  Arg [-SOURCE] :
    string - the name of the source where the variation comes from

  Arg [-TYPE] :
     string - the class of structural variation e.g. 'CNV'

  Arg [-VARIATION_ID] :
    int - the internal id of the variation object associated with this
    identifier. This may be provided instead of a variation object so that
    the variation may be lazy-loaded from the database on demand.

  Example    :
    $svf = Bio::EnsEMBL::Variation::StructuralVariationFeature->new
       (-start   => 100,
        -end     => 200,
        -strand  => 1,
        -slice   => $slice,
        -variation_name => 'esv25480',
        -map_weight  => 1,
        -variation => $v,
		-class => 'CNV');

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
  my ($var_name, $map_weight, $variation, $variation_id, $source, $sv_class, $bound_start, $bound_end) =
    rearrange([qw(VARIATION_NAME MAP_WEIGHT VARIATION _VARIATION_ID
				  SOURCE TYPE BOUND_START BOUND_END)], @_);

  $self->{'variation_name'}   = $var_name;
  $self->{'map_weight'}       = $map_weight;
  $self->{'variation'}        = $variation;
  $self->{'_variation_id'}    = $variation_id;
  $self->{'source'}           = $source;
  $self->{'class'}  = $sv_class;
  $self->{'bound_start'} = $bound_start;
  $self->{'bound_end'} = $bound_end;
 
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
               StructuralVariationFeatures this is simply the name of the variation
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



=head2 map_weight

  Arg [1]    : int $newval (optional) 
               The new value to set the map_weight attribute to
  Example    : $map_weight = $obj->map_weight()
  Description: Getter/Setter for the map_weight attribute. The map_weight
               is the number of times this features variation was mapped to
               the genome.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub map_weight{
  my $self = shift;
  return $self->{'map_weight'} = shift if(@_);
  return $self->{'map_weight'};
}

=head2 get_nearest_Gene

  Example     : $svf->get_nearest_Gene($flanking_size);
  Description : Getter a Gene which is associated to or nearest to the StructuralVariationFeature
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


=head2 variation

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Variation $variation
  Example    : $v = $svf->variation();
  Description: Getter/Setter for the variation associated with this feature.
               If not set, and this StructuralVariationFeature has an associated adaptor
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


=head2 class

    Args         : None
    Example      : my $sv_class = $svf->class()
    Description  : Getter/setter for the class of structural variation
    ReturnType   : String $sv_class
    Exceptions   : none
    Caller       : General
    Status       : At Risk

=cut

sub class{
  my $self = shift;
  
  $self->{'class'} = shift if @_;
  
  return $self->{'class'};
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


=head2 get_all_sources

    Args        : none
    Example     : my @sources = @{$svf->get_all_sources()};
    Description : returns a list of all the sources for this
                  StructuralVariationFeature
    ReturnType  : reference to list of strings
    Exceptions  : none
    Caller      : general
    Status      : At Risk
                : Variation database is under development.
=cut

sub get_all_sources{
    my $self = shift;
   
    my @sources;
    my %sources;
    if ($self->{'adaptor'}){
	map {$sources{$_}++} @{$self->{'adaptor'}->get_all_synonym_sources($self)};
	$sources{$self->source}++;
	@sources = keys %sources;
	return \@sources;
    }
    return \@sources;
}


=head2 bound_start

	Arg [1]    : int $bound_start (optional)
				The new value to set the bound_start attribute to
    Example     : my $bound_start = $svf->bound_start();
    Description : Getter/setter for the 5'-most coordinate defined for this StructuralVariationFeature
    ReturnType  : int
    Exceptions  : none
    Caller      : general
    Status      : At Risk
=cut

sub bound_start{
  my $self = shift;
  return $self->{'bound_start'} = shift if(@_);
  return $self->{'bound_start'};
}


=head2 bound_end

	Arg [1]    : int $bound_end (optional)
				The new value to set the bound_end attribute to
    Example     : my $bound_end = $svf->bound_end();
    Description : Getter/setter for the 3'-most coordinate defined for this StructuralVariationFeature
    ReturnType  : int
    Exceptions  : none
    Caller      : general
    Status      : At Risk
=cut

sub bound_end{
  my $self = shift;
  return $self->{'bound_end'} = shift if(@_);
  return $self->{'bound_end'};
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
	
	# fit the bound_start and bound_end coords to the new coords
	if(defined $self->{'bound_start'}) {
	  $transformed->{'bound_start'} = $transformed->start - ($self->start - $self->{'bound_start'});
	}
	
	if(defined $self->{'bound_end'}) {
	  $transformed->{'bound_end'} = $transformed->end + ($self->{'bound_end'} - $self->end);
	}
  }
  
  return $transformed;
}


sub transfer {
  my $self = shift;
  
  # run the transfer method from the parent class
  my $transferred = $self->SUPER::transfer(@_);
  
  if(defined $transferred) {
	
	# fit the bound_start and bound_end coords to the new coords
	if(defined $self->{'bound_start'}) {
	  $transferred->{'bound_start'} = $transferred->start - ($self->start - $self->{'bound_start'});
	}
	
	if(defined $self->{'bound_end'}) {
	  $transferred->{'bound_end'} = $transferred->end + ($self->{'bound_end'} - $self->end);
	}
  }
  
  return $transferred;
}

1;
