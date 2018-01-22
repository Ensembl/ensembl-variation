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

# Ensembl module for Bio::EnsEMBL::Variation::StructuralVariationFeature
#
#


=head1 NAME

Bio::EnsEMBL::Variation::StructuralVariationFeature - A genomic position for a structural variation.

=head1 SYNOPSIS

    # Structural variation feature representing a CNV
    $svf = Bio::EnsEMBL::Variation::StructuralVariationFeature->new
       (-start   => 100,
        -end     => 200,
        -strand  => 1,
        -slice   => $slice, # a Slice object (Core API) 
        -variation_name => 'esv1001',
        -class_so_term => 'copy_number_variation',
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

use Scalar::Util qw(weaken isweak);

use Bio::EnsEMBL::Variation::BaseVariationFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Variation::Utils::Constants qw($DEFAULT_OVERLAP_CONSEQUENCE %VARIATION_CLASSES); 
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT);
use Bio::EnsEMBL::Variation::StructuralVariationOverlap;
use Bio::EnsEMBL::Variation::TranscriptStructuralVariation;
use Bio::EnsEMBL::Variation::IntergenicStructuralVariation;

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
    object ref - the source object describing where the structural variant comes from.
   
  Arg [-STUDY] :
    object ref - the study object describing where the structural variant comes from.
  
  Arg [-IS_SOMATIC] :
    int - flag to inform whether the structural variant is a somatic (1) or germline (0).

  Arg [-BREAKPOINT_ORDER] :
    int - For a structural variant with multiple breakpoints, this gives the predicted order of the breakpoint event.
  
  Arg [-LENGTH] :
    int - Length of the structural variant. Useful when the structural variant is an insertion with given length. 

  Arg [-STRUCTURAL_VARIATION] :
    object ref - the structural variation object for this feature

  Example    :
    $svf = Bio::EnsEMBL::Variation::StructuralVariationFeature->new
       (-start   => 100,
        -end     => 200,
        -strand  => 1,
        -slice   => $slice,
        -variation_name => 'esv25480',
        -class_so_term => 'copy_number_variation'
       );

  Description: Constructor. Instantiates a new StructuralVariationFeature object.
  Returntype : Bio::EnsEMBL::Variation::StructuralVariationFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);
  
  my (
    $var_name, 
    $structural_variation_id,
    $source_id,
    $study_id, 
    $class_so_term, 
    $inner_start, 
    $inner_end,
    $outer_start,
    $outer_end, 
    $allele_string,
    $is_somatic,
    $breakpoint_order,
    $length,
    $source,
    $study,
    $structural_variation
  ) = rearrange([qw(
  VARIATION_NAME 
  _STRUCTURAL_VARIATION_ID
  _SOURCE_ID 
  _STUDY_ID
  CLASS_SO_TERM
  INNER_START 
  INNER_END 
  OUTER_START
  OUTER_END
  ALLELE_STRING
  IS_SOMATIC
  BREAKPOINT_ORDER
  LENGTH
  SOURCE
  STUDY
  STRUCTURAL_VARIATION
  )], @_);


  $self->{'variation_name'}           = $var_name;
  $self->{'_structural_variation_id'} = $structural_variation_id;
  $self->{'_source_id'}               = $source_id;
  $self->{'_study_id'}                = $study_id;
  $self->{'class_SO_term'}            = $class_so_term;
  $self->{'inner_start'}              = $inner_start;
  $self->{'inner_end'}                = $inner_end;
  $self->{'outer_start'}              = $outer_start;
  $self->{'outer_end'}                = $outer_end;
  $self->{'allele_string'}            = $allele_string;
  $self->{'is_somatic'}               = $is_somatic || 0;
  $self->{'breakpoint_order'}         = $breakpoint_order;
  $self->{'length'}                   = $length;
  $self->{'source'}                   = $source;
  $self->{'study'}                    = $study;
  $self->{'structural_variation'}     = $structural_variation;
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
  Status     : Stable

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
        defined($self->{'_structural_variation_id'})) {
    # lazy-load from database on demand
    my $sva = $self->{'adaptor'}->db()->get_StructuralVariationAdaptor();
    $self->{'structural_variation'} = $sva->fetch_by_dbID($self->{'_structural_variation_id'});
    if (!defined($self->{'structural_variation'})) {
      $sva = $self->{'adaptor'}->db()->get_SupportingStructuralVariationAdaptor();
      $self->{'structural_variation'} = $sva->fetch_by_dbID($self->{'_structural_variation_id'});
    }
  }

  return $self->{'structural_variation'};
}


=head2 length

  Example    : $length = $obj->length()
  Description: Getter for the length of the structural variant, if possible.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub length{
  my $self = shift;
  return $self->{'length'} if (defined($self->{'length'}));
  my $SO_term = $self->class_SO_term;
  if ($SO_term =~ /copy|deletion|duplication|inversion/i && !$self->breakpoint_order) {
    return ($self->seq_region_end-$self->seq_region_start)+1;
  }
}


=head2 get_all_VariationSets

    Args        : none
    Example     : my @vs = @{$svf->get_all_VariationSets()};
    Description : returns a reference to a list of all the VariationSets this
                  StructuralVariationFeature is a member of
    ReturnType  : reference to list of Bio::EnsEMBL::Variation::VariationSet objects
    Exceptions  : if no adaptor is attached to this object
    Caller      : general
    Status      : Stable
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

=head2 get_all_StructuralVariationOverlaps

  Description : Get all the StructuralVariationOverlaps associated with this StructuralVariation, this
                includes TranscriptStructuralVariations and regulatory feature overlap object.
  Returntype  : listref of Bio::EnsEMBL::Variation::StructuralVariationOverlap objects
  Exceptions  : none
  Status      : At Risk

=cut

sub get_all_StructuralVariationOverlaps {
  my $self = shift;
  
  my $vfos =  [
  @{ $self->get_all_TranscriptStructuralVariations },
  @{ $self->get_all_RegulatoryFeatureStructuralVariations },
  @{ $self->get_all_MotifFeatureStructuralVariations },
  ];

  if (my $iv = $self->get_IntergenicStructuralVariation) {
  push @$vfos, $iv;
  }

  return $vfos;
}

=head2 get_all_TranscriptStructuralVariations

  Arg [1]     : (optional) listref of Bio::EnsEMBL::Transcript objects
  Example     : $svf->get_all_TranscriptStructuralVariations;
  Description : Get all the TranscriptStructuralVariations associated with this
                StructuralVariationFeature. If the optional list of Transcripts
        is supplied, get only TranscriptStructuralVariations
            associated with those Transcripts.
  Returntype  : listref of Bio::EnsEMBL::Variation::TranscriptVariation objects
  Exceptions  : Thrown on wrong argument type
  Caller      : general
  Status      : At Risk

=cut

sub get_all_TranscriptStructuralVariations {
  my ($self, $transcripts) = @_;
  
  if ($transcripts) {
  assert_ref($transcripts, 'ARRAY');
  map { assert_ref($_, 'Bio::EnsEMBL::Transcript') } @$transcripts;
  }
  
  elsif (not defined $self->{transcript_structural_variations}) {
    # this VariationFeature is not in the database so we have to build the 
    # TranscriptVariations ourselves
    
    unless ($transcripts) {
      # if the caller didn't supply some transcripts fetch those around this VariationFeature
      # get a slice around this transcript including the maximum distance up and down-stream
      # that we still call consequences for
      my $slice = $self->feature_Slice->expand(
      MAX_DISTANCE_FROM_TRANSCRIPT, 
      MAX_DISTANCE_FROM_TRANSCRIPT
      );
      
      # fetch all transcripts on this slice 
      $transcripts = $slice->get_all_Transcripts(1);
    }
    
    my @unfetched_transcripts = grep { 
      not exists $self->{transcript_structural_variations}->{$self->_get_transcript_key($_)}
      and not exists $self->{transcript_structural_variations}->{$_->stable_id}
    } @$transcripts;
    
    for my $transcript (@unfetched_transcripts) {
      $self->add_TranscriptStructuralVariation(
        Bio::EnsEMBL::Variation::TranscriptStructuralVariation->new(
          -structural_variation_feature  => $self,
          -transcript                    => $transcript,
          -adaptor                       => undef,
        )
      );
    }
  }
  
  if ($transcripts) {
    # just return TranscriptVariations for the requested Transcripts
    return [
      map {
        $self->{transcript_structural_variations}->{$self->_get_transcript_key($_)} ||
        $self->{transcript_structural_variations}->{$_->stable_id}
      } @$transcripts
    ];
  }
  else {
    # return all TranscriptVariations
    return [ map {$self->{transcript_structural_variations}->{$_}} sort keys %{$self->{transcript_structural_variations}} ];
  }
}

=head2 get_all_RegulatoryFeatureStructuralVariations

  Description : Get all the RegulatoryFeatureStructuralVariations associated with this VariationFeature.
  Returntype  : listref of Bio::EnsEMBL::Variation::StructuralVariationOverlap objects
  Exceptions  : none
  Status      : At Risk

=cut

sub get_all_RegulatoryFeatureStructuralVariations {
  my $self = shift;
  return $self->_get_all_RegulationStructuralVariations('RegulatoryFeature', @_);
}

=head2 get_all_MotifFeatureStructuralVariations

  Description : Get all the MotifFeatureStructuralVariations associated with this VariationFeature.
  Returntype  : listref of Bio::EnsEMBL::Variation::StructuralVariationOverlap objects
  Exceptions  : none
  Status      : At Risk

=cut

sub get_all_MotifFeatureStructuralVariations {
  my $self = shift;
  return $self->_get_all_RegulationStructuralVariations('MotifFeature', @_);
}

=head2 get_all_ExternalFeatureStructuralVariations

  Description : Get all the ExternalFeatureStructuralVariations associated with this VariationFeature.
  Returntype  : listref of Bio::EnsEMBL::Variation::StructuralVariationOverlap objects
  Exceptions  : none
  Status      : At Risk

=cut

sub get_all_ExternalFeatureStructuralVariations {
  my $self = shift;
  return $self->_get_all_RegulationStructuralVariations('ExternalFeature', @_);
}

sub _get_all_RegulationStructuralVariations {
  my ($self, $type) = @_;
  
  unless ($type && ($type eq 'RegulatoryFeature' || $type eq 'MotifFeature' || $type eq 'ExternalFeature')) {
  throw("Invalid Ensembl Regulation type '$type'");
  }
  
  unless ($self->{regulation_structural_variations}->{$type}) {
  my $fg_adaptor;
  
  if (my $adap = $self->adaptor) {
    if(my $db = $adap->db) {
    $fg_adaptor = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(
      -species  => $adap->db->species, 
      -type     => $type,
    );      
    }
    
    unless ($fg_adaptor) {
    warning("Failed to get adaptor for $type");
    return [];
    }
  }
  else {
    warning('Cannot get variation features without attached adaptor');
    return [];
  }
  
  my $slice = $self->feature_Slice;
  
  my $constructor = 'Bio::EnsEMBL::Variation::StructuralVariationOverlap';
  
  eval {
    $self->{regulation_structural_variations}->{$type} = [ 
    map {  
      $constructor->new(
      -structural_variation_feature  => $self,
      -feature                       => $_,
      );
    } map { $_->transfer($self->slice) } @{ $fg_adaptor->fetch_all_by_Slice($slice) } 
    ];
  };
  
  $self->{regulation_structural_variations}->{$type} ||= [];
  }
  
  return $self->{regulation_structural_variations}->{$type};
}


sub get_IntergenicStructuralVariation {
  my $self = shift;
  my $no_ref_check = shift;
  
  unless (exists $self->{intergenic_structural_variation}) {
  if (scalar(@{ $self->get_all_TranscriptStructuralVariations }) == 0) {
    $self->{intergenic_structural_variation} = Bio::EnsEMBL::Variation::IntergenicStructuralVariation->new(
    -structural_variation_feature  => $self,
    -no_ref_check                  => $no_ref_check,
    );
  }
  else {
    $self->{intergenic_structural_variation} = undef;
  }
  }
  
  return $self->{intergenic_structural_variation};
}



=head2 TranscriptStructuralVariation

  Arg [1]     : Bio::EnsEMBL::Variation::TranscriptStructuralVariation
  Example     : $vf->add_TranscriptStructuralVariation($tsv);
  Description : Adds a TranscriptStructuralVariation to the structural variation
                feature object.
  Exceptions  : thrown on bad argument
  Caller      : Bio::EnsEMBL::Variation::StructuralVariationFeature,
                Bio::EnsEMBL::Varaition::Utils::VEP
  Status      : Stable

=cut

sub add_TranscriptStructuralVariation {
  my ($self, $tsv) = @_;
  assert_ref($tsv, 'Bio::EnsEMBL::Variation::TranscriptStructuralVariation');
  # we need to weaken the reference back to us to avoid a circular reference
  weaken($tsv->{base_variation_feature}) unless isweak($tsv->{base_variation_feature});

  # use a different method for speed
  my $tr_stable_id;

  # best for cache/VEP
  if(my $tr = $tsv->transcript) {
    $tr_stable_id = $self->_get_transcript_key($tr);
  }

  # best for API/DB
  else {
    $tr_stable_id = $tsv->transcript_stable_id;
  }

  $self->{transcript_structural_variations}->{$tr_stable_id} = $tsv;
}


=head2 get_all_OverlapConsequences

  Description: Get a list of all the unique OverlapConsequences of this StructuralVariationFeature, 
               calculating them on the fly from the StructuralTranscriptVariations if necessary
  Returntype : listref of Bio::EnsEMBL::Variation::OverlapConsequence objects
  Exceptions : none
  Status     : Stable

=cut

sub get_all_OverlapConsequences {
    my $self = shift;

    unless ($self->{overlap_consequences}) {
        
        # work them out and store them in a hash keyed by SO_term as we don't 
        # want duplicates from different VFOs

        my %overlap_cons;

        for my $vfo (@{ $self->get_all_TranscriptStructuralVariations }) {
            for my $allele (@{ $vfo->get_all_alternate_TranscriptStructuralVariationAlleles }) {
                for my $cons (@{ $allele->get_all_OverlapConsequences }) {
                    $overlap_cons{$cons->SO_term} = $cons;
                }
            }
        }

        # if we don't have any consequences we use a default from Constants.pm 
        # (currently set to the intergenic consequence)

        $self->{overlap_consequences} = [ 
            %overlap_cons ? values %overlap_cons : $DEFAULT_OVERLAP_CONSEQUENCE
        ];
    }

    return $self->{overlap_consequences};
}


=head2 get_all_supporting_evidence_classes

  Example     : $svf->_get_all_supporting_evidence_classes
  Description : Retrieves the classes (SO term) of the supporting evidence associated 
                with this structural variation feature (mainly, for web purposes).
                Return empty list if there are none.
  Returntype  : reference to list of string
  Exceptions  : None
  Caller      : web
  Status      : Stable

=cut

sub get_all_supporting_evidence_classes {
  my $self = shift;
  my $ssva = $self->adaptor->db->get_SupportingStructuralVariationAdaptor();
        
        my $ssv_SO_list = $ssva->fetch_all_SO_term_by_structural_variation_dbID($self->{_structural_variation_id});
  
        return $ssv_SO_list;
}


=head2 var_class

    Args        : None
    Example     : my $sv_class = $svf->var_class()
    Description : Getter for the class of structural variation
    ReturnType  : String
    Exceptions  : none
    Caller      : General
    Status      : Stable

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

    Args        : None
    Example     : my $sv_so_term = $svf->class_SO_term()
    Description : Getter for the class of structural variation, returning the SO term
    ReturnType  : String
    Exceptions  : none
    Caller      : General
    Status      : Stable

=cut

sub class_SO_term {
  my $self = shift;

  return $self->{class_SO_term};
}


=head2 source

  Arg [1]    : Bio::EnsEMBL::Variation::Source $src (optional)
               The new value to set the source attribute to
  Example    : $source = $svf->source()
  Description: Getter/Setter for the source object attribute
  Returntype : Bio::EnsEMBL::Variation::Source
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source{
  my $self = shift;
  
  # set
  if(@_) {
    if(!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::Source')) {
      throw("Bio::EnsEMBL::Variation::Source argument expected");
    }
    $self->{'source'} = shift;
  }
  # get
  elsif(!defined($self->{'source'}) && $self->adaptor() && defined($self->{'_source_id'})) {
    # lazy-load from database on demand
    my $sa = $self->adaptor->db()->get_SourceAdaptor();
    $self->{'source'} = $sa->fetch_by_dbID($self->{'_source_id'});
  }
  
  return $self->{'source'};
}


=head2 source_name

  Arg [1]    : string $source_name (optional)
               The new value to set the source name attribute to
  Example    : $source_name = $svf->source_name()
  Description: Getter/Setter for the source name attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source_name{
  my $self = shift;
  my $source = $self->source;
  return unless defined $source;
  
  $source->name(@_) if(@_);
  return $source->name;
}


=head2 source_version

  Arg [1]    : string $source_version (optional)
               The new value to set the source version attribute to
  Example    : $source_version = $svf->source_version()
  Description: Getter/Setter for the source version attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source_version{
  my $self = shift;
  my $source = $self->source;
  return unless defined $source;
  
  $source->version(@_) if(@_);
  return $source->version;
}


=head2 source_description

  Arg [1]    : string $source_description (optional)
               The new value to set the source description attribute to
  Example    : $source_description = $svf->source_description()
  Description: Getter/Setter for the source description attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source_description{
  my $self = shift;
  my $source = $self->source;
  return unless defined $source;
  
  $source->description(@_) if(@_);
  return $source->description;
}


=head2 bound_start

    Args        : None
    Example     : my $bound_start = $svf->bound_start();
    Description : Getter/setter for the 5'-most coordinate defined for this StructuralVariationFeature (outer_start or start)
    ReturnType  : int
    Exceptions  : none
    Caller      : general
    Status      : Stable
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
    Status      : Stable
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
    Status      : Stable
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
    Status      : Stable
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


=head2 study

  Arg [1]    : Bio::EnsEMBL::Variation::Study (optional)
  Example    : $study = $svf->study()
  Description: Getter/Setter for the study object
  Returntype : Bio::EnsEMBL::Variation::Study
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub study {
  my $self = shift;
  
  # set
 if(@_) {
    if(!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::Study')) {
      throw("Bio::EnsEMBL::Variation::Study argument expected");
    }
    $self->{'study'} = shift;
  }
  # get
  elsif(!defined($self->{'study'}) && $self->adaptor() && defined($self->{'_study_id'})) {
    # lazy-load from database on demand
    my $sa = $self->adaptor->db()->get_StudyAdaptor();
    $self->{'study'} = $sa->fetch_by_dbID($self->{'_study_id'});
  }
  
  return $self->{'study'};
}


=head2 get_reference_sequence

    Args        : none
    Example     : my $seq = $svf->get_reference_sequence
    Description : returns a string containing the reference sequence for the region
                  covered by this StructuralVariationFeature
    ReturnType  : String
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


=head2 to_VCF_record

  Example    : $vcf_arrayref = $svf->to_VCF_record();
  Description: Converts this StructuralVariationFeature object to an arrayref
               representing the columns of a VCF line.
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : VEP
  Status     : Stable

=cut

sub to_VCF_record {
  my $self = shift;

  # convert to SO term
  my %terms = (
    insertion => 'INS',
    deletion => 'DEL',
    tandem_duplication => 'TDUP',
    duplication => 'DUP'
  );

  my $alt = '<'.($terms{$self->class_SO_term} || $self->class_SO_term).'>';

  return [
    $self->{chr} || $self->seq_region_name,
    $self->start - 1,
    $self->variation_name || '.',
    $self->_get_prev_base(),
    $alt,
    '.', '.',
    'END='.$self->end
  ];
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

# used by VEP - fills out the consequence type hash keys
# so that the API doesn't try to go to the DB to fill them
sub _finish_annotation {
  my $self = shift;
  $self->{$_.'_structural_variations'} ||= {} for qw(transcript regulation);
  $self->{regulation_structural_variations}->{$_} ||= [] for qw(ExternalFeature MotifFeature RegulatoryFeature);
  $self->get_IntergenicStructuralVariation(1);
}

1;
