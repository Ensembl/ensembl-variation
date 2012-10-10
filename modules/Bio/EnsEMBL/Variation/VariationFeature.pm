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

# Ensembl module for Bio::EnsEMBL::Variation::VariationFeature
#
# Copyright (c) 2004 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::VariationFeature - A genomic position for a nucleotide variation.

=head1 SYNOPSIS

    # Variation feature representing a single nucleotide polymorphism
    $vf = Bio::EnsEMBL::Variation::VariationFeature->new
       (-start   => 100,
        -end     => 100,
        -strand  => 1,
        -slice   => $slice,
        -allele_string => 'A/T',
        -variation_name => 'rs635421',
        -map_weight  => 1,
        -variation => $v);

    # Variation feature representing a 2bp insertion
    $vf = Bio::EnsEMBL::Variation::VariationFeature->new
       (-start   => 1522,
        -end     => 1521, # end = start-1 for insert
        -strand  => -1,
        -slice   => $slice,
        -allele_string => '-/AA',
        -variation_name => 'rs12111',
        -map_weight  => 1,
        -variation => $v2);

    ...

    # a variation feature is like any other ensembl feature, can be
    # transformed etc.
    $vf = $vf->transform('supercontig');

    print $vf->start(), "-", $vf->end(), '(', $vf->strand(), ')', "\n";

    print $vf->name(), ":", $vf->allele_string();

    # Get the Variation object which this feature represents the genomic
    # position of. If not already retrieved from the DB, this will be
    # transparently lazy-loaded
    my $v = $vf->variation();

=head1 DESCRIPTION

This is a class representing the genomic position of a nucleotide variation
from the ensembl-variation database.  The actual variation information is
represented by an associated Bio::EnsEMBL::Variation::Variation object. Some
of the information has been denormalized and is available on the feature for
speed purposes.  A VariationFeature behaves as any other Ensembl feature.
See B<Bio::EnsEMBL::Feature> and B<Bio::EnsEMBL::Variation::Variation>.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::VariationFeature;

use Scalar::Util qw(weaken isweak);

use Bio::EnsEMBL::Variation::BaseVariationFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp); 
use Bio::EnsEMBL::Variation::Utils::Sequence qw(ambiguity_code hgvs_variant_notation SO_variation_class format_hgvs_string);
use Bio::EnsEMBL::Variation::Utils::Sequence;
use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT);
use Bio::EnsEMBL::Variation::Utils::Constants qw($DEFAULT_OVERLAP_CONSEQUENCE %VARIATION_CLASSES); 
use Bio::EnsEMBL::Variation::RegulatoryFeatureVariation;
use Bio::EnsEMBL::Variation::MotifFeatureVariation;
use Bio::EnsEMBL::Variation::ExternalFeatureVariation;
use Bio::EnsEMBL::Variation::IntergenicVariation;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;
use Bio::PrimarySeq;
use Bio::SeqUtils;

our @ISA = ('Bio::EnsEMBL::Variation::BaseVariationFeature');

our $DEBUG = 0;
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
    string - the name of the source where the SNP comes from
  
  Arg [-SOURCE_VERSION] :
    number - the version of the source where the SNP comes from

  Arg [-VALIDATION_CODE] :
     reference to list of strings

  Arg [-OVERLAP_CONSEQUENCES] :
     listref of Bio::EnsEMBL::Variation::OverlapConsequences - all the consequences of this VariationFeature

  Arg [-VARIATION_ID] :
    int - the internal id of the variation object associated with this
    identifier. This may be provided instead of a variation object so that
    the variation may be lazy-loaded from the database on demand.

  Example    :
    $vf = Bio::EnsEMBL::Variation::VariationFeature->new(
        -start   => 100,
        -end     => 100,
        -strand  => 1,
        -slice   => $slice,
        -allele_string => 'A/T',
        -variation_name => 'rs635421',
        -map_weight  => 1,
	    -source  => 'dbSNP',
	    -validation_code => ['cluster','doublehit'],
        -variation => $v
    );

  Description: Constructor. Instantiates a new VariationFeature object.
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
    
  my $self = $class->SUPER::new(@_);

  my (
      $allele_str, 
      $var_name, 
      $map_weight, 
      $variation,
      $variation_id, 
      $source, 
      $source_version, 
      $is_somatic, 
      $validation_code, 
      $overlap_consequences,
      $class_so_term,
      $minor_allele,
      $minor_allele_freq,
      $minor_allele_count
  ) = rearrange([qw(
          ALLELE_STRING 
          VARIATION_NAME 
          MAP_WEIGHT 
          VARIATION 
          _VARIATION_ID 
          SOURCE 
          SOURCE_VERSION
          IS_SOMATIC 
          VALIDATION_CODE 
		  OVERLAP_CONSEQUENCES 
          CLASS_SO_TERM
          MINOR_ALLELE
          MINOR_ALLELE_FREQUENCY
          MINOR_ALLELE_COUNT
        )], @_);

  $self->{'allele_string'}          = $allele_str;
  $self->{'variation_name'}         = $var_name;
  $self->{'map_weight'}             = $map_weight;
  $self->{'variation'}              = $variation;
  $self->{'_variation_id'}          = $variation_id;
  $self->{'source'}                 = $source;
  $self->{'source_version'}         = $source_version;
  $self->{'is_somatic'}             = $is_somatic;
  $self->{'validation_code'}        = $validation_code;
  $self->{'overlap_consequences'}   = $overlap_consequences;
  $self->{'class_SO_term'}          = $class_so_term;
  $self->{'minor_allele'}           = $minor_allele;
  $self->{'minor_allele_frequency'} = $minor_allele_freq;
  $self->{'minor_allele_count'}     = $minor_allele_count;
  
  return $self;
}



sub new_fast {

  my $class = shift;
  my $hashref = shift;
  my $self = bless $hashref, $class;
  weaken($self->{'adaptor'})  if ( ! isweak($self->{'adaptor'}) );
  return $self;

}


=head2 allele_string

  Arg [1]    : string $newval (optional)
               The new value to set the allele_string attribute to
  Example    : $allele_string = $obj->allele_string()
  Description: Getter/Setter for the allele_string attribute.
               The allele_string is a '/' demimited string representing the
               alleles associated with this features variation.
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

=head2 display_id

  Arg [1]    : none
  Example    : print $vf->display_id(), "\n";
  Description: Returns the 'display' identifier for this feature. For
               VariationFeatures this is simply the name of the variation
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

=head2 minor_allele

  Arg [1]    : string $minor_allele (optional)
               The new minor allele string
  Example    : $ma = $obj->minor_allele()
  Description: Get/set the minor allele of this variation, as reported by dbSNP
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub minor_allele {
    my ($self, $minor_allele) = @_;
    $self->{minor_allele} = $minor_allele if defined $minor_allele;
    return $self->{minor_allele}
}

=head2 minor_allele_frequency

  Arg [1]    : float $minor_allele_frequency (optional)
               The new minor allele frequency
  Example    : $maf = $obj->minor_allele_frequency()
  Description: Get/set the frequency of the minor allele of this variation, as reported by dbSNP
  Returntype : float
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub minor_allele_frequency {
    my ($self, $minor_allele_frequency) = @_;
    $self->{minor_allele_frequency} = $minor_allele_frequency if defined $minor_allele_frequency;
    return $self->{minor_allele_frequency}
}

=head2 minor_allele_count

  Arg [1]    : int $minor_allele_count (optional)
               The new minor allele count
  Example    : $maf_count = $obj->minor_allele_count()
  Description: Get/set the sample count of the minor allele of this variation, as reported by dbSNP
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub minor_allele_count {
    my ($self, $minor_allele_count) = @_;
    $self->{minor_allele_count} = $minor_allele_count if defined $minor_allele_count;
    return $self->{minor_allele_count}
}



=head2 get_all_TranscriptVariations

  Arg [1]     : (optional) listref of Bio::EnsEMBL::Transcript objects
  Example     : $vf->get_all_TranscriptVariations;
  Description : Get all the TranscriptVariations associated with this VariationFeature.
                If the optional list of Transcripts is supplied, get only TranscriptVariations
		        associated with those Transcripts.
  Returntype  : listref of Bio::EnsEMBL::Variation::TranscriptVariation objects
  Exceptions  : Thrown on wrong argument type
  Caller      : general
  Status      : At Risk

=cut

sub get_all_TranscriptVariations {
    
    my ($self, $transcripts) = @_;

    if ($transcripts) {
        assert_ref($transcripts, 'ARRAY');
        map { assert_ref($_, 'Bio::EnsEMBL::Transcript') } @$transcripts;
    }

    #die unless $self->{transcript_variations};

    if ($self->dbID && not defined $self->{transcript_variations}) {
        # this VariationFeature is from the database, so we can just fetch the 
        # TranscriptVariations from the database as well

        if (my $db = $self->adaptor->db) {
            my $tva = $db->get_TranscriptVariationAdaptor;

            # just fetch TVs for all Transcripts because that's more efficient,
            # we'll only return those asked for later on

            my $tvs = $tva->fetch_all_by_VariationFeatures([$self]);

            map { $self->add_TranscriptVariation($_) } @$tvs;
        }
    }
    elsif (not defined $self->{transcript_variations}) {
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
            not exists $self->{transcript_variations}->{$_->stable_id} 
        } @$transcripts;

        for my $transcript (@unfetched_transcripts) {
            $self->add_TranscriptVariation(
                Bio::EnsEMBL::Variation::TranscriptVariation->new(
                    -variation_feature  => $self,
                    -transcript         => $transcript,
                    -adaptor            => ($self->adaptor->db ? $self->adaptor->db->get_TranscriptVariationAdaptor : undef),
                )
            );
        }
    }

    if ($transcripts) {
        # just return TranscriptVariations for the requested Transcripts
        return [ map { $self->{transcript_variations}->{$_->stable_id} } @$transcripts ];
    }
    else {
        # return all TranscriptVariations
        return [ values %{ $self->{transcript_variations} } ];
    }
}

=head2 get_all_RegulatoryFeatureVariations

  Description : Get all the RegulatoryFeatureVariations associated with this VariationFeature.
  Returntype  : listref of Bio::EnsEMBL::Variation::RegulatoryFeatureVariation objects
  Exceptions  : none
  Status      : At Risk

=cut

sub get_all_RegulatoryFeatureVariations {
    my $self = shift;
    return $self->_get_all_RegulationVariations('RegulatoryFeature', @_);
}

=head2 get_all_MotifFeatureVariations

  Description : Get all the MotifFeatureVariations associated with this VariationFeature.
  Returntype  : listref of Bio::EnsEMBL::Variation::MotifFeatureVariation objects
  Exceptions  : none
  Status      : At Risk

=cut

sub get_all_MotifFeatureVariations {
    my $self = shift;
    return $self->_get_all_RegulationVariations('MotifFeature', @_);
}

=head2 get_all_ExternalFeatureVariations

  Description : Get all the ExternalFeatureVariations associated with this VariationFeature.
  Returntype  : listref of Bio::EnsEMBL::Variation::ExternalFeatureVariation objects
  Exceptions  : none
  Status      : At Risk

=cut

sub get_all_ExternalFeatureVariations {
    my $self = shift;
    return $self->_get_all_RegulationVariations('ExternalFeature', @_);
}

sub _get_all_RegulationVariations {
    my ($self, $type) = @_;

    unless ($type && ($type eq 'RegulatoryFeature' || $type eq 'MotifFeature' || $type eq 'ExternalFeature')) {
        throw("Invalid Ensembl Regulation type '$type'");
    }

    unless ($self->{regulation_variations}->{$type}) {
    
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
                
        my $constructor = 'Bio::EnsEMBL::Variation::'.$type.'Variation';

		eval {
		  $self->{regulation_variations}->{$type} = [ 
			  map {  
				  $constructor->new(
					  -variation_feature  => $self,
					  -feature            => $_,
				  );
			  } map { $_->transfer($self->slice) } @{ $fg_adaptor->fetch_all_by_Slice($slice) } 
		  ];
		};
		
		$self->{regulation_variations}->{$type} ||= [];
    }

    return $self->{regulation_variations}->{$type};
}

sub get_IntergenicVariation {
    my $self = shift;
    my $no_ref_check = shift;

    unless (exists $self->{intergenic_variation}) {
        if (scalar(@{ $self->get_all_TranscriptVariations }) == 0) {
            $self->{intergenic_variation} = Bio::EnsEMBL::Variation::IntergenicVariation->new(
                -variation_feature  => $self,
                -no_ref_check       => $no_ref_check,
            );
        }
        else {
            $self->{intergenic_variation} = undef;
        }
    }

    return $self->{intergenic_variation};
}

=head2 get_all_VariationFeatureOverlaps

  Description : Get all the VariationFeatureOverlaps associated with this VariationFeature, this
                includes TranscriptVariations and regulatory feature overlap object.
  Returntype  : listref of Bio::EnsEMBL::Variation::VariationFeatureOverlap objects
  Exceptions  : none
  Status      : At Risk

=cut

sub get_all_VariationFeatureOverlaps {
    my $self = shift;
    
    my $vfos =  [
        @{ $self->get_all_TranscriptVariations },
        @{ $self->get_all_RegulatoryFeatureVariations },
        @{ $self->get_all_MotifFeatureVariations },
        @{ $self->get_all_ExternalFeatureVariations },
    ];

    if (my $iv = $self->get_IntergenicVariation) {
        push @$vfos, $iv;
    }

    return $vfos;
}

=head2 add_TranscriptVariation

   Arg [1]     : Bio::EnsEMBL::Variation::TranscriptVariation
   Example     : $vf->add_TranscriptVariation($tv);
   Description : Adds a TranscriptVariation to the variation feature object.
   Exceptions  : thrown on bad argument
   Caller      : Bio::EnsEMBL::Variation::TranscriptVariationAdaptor
   Status      : At Risk

=cut

sub add_TranscriptVariation {
    my ($self, $tv) = @_;
    assert_ref($tv, 'Bio::EnsEMBL::Variation::TranscriptVariation');
    # we need to weaken the reference back to us to avoid a circular reference
    weaken($tv->{base_variation_feature});
    $self->{transcript_variations}->{$tv->transcript_stable_id} = $tv;
}

=head2 variation

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Variation $variation
  Example    : $v = $vf->variation();
  Description: Getter/Setter for the variation associated with this feature.
               If not set, and this VariationFeature has an associated adaptor
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
  elsif(!defined($self->{'variation'}) && $self->adaptor() &&
        defined($self->{'_variation_id'})) {
    # lazy-load from database on demand
    my $va = $self->adaptor->db()->get_VariationAdaptor();
    $self->{'variation'} = $va->fetch_by_dbID($self->{'_variation_id'});
  }

  return $self->{'variation'};
}

=head2 consequence_type

  Arg [1]    : (optional) String $term_type
  Description: Get a list of all the unique consequence terms of this 
               VariationFeature. By default returns Ensembl display terms
               (e.g. 'NON_SYNONYMOUS_CODING'). $term_type can also be 'label'
               (e.g. 'Non-synonymous coding'), 'SO' (Sequence Ontology, e.g.
               'non_synonymous_codon') or 'NCBI' (e.g. 'missense')
  Returntype : listref of strings
  Exceptions : none
  Status     : At Risk

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

=head2 get_all_OverlapConsequences

  Description: Get a list of all the unique OverlapConsequences of this VariationFeature, 
               calculating them on the fly from the TranscriptVariations if necessary
  Returntype : listref of Bio::EnsEMBL::Variation::OverlapConsequence objects
  Exceptions : none
  Status     : At Risk

=cut

sub get_all_OverlapConsequences {
    my $self = shift;

    unless ($self->{overlap_consequences}) {
        
        # work them out and store them in a hash keyed by SO_term as we don't 
        # want duplicates from different VFOs

        my %overlap_cons;

        for my $vfo (@{ $self->get_all_TranscriptVariations }) {
            for my $allele (@{ $vfo->get_all_alternate_VariationFeatureOverlapAlleles }) {
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

=head2 add_OverlapConsequence

  Arg [1]    : Bio::EnsEMBL::Variation::OverlapConsequence instance
  Description: Add an OverlapConsequence to this VariationFeature's list 
  Returntype : none
  Exceptions : throws if the argument is the wrong type
  Status     : At Risk

=cut

sub add_OverlapConsequence {
    my ($self, $oc) = @_;
    assert_ref($oc, 'Bio::EnsEMBL::Variation::OverlapConsequence');
    push @{ $self->{overlap_consequences} ||= [] }, $oc;
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

=head2 display_consequence

  Arg [1]    : (optional) String $term_type
  Description: Get the term for the most severe consequence of this 
               VariationFeature. By default returns Ensembl display terms
               (e.g. 'NON_SYNONYMOUS_CODING'). $term_type can also be 'label'
               (e.g. 'Non-synonymous coding'), 'SO' (Sequence Ontology, e.g.
               'non_synonymous_codon') or 'NCBI' (e.g. 'missense')
  Returntype : string
  Exceptions : none
  Status     : At Risk

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

=head2 add_consequence_type

    Status : Deprecated, use add_OverlapConsequence instead

=cut

sub add_consequence_type{
    my $self = shift;
    warning('Deprecated method, use add_OverlapConsequence instead');
    return $self->add_OverlapConsequence(@_);
}

=head2 get_consequence_type

    Status : Deprecated, use consequence_type instead

=cut

sub get_consequence_type {
    my $self = shift;
    warning('Deprecated method, use consequence_type instead');
    return $self->consequence_type;
}

=head2 ambig_code

    Args         : None
    Example      : my $ambiguity_code = $vf->ambig_code()
    Description  : Returns the ambigutiy code for the alleles in the VariationFeature
    ReturnType   : String $ambiguity_code
    Exceptions   : none    
    Caller       : General
    Status       : At Risk

=cut 

sub ambig_code{
    my $self = shift;
    
    return &ambiguity_code($self->allele_string());
}

=head2 var_class

    Args[1]      : (optional) no_db - don't use the term from the database, always calculate it from the allele string 
                   (used by the ensembl variation pipeline)
    Example      : my $variation_class = $vf->var_class
    Description  : returns the Ensembl term for the class of this variation
    ReturnType   : string
    Exceptions   : throws if we can't find a corresponding display term for an SO term
    Caller       : General
    Status       : At Risk

=cut

sub var_class {

    my $self    = shift;
    my $no_db   = shift;
    
    unless ($self->{class_display_term}) {
        
        my $so_term = $self->class_SO_term(undef, $no_db);

        # convert the SO term to the ensembl display term
       
        $self->{class_display_term} = $self->is_somatic ? 
            $VARIATION_CLASSES{$so_term}->{somatic_display_term} : 
            $VARIATION_CLASSES{$so_term}->{display_term};
    }
    
    return $self->{class_display_term};
}

=head2 class_SO_term

    Args[1]      : (optional) class_SO_term - the SO term for the class of this variation feature
    Args[2]      : (optional) no_db - don't use the term from the database, always calculate it from the allele string 
                   (used by the ensembl variation pipeline)
    Example      : my $SO_variation_class = $vf->class_SO_term()
    Description  : Get/set the SO term for the class of this variation
    ReturnType   : string
    Exceptions   : none
    Caller       : General
    Status       : At Risk

=cut

sub class_SO_term {
    my ($self, $class_SO_term, $no_db) = @_;
   
    $self->{class_SO_term} = $class_SO_term if $class_SO_term;

    if ($no_db || !$self->{class_SO_term}) {
        $self->{class_SO_term} = SO_variation_class($self->allele_string);
    }

    return $self->{class_SO_term};
}

=head2 get_all_validation_states

  Arg [1]    : none
  Example    : my @vstates = @{$vf->get_all_validation_states()};
  Description: Retrieves all validation states for this variationFeature.  Current
               possible validation statuses are 'cluster','freq','submitter',
               'doublehit', 'hapmap'
  Returntype : reference to list of strings
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub get_all_validation_states {
    my $self = shift;
    return Bio::EnsEMBL::Variation::Utils::Sequence::get_all_validation_states($self->{'validation_code'});
}


=head2 add_validation_state

  Arg [1]    : string $state
  Example    : $vf->add_validation_state('cluster');
  Description: Adds a validation state to this variation.
  Returntype : none
  Exceptions : warning if validation state is not a recognised type
  Caller     : general
  Status     : At Risk

=cut

sub add_validation_state {
    Bio::EnsEMBL::Variation::Utils::Sequence::add_validation_state(@_);
}

=head2 source

  Arg [1]    : string $source_name (optional) - the new value to set the source attribute to
  Example    : $source = $vf->source;
  Description: Getter/Setter for the source attribute
  Returntype : the source name as a string, 
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub source {
  my ($self, $source) = @_;
  $self->{source} = $source if $source;
  return $self->{source};
}

=head2 source_version

  Arg [1]    : number $source_version (optional) - the new value to set the source version attribute to
  Example    : $source_version = $vf->source_version;
  Description: Getter/Setter for the source version attribute
  Returntype : the source version as a number 
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub source_version {
  my ($self, $source_version) = @_;
  $self->{source_version} = $source_version if $source_version;
  return $self->{source_version};
}

=head2 is_somatic

  Arg [1]    : boolean $is_somatic (optional)
               The new value to set the is_somatic flag to
  Example    : $is_somatic = $vf->is_somatic
  Description: Getter/Setter for the is_somatic flag, which identifies this variation feature as either somatic or germline
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

=head2 is_tagged

  Args        : None
  Example     : my $populations = $vf->is_tagged();
  Description : If the variation is tagged in any population, returns an array with the populations where the variation_feature
                is tagged (using a criteria of r2 > 0.99). Otherwise, returns null
  ReturnType  : list of Bio::EnsEMBL::Variation::Population
  Exceptions  : none
  Caller      : general
  Status      : At Risk
  
=cut

sub is_tagged{
    my $self = shift;
    
    if ($self->adaptor()){
	my $population_adaptor = $self->adaptor()->db()->get_PopulationAdaptor();
	return $population_adaptor->fetch_tagged_Population($self);
    }
}

=head2 is_tag

  Args        : None
  Example     : my $populations = $vf->is_tag();
  Description : Returns an array of populations in which this variation feature
                is a tag SNP.
  ReturnType  : list of Bio::EnsEMBL::Variation::Population
  Exceptions  : none
  Caller      : general
  Status      : At Risk
  
=cut

sub is_tag{
    my $self = shift;
    
    if ($self->adaptor()){
	my $population_adaptor = $self->adaptor()->db()->get_PopulationAdaptor();
	return $population_adaptor->fetch_tag_Population($self);
    }
}

=head2 get_all_tagged_VariationFeatures

  Args        : Bio::EnsEMBL::Variation::Population $pop (optional)
  Example     : my $vfs = $vf->get_all_tagged_VariationFeatures();
  Description : Returns an arrayref of variation features that are tagged by
                this variation feature, in the population $pop if specified.
  ReturnType  : list of Bio::EnsEMBL::Variation::VariationFeature
  Exceptions  : none
  Caller      : general
  Status      : At Risk
  
=cut

sub get_all_tagged_VariationFeatures {
  return $_[0]->adaptor->fetch_all_tagged_by_VariationFeature(@_);
}

=head2 get_all_tag_VariationFeatures

  Args        : Bio::EnsEMBL::Variation::Population $pop (optional)
  Example     : my $vfs = $vf->get_all_tag_VariationFeatures();
  Description : Returns an arrayref of variation features that tag this
                variation feature, in the population $pop if specified.
  ReturnType  : list of Bio::EnsEMBL::Variation::VariationFeature
  Exceptions  : none
  Caller      : general
  Status      : At Risk
  
=cut

sub get_all_tag_VariationFeatures {
  return $_[0]->adaptor->fetch_all_tags_by_VariationFeature(@_);
}

=head2 get_all_tag_and_tagged_VariationFeatures

  Args        : Bio::EnsEMBL::Variation::Population $pop (optional)
  Example     : my $vfs = $vf->get_all_tag_and_tagged_VariationFeatures();
  Description : Returns an arrayref of variation features that either tag or are
                tagged by this variation feature, in the population $pop if
				specified.
  ReturnType  : list of Bio::EnsEMBL::Variation::VariationFeature
  Exceptions  : none
  Caller      : general
  Status      : At Risk
  
=cut

sub get_all_tag_and_tagged_VariationFeatures {
  return $_[0]->adaptor->fetch_all_tags_and_tagged_by_VariationFeature(@_);
}



=head2 is_reference
  Arg        : none
  Example    : my $reference = $vf->is_reference()
  Description: Returns 1 if VF's slice is a reference slice else 0
  Returntype : int
  Caller     : general
  Status     : At Risk

=cut

sub is_reference {
  my ($self) = @_;
  my $slice = $self->slice;

  if ( !defined( $self->{'is_reference'} ) ) {
    $self->{'is_reference'} = $slice->is_reference();
  }

  return $self->{'is_reference'};
}

=head2 convert_to_SNP

  Args        : None
  Example     : my $snp = $vf->convert_to_SNP()
  Description : Creates a Bio::EnsEMBL::SNP object from Bio::EnsEMBL::VariationFeature. Mainly used for
                backwards compatibility
  ReturnType  : Bio::EnsEMBL::SNP
  Exceptions  : None
  Caller      : general      
  Status      : At Risk

=cut

sub convert_to_SNP{
    my $self = shift;

    require Bio::EnsEMBL::SNP;  #for backwards compatibility. It will only be loaded if the function is called

    my $snp = Bio::EnsEMBL::SNP->new_fast({
	        'dbID'       => $self->variation()->dbID(),
		'_gsf_start'  => $self->start,
		'_gsf_end'    => $self->end,
		'_snp_strand' => $self->strand,
		'_gsf_score'  => 1,
		'_type'       => $self->var_class,
		'_validated'  => $self->get_all_validation_states(),
		'alleles'    => $self->allele_string,
		'_ambiguity_code' => $self->ambig_code,
		'_mapweight'  => $self->map_weight,
		'_source' => $self->source
		});
    return $snp;
}

=head2 get_all_LD_values

    Args        : none
    Description : returns all LD values for this variation feature. This function will only work correctly if the variation
                  database has been attached to the core database. 
    ReturnType  : Bio::EnsEMBL::Variation::LDFeatureContainer
    Exceptions  : none
    Caller      : snpview
    Status      : At Risk
                : Variation database is under development.

=cut

sub get_all_LD_values{
    my $self = shift;
    
    if ($self->adaptor()){
	my $ld_adaptor = $self->adaptor()->db()->get_LDFeatureContainerAdaptor();
	return $ld_adaptor->fetch_by_VariationFeature($self);
    }
    return {};
}

=head2 get_all_LD_Populations

    Args        : none
    Description : returns a list of populations that could produces LD values
	              for this VariationFeature
    ReturnType  : listref of Bio::EnsEMBL::Variation::Population objects
    Exceptions  : none
    Caller      : snpview
    Status      : At Risk

=cut

sub get_all_LD_Populations{
    my $self = shift;
    
	my $pa = $self->adaptor->db->get_PopulationAdaptor;
	return [] unless $pa;
	
	my $ld_pops = $pa->fetch_all_LD_Populations;
	return [] unless $ld_pops;
	
	my $sth = $self->adaptor->dbc->prepare(qq{
	  SELECT ip.population_sample_id, c.seq_region_start, c.genotypes
	  FROM compressed_genotype_region c, individual_population ip
	  WHERE c.sample_id = ip.individual_sample_id
	  AND c.seq_region_id = ?
	  AND c.seq_region_start < ?
	  AND c.seq_region_end > ?
	});
	
	my $this_vf_start = $self->seq_region_start;
	
	$sth->bind_param(1, $self->feature_Slice->get_seq_region_id);
	$sth->bind_param(2, $self->seq_region_end);
	$sth->bind_param(3, $this_vf_start);
	
	$sth->execute;
	
	my ($sample_id, $seq_region_start, $genotypes);
	$sth->bind_columns(\$sample_id, \$seq_region_start, \$genotypes);
	
	my %have_genotypes = ();
	
	while($sth->fetch()) {
	  
	  next if $have_genotypes{$sample_id};
	  
	  if($seq_region_start == $this_vf_start) {
		$have_genotypes{$sample_id} = 1;
		next;
	  }
	  
	  my @genotypes = unpack '(www)*', $genotypes;
	  my $gt_start = $seq_region_start;
	  
	  while(my( $var_id, $gt_code, $gap ) = splice @genotypes, 0, 3 ) {
		if($gt_start == $this_vf_start) {
		  $have_genotypes{$sample_id} = 1;
		  last;
		}
		$gt_start += $gap + 1 if defined $gap;
	  }
	}
	
	my @final_list = grep {$have_genotypes{$_->dbID}} @$ld_pops;
	
	return \@final_list;
}

=head2 get_all_sources

    Args        : none
    Example     : my @sources = @{$vf->get_all_sources()};
    Description : returns a list of all the sources for this
                  VariationFeature
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
    if ($self->adaptor()){
	map {$sources{$_}++} @{$self->adaptor()->get_all_synonym_sources($self)};
	$sources{$self->source}++;
	@sources = keys %sources;
	return \@sources;
    }
    return \@sources;
}

=head2 ref_allele_string

  Args       : none
  Example    : $reference_allele_string = $self->ref_allele_string()
  Description: Getter for the reference allele_string, always the first.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub ref_allele_string{
    my $self = shift;

    my @alleles = split /[\|\\\/]/,$self->allele_string;
    return $alleles[0];
}


=head2 get_all_VariationSets

    Args        : none
    Example     : my @vs = @{$vf->get_all_VariationSets()};
    Description : returns a reference to a list of all the VariationSets this
                  VariationFeature is a member of
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
    my $variation_sets = $vs_adaptor->fetch_all_by_Variation($self->variation());
    
    return $variation_sets;
}


=head2 get_all_Alleles

  Args       : none
  Example    : @alleles = @{$vf->get_all_Alleles}
  Description: Gets all Allele objects from the underlying variation object,
			   with reference alleles first.
  Returntype : listref of Bio::EnsEMBL::Variation::Allele objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Alleles{
    my $self = shift;
	
	my @alleles = @{$self->variation->get_all_Alleles};
	
	# put all alleles in a hash
	my %order = ();
	foreach my $allele(@alleles) {
	  $order{$allele->allele} = 1;
	}
	
	$order{$self->ref_allele_string} = 2;
	
	# now sort them by population, submitter, allele
	my @new_alleles = sort {
	  ($a->population ? $a->population->name : "") cmp ($b->population ? $b->population->name : "") ||
	  ($a->subsnp ? $a->subsnp : "") cmp ($b->subsnp ? $b->subsnp : "") ||
	  $order{$b->allele} <=> $order{$a->allele}
	} @alleles;
	
	return \@new_alleles;
}


=head2 get_all_PopulationGenotypes

  Args       : none
  Example    : @pop_gens = @{$vf->get_all_PopulationGenotypes}
  Description: Gets all PopulationGenotype objects from the underlying variation
			   object, with reference genotypes first.
  Returntype : listref of Bio::EnsEMBL::Variation::PopulationGenotype objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_PopulationGenotypes{
    my $self = shift;
	
	my @gens = @{$self->variation->get_all_PopulationGenotypes};
	
	# put all alleles in a hash
	my %order = ();
	foreach my $gen(@gens) {
	  # homs low priority, hets higher
	  $order{$gen->allele1.$gen->allele2} = ($gen->allele1 eq $gen->allele2 ? 1 : 2);
	}
	
	# ref hom highest priority
	$order{$self->ref_allele_string x 2} = 3;
	
	# now sort them by population, submitter, genotype
	my @new_gens = sort {
	  ($a->population ? $a->population->name : "") cmp ($b->population ? $b->population->name : "") ||
	  ($a->subsnp ? $a->subsnp : "") cmp ($b->subsnp ? $b->subsnp : "") ||
	  $order{$b->allele1.$b->allele2} <=> $order{$a->allele1.$a->allele2}
	} @gens;
	
	return \@new_gens;
}

=head2 get_all_hgvs_notations

  Arg [1]    : Bio::EnsEMBL::Feature $ref_feature (optional)
               Get the HGVS notation of this VariationFeature relative to the slice it is on. If an optional reference feature is supplied, returns the coordinates
	       relative to this feature.
  Arg [2]    : string (Optional)
	       Indicate whether the HGVS notation should be reported in genomic coordinates or cDNA coordinates.
	       'g' -> Genomic position numbering
	       'c' -> cDNA position numbering
	       'p' -> protein position numbering
  Arg [3]    : string (Optional)
               A name to use for the reference can be supplied. By default the name returned by the display_id() method of the reference feature will be used.
  Arg [4]    : string (Optional)
               Return just the HGVS notation corresponding to this allele
	       
  Example    : my $vf = $variation_feature_adaptor->fetch_by_dbID(565770);
	       my $tr = $transcript_adaptor->fetch_by_stable_id('ENST00000335295');
	       my $hgvs = $vf->get_all_hgvs_notations($tr,'p');
	       while (my ($allele,$hgvs_str) = each(%{$hgvs})) {
		print "Allele $allele :\t$hgvs_str\n"; # Will print 'Allele - : ENSP00000333994.3:p.Val34_Tyr36delinsAsp'
	       }
	       
  Description: Returns a reference to a hash with the allele as key and a string with the HGVS notation of this VariationFeature as value. By default uses the
               slice it is plcaed on as reference but a different reference feature can be supplied.
  Returntype : Hash reference
  Exceptions : Throws exception if VariationFeature can not be described relative to the feature_Slice of the supplied reference feature
  Caller     : general
  Status     : Experimental

=cut
sub get_all_hgvs_notations {
    
    my $self                 = shift;
    my $ref_feature          = shift;
    my $numbering            = shift;    ## HGVS system g=genomic, c=coding, p=protein
    my $reference_name       = shift;    ## If the ref_feature is a slice, this is over-written
    my $use_allele           = shift;    ## optional single allele to check
    my $transcript_variation = shift;    ## optional transcript variation - looked up for c|p if not supplied
    
    my %hgvs;

    ##### don't get them for HGMD mutations or CNV probes
    return {} if ($self->allele_string =~ /INS|DEL|HGMD|CNV/ig || $self->var_class() =~ /microsat/i);
    ##### By default, use genomic position numbering
    $numbering ||= 'g';
      
    # If no reference feature is supplied, set it to the slice underlying this VariationFeature    
    $ref_feature  ||= $self->slice();
  	    
    # Special parsing for LRG
    if (defined $reference_name && $reference_name =~ /^LRG_/) {
		# Remove version
	if ($reference_name =~ /(.+)\.\d+$/) {
	    $reference_name = $1;
	}
    }    

    ### Check/get transcript variation available for protein & coding 
    if ($ref_feature->isa('Bio::EnsEMBL::Transcript')) {	
	
	# Get a TranscriptVariation object for this VariationFeature and the supplied Transcript if it wasn't passed in the call
	$transcript_variation = $self->get_all_TranscriptVariations([$ref_feature])->[0] if (!defined($transcript_variation));	
	
	##### call new TranscriptVariationAllele method for each allele
    }
    
      
    if ($numbering eq 'p') {
	
	#### If there is no transcript variation supplied and the variant
	#### is not in the translated region there is no protein change
	return {} if (!defined($transcript_variation) || 
		      !defined($transcript_variation->translation_start()) || 
		      !defined($transcript_variation->translation_end()));
	
	##### call TranscriptVariationAllele method for each allele
	foreach my $transcriptVariationAllele (@{$transcript_variation->get_all_alternate_TranscriptVariationAlleles()} ){

	    my $allele_string    = $transcriptVariationAllele->feature_seq();
	    my $hgvs_full_string = $transcriptVariationAllele->hgvs_protein();

	    if($allele_string ne (split/\//,$self->allele_string())[1] ){
		reverse_comp(\$allele_string);    ### hash returned relative to input variation feature strand regardless of transcript strand
	    }
	    $hgvs{$allele_string} = $hgvs_full_string ;
	} 
	return \%hgvs;
    }
    
    elsif ( $numbering =~ m/c|n/) { ### coding or non- coding transcript
	
	return {} if (!defined $transcript_variation);
	
	foreach my $transcriptVariationAllele (@{$transcript_variation->get_all_alternate_TranscriptVariationAlleles()} ){

	    my $allele_string    = $transcriptVariationAllele->feature_seq();
	    my $hgvs_full_string = $transcriptVariationAllele->hgvs_transcript();

	    if($allele_string ne (split/\//,$self->allele_string())[1] ){
		 reverse_comp(\$allele_string);    ### hash returned relative to input variation feature strand regardless of transcript strand
	    }
	    $hgvs{$allele_string} = $hgvs_full_string ;
	} 
	return \%hgvs;
    }
    
    elsif( $numbering =~ m/g/ ) {
	#### handling both alleles together locally for genomic class
	my $hgvs = $self->hgvs_genomic($ref_feature, $reference_name, $use_allele );	
	return $hgvs;
    }
    else{
	warn("HGVS notation is not available for coordinate system: $numbering");
	return {};
    }
}

### HGVS short flanking sequence required to check if HGVS variant class should be dup rather than ins
sub _get_flank_seq{

    my $self = shift;

    # Get the underlying slice and sequence
    my $ref_slice = $self->slice();
    
    my @allele = split(/\//,$self->allele_string());
    #### add flank at least as long as longest allele to allow checking 
    my $add_length = 0;

    foreach my $al(@allele){ ## may have >2 alleles
	if(length($al) > $add_length){
	    $add_length = length $al ;
	}
    }
    $add_length++;
    
    my $ref_start = $add_length ;
    my $ref_end   = $add_length + ($self->end() - $self->start());
    my $seq_start = ($self->start() - $ref_start);
    
    # Should we be at the beginning of the sequence, adjust the coordinates to not cause an exception
    if ($seq_start < 0) {
	$ref_start += $seq_start;
	$ref_end   += $seq_start;
	$seq_start  = 0;
    }
  	
    my $flank_seq   = $ref_slice->subseq($seq_start +1 , $seq_start + $ref_end, 1);
    
    
    return ($flank_seq, $ref_start, $ref_end );
}

#### format HGVS hash for genomic reference sequence
### Input:  Variation feature
### Output: $hgvs_notation hash



=head2  hgvs_genomic

  Arg [1]    : Bio::EnsEMBL::Feature $ref_feature (optional)
               Get the HGVS notation of this VariationFeature relative to the slice it is on. If an optional reference feature is supplied, returns the coordinates
	       relative to this feature.
  Arg [2]    : string (Optional)
               A name to use for the reference can be supplied. By default the name returned by the display_id() method of the reference feature will be used.
  Arg [4]    : string (Optional)
               Return just the HGVS notation corresponding to this allele
	       
  
	       
  Description: Returns a reference to a hash with the allele as key and a string with the genomic HGVS notation of this VariationFeature as value. By default uses the
               slice it is placed on as reference but a different reference feature can be supplied.
  Returntype : Hash reference
  Exceptions : Throws exception if VariationFeature can not be described relative to the feature_Slice of the supplied reference feature
  Caller     : general
  Status     : Experimental

=cut
sub hgvs_genomic{
    
    my $self             = shift;
    my $ref_feature      = shift;    ## can be a transcript
    my $reference_name   = shift;    ## If the ref_feature is a slice, this is over-written
    my $use_allele       = shift;    ## optional single allele to check

    my %hgvs;
    ########set up sequence reference
    my $ref_slice;  #### need this for genomic notation
    
    if ($ref_feature->isa('Bio::EnsEMBL::Slice')) {
	$ref_slice = $ref_feature;
    }
    else {	
	# get slice if not supplied
	$ref_slice = $ref_feature->feature_Slice;	    
    }         
    
    # Create new VariationFeature on the slice of the reference feature (unless the reference feature is the slice the VF is on)
    my $tr_vf;
    if ( $self->slice->coord_system->name() eq "chromosome") {
	$tr_vf = $self;
    }
    else {	
	$tr_vf = $self->transfer($ref_slice);
    }
    # Return undef if this VariationFeature could not be transferred
    return {} if (!defined($tr_vf));
	
	
    # Return undef if this VariationFeature does not fall within the supplied feature.
    return {} if ($tr_vf->start < 1 || 
		  $tr_vf->end   < 1 || 
		  $tr_vf->start > ($ref_feature->end - $ref_feature->start + 1) || 
		  $tr_vf->end   > ($ref_feature->end - $ref_feature->start + 1));
    
    #########   define reference sequence name ###################################
        
    # If the reference is a slice, use the seq_region_name as identifier
    $reference_name ||= $ref_feature->seq_region_name if ($ref_feature->isa('Bio::EnsEMBL::Slice'));
      
    # Use the feature's display id as reference name unless specified otherwise. 
    # If the feature is a transcript or translation, append the version number as well
    $reference_name ||= $ref_feature->display_id() . ($ref_feature->isa('Bio::EnsEMBL::Transcript') && 
						      $ref_feature->display_id !~ /\.\d+$/ ? '.' . $ref_feature->version() : '');
    
  
    ##### get short flank sequence for duplication checking & adjusted variation coordinates
    my ($ref_seq, $ref_start, $ref_end) = _get_flank_seq($tr_vf);

    my @all_alleles = split(/\//,$tr_vf->allele_string());    
    shift @all_alleles;  ## remove reference allele - not useful for HGVS
    foreach my $allele ( @all_alleles ) {
	
	## If a particular allele was requested, ignore others
	next if  (defined($use_allele) && $allele ne $use_allele);

	# Skip if the allele contains weird characters
	next if $allele =~ m/[^ACGT\-]/ig;   
	
	##### vf strand is relative to slice - if transcript feature slice, may need complimenting
	my $check_allele = $allele;
	if( $self->strand() <0 && $ref_slice->strand >0 ||
	    $self->strand() >0 && $ref_slice->strand < 0
	    ){	    
	    reverse_comp(\$check_allele);
	    if($DEBUG ==1){print "***************Flipping alt allele $allele => $check_allele to match coding strand\n";}
	}

	## work out chrom coord for hgvs string if transcript slice supplied
	my ($chr_start,$chr_end);  
	if ( $tr_vf->slice->is_toplevel() == 1) {
	    $chr_start = $tr_vf->seq_region_start();
	    $chr_end   = $tr_vf->seq_region_end();
	}
	else{
	    ## add feature start to start of var-in-feature
	    if( $tr_vf->seq_region_start() <  $tr_vf->seq_region_end()){
		$chr_start =  $tr_vf->start() + $tr_vf->seq_region_start() ; 
		$chr_end   =  $tr_vf->end()   + $tr_vf->seq_region_start() ;
	    }
	    else{
		$chr_start =  $tr_vf->seq_region_start() - $tr_vf->start()  ;
		$chr_end   =  $tr_vf->seq_region_start() - $tr_vf->end()  ;
	    }
	}
	   

	my $hgvs_notation = hgvs_variant_notation( $check_allele,          ## alt allele in refseq strand orientation
						   $ref_seq,               ## substring of slice for ref allele extraction
						   $ref_start,             ## start on substring of slice for ref allele extraction
						   $ref_end,
						   $chr_start,             ## start wrt seq region slice is on (eg. chrom)
						   $chr_end,   
						   $self->variation_name()); ## for error message 
      

	# Skip if e.g. allele is identical to the reference slice
	next if (!defined($hgvs_notation));
	
	# Add the name of the reference
	$hgvs_notation->{'ref_name'} = $reference_name;
	# Add the position_numbering scheme
	$hgvs_notation->{'numbering'} = 'g';     
	
	# Construct the HGVS notation from the data in the hash
	$hgvs_notation->{'hgvs'} = format_hgvs_string( $hgvs_notation);
	
	$hgvs{$allele} = $hgvs_notation->{'hgvs'};
    }
    return \%hgvs;

}



sub length {
  my $self = shift;
  return $self->{'end'} - $self->{'start'} + 1;
}

=head2 summary_as_hash

  Example       : $feature_summary = $feature->summary_as_hash();
  Description   : Extends Feature::summary_as_hash
                  Retrieves a summary of this VariationFeature object.
					                        
  Returns       : hashref of descriptive strings

=cut

sub summary_as_hash {
  my $self = shift;
  my $summary_ref = $self->SUPER::summary_as_hash;
  $summary_ref->{'consequence_type'} = $self->display_consequence;
  my @allele_list = split(/\//,$self->allele_string);
  $summary_ref->{'alt_alleles'} = \@allele_list;
  return $summary_ref;
}

=head2 flank_match

  Arg [1]    : int $newval (optional)
               The new value to set the flank_match attribute to
  Example    : $flank_match = $obj->flank_match()
  Description: Getter/Setter for the flank_match attribute.
               Return values:
               1 = submitted flank has perfect match to genomic reference (allowing for neighbouring variants)
               0 = imperfect match
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub flank_match{
  my $self = shift;
  return $self->{'flank_match'} = shift if(@_);
  return $self->{'flank_match'};
}



1;
