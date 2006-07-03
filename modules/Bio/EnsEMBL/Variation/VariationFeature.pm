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

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::VariationFeature;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(ambiguity_code variation_class);
use Bio::EnsEMBL::Variation::ConsequenceType;
use Bio::EnsEMBL::Variation::Variation;

our @ISA = ('Bio::EnsEMBL::Feature');

my %CONSEQUENCE_TYPES = %Bio::EnsEMBL::Variation::ConsequenceType::CONSEQUENCE_TYPES;

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

  Arg [-VALIDATION_CODE] :
     reference to list of strings

  Arg [-CONSEQUENCE_TYPE] :
     string - highest consequence type for the transcripts of the VariationFeature

  Arg [-VARIATION_ID] :
    int - the internal id of the variation object associated with this
    identifier. This may be provided instead of a variation object so that
    the variation may be lazy-loaded from the database on demand.

  Example    :
    $vf = Bio::EnsEMBL::Variation::VariationFeature->new
       (-start   => 100,
        -end     => 100,
        -strand  => 1,
        -slice   => $slice,
        -allele_string => 'A/T',
        -variation_name => 'rs635421',
        -map_weight  => 1,
	-source  => 'dbSNP',
	-validation_code => ['cluster','doublehit'],
	-consequence_type => 'INTRONIC',
        -variation => $v);

  Description: Constructor. Instantiates a new VariationFeature object.
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : none
  Caller     : general

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);
  my ($allele_str, $var_name, $map_weight, $variation, $variation_id, $source, $validation_code, $consequence_type) =
    rearrange([qw(ALLELE_STRING VARIATION_NAME 
                  MAP_WEIGHT VARIATION VARIATION_ID SOURCE VALIDATION_CODE 
		  CONSEQUENCE_TYPE)], @_);

  $self->{'allele_string'}    = $allele_str;
  $self->{'variation_name'}   = $var_name;
  $self->{'map_weight'}       = $map_weight;
  $self->{'variation'}        = $variation;
  $self->{'_variation_id'}    = $variation_id;
  $self->{'source'}           = $source;
  $self->{'validation_code'}  = $validation_code;
  $self->{'consequence_type'} = $consequence_type || 'INTERGENIC';
 
  return $self;
}



sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
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

=cut

sub map_weight{
  my $self = shift;
  return $self->{'map_weight'} = shift if(@_);
  return $self->{'map_weight'};
}


=head2 get_all_TranscriptVariations

  Example     : $vf->get_all_TranscriptVariations;
  Description : Getter a list with all the TranscriptVariations associated associated to the VariationFeature
  Returntype  : ref to Bio::EnsEMBL::Variation::VariationFeature
  Exceptions  : None
  Caller      : general

=cut

sub get_all_TranscriptVariations{
    my $self = shift;
    
    if(!defined($self->{'transcriptVariations'}) && $self->{'adaptor'})    {
	#lazy-load from database on demand
	my $tva = $self->{'adaptor'}->db()->get_TranscriptVariationAdaptor();
	$tva->fetch_all_by_VariationFeatures([$self]);
	$self->{'transcriptVariations'} ||= [];
    }
    return $self->{'transcriptVariations'};
}

=head2

   Arg [1]     : Bio::EnsEMBL::Variation::TranscriptVariation
   Example     : $vf->add_TranscriptVariation($tv);
   Description : Adds another Transcript variation to the variation feature object
   Exceptions  : thrown on bad argument
   Caller      : Bio::EnsEMBL::Variation::TranscriptVariationAdaptor

=cut

sub add_TranscriptVariation{
    my $self= shift;
    if (@_){
	if(!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::TranscriptVariation')) {
	    throw("Bio::EnsEMBL::Variation::TranscriptVariation argument expected");
	}
	#a variation feature can have multiple transcript Variations
	push @{$self->{'transcriptVariations'}},shift;
    }
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


=head2 display_consequence

  Arg [1]    : (optional) string $consequence_type
  Example    : $display_consequence = $vf->display_consequence();
  Description: Getter/Setter for the consequence type to display,
               when more than one
  Returntype : string
  Exceptions : throw on incorrect argument
  Caller     : webteam

=cut

sub display_consequence{
    my $self = shift;
    
    if (@_){
	my $display_consequence = shift;
	if (! defined $CONSEQUENCE_TYPES{$display_consequence}){
	    my $valid = join(',',map({"'$_'"} keys(%CONSEQUENCE_TYPES)));
	    throw("Type argument must be one of: $valid");
	}
	return $display_consequence;
    }
    else{
	#get the value to display from the consequence_type attribute
	my $highest_ct = 'INTERGENIC';
	foreach my $ct (@{$self->get_consequence_type}){
	    if ($CONSEQUENCE_TYPES{$ct} < $CONSEQUENCE_TYPES{$highest_ct}){
		$highest_ct = $ct;
	    }
	}
	return $highest_ct;
    }
}
=head2 add_consequence_type

    Arg [1]     : string $consequence_type
    Example     : $vf->add_consequence_type("UPSTREAM")
    Description : Setter for the consequence type of this VariationFeature
                  Allowed values are: 'ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','FRAMESHIFT_CODING',
		  'NON_SYNONYMOUS_CODING','SPLICE_SITE','SYNONYMOUS_CODING','REGULATORY_REGION',
		  '5PRIME_UTR','3PRIME_UTR','INTRONIC','UPSTREAM','DOWNSTREAM','INTERGENIC'
    ReturnType  : string
    Exceptions  : none
    Caller      : general

=cut

sub add_consequence_type{
    my $self = shift;
    my $consequence_type = shift;

    if ($CONSEQUENCE_TYPES{$consequence_type}){
	push @{$self->{'consequence_type'}}, $consequence_type;
	return $self->{'consequence_type'};
    }
    warning("You are trying to set the consequence type to a non-allowed type. The allowed types are: ", keys %CONSEQUENCE_TYPES);
    return '';
}

=head2 get_consequence_type

   Arg[1]      : (optional) Bio::EnsEMBL::Gene $g
   Example     : if($vf->get_consequence_type eq 'INTRONIC'){do_something();}
   Description : Getter for the consequence type of this variation, which is the highest of the transcripts that has.
                 If an argument provided, gets the highest of the transcripts where the gene appears
                 Allowed values are:'ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','FRAMESHIFT_CODING',
		  'NON_SYNONYMOUS_CODING','SPLICE_SITE','SYNONYMOUS_CODING','REGULATORY_REGION',
		  '5PRIME_UTR','3PRIME_UTR','INTRONIC','UPSTREAM','DOWNSTREAM','INTERGENIC'
   Returntype : ref to array of strings
   Exceptions : throw if provided argument not a gene
   Caller     : general

=cut

sub get_consequence_type {
  my $self = shift;
  my $gene = shift;
    
  if(!defined $gene){
    return $self->{'consequence_type'};
  } else{
    my $highest_priority;
    #first, get all the transcripts, if any
    my $transcript_variations = $self->get_all_TranscriptVariations();
    #if no transcripts, return INTERGENIC type
    if (!defined $transcript_variations){
      return 'INTERGENIC';
    }
    if (!ref $gene || !$gene->isa("Bio::EnsEMBL::Gene")){
      throw("$gene is not a Bio::EnsEMBL::Gene type!");
    }
    my $transcripts = $gene->get_all_Transcripts();
	my %transcripts_genes;
	my @new_transcripts;
	map {$transcripts_genes{$_->dbID()}++} @{$transcripts};
	foreach my $transcript_variation (@{$transcript_variations}){
	    if (exists $transcripts_genes{$transcript_variation->transcript->dbID()}){
		push @new_transcripts,$transcript_variation;
	    }
	}
	$highest_priority = $self->_highest_priority(\@new_transcripts);	
	return $highest_priority;
    }

}


#for a list of transcript variations, gets the one with highest priority
sub _highest_priority{
    my $self= shift;
    my $transcript_variations = shift;
    my $highest_type = 'INTERGENIC';
    my $highest_priority; #ref to array with highest types
    foreach my $tv (@{$transcript_variations}){
 	#with a frameshift coding, return, is the highest value
	my $consequences = $tv->consequence_type; #returns a ref to array
	foreach my $consequence_type (@{$consequences}){
	    if (defined $CONSEQUENCE_TYPES{$consequence_type} && $CONSEQUENCE_TYPES{$consequence_type} < $CONSEQUENCE_TYPES{$highest_type}){
		$highest_type = $consequence_type;
		$highest_priority = $tv->consequence_type; 
	    }
	}
    }    
    return $highest_priority;
}

=head2 ambig_code

    Args         : None
    Example      : my $ambiguity_code = $vf->ambig_code()
    Description  : Returns the ambigutiy code for the alleles in the VariationFeature
    ReturnType   : String $ambiguity_code
    Exceptions   : none    
    Caller       : General

=cut 

sub ambig_code{
    my $self = shift;
    
    return &ambiguity_code($self->allele_string());
}

=head2 var_class

    Args         : None
    Example      : my $variation_class = $vf->var_class()
    Description  : returns the class for the variation, according to dbSNP classification
    ReturnType   : String $variation_class
    Exceptions   : none
    Caller       : General
=cut

sub var_class{
    my $self = shift;
    return &variation_class($self->allele_string());
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

=cut

sub get_all_validation_states {
  my $self = shift;

  my @VSTATES = @Bio::EnsEMBL::Variation::Variation::VSTATES;

  my $code = $self->{'validation_code'};
  # convert the validation state strings into a bit field
  # this preserves the same order and representation as in the database
  # and filters out invalid states

  my %VSTATE2BIT = %Bio::EnsEMBL::Variation::Variation::VSTATE2BIT;
  my $vcode = 0;
  $code ||= [];
  foreach my $vstate (@$code) {
    $vcode |= $VSTATE2BIT{lc($vstate)} || 0;
  }

  # convert the bit field into an ordered array
  my @states;
  for(my $i = 0; $i < @VSTATES; $i++) {
    push @states, $VSTATES[$i] if((1 << $i) & $vcode);
  }
  return \@states;
}




=head2 add_validation_state

  Arg [1]    : string $state
  Example    : $vf->add_validation_state('cluster');
  Description: Adds a validation state to this variation.
  Returntype : none
  Exceptions : warning if validation state is not a recognised type
  Caller     : general

=cut

sub add_validation_state {
  my $self  = shift;
  my $state = shift;

  my %VSTATE2BIT = %Bio::EnsEMBL::Variation::Variation::VSTATE2BIT;
  my @VSTATES = @Bio::EnsEMBL::Variation::Variation::VSTATES;
  # convert string to bit value and add it to the existing bitfield
  my $bitval = $VSTATE2BIT{lc($state)};

  if(!$bitval) {
    warning("$state is not a recognised validation status. Recognised " .
            "validation states are: @VSTATES");
    return;
  }

  $self->{'validation_code'} |= $bitval;

  return;
}



=head2 source

  Arg [1]    : string $source (optional)
               The new value to set the source attribute to
  Example    : $source = $vf->source()
  Description: Getter/Setter for the source attribute
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub source{
  my $self = shift;
  return $self->{'source'} = shift if(@_);
  return $self->{'source'};
}

=head2 is_tagged

  Args        : None
  Example     : my $populations = $vf->is_tagged();
  Description : If the variation is tagged in any population, returns an array with the populations where the variation_feature
                is tagged (using a criteria of r2 > 0.99). Otherwise, returns null
  ReturnType  : list of Bio::EnsEMBL::Variation::Population
  Exceptions  : none
  Caller      : general
  
=cut

sub is_tagged{
    my $self = shift;
    
    if ($self->{'adaptor'}){
	my $population_adaptor = $self->{'adaptor'}->db()->get_PopulationAdaptor();
	return $population_adaptor->fetch_tagged_Population($self);
    }
}

=head2 convert_to_SNP

  Args        : None
  Example     : my $snp = $vf->convert_to_SNP()
  Description : Creates a Bio::EnsEMBL::SNP object from Bio::EnsEMBL::VariationFeature. Mainly used for
                backwards comnpatibility
  ReturnType  : Bio::EnsEMBL::SNP
  Exceptions  : None
  Caller      : general      

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
    
    if ($self->{'adaptor'}){
	my $ld_adaptor = $self->{'adaptor'}->db()->get_LDFeatureContainerAdaptor();
	return $ld_adaptor->fetch_by_VariationFeature($self);
    }
    return {};
}

1;
