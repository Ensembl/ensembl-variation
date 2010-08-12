# Ensembl module for Bio::EnsEMBL::Variation::Variation
#
# Copyright (c) 2004 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::Variation - Ensembl representation of a nucleotide variation.

=head1 SYNOPSIS

    $v = Bio::EnsEMBL::Variation::Variation->new(-name   => 'rs123',
                                                 -source => 'dbSNP');

    # add additional synonyms for the same SNP
    $v->add_synonym('dbSNP', 'ss3242');
    $v->add_synonym('TSC', '53253');

    # add some validation states for this SNP
    $v->add_validation_status('freq');
    $v->add_validation_status('cluster');

    # add alleles associated with this SNP
    $a1 = Bio::EnsEMBL::Allele->new(...);
    $a2 = Bio::EnsEMBL::Allele->new(...);
    $v->add_Allele($a1);
    $v->add_Allele($a2);

    # set the flanking sequences
    $v->five_prime_flanking_seq($seq);
    $v->three_prime_flanking_seq($seq);


    ...

    # print out the default name and source of the variation and the version
    print $v->source(), ':',$v->name(), ".",$v->version,"\n";

    # print out every synonym associated with this variation
    @synonyms = @{$v->get_all_synonyms()};
    print "@synonyms\n";

    # print out synonyms and their database associations
    my $sources = $v->get_all_synonym_sources();
    foreach my $src (@$sources) {
      @synonyms = $v->get_all_synonyms($src);
      print "$src: @synonyms\n";
    }


    # print out validation states
    my @vstates = @{$v->get_all_validation_states()};
    print "@validation_states\n";

    # print out flanking sequences
    print "5' flanking: ", $v->five_prime_flanking_seq(), "\n";
    print "3' flanking: ", $v->three_prime_flanking_seq(), "\n";


=head1 DESCRIPTION

This is a class representing a nucleotide variation from the
ensembl-variation database. A variation may be a SNP a multi-base substitution
or an insertion/deletion.  The objects Alleles associated with a Variation
object describe the nucleotide change that Variation represents.

A Variation object has an associated identifier and 0 or more additional
synonyms.  The position of a Variation object on the Genome is represented
by the B<Bio::EnsEMBL::Variation::VariationFeature> class.

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::Variation;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(ambiguity_code variation_class);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);


# List of validation states. Order must match that of set in database
our @VSTATES = ('cluster','freq','submitter','doublehit','hapmap', 'failed', 'non-polymorphic', 'observed');

# Conversion of validation state to bit value
our %VSTATE2BIT = ('cluster'   => 1,   # 00000001
                   'freq'      => 2,   # 00000010
                   'submitter' => 4,   # 00000100
                   'doublehit' => 8,   # 00001000
                   'hapmap'    => 16,  # 00010000
		   'failed'    => 32,
		   'non-polymorphic' => 64, # 00100000
		   'observed'  => 128   # 01000000
		   
		   ); 


=head2 new

  Arg [-dbID] :
    int - unique internal identifier for snp

  Arg [-ADAPTOR] :
    Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor
    Adaptor which provides database connectivity for this Variation object

  Arg [-NAME] :
    string - the name of this SNP

  Arg [-SOURCE] :
    string - the source of this SNP

  Arg [-SOURCE_DESCRIPTION] :
    string - description of the SNP source

  Arg [-SYNONYMS] :
    reference to hash with list reference values -  keys are source
    names and values are lists of identifiers from that db.
    e.g.: {'dbSNP' => ['ss1231', '1231'], 'TSC' => ['1452']}

  Arg [-ANCESTRAL_ALLELES] :
    string - the ancestral allele of this SNP

  Arg [-ALLELES] :
    reference to list of Bio::EnsEMBL::Variation::Allele objects

  Arg [-VALIDATION_STATES] :
    reference to list of strings

  Arg [-MOLTYPE] :
    string - the moltype of this SNP

  Arg [-FIVE_PRIME_FLANKING_SEQ] :
    string - the five prime flanking nucleotide sequence

  Arg [-THREE_PRIME_FLANKING_SEQ] :
    string - the three prime flanking nucleotide sequence

  Arg [-FAILED_DESCRIPTION] :
    string - description why the variation failed Ensembl pipeline

  Example    : $v = Bio::EnsEMBL::Variation::Variation->new
                    (-name   => 'rs123',
                     -source => 'dbSNP');

  Description: Constructor. Instantiates a new Variation object.
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut


sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($dbID, $adaptor, $name, $src, $src_desc, $src_url, $is_somatic, $syns, $ancestral_allele,
      $alleles, $valid_states, $moltype, $five_seq, $three_seq, $failed_description, $flank_flag) =
        rearrange([qw(dbID ADAPTOR NAME SOURCE SOURCE_DESCRIPTION SOURCE_URL IS_SOMATIC SYNONYMS ANCESTRAL_ALLELE ALLELES
                      VALIDATION_STATES MOLTYPE FIVE_PRIME_FLANKING_SEQ
                      THREE_PRIME_FLANKING_SEQ FAILED_DESCRIPTION FLANK_FLAG)],@_);


  # convert the validation state strings into a bit field
  # this preserves the same order and representation as in the database
  # and filters out invalid states
  my $vcode = 0;
  $valid_states ||= [];
  foreach my $vstate (@$valid_states) {
    $vcode |= $VSTATE2BIT{lc($vstate)} || 0;
  }

  return bless {'dbID' => $dbID,
                'adaptor' => $adaptor,
                'name'   => $name,
                'source' => $src,
				'source_description' => $src_desc,
				'source_url' => $src_url,
				'is_somatic' => $is_somatic,
                'synonyms' => $syns || {},
		        'ancestral_allele' => $ancestral_allele,
                'alleles' => $alleles || [],
                'validation_code' => $vcode,
		        'moltype' => $moltype,
                'five_prime_flanking_seq' => $five_seq,
                'three_prime_flanking_seq' => $three_seq,
	            'failed_description' => $failed_description,
			    'flank_flag' => $flank_flag}, $class;
}



=head2 name

  Arg [1]    : string $newval (optional)
               The new value to set the name attribute to
  Example    : $name = $obj->name()
  Description: Getter/Setter for the name attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub name{
  my $self = shift;
  return $self->{'name'} = shift if(@_);
  return $self->{'name'};
}


=head2 get_all_Genes

  Args        : None
  Example     : $genes = $v->get_all_genes();
  Description : Retrieves all the genes where this Variation
                has a consequence.
  ReturnType  : reference to list of Bio::EnsEMBL::Gene
  Exceptions  : None
  Caller      : general
  Status      : At Risk

=cut

sub get_all_Genes{
    my $self = shift;
    my $genes;
    if (defined $self->{'adaptor'}){
	my $UPSTREAM = 5000;
	my $DOWNSTREAM = 5000;
	my $vf_adaptor = $self->adaptor()->db()->get_VariationFeatureAdaptor();
	my $vf_list = $vf_adaptor->fetch_all_by_Variation($self);
	#foreach vf, get the slice is on, us ethe USTREAM and DOWNSTREAM limits to get all the genes, and see if SNP is within the gene
	my $new_slice;
	my $gene_list;
	my $gene_hash;

	foreach my $vf (@{$vf_list}){
	    #expand the slice UPSTREAM and DOWNSTREAM
	    $new_slice = $vf->feature_Slice()->expand($UPSTREAM,$DOWNSTREAM);
	    #get the genes in the new slice
	    $gene_list = $new_slice->get_all_Genes();
	    foreach my $gene (@{$gene_list}){
		if (($vf->start >= $gene->seq_region_start - $UPSTREAM) && ($vf->start <= $gene->seq_region_end + $DOWNSTREAM) && ($vf->end <= $gene->seq_region_end + $DOWNSTREAM)){
		    #the vf is affecting the gene, add to the hash if not present already
		    if (!exists $gene_hash->{$gene->dbID}){
			$gene_hash->{$gene->dbID} = $gene;
		    }
		}
	    }
	}
	#and return all the genes
	push @{$genes}, values %{$gene_hash};
    }
    return $genes;
}




=head2 get_all_VariationFeatures

  Args        : None
  Example     : $vfs = $v->get_all_VariationFeatures();
  Description : Retrieves all VariationFeatures for this Variation
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::VariationFeature
  Exceptions  : None
  Caller      : general
  Status      : At Risk

=cut

sub get_all_VariationFeatures{
  my $self = shift;
  
  if(defined $self->adaptor) {
	
	# get variation feature adaptor
	my $vf_adaptor = $self->adaptor()->db()->get_VariationFeatureAdaptor();
	
	return $vf_adaptor->fetch_all_by_Variation($self);
  }
  
  else {
	warn("No variation database attached");
	return [];
  }
}

=head2 get_VariationFeature_by_dbID

  Args        : None
  Example     : $vf = $v->get_VariationFeature_by_dbID();
  Description : Retrieves a VariationFeature for this Variation by it's internal
				database identifier
  ReturnType  : Bio::EnsEMBL::Variation::VariationFeature
  Exceptions  : None
  Caller      : general
  Status      : At Risk

=cut

sub get_VariationFeature_by_dbID{
  my $self = shift;
  my $dbID = shift;
  
  throw("No dbID defined") unless defined $dbID;
  
  if(defined $self->adaptor) {
	
	# get variation feature adaptor
	my $vf_adaptor = $self->adaptor()->db()->get_VariationFeatureAdaptor();
	
	my $vf = $vf_adaptor->fetch_by_dbID($dbID);
	
	# check defined
	if(defined($vf)) {
	  
	  # check it is the same variation ID
	  if($vf->{_variation_id} == $self->dbID) {
		return $vf;
	  }
	  
	  else {
		warn("Variation dbID for Variation Feature does not match this Variation's dbID");
		return undef;
	  }
	}
	
	else {
	  return undef;
	}
  }
  
  else {
	warn("No variation database attached");
	return undef;
  }  
}



=head2 get_all_synonyms

  Arg [1]    : (optional) string $source - the source of the synonyms to
               return.
  Example    : @dbsnp_syns = @{$v->get_all_synonyms('dbSNP')};
               @all_syns = @{$v->get_all_synonyms()};
  Description: Retrieves synonyms for this Variation. If a source argument
               is provided all synonyms from that source are returned,
               otherwise all synonyms are returned.
  Returntype : reference to list of strings
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub get_all_synonyms {
  my $self = shift;
  my $source = shift;

  if($source) {
    return $self->{'synonyms'}->{$source} || []
  }

  my @synonyms = map {@$_} values %{$self->{'synonyms'}};

  return \@synonyms;
}



=head2 get_all_synonym_sources

  Arg [1]    : none
  Example    : my @sources = @{$v->get_all_synonym_sources()};
  Description: Retrieves a list of all the sources for synonyms of this
               Variation.
  Returntype : reference to a list of strings
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub get_all_synonym_sources {
  my $self = shift;
  my @sources = keys %{$self->{'synonyms'}};
  return \@sources;
}



=head2 add_synonym

  Arg [1]    : string $source
  Arg [2]    : string $syn
  Example    : $v->add_synonym('dbSNP', 'ss55331');
  Description: Adds a synonym to this variation.
  Returntype : none
  Exceptions : throw if $source argument is not provided
               throw if $syn argument is not provided
  Caller     : general
  Status     : At Risk

=cut

sub add_synonym {
  my $self   = shift;
  my $source = shift;
  my $syn    = shift;

  throw("source argument is required") if(!$source);
  throw("syn argument is required") if(!$syn);

  $self->{'synonyms'}->{$source} ||= [];

  push @{$self->{'synonyms'}->{$source}}, $syn;

  return;
}



=head2 get_all_validation_states

  Arg [1]    : none
  Example    : my @vstates = @{$v->get_all_validation_states()};
  Description: Retrieves all validation states for this variation.  Current
               possible validation statuses are 'cluster','freq','submitter',
               'doublehit', 'hapmap'
  Returntype : reference to list of strings
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub get_all_validation_states {
  my $self = shift;


  my $code = $self->{'validation_code'};

  # convert the bit field into an ordered array
  my @states;
  for(my $i = 0; $i < @VSTATES; $i++) {
    push @states, $VSTATES[$i] if((1 << $i) & $code);
  }

  return \@states;
}




=head2 add_validation_state

  Arg [1]    : string $state
  Example    : $v->add_validation_state('cluster');
  Description: Adds a validation state to this variation.
  Returntype : none
  Exceptions : warning if validation state is not a recognised type
  Caller     : general
  Status     : At Risk

=cut

sub add_validation_state {
  my $self  = shift;
  my $state = shift;

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
  Example    : $source = $v->source()
  Description: Getter/Setter for the source attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source{
  my $self = shift;
  return $self->{'source'} = shift if(@_);
  return $self->{'source'};
}



=head2 source_description

  Arg [1]    : string $source_description (optional)
               The new value to set the source description attribute to
  Example    : $source_description = $v->source_description()
  Description: Getter/Setter for the source description attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source_description{
  my $self = shift;
  return $self->{'source_description'} = shift if(@_);
  return $self->{'source_description'};
}



=head2 source_url

  Arg [1]    : string $source_url (optional)
               The new value to set the source URL attribute to
  Example    : $source_url = $v->source_url()
  Description: Getter/Setter for the source URL attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source_url{
  my $self = shift;
  return $self->{'source_url'} = shift if(@_);
  return $self->{'source_url'};
}

=head2 is_somatic

  Arg [1]    : boolean $is_somatic (optional)
               The new value to set the is_somatic flag to
  Example    : $is_somatic = $v->is_somatic
  Description: Getter/Setter for the is_somatic flag, which identifies this variation as either somatic or germline
  Returntype : boolean
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub is_somatic {
  my ($self, $is_somatic) = @_;
  $self->{is_somatic} = $is_somatic if defined $is_somatic;
  return $self->{is_somatic};
}

=head2 get_all_Alleles

  Arg [1]    : none
  Example    : @alleles = @{$v->get_all_Alleles()};
  Description: Retrieves all Alleles associated with this variation
  Returntype : reference to list of Bio::EnsEMBL::Variation::Allele objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Alleles {
  my $self = shift;
  return $self->{'alleles'};
}



=head2 add_Allele

  Arg [1]    : Bio::EnsEMBL::Variation::Allele $allele
  Example    : $v->add_Allele(Bio::EnsEMBL::Variation::Alelele->new(...));
  Description: Associates an Allele with this variation
  Returntype : none
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk

=cut

sub add_Allele {
  my $self = shift;
  my $allele = shift;

  if(!ref($allele) || !$allele->isa('Bio::EnsEMBL::Variation::Allele')) {
    throw("Bio::EnsEMBL::Variation::Allele argument expected");
  }

  push @{$self->{'alleles'}}, $allele;
}

=head2 ancestral_allele

  Arg [1]    : string $ancestral_allele (optional)
  Example    : $ancestral_allele = v->ancestral_allele();
  Description: Getter/Setter ancestral allele associated with this variation
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub ancestral_allele {
  my $self = shift;
  return $self->{'ancestral_allele'} = shift if(@_);
  return $self->{'ancestral_allele'};
}

=head2 moltype

  Arg [1]    : string $moltype (optional)
               The new value to set the moltype attribute to
  Example    : $moltype = v->moltype();
  Description: Getter/Setter moltype associated with this variation
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub moltype {
  my $self = shift;
  return $self->{'moltype'} = shift if(@_);
  return $self->{'moltype'};
}


=head2 five_prime_flanking_seq

  Arg [1]    : string $newval (optional) 
               The new value to set the five_prime_flanking_seq attribute to
  Example    : $five_prime_flanking_seq = $obj->five_prime_flanking_seq()
  Description: Getter/Setter for the five_prime_flanking_seq attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub five_prime_flanking_seq{
  my $self = shift;

  #setter of the flanking sequence
  return $self->{'five_prime_flanking_seq'} = shift if(@_);
  #lazy-load the flanking sequence from the database
  if (!defined $self->{'five_prime_flanking_seq'} && $self->{'adaptor'}){
      my $variation_adaptor = $self->adaptor()->db()->get_VariationAdaptor();
      ($self->{'three_prime_flanking_seq'},$self->{'five_prime_flanking_seq'}) = @{$variation_adaptor->get_flanking_sequence($self->{'dbID'})};
  }
  return $self->{'five_prime_flanking_seq'};
}




=head2 three_prime_flanking_seq

  Arg [1]    : string $newval (optional) 
               The new value to set the three_prime_flanking_seq attribute to
  Example    : $three_prime_flanking_seq = $obj->three_prime_flanking_seq()
  Description: Getter/Setter for the three_prime_flanking_seq attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub three_prime_flanking_seq{
  my $self = shift;

  #setter of the flanking sequence
  return $self->{'three_prime_flanking_seq'} = shift if(@_);
  #lazy-load the flanking sequence from the database
  if (!defined $self->{'three_prime_flanking_seq'} && $self->{'adaptor'}){
      my $variation_adaptor = $self->adaptor()->db()->get_VariationAdaptor();
      ($self->{'three_prime_flanking_seq'},$self->{'five_prime_flanking_seq'}) = @{$variation_adaptor->get_flanking_sequence($self->{'dbID'})};
  }
  return $self->{'three_prime_flanking_seq'};
}


=head2 get_all_IndividualGenotypes

  Args       : none
  Example    : $ind_genotypes = $var->get_all_IndividualGenotypes()
  Description: Getter for IndividualGenotypes for this Variation, returns empty list if 
               there are none 
  Returntype : listref of IndividualGenotypes
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub get_all_IndividualGenotypes {
    my $self = shift;
    if (defined ($self->{'adaptor'})){
	my $igtya = $self->{'adaptor'}->db()->get_IndividualGenotypeAdaptor();
	
	return $igtya->fetch_all_by_Variation($self);
    }
    return [];
}

=head2 get_all_PopulationGenotypes

  Args       : none
  Example    : $pop_genotypes = $var->get_all_PopulationGenotypes()
  Description: Getter for PopulationGenotypes for this Variation, returns empty list if 
               there are none. 
  Returntype : listref of PopulationGenotypes
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub get_all_PopulationGenotypes {
    my $self = shift;

    #simulate a lazy-load on demand situation, used by the Glovar team
    if (!defined($self->{'populationGenotypes'}) && defined ($self->{'adaptor'})){
	my $pgtya = $self->{'adaptor'}->db()->get_PopulationGenotypeAdaptor();
	
	return $pgtya->fetch_all_by_Variation($self);
    }
    return $self->{'populationGenotypes'};

}


=head2 add_PopulationGenotype

    Arg [1]     : Bio::EnsEMBL::Variation::PopulationGenotype
    Example     : $v->add_PopulationGenotype($pop_genotype)
    Description : Adds another PopulationGenotype to the Variation object
    Exceptions  : thrown on bad argument
    Caller      : general
    Status      : At Risk

=cut

sub add_PopulationGenotype{
    my $self = shift;

    if (@_){
	if(!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::PopulationGenotype')) {
	    throw("Bio::EnsEMBL::Variation::PopulationGenotype argument expected");
	}
	#a variation can have multiple PopulationGenotypes
	push @{$self->{'populationGenotypes'}},shift;
    }

}


=head2 ambig_code

    Args         : None
    Example      : my $ambiguity_code = $v->ambig_code()
    Description  : Returns the ambigutiy code for the alleles in the Variation
    ReturnType   : String $ambiguity_code
    Exceptions   : none    
    Caller       : General
    Status       : At Risk

=cut 

sub ambig_code{
    my $self = shift;
    my $alleles = $self->get_all_Alleles(); #get all Allele objects
    my %alleles; #to get all the different alleles in the Variation
    map {$alleles{$_->allele}++} @{$alleles};
    my $allele_string = join "|",keys %alleles;
    return &ambiguity_code($allele_string);
}

=head2 var_class

    Args         : None
    Example      : my $variation_class = $vf->var_class()
    Description  : returns the class for the variation, according to dbSNP classification
    ReturnType   : String $variation_class
    Exceptions   : none
    Caller       : General
    Status       : At Risk

=cut

sub var_class{
    my $self = shift;
    my $alleles = $self->get_all_Alleles(); #get all Allele objects
    my %alleles; #to get all the different alleles in the Variation
    map {$alleles{$_->allele}++} @{$alleles};
    my $allele_string = join "|",keys %alleles;
    return &variation_class($allele_string, $self->is_somatic);
}

=head2 failed_description

  Arg [1]    : string $failed_description (optional)
               The new value to set the failed_description attribute to
  Example    : $failed_description = $v->failed_description()
  Description: Getter/Setter for the failed_description attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub failed_description{
  my $self = shift;
  return $self->{'failed_description'} = shift if(@_);
  return $self->{'failed_description'};
}

=head2 derived_allele_frequency

  Arg[1]     : Bio::EnsEMBL::Variation::Population  $population 
  Example    : $daf = $variation->derived_allele_frequency($population);
  Description: Gets the derived allele frequency for the population. 
               The DAF is the frequency of the reference allele that is 
               different from the allele in Chimp. If none of the alleles
               is the same as the ancestral, will return reference allele
               frequency
  Returntype : float
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub derived_allele_frequency{
  my $self = shift;
  my $population = shift;
  my $daf;

  if(!ref($population) || !$population->isa('Bio::EnsEMBL::Variation::Population')) {
      throw('Bio::EnsEMBL::Variation::Population argument expected.');
  }
  my $ancestral_allele = $self->ancestral_allele();
  if (defined $ancestral_allele){
	#get reference allele
	my $vf_adaptor = $self->adaptor->db->get_VariationFeatureAdaptor();
	my $vf = shift @{$vf_adaptor->fetch_all_by_Variation($self)};
	my $ref_freq;
	#get allele in population
	my $alleles = $self->get_all_Alleles();
	
	foreach my $allele (@{$alleles}){
	  next unless defined $allele->population;
	  
	  if (($allele->allele eq $vf->ref_allele_string) and ($allele->population->name eq $population->name)){
		$ref_freq = $allele->frequency;
	  }
	}
	
	if(defined $ref_freq) {
	  if ($ancestral_allele eq $vf->ref_allele_string){
		$daf = 1 - $ref_freq
	  }
	  elsif ($ancestral_allele ne $vf->ref_allele_string){
		$daf = $ref_freq;
	  }
	}
  }
  
  return $daf;
}

=head2 derived_allele

  Arg[1]     : Bio::EnsEMBL::Variation::Population  $population 
  Example    : $da = $variation->derived_allele($population);
  Description: Gets the derived allele for the population. 
  Returntype : float
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub derived_allele {
     my $self = shift();
     my $population = shift();

     my $population_dbID = $population->dbID();
     my $ancestral_allele_str = $self->ancestral_allele();

     if (not defined($ancestral_allele_str)) {
         return;
     }

     my $alleles = $self->get_all_Alleles();

     my $derived_allele_str;

     foreach my $allele (@{$alleles}) {
         my $allele_population = $allele->population();

         if (defined($allele_population) and
             $allele_population->dbID() == $population_dbID)
         {
             my $allele_str = $allele->allele();

             if ($ancestral_allele_str ne $allele_str) {
                 if (defined($derived_allele_str)) {
                     return;
                 } else {
                     $derived_allele_str = $allele_str;
                 }
             }
         }
     }
     return $derived_allele_str;
}

1;
