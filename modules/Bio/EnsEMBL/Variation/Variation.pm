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



=head1 NAME

Bio::EnsEMBL::Variation::Variation - Ensembl representation of a nucleotide variation.

=head1 SYNOPSIS

    $v = Bio::EnsEMBL::Variation::Variation->new(-name   => 'rs123',
                                                 -source => $source);

    # add additional synonyms for the same SNP
    $v->add_synonym('dbSNP', 'ss3242');
    $v->add_synonym('TSC', '53253');

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
    print $v->source_name(), ':',$v->name(), ".",$v->source_version(),"\n";

    # print out every synonym associated with this variation
    @synonyms = @{$v->get_all_synonyms()};
    print "@synonyms\n";

    # print out synonyms and their database associations
    my $sources = $v->get_all_synonym_sources();
    foreach my $src (@$sources) {
      @synonyms = $v->get_all_synonyms($src);
      print "$src: @synonyms\n";
    }

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

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::Variation;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref wrap_array);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(ambiguity_code SO_variation_class);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);
use Bio::EnsEMBL::Variation::Utils::Sequence;
use Bio::EnsEMBL::Variation::Utils::Constants qw(%VARIATION_CLASSES); 
use Bio::EnsEMBL::Variation::Utils::Sequence  qw(%EVIDENCE_VALUES); 
use Bio::EnsEMBL::Variation::Failable;
use vars qw(@ISA);
use Scalar::Util qw(weaken);

@ISA = qw(Bio::EnsEMBL::Storable Bio::EnsEMBL::Variation::Failable);

=head2 new

  Arg [-dbID] :
    int - unique internal identifier for snp

  Arg [-ADAPTOR] :
    Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor
    Adaptor which provides database connectivity for this Variation object

  Arg [-NAME] :
    string - the name of this variant

  Arg [-SOURCE] :
    object ref - the source object describing where the variant comes from. 

  Arg [-SYNONYMS] :
    reference to hash with list reference values -  keys are source
    names and values are lists of identifiers from that db.
    e.g.: {'dbSNP' => ['ss1231', '1231'], 'TSC' => ['1452']}

  Arg [-ANCESTRAL_ALLELE] :
    string - the ancestral allele of this SNP

  Arg [-ALLELES] :
    Listref of Bio::EnsEMBL::Variation::Allele objects

  Arg [-EVIDENCE] :
    Listref of strings

  Arg [-FIVE_PRIME_FLANKING_SEQ] :
    string - the five prime flanking nucleotide sequence

  Arg [-THREE_PRIME_FLANKING_SEQ] :
    string - the three prime flanking nucleotide sequence

  Example    : $v = Bio::EnsEMBL::Variation::Variation->new
                    (-name   => 'rs123',
                     -source => $source);

  Description: Constructor. Instantiates a new Variation object.
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($dbID, $adaptor, $name, $class_so_term, $source_id, $source, $is_somatic, $flipped, $syns,
      $ancestral_allele, $alleles, $five_seq, $three_seq, $flank_flag, $minor_allele, $minor_allele_frequency,
      $minor_allele_count, $clinical_significance, $evidence ) =
        rearrange([qw(dbID ADAPTOR NAME CLASS_SO_TERM _SOURCE_ID SOURCE IS_SOMATIC FLIPPED SYNONYMS ANCESTRAL_ALLELE ALLELES
                      FIVE_PRIME_FLANKING_SEQ THREE_PRIME_FLANKING_SEQ FLANK_FLAG 
                      MINOR_ALLELE MINOR_ALLELE_FREQUENCY MINOR_ALLELE_COUNT CLINICAL_SIGNIFICANCE EVIDENCE)],@_);

  my $self = bless {
    'dbID' => $dbID,
    'adaptor' => $adaptor,
    'name'   => $name,
    'class_SO_term' => $class_so_term,
    '_source_id' => $source_id,
    'source' => $source,
    'is_somatic' => $is_somatic,
    'flipped' => $flipped,
    'synonyms' => $syns || {},
    'ancestral_allele' => $ancestral_allele,
    'five_prime_flanking_seq' => $five_seq,
    'three_prime_flanking_seq' => $three_seq,
    'flank_flag' => $flank_flag,
    'minor_allele' => $minor_allele,
    'minor_allele_frequency' => $minor_allele_frequency,
    'minor_allele_count' => $minor_allele_count,
    'clinical_significance' => $clinical_significance,
    'evidence' => $evidence
  }, $class;
  
  $self->add_Allele($alleles) if defined($alleles);
  
  return $self;
}

sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}


=head2 has_failed_alleles

  Example    : print "Variation '" . $var->name() . "' has " . ($var->has_failed_alleles() ? "" : "no ") . " failed alleles\n";
  Description: Returns true if this variation has alleles that are flagged as failed
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub has_failed_alleles {
    my $self = shift;
  
    map {return 1 if ($_->is_failed())} @{$self->get_all_Alleles()};
    return 0;
}


=head2 add_Allele

  Arg [1]    : Bio::EnsEMBL::Variation::Allele $allele 
  Example    : $v->add_allele(Bio::EnsEMBL::Variation::Allele->new(...));
  Description: Add an Allele to this variation.
  Returntype : none
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk

=cut

sub add_Allele {
    my $self = shift;
    my $allele = shift;

    # This method also accepts a list of alleles so wrap the argument in an array and treat as such
    $allele = wrap_array($allele);  
    map {assert_ref($_,'Bio::EnsEMBL::Variation::Allele')} @{$allele};
    
    # Store the allele in the private hash
    $self->{alleles} = [] unless (exists($self->{alleles}));
    push(@{$self->{alleles}},@{$allele});
    
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

sub stable_id {
  my $self = shift;
  return $self->name(@_);
}


=head2 get_all_Genes

  Args        : None
  Example     : $genes = $v->get_all_genes();
  Description : Retrieves all the genes where this Variation
                has a consequence.
  ReturnType  : Listref of Bio::EnsEMBL::Gene
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
  ReturnType  : Listref of Bio::EnsEMBL::Variation::VariationFeature
  Exceptions  : None
  Caller      : general
  Status      : Stable

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
  Returntype : Listref of strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_synonyms {
    my $self = shift;
    my $source = shift;

    if ($source) {
        $source = [$source];
    }
    else {
        $source = $self->get_all_synonym_sources();
    }
    
    my @synonyms;
    map {push(@synonyms,keys(%{$self->{synonyms}{$_}}))} @{$source};

    return \@synonyms;
}



=head2 get_all_synonym_sources

  Arg [1]    : none
  Example    : my @sources = @{$v->get_all_synonym_sources()};
  Description: Retrieves a list of all the sources for synonyms of this
               Variation.
  Returntype : Listref of strings
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

  $self->{'synonyms'}{$source}{$syn}++;

  return;
}


=head2 get_all_evidence_values

  Arg [1]    : none
  Example    : my @evidence = @{$v->get_all_evidence_values()};
  Description: Retrieves all evidence values for this variation. Current
               possible evidence values are 'Multiple_observations', 'Frequency',
               'HapMap', '1000Genomes', 'ESP', 'Cited', 'Phenotype_or_Disease', 'ExAC'
  Returntype : Listref of strings
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub get_all_evidence_values {
    my $self = shift;
    return $self->{'evidence'} || [];
}


=head2 add_evidence_value

  Arg [1]    : string $state
  Example    : $v->add_evidence_value('Frequency');
  Description: Adds an evidence value  to this variation.
  Returntype : none
  Exceptions : 
  Caller     : general
  Status     : At Risk

=cut

sub add_evidence_value {
    
    my $self = shift;
    my $add_ev = shift if(@_);

    ## do not add evidence value unless it is in the list of permitted values
    return $self->{'evidence'} unless $EVIDENCE_VALUES{$add_ev};

    push @{$self->{'evidence'}}, $add_ev;
    my %unique = map { $_ => 1 } @{$self->{'evidence'}};
    @{$self->{'evidence'}} = keys %unique;

    return $self->{'evidence'};    
}


=head2 source

  Arg [1]    : Bio::EnsEMBL::Variation::Source $src (optional)
               The new value to set the source attribute to
  Example    : $source = $v->source()
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
  Example    : $source_name = $v->source_name()
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


=head2 source_type

  Arg [1]    : string $source_type (optional)
               The new value to set the source type attribute to
  Example    : $source_type = $v->source_type()
  Description: Getter/Setter for the source type attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub source_type{
  my $self = shift;
  my $source = $self->source;
  return unless defined $source;
  
  $source->type(@_) if(@_);
  return $source->type;
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
  my $source = $self->source;
  return unless defined $source;
  
  $source->description(@_) if(@_);
  return $source->description;
}



=head2 source_version

  Arg [1]    : string $source_version (optional)
               The new value to set the source version attribute to
  Example    : $source_version = $v->source_version()
  Description: Getter/Setter for the source version attribute
  Returntype : int
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
  my $source = $self->source;
  return unless defined $source;
  
  $source->url(@_) if(@_);
  return $source->url;
}

=head2 source_somatic_status

  Arg [1]    : string $source_somatic_status (optional)
               The new value to set the source status attribute to
  Example    : $source_somatic_status = $v->source_somatic_status()
  Description: Getter/Setter for the source somatic status attribute, which identifies 
               the source of this variation as somatic, germline or mixed
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub source_somatic_status{
  my $self = shift;
  my $source = $self->source;
  return unless defined $source;
  
  $source->somatic_status(@_) if(@_);
  return $source->somatic_status;
}

=head2 has_somatic_source

  Arg [1]    : boolean $has_somatic_source (optional)
               The new value to set the has_somatic_source flag to
  Example    : $has_somatic_source = $v->has_somatic_source
  Description: Getter/Setter for the has_somatic_source flag, which identifies if this variation 
               comes from a somatic or a germline/mixed source
  Returntype : boolean
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub has_somatic_source {
  my ($self, $has_somatic_source) = @_;
  $self->{has_somatic_source} = (defined $has_somatic_source) ? $has_somatic_source : ($self->source_somatic_status && $self->source_somatic_status eq 'somatic' ? 1 : 0);
  return $self->{has_somatic_source};
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

=head2 flipped

  Arg [1]    : boolean $flipped (optional)
               The new value to set the flipped flag to
  Example    : $flipped = $v->flipped
  Description: Getter/Setter for the flipped flag, which identifies if this
               variation's strand has been flipped during the import process
  Returntype : boolean
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub flipped {
  my ($self, $flipped) = @_;
  $self->{flipped} = $flipped if defined $flipped;
  return $self->{flipped};
}

=head2 get_all_Alleles

  Arg [1]    : Bio::EnsEMBL::Variation::Population $pop (optional)
  Example    : @alleles = @{$v->get_all_Alleles()};
  Description: Retrieves all Alleles associated with this variation
  Returntype : Listref of Bio::EnsEMBL::Variation::Allele objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Alleles {
  my $self = shift;
  my $pop  = shift;
  
  # If the private hash key 'alleles' does not exist, no attempt has been made to load them, so do that
  unless (exists($self->{alleles})) {
      
    # Get an AlleleAdaptor
    assert_ref($self->adaptor(),'Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');
    my $allele_adaptor = $self->adaptor->db->get_AlleleAdaptor();
      
    $self->add_Allele($allele_adaptor->fetch_all_by_Variation($self, $pop));
  }
  
  if($pop && assert_ref($pop, 'Bio::EnsEMBL::Variation::Population')) {
    return [grep {$_->population && $_->population->dbID == $pop->dbID} @{$self->{alleles}}];
  }

  return $self->{alleles};
}



=head2 ancestral_allele

  Arg [1]    : string $ancestral_allele (optional)
  Example    : $ancestral_allele = v->ancestral_allele();
  Description: Getter/Setter ancestral allele associated with this variation
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub ancestral_allele {
  my $self = shift;
  return $self->{'ancestral_allele'} = shift if(@_);
  return $self->{'ancestral_allele'};
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


=head2 get_all_SampleGenotypes

  Args       : none
  Example    : $sample_genotypes = $var->get_all_SampleGenotypes()
  Description: Getter for SampleGenotypes for this Variation, returns empty list if 
               there are none 
  Returntype : Listref of SampleGenotypes
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_SampleGenotypes {
  my $self = shift;
  my $individual = shift;
  
  if (defined ($self->{'adaptor'})){
    my $sgtya = $self->{'adaptor'}->db()->get_SampleGenotypeAdaptor();
    
    return $sgtya->fetch_all_by_Variation($self, $individual);
  }
  
  return [];
}

=head2 get_all_PopulationGenotypes

  Arg [1]    : Bio::EnsEMBL::Variation::Population $pop (optional)
  Example    : $pop_genotypes = $var->get_all_PopulationGenotypes()
  Description: Getter for PopulationGenotypes for this Variation, returns empty list if 
               there are none. 
  Returntype : Listref of PopulationGenotypes
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_PopulationGenotypes {
  my $self = shift;
  my $pop  = shift;

  #simulate a lazy-load on demand situation, used by the Glovar team
  if (!defined($self->{populationGenotypes}) && defined ($self->{'adaptor'})){
    my $pgtya = $self->adaptor->db()->get_PopulationGenotypeAdaptor();
  
    $self->{populationGenotypes} = $pgtya->fetch_all_by_Variation($self, $pop);
  }
  
  if($pop && assert_ref($pop, 'Bio::EnsEMBL::Variation::Population')) {
    return [grep {$_->population && $_->population->dbID == $pop->dbID} @{$self->{populationGenotypes}}];
  }
  
  return $self->{populationGenotypes};
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

    Args        : None
    Example     : my $ambiguity_code = $v->ambig_code()
    Description : Returns the ambigutiy code for the alleles in the Variation
    ReturnType  : String $ambiguity_code
    Exceptions  : none
    Caller      : General
    Status      : Stable

=cut 

sub ambig_code{
    my $self = shift;
  
  my $code;
  
  # first try via VF
  if(my @vfs = @{$self->get_all_VariationFeatures}) {
    if(scalar @vfs) {
    $code = $vfs[0]->ambig_code;
    }
  }
  
  # otherwise get it via alleles attatched to this object already
  if(!defined($code)) {
    my $alleles = $self->get_all_Alleles(); #get all Allele objects
    my %alleles; #to get all the different alleles in the Variation
    map {$alleles{$_->allele}++} @{$alleles};
    my $allele_string = join "|",keys %alleles;
    $code = &ambiguity_code($allele_string);
  }
  
  return $code;
}

=head2 var_class

    Args        : None
    Example     : my $variation_class = $vf->var_class()
    Description : returns the class for the variation, according to dbSNP classification
    ReturnType  : String $variation_class
    Exceptions  : none
    Caller      : General
    Status      : Stable

=cut

sub var_class{
    my $self = shift;
    
    unless ($self->{class_display_term}) {
       
        unless ($self->{class_SO_term}) {
            # work out the term from the alleles
            
            my $alleles = $self->get_all_Alleles(); #get all Allele objects
            my %alleles; #to get all the different alleles in the Variation
            map {$alleles{$_->allele}++} @{$alleles};
            my $allele_string = join '/',keys %alleles;

            # if no allele objects, try to fetch from variation_feature
            unless($allele_string || @$alleles) {
              my ($vf) = @{$self->get_all_VariationFeatures};
              $allele_string = $vf->allele_string() if $vf;
            }

            $self->{class_SO_term} = SO_variation_class($allele_string);
        }

        # convert the SO term to the ensembl display term

        $self->{class_display_term} = $self->is_somatic ? 
            $VARIATION_CLASSES{$self->{class_SO_term}}->{somatic_display_term} : 
            $VARIATION_CLASSES{$self->{class_SO_term}}->{display_term};
    }
    
    return $self->{class_display_term};
}

=head2 derived_allele_frequency

  Arg[1]     : Bio::EnsEMBL::Variation::Population  $population 
  Example    : $daf = $variation->derived_allele_frequency($population);
  Description: Gets the derived allele frequency for the population. 
               The DAF is the frequency of the reference allele that is 
               different from the allele in Chimp. If none of the alleles
               is the same as the ancestral, will return reference allele
               frequency
  Returntype : Float
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
    #$self->_fill_in_maf_data unless defined($self->{minor_allele});
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
    $self->_fill_in_maf_data unless defined($self->{minor_allele_frequency});
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
    $self->_fill_in_maf_data unless defined($self->{minor_allele_count});
    return $self->{minor_allele_count};
}


##############################################
## HACK FIX FOR RELEASE 90 MISSING MAF DATA ##
##############################################
sub _fill_in_maf_data {
  my $self = shift;

  # only do this for vars with 1KG evidence code
  return unless grep {$_ eq '1000Genomes'} @{$self->get_all_evidence_values || []};

  my (%counts, $total);

  foreach my $allele(grep {$_->population->name =~ /^1000GENOMES:[A-Z]{3}$/} @{$self->get_all_Alleles}) {
    my $count = $allele->count;
    $counts{$allele->allele} += $count;
    $total += $count;
  }

  return unless $total;
  return unless scalar keys %counts > 1;

  # find the second most frequent allele
  my @sorted = sort {$counts{$b} <=> $counts{$a}} keys %counts;

  $self->{minor_allele} = $sorted[1];
  $self->{minor_allele_count} = $counts{$sorted[1]};
  $self->{minor_allele_frequency} = sprintf("%.4g", $counts{$sorted[1]} / $total);
}


=head2 get_all_clinical_significance_states

  Arg [1]    : none
  Example    : my @csstates = @{$v->get_all_clinical_significance_states()};
  Description: Retrieves all clinical_significance states for this variation, as reported by dbSNP.
               When available, this will contain one or more of the following strings:
                unknown 
                untested
                non-pathogenic
                probable-non-pathogenic
                probable-pathogenic
                pathogenic
                drug-response
                histocompatibility
                other 

  Returntype : Listref of strings
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub get_all_clinical_significance_states {
    my $self = shift;

    return $self->{clinical_significance} || [];
}

    
=head2 get_all_PhenotypeFeatures

  Args       : none
  Example    : my $pfs = $var->get_all_PhenotypeFeatures()
  Description: Getter for PhenotypeFeatures for this Variation, returns empty list if 
               there are none. 
  Returntype : Listref of PhenotypeFeatures
  Exceptions : none
  Caller     : general

=cut

sub get_all_PhenotypeFeatures {
    my $self = shift;

    #Assert the adaptor reference
    assert_ref($self->adaptor(),'Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');
    
    # Get the annotations from the database
    return $self->adaptor->db->get_PhenotypeFeatureAdaptor()->fetch_all_by_Variation($self);

}

=head2 display_consequence

  Arg [1]    : (optional) String $term_type
  Description: Get the term for the most severe consequence of this 
               VariationFeature. By default returns the Sequence Ontology term
               (e.g. 'missense_variant'). $term_type can also be 'label'
               (e.g. 'Missense variant'), 'Ensembl' ((e.g. 'NON_SYNONYMOUS_CODING') 
               or 'NCBI' (e.g. 'missense')
  Returntype : string
  Exceptions : none
  Status     : Stable

=cut
sub display_consequence {
    my $self = shift;
    my $term_type = shift;
        
    my @ocs = map {@{$_->get_all_OverlapConsequences}} @{$self->get_all_VariationFeatures};
    
    ## there is no consequence if there is no genomic location
    return undef if scalar @ocs == 0;

    my $highest;
    
    for my $cons (@ocs) {
        $highest ||= $cons;
        if ($cons->rank < $highest->rank) {
            $highest = $cons;
        }
    }
	
	  my $method_name;
    if(defined($term_type)) {
        $method_name = $term_type.($term_type eq 'label' ? '' : '_term');
        $method_name = 'SO_term' unless @ocs && $ocs[0]->can($method_name);
    }

    $method_name ||= 'SO_term';
	
    return $highest->$method_name;
}

    
=head2 get_all_Publications

  Args       : none
  Example    : my $pubs = $var->get_all_Publications()
  Description: Getter for Publications citing this Variation, returns empty list if 
               there are none. 
  Returntype : Listref of Publications
  Exceptions : none
  Caller     : general

=cut
sub get_all_Publications{
    my $self = shift;

    #Assert the adaptor reference
    assert_ref($self->adaptor(),'Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');
    
    # Get the citations from the database
    return $self->adaptor->db->get_PublicationAdaptor()->fetch_all_by_Variation($self);

}

=head2 get_all_attributes

Example    : %attribs = %{$obj->get_all_attributes}
Description: Retrieves attributes of this Variation as a hash reference
             containing key-value pairs e.g. "A_homoeologue_locus" => "EPlTAEV01066791"
             Here is a list of the keys used: A_homoeologue_locus, B_homoeologue_locus, D_homoeologue_locus
             This is used to cross-reference SNP loci in the context of polyploid genomes
Returntype : hashref
Exceptions : none
Caller     : general
Status     : Stable
                                                        
=cut

sub get_all_attributes {
    my $self = shift;

    if(!defined($self->{attribs})) {
        $self->{attribs} = $self->adaptor->_fetch_attribs_by_dbID($self->dbID);
    }

    return $self->{attribs};
}

=head2 update_attributes

Example    : %attribs = $obj->update_attributes(\%attribs);
Description: Updates the attributes of this Variation 
             May add new types to pre-existing data; does not over-write
Returntype : nono
Exceptions : none
Caller     : general
Status     : Experimental
                                                        
=cut
sub update_attributes{
    my $self       = shift;
    my $attributes = shift;

    foreach my $attrib_type (keys %$attributes){
        $self->{attribs}->{$attrib_type} = $attributes->{$attrib_type};
    }
}

1;
