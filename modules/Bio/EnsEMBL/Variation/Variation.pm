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
    $v->five_prime_flanking($seq);
    $v->three_prime_flanking($seq);


    ...

    # print out the default name and source of the variation
    print $v->source(), ':',$v->name(), "\n";

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
    print "5' flanking: ", $v->five_prime_flanking(), "\n";
    print "3' flanking: ", $v->five_prime_flanking(), "\n";


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
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);


# List of validation states. Order must match that of set in database
our @VSTATES = ['cluster','freq','submitter','doublehit','hapmap'];

# Conversion of validation state to bit value
our %VSTATE2BIT = ('cluster'   => 1,   # 00000001
                   'freq'      => 2,   # 00000010
                   'submitter' => 4,   # 00000100
                   'doublehit' => 8,   # 00001000
                   'hapmap'    => 16); # 00010000


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

  Arg [-SYNONYMS] :
    reference to hash with list reference values -  keys are source
    names and values are lists of identifiers from that db.
    e.g.: {'dbSNP' => ['ss1231', '1231'], 'TSC' => ['1452']}

  Arg [-ALLELES] :
    reference to list of Bio::EnsEMBL::Variation::Allele objects

  Arg [-VALIDATION_STATES] :
    reference to list of strings

  Arg [-FIVE_PRIME_FLANKING_SEQ] :
    string - the five prime flanking nucleotide sequence

  Arg [-THREE_PRIME_FLANKING_SEQ] :
    string - the three prime flanking nucleotide sequence

  Example    : $v = Bio::EnsEMBL::Variation::Variation->new
                    (-name   => 'rs123',
                     -source => 'dbSNP');

  Description: Constructor. Instantiates a new Variation object.
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : none
  Caller     : general

=cut


sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($dbID, $adaptor, $name, $src, $syns,
      $alleles, $valid_states, $five_seq, $three_seq) =
        rearrange([qw(dbID ADAPTOR NAME SOURCE SYNONYMS ALLELES
                      VALIDATION_STATES FIVE_PRIME_FLANKING_SEQ
                      THREE_PRIME_FLANKING_SEQ)],@_);


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
                'source' => $src,
                'synonyms' => $syns || {},
                'alleles' => $alleles || [],
                'validation_code' => $vcode,
                'five_prime_flanking_seq' => $five_seq,
                'three_prime_flanking_seq' => $three_seq}, $class;
}



=head2 name

  Arg [1]    : string $newval (optional)
               The new value to set the name attribute to
  Example    : $name = $obj->name()
  Description: Getter/Setter for the name attribute
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub name{
  my $self = shift;
  return $self->{'name'} = shift if(@_);
  return $self->{'name'};
}



=head2 source

  Arg [1]    : string $newval (optional) 
               The new value to set the source attribute to
  Example    : $source = $obj->source()
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

=cut

sub get_all_synonyms {
  my $self = shift;
  my $source = shift;

  if($source) {
    return $self->{'synonyms'}->{$source} || []
  }

  my @synonyms = values %{$self->{'synonyms'}};

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

=cut

sub get_all_validation_states {
  my $self = shift;


  my $code = $self->{'validation_code'};

  # convert the bit field into an ordered array
  my @states;
  for(my $i = 0; $i < @VSTATES; $i++) {
    push @states, $VSTATES[$i] if((1 << $i) | $code);
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

  Arg [1]    : string $newval (optional) 
               The new value to set the source attribute to
  Example    : $source = $obj->source()
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




=head2 get_all_Alleles

  Arg [1]    : none
  Example    : @alleles = @{v->get_all_Alleles()};
  Description: Retrieves all Alleles associated with this variation
  Returntype : reference to list of Bio::EnsEMBL::Variation::Allele objects
  Exceptions : none
  Caller     : general

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

=cut

sub add_Allele {
  my $self = shift;
  my $allele = shift;

  if(!ref($allele) || !$allele->isa('Bio::EnsEMBL::Variation::Allele')) {
    throw("Bio::EnsEMBL::Variation::Allele argument expected");
  }

  push @{$self->{'alleles'}}, $allele;
}




=head2 five_prime_flanking_seq

  Arg [1]    : string $newval (optional) 
               The new value to set the five_prime_flanking_seq attribute to
  Example    : $five_prime_flanking_seq = $obj->five_prime_flanking_seq()
  Description: Getter/Setter for the five_prime_flanking_seq attribute
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub five_prime_flanking_seq{
  my $self = shift;

  ### TODO: this should maybe be lazy-loaded from database

  return $self->{'five_prime_flanking_seq'} = shift if(@_);
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

=cut

sub three_prime_flanking_seq{
  my $self = shift;

  ### TODO: this should maybe be lazy-loaded from database

  return $self->{'three_prime_flanking_seq'} = shift if(@_);
  return $self->{'three_prime_flanking_seq'};
}

1;
