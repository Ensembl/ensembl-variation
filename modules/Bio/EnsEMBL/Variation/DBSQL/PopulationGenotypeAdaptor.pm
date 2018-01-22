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

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::PopulationGenotypeAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::PopulationGenotypeAdaptor

=head1 SYNOPSIS
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous');

  my $pga = $registry->get_adaptor('human', 'variation', 'populationgenotype');
  my $va = $registry->get_adaptor('human', 'variation', 'variation');

  # Get a PopulationGenotype by its internal identifier
  my $pgtype = $pga->fetch_by_dbID(145);
  print join(' ', $pgtype->population()->name(), $pgtype->allele1(), $pgtype->allele2(), $pgtype->frequency()), "\n";

  # Get all PopulationGenotypes for a Variation
  my $variation = $va->fetch_by_name('rs1121');
  foreach $pgtype (@{ $pga->fetch_all_by_Variation($variation) }) {
    print join(' ', $pgtype->population()->name(), $pgtype->allele1(), $pgtype->allele2(), $pgtype->frequency()), "\n";
  }

=head1 DESCRIPTION

This adaptor provides database connectivity for PopulationGenotype objects.
PopulationGenotypes may be retrieved from the Ensembl variation database by
several means using this module.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::PopulationGenotypeAdaptor;

use Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);

use Bio::EnsEMBL::Variation::PopulationGenotype;

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor');

our $CACHE_SIZE = 5;




sub store {
  my ($self, $popgt) = @_;
  
  my $dbh = $self->dbc->db_handle;
  
  # get genotype code
  my $gt_code = $self->_genotype_code($popgt->genotype);
  
  my $sth = $dbh->prepare_cached(q{
	INSERT DELAYED INTO population_genotype (
	  variation_id,
	  subsnp_id,
	  genotype_code_id,
	  frequency,
	  population_id,
	  count			
	) VALUES (?,?,?,?,?,?)
  });
  
  $sth->execute(
	$popgt->{_variation_id} || $popgt->variation->dbID,
	$popgt->{subsnp},
	$gt_code,
	$popgt->frequency,
	$popgt->population ? $popgt->population->dbID : undef,
	$popgt->count
  );
  
  # reset cache
  delete $self->{_cache} if $self->{_cache};
  
  $sth->finish;
}



sub store_multiple {
  my ($self, $popgts) = @_;
  
  my $dbh = $self->dbc->db_handle;
  
  my $q_string = join ",", map {'(?,?,?,?,?,?)'} @$popgts;
  
  my @args = map {
	$_->{_variation_id} || $_->variation->dbID,
	$_->{subsnp},
	$self->_genotype_code($_->genotype),
	$_->frequency,
	$_->population ? $_->population->dbID : undef,
	$_->count
  } @$popgts;
  
  my $sth = $dbh->prepare_cached(qq{
	INSERT INTO population_genotype (
	  variation_id,
	  subsnp_id,
	  genotype_code_id,
	  frequency,
	  population_id,
	  count				
	) VALUES $q_string
  });
  
  $sth->execute(@args);
  
  # reset cache
  delete $self->{_cache} if $self->{_cache};
  
  $sth->finish;
}

sub store_to_file_handle {
  my ($self, $popgt, $file_handle) = @_;
  
  print $file_handle join("\t",
	$popgt->{_variation_id} || $popgt->variation->dbID || '\N',
	$popgt->{subsnp} || '\N',
	$self->_genotype_code($popgt->genotype),
	defined($popgt->frequency) ? $popgt->frequency :  '\N',
	$popgt->population ? $popgt->population->dbID : '\N',
	defined($popgt->count) ? $popgt->count : '\N',
  )."\n";
}

=head2 fetch_by_dbID

  Arg [1]    : int $dbID
  Example    : $pgtype = $pgtype_adaptor->fetch_by_dbID(15767);
  Description: Retrieves a population genotype via its unique internal
               identifier.  undef is returned if no such population genotype
               exists.
  Returntype : Bio::EnsEMBL::Variation::Variation::PopulationGenotype or undef
  Exceptions : throw if no dbID argument is provided
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;
  
  if (! $dbID){
    throw('no dbID argument provided');
  }
  return shift @{$self->generic_fetch("pg.population_genotype_id = " . $dbID)};

}




=head2 fetch_all_by_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population
  Example    : $pop = $pop_adaptor->fetch_by_dbID(1345);
               @gtys = $pgty_adaptor->fetch_all_by_Population($pop);
  Description: Retrieves all genotypes which are stored for a specified
               population.
  Returntype : Bio::EnsEMBL::Variation::PopulationGenotype
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Population {
  my $self = shift;
  my $pop = shift;

  if(!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population')) {
    throw('Bio::EnsEMBL::Variation::Population argument expected');
  }

  if(!defined($pop->dbID())) {
    warning("Cannot retrieve genotypes for population without set dbID");
    return [];
  }

  my $constraint = "pg.population_id = " . $pop->dbID();
  
  # Add the constraint for failed variations
  $constraint .= " AND " . $self->db->_exclude_failed_variations_constraint()
      unless $self->db->include_failed_variations() ;
  
  return $self->generic_fetch($constraint);
}



=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL::Variation $variation
  Arg [2]    : Bio::EnsEMBL::Variation::Population (optional)
  Example    : my $var = $variation_adaptor->fetch_by_name( "rs1121" )
               $poptypes = $poptype_adaptor->fetch_all_by_Variation( $var )
  Description: Retrieves a list of population genotypes for the given Variation.
               If none are available an empty listref is returned.
  Returntype : listref Bio::EnsEMBL::Variation::PopulationGenotype 
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut


sub fetch_all_by_Variation {
  my $self = shift;
  my $variation = shift;
  my $population = shift;
  
  # Make sure that we are passed a Variation object
  assert_ref($variation,'Bio::EnsEMBL::Variation::Variation');
  
  # If we got a population argument, make sure that it is a Population object
  assert_ref($population,'Bio::EnsEMBL::Variation::Population') if (defined($population));
  
  my $variation_id = $variation->dbID();
  
	my $cached;
  my $return = [];

  if(defined($self->{_cache})) {
    foreach my $stored(@{$self->{_cache}}) {
      my @keys = keys %{$stored};
      $cached = $stored->{$keys[0]} if $keys[0] eq $variation_id;
      last if defined($cached);
    }
  }

  if(!defined($cached)) {
    
    # Add a constraint on the variation_id column and pass to generic fetch
    my $constraint = qq{ pg.variation_id = $variation_id };
    
    # If required, add a constraint on the population id
    if (defined($population)) {
      my $population_id = $population->dbID();
      $constraint .= qq{ AND pg.population_id = $population_id };
    }
  
    $cached = $self->generic_fetch($constraint);
    
    # If a population was specified, attach the population to the object
    map {$_->population($population)} @{$cached} if (defined($population));

    # add freqs from genotypes for human (1KG data)
    push @$cached, @{$self->_fetch_all_by_Variation_from_Genotypes($variation, $population)};
		
    # don't store if population specified
    return $cached if defined($population);
    
    # add genotypes for this variant to the cache
    push @{$self->{_cache}}, {$variation_id => $cached};
	
    # shift off first element to keep cache within size limit
    shift @{$self->{_cache}} if scalar @{$self->{_cache}} > $CACHE_SIZE;
  }
  
  if(defined($population)) {
		@$return = grep {$_->dbID eq $population->dbID} @{$cached};
  }
  else {
    $return = $cached;
  }
	
  return $return;
}

sub _fetch_all_by_Variation_from_Genotypes {
  my $self = shift;
  return $self->_generic_fetch_all_by_Variation_from_Genotypes(@_, 'PopulationGenotype');
}

=head2 fetch_all

  Description: Retrieves a list of all population genotypes.
  Returntype : listref Bio::EnsEMBL::Variation::PopulationGenotype 
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut


sub fetch_all {
  my $self = shift;
  
  my $constraint;
  # Add the constraint for failed variations if required
  $constraint = $self->db->_exclude_failed_variations_constraint()  unless $self->db->include_failed_variations() ;
  
  return $self->generic_fetch($constraint);
}

sub _tables{
    my $self = shift;

    my @tables = (
        ['population_genotype','pg']
        );
    # join to variation table if fails are not to be included    
    push @tables, ['variation','v']   unless $self->db->include_failed_variations() ;

    return @tables;
}

#Add a left join to the variation table if fails are not to be included
sub _left_join { 
  my $self = shift;

  return undef if $self->db->include_failed_variations() ;
  return ([ 'variation', 'v.variation_id = pg.variation_id']) ; 

}

sub _columns{
  return qw(pg.population_genotype_id pg.variation_id pg.subsnp_id pg.population_id pg.genotype_code_id pg.frequency pg.count)
}

sub _write_columns {
  return qw(variation_id subsnp_id genotype_code_id frequency population_id count);
}

sub _objs_from_sth{
  
  my $self = shift;
  my $sth = shift;
  
  my ($dbID, $variation_id, $subsnp_id, $population_id, $gt_code, $freq, $count);
  
  $sth->bind_columns(\$dbID, \$variation_id, \$subsnp_id, \$population_id, \$gt_code, \$freq, \$count);
  
  my (%pop_hash, %gt_code_hash, @results);
  
  while($sth->fetch) {
  
	my $pgtype  = Bio::EnsEMBL::Variation::PopulationGenotype->new_fast({
	  _variation_id => $variation_id,
	  subsnp        => $subsnp_id,
	  adaptor       => $self,
	  frequency     => $freq,
	  count         => $count,
      dbID          => $dbID,
	});
  
	$pop_hash{$population_id} ||= [];
	push @{$pop_hash{$population_id}}, $pgtype;
	
	$gt_code_hash{$gt_code} ||= [];
	push @{$gt_code_hash{$gt_code}}, $pgtype;
	
	push @results, $pgtype;
  }
  
  # fetch populations
  my $pa = $self->db()->get_PopulationAdaptor();
  my $pops = $pa->fetch_all_by_dbID_list([keys %pop_hash]);
  
  foreach my $p (@$pops) {
	foreach my $pgty (@{$pop_hash{$p->dbID()}}) {
	  $pgty->{population} = $p;
	}
  }
  
  # get all genotypes from codes
  my $gtca = $self->db->get_GenotypeCodeAdaptor();
  my $gtcs = $gtca->fetch_all_by_dbID_list([keys %gt_code_hash]);
  
  foreach my $gtc(@$gtcs) {
	foreach my $pgty(@{$gt_code_hash{$gtc->dbID}}) {
	  $pgty->{genotype} = $gtc->genotype;
	}
  }
  
  return \@results;
}

1;
