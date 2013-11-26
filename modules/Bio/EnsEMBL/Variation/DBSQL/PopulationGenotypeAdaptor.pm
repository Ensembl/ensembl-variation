=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::PopulationGenotypeAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::PopulationGenotypeAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $pga = $reg->get_adaptor("human","variation","populationgenotype");  
  $pa = $reg->get_adaptor("human","variation","population");

  # Get a PopulationGenotype by its internal identifier
  $pgtype = $ia->fetch_by_dbID(145);

  print $pgtype->population->name(), " ",
        $pgtype->allele1(), ' ', $pgtype->allele2(), ' ', $pgtype->frequency();

  # Get all population genotypes for an population
  $pop = $pa->fetch_by_dbID(1219);

  foreach $pgtype (@{$pga->fetch_all_by_Population($pop)}) {
    print $pgtype->variation()->name(),  ' ',
          $pgtype->frequency();
          $pgtype->allele1(), '/', $pgtype->allele2(), "\n";
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

use Scalar::Util qw(weaken);

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor');




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
  $constraint .= " AND " . $self->db->_exclude_failed_variations_constraint();
  
  return $self->generic_fetch($constraint);
}



=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL::Variation $variation
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
  
  if(!ref($variation) || !$variation->isa('Bio::EnsEMBL::Variation::Variation')) {
	throw('Bio::EnsEMBL::Variation::Variation argument expected');
  }
  
  if(!defined($variation->dbID())) {
	warning("Cannot retrieve genotypes for variation without set dbID");
	return [];
  }
  
  my $pgs = $self->generic_fetch("pg.variation_id = " . $variation->dbID());
  
  # fetch pop GTs from ind GTs for human (1KG data)
  push @$pgs, @{$self->_fetch_all_by_Variation_from_Genotypes($variation)};
  
  return $pgs;
}



sub _fetch_all_by_Variation_from_Genotypes {
  my $self = shift;
  my $variation = shift;
  my $population = shift;
  
  # Make sure that we are passed a Variation object
  assert_ref($variation,'Bio::EnsEMBL::Variation::Variation');
  
  # If we got a population argument, make sure that it is a Population object
  assert_ref($population,'Bio::EnsEMBL::Variation::Population') if (defined($population));
  
  # fetch all genotypes
  my $genotypes = $variation->get_all_IndividualGenotypes();
  
  return [] unless scalar @$genotypes;
  
  # get populations for individuals
  my (@pop_list, %pop_hash);
  
  if(defined($population)) {
	@pop_list = ($population);
	map {$pop_hash{$population->dbID}{$_->dbID} = 1} @{$population->get_all_Individuals};
  }
  else {
	my $pa = $self->db->get_PopulationAdaptor();
	%pop_hash = %{$pa->_get_individual_population_hash([map {$_->individual->dbID} @$genotypes])};
	return [] unless %pop_hash;
	
	@pop_list = @{$pa->fetch_all_by_dbID_list([keys %pop_hash])};
  }
  
  return [] unless @pop_list and %pop_hash;
  
  my %ss_list = map {$_->subsnp || '' => 1} @$genotypes;
  
  my @objs;
  
  foreach my $pop(@pop_list) {
	
	next unless $pop->_freqs_from_gts;
	
	foreach my $ss(keys %ss_list) {
	  my (%counts, $total, @freqs);
	  map {$counts{$_->genotype_string(1)}++}
		grep {$pop_hash{$pop->dbID}{$_->individual->dbID}}
		grep {$_->subsnp || '' eq $ss}
		@$genotypes;
	  
	  next unless %counts;
	  
	  my @alleles = keys %counts;
	  $total += $_ for values %counts;
	  next unless $total;
	  
	  @freqs = map {defined($counts{$_}) ? ($counts{$_} / $total) : 0} @alleles;
	  
	  for my $i(0..$#alleles) {
		push @objs, Bio::EnsEMBL::Variation::PopulationGenotype->new_fast({
		  genotype   => [split /\|/, $alleles[$i]],
		  count      => scalar keys %counts ? ($counts{$alleles[$i]} || 0) : undef,
		  frequency  => @freqs ? $freqs[$i] : undef,
		  population => $pop,
		  variation  => $variation,
		  adaptor    => $self,
		  subsnp     => $ss eq '' ? undef : $ss,
		});

                weaken($objs[-1]->{'variation'});
	  }
	}
  }
  
  return \@objs;
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
  
  # Add the constraint for failed variations
  my $constraint = $self->db->_exclude_failed_variations_constraint();
  
  return $self->generic_fetch($constraint);
}

sub _tables{return (
  ['population_genotype','pg'],
  ['failed_variation','fv']
)}

#Add a left join to the failed_variation table
sub _left_join { return ([ 'failed_variation', 'fv.variation_id = pg.variation_id']); }

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
