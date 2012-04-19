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

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::PopulationGenotypeAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
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

use Bio::EnsEMBL::Variation::PopulationGenotype;

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
	  sample_id,
	  frequency,
	  count			
	) VALUES (?,?,?,?,?,?)
  });
  
  $sth->execute(
	$popgt->{_variation_id} || $popgt->variation->dbID,
	$popgt->{subsnp},
	$gt_code,
	$popgt->population ? $popgt->population->dbID : undef,
	$popgt->frequency,
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
	$_->population ? $_->population->dbID : undef,
	$_->frequency,
	$_->count
  } @$popgts;
  
  my $sth = $dbh->prepare_cached(qq{
	INSERT INTO population_genotype (
	  variation_id,
	  subsnp_id,
	  genotype_code_id,
	  sample_id,
	  frequency,
	  count				
	) VALUES $q_string
  });
  
  $sth->execute(@args);
  
  $sth->finish;
}

sub store_to_file_handle {
	my ($self, $popgt, $file_handle) = @_;
	
	my $dbh = $self->dbc->db_handle;
	
	print $file_handle join("\t",
		$popgt->{_variation_id} || $popgt->variation->dbID || '\N',
		$popgt->{subsnp} || '\N',
		$self->_genotype_code($popgt->genotype),
		$popgt->population ? $popgt->population->dbID : '\N',
		defined($popgt->frequency) ? $popgt->frequency :  '\N',
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
  Status     : At Risk

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
  Status     : At Risk

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

  my $constraint = "pg.sample_id = " . $pop->dbID();
  
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
  Status     : At Risk

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
  
  return $self->generic_fetch("pg.variation_id = " . $variation->dbID());
}

=head2 fetch_all

  Description: Retrieves a list of all population genotypes.
  Returntype : listref Bio::EnsEMBL::Variation::PopulationGenotype 
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

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

#ÊAdd a left join to the failed_variation table
sub _left_join { return ([ 'failed_variation', 'fv.variation_id = pg.variation_id']); }

sub _columns{
  return qw(pg.population_genotype_id pg.variation_id pg.subsnp_id pg.sample_id pg.genotype_code_id pg.frequency pg.count)
}

sub _write_columns {
  return qw(variation_id subsnp_id genotype_code_id sample_id frequency count);
}

sub _objs_from_sth{
  
  my $self = shift;
  my $sth = shift;
  
  my ($dbID, $variation_id, $subsnp_id, $sample_id, $gt_code, $freq, $count);
  
  $sth->bind_columns(\$dbID, \$variation_id, \$subsnp_id, \$sample_id, \$gt_code, \$freq, \$count);
  
  my (%pop_hash, %gt_code_hash, @results);
  
  while($sth->fetch) {
  
	my $pgtype  = Bio::EnsEMBL::Variation::PopulationGenotype->new_fast({
	  _variation_id => $variation_id,
	  subsnp        => $subsnp_id,
	  adaptor       => $self,
	  frequency     => $freq,
	  count         => $count
	});
  
	$pop_hash{$sample_id} ||= [];
	push @{$pop_hash{$sample_id}}, $pgtype;
	
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
