=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
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
	['genotype_code','gc'],
	['allele_code','ac'],
	['failed_variation','fv']
)}

#ÊAdd a left join to the failed_variation table
sub _left_join { return ([ 'failed_variation', 'fv.variation_id = pg.variation_id']); }

sub _columns{
    return qw(pg.population_genotype_id pg.variation_id pg.subsnp_id pg.sample_id ac.allele gc.haplotype_id pg.frequency pg.count)
}

sub _default_where_clause  {
	return 'pg.genotype_code_id = gc.genotype_code_id AND gc.allele_code_id = ac.allele_code_id';
}

sub _objs_from_sth{
  
	my $self = shift;
	my $sth = shift;
	
	my ($dbID, $variation_id, $subsnp_id, $sample_id, $allele, $haplotype_id, $freq, $count);
	
	$sth->bind_columns(\$dbID, \$variation_id, \$subsnp_id, \$sample_id, \$allele, \$haplotype_id, \$freq, \$count);
	
	my (@results, %genotypes, %variation_ids, %samples, %freqs, %counts);
	
	while($sth->fetch) {
		# some rows won't have ssIDs, make one up
		$subsnp_id ||= "s".$variation_id;
		
		$genotypes{$subsnp_id}{$sample_id}{$haplotype_id} = $allele;
		$freqs{$subsnp_id}{$sample_id}                    = $freq;
		$counts{$subsnp_id}{$sample_id}                   = $count;
		$variation_ids{$subsnp_id}                        = $variation_id;
		$samples{$sample_id}                              = undef;
	}
	
	# fetch populations
	my $pa = $self->db->get_PopulationAdaptor;
	$samples{$_->dbID} = $_ for @{$pa->fetch_all_by_dbID_list([keys %samples])};
	
	my $ploidy = $self->ploidy;
	
	# create genotype objects from hash
	foreach $subsnp_id(keys %genotypes) {
		foreach $sample_id(keys %{$genotypes{$subsnp_id}}) {
			my @gt = map {$genotypes{$subsnp_id}{$sample_id}{$_}} sort {$a <=> $b} keys %{$genotypes{$subsnp_id}{$sample_id}};
			
			# splice it down to ploidy size
			@gt = splice @gt, 0, $ploidy;
			
			push @results, Bio::EnsEMBL::Variation::PopulationGenotype->new_fast({
				_variation_id => $variation_ids{$subsnp_id},
				subsnp        => ($subsnp_id =~ /^s/ ? undef : $subsnp_id),
				population    => $samples{$sample_id},
				genotype      => \@gt,
				adaptor       => $self,
				frequency     => $freqs{$subsnp_id}{$sample_id},
				count         => $counts{$subsnp_id}{$sample_id},
			});
		}
	}
	
	return \@results;
}
1;
