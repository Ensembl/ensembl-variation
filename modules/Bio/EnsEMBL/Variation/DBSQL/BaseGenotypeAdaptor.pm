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

Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor

=head1 DESCRIPTION

Abstract adaptor class for fetching genotypes. Should not be invoked directly.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor;

use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);


our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');

=head2 get_subsnp_handle

  Arg[1]      : Bio::EnsEMBL::Variation::Allele
                The allele object to get the subsnp handle for
  Example     : my $handle = $adaptor->get_subsnp_handle($allele);
                print "The allele '" . $allele->allele() . "' of subsnp 'ss" . $allele->subsnp_id() . "' was submitted by '$handle'\n";
		
  Description : Gets the submitter handle for the specified genotype
  ReturnType  : string
  Exceptions  : thrown on incorrect argument
  Caller      : general
  Status      : At Risk

=cut

sub get_subsnp_handle {
    my $self = shift;
    my $gt = shift;
    
    # Assert that the object passed is a Genotype
    assert_ref($gt,'Bio::EnsEMBL::Variation::Genotype');
    
    # Get the subsnp id and get rid of any 'ss' prefix
    my $ssid = $gt->subsnp() || "";
    $ssid =~ s/^ss//;
    
    my $stmt = qq{
        SELECT
            handle
        FROM
            subsnp_handle
        WHERE
            subsnp_id = ?
        LIMIT 1
    };
    my $sth = $self->prepare($stmt);
    $sth->execute($ssid);
    my $handle;
    $sth->bind_columns(\$handle);
    $sth->fetch();
    
    return $handle;
}

# caches genotype codes in a hash on the adaptor object
sub _cache_genotype_codes {
	my $self = shift;
	
	if(!defined($self->{_genotype_code_adaptor})) {
		$self->{_genotype_code_adaptor} = $self->db->get_GenotypeCodeAdaptor;
	}
	
	my %gt_codes = map {(join "|", @{$_->genotype}) => $_->dbID} @{$self->{_genotype_code_adaptor}->fetch_all()};
	
	$self->db->{_genotype_codes} = \%gt_codes;
	
	return $self->db->{_genotype_codes};
}

# get or (if not yet in DB) add a new GT code
sub _genotype_code {
	my ($self, $genotype) = @_;
	
	# check if cache is loaded
	my $just_loaded = 0;
	
	if(!exists($self->db->{_genotype_codes})) {
		$self->_cache_genotype_codes;
		$just_loaded = 1;
	}
	
	my $gt_string = join "|", @$genotype;
	
	if(!exists($self->db->{_genotype_codes}->{$gt_string})) {
		
		# check another process hasn't created it by reloading the cache
		$self->_cache_genotype_codes unless $just_loaded;
		
		# if it still doesn't exist
		if(!exists($self->db->{_genotype_codes}->{$gt_string})) {
			
			# get allele codes
			if(!defined($self->{_allele_adaptor})) {
				$self->{_allele_adaptor} = $self->db->get_AlleleAdaptor;
			}
			
			my %allele_codes = map {$_ => $self->{_allele_adaptor}->_allele_code($_)} @$genotype;
			
			my $dbh = $self->dbc->db_handle;
			
			my $sth = $dbh->prepare(q{
				SELECT max(genotype_code_id) FROM genotype_code
			});
			$sth->execute();
			
			my $max_id;
			$sth->bind_columns(\$max_id);
			$sth->fetch;
			$sth->finish;
			$max_id ||= 0;
			
			my $gt_code = $max_id + 1;
			
			$sth = $dbh->prepare(q{
				INSERT INTO genotype_code (
					genotype_code_id, allele_code_id, haplotype_id
				)
				VALUES (?,?,?)
			});
			
			for my $hap_id(1..(scalar @$genotype)) {
				$sth->execute($gt_code, $allele_codes{$genotype->[$hap_id-1]}, $hap_id);
			}
			
			$sth->finish;
			
			$self->db->{_genotype_codes}->{$gt_string} = $gt_code;
		}
	}
	
	return $self->db->{_genotype_codes}->{$gt_string};
}


1;
