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
  ReturnType  : String
  Exceptions  : thrown on incorrect argument
  Caller      : general
  Status      : Stable

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
  
  # only fetch codes with dbID higher than what we have already
  $self->db->{_max_genotype_code} ||= 0;
	
	my %gt_codes = map {(join "|", (@{$_->genotype}, defined($_->phased) ? $_->phased : "NULL")) => $_->dbID} sort {$b->dbID <=> $a->dbID} @{$self->{_genotype_code_adaptor}->generic_fetch('gc.genotype_code_id > '.$self->db->{_max_genotype_code})};
  
  # store highest code
  $self->db->{_max_genotype_code} = (sort {$a <=> $b} values %gt_codes)[-1];
	
  # add new codes to hash
	$self->db->{_genotype_codes}->{$_} = $gt_codes{$_} for keys %gt_codes;
	
	return $self->db->{_genotype_codes};
}

# get or (if not yet in DB) add a new GT code
sub _genotype_code {
	my ($self, $genotype, $phased) = @_;
	
	# check if cache is loaded
	my $just_loaded = 0;
	
	if(!exists($self->db->{_genotype_codes})) {
		$self->_cache_genotype_codes;
		$just_loaded = 1;
	}
	
	# include phased status in $gt_string to make unique
	my $gt_string = join "|", (@$genotype, defined($phased) ? $phased : 'NULL');
	
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
					genotype_code_id, allele_code_id, haplotype_id, phased
				)
				VALUES (?,?,?,?)
			});
			
			$phased = undef if defined($phased) && $phased eq 'NULL';
			
			for my $hap_id(1..(scalar @$genotype)) {
				$sth->execute($gt_code, $allele_codes{$genotype->[$hap_id-1]}, $hap_id, $phased);
			}
			
			$sth->finish;
			
			$self->db->{_genotype_codes}->{$gt_string} = $gt_code;
		}
	}
	
	return $self->db->{_genotype_codes}->{$gt_string};
}


1;
