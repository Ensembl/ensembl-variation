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

1;
