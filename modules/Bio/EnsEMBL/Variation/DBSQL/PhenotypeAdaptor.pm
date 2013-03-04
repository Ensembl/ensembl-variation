=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::PhenotypeAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::PhenotypeAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $pa = $reg->get_adaptor("human","variation","phenotype");

  # Get a list of all phenotypes.
  $phenotypes = $pa->fetch_all();

=head1 DESCRIPTION

This adaptor provides database connectivity for Phenotype objects.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::PhenotypeAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Variation::Phenotype;

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');

sub fetch_by_description {
    my $self = shift;
    my $desc = shift;
    return $self->generic_fetch("p.description = '$desc'");
}

sub _tables {
    return (['phenotype', 'p']);
}

sub _columns {
    return qw(p.phenotype_id p.name p.description);
}

sub _objs_from_sth {
    my ($self, $sth) = @_;
    
    my ($phenotype_id, $name, $description, @result);
    
    $sth->bind_columns(\$phenotype_id, \$name, \$description);
    
    push @result, Bio::EnsEMBL::Variation::Phenotype->new_fast({
        dbID        => $phenotype_id,
        name        => $name,
        description => $description,
        adaptor     => $self,
    }) while $sth->fetch;
    
    $sth->finish;
    
    return \@result;
}

1;