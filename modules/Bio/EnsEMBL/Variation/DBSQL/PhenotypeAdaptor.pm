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
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::PhenotypeAdaptor
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

sub store{
   my ($self, $pheno) = @_;

    my $dbh = $self->dbc->db_handle;

    my $sth = $dbh->prepare(qq{
        INSERT INTO phenotype (
             name,
             description
        ) VALUES (?,?)
    });

    $sth->execute(        
        $pheno->{name},
        $pheno->{description}        
    );

    $sth->finish;

    # get dbID
    my $dbID = $dbh->last_insert_id(undef, undef, 'phenotype', 'phenotype_id');
    $pheno->{dbID}    = $dbID;
    $pheno->{adaptor} = $self;

}



1;
