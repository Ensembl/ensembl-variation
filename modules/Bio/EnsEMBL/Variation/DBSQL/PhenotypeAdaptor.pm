=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw);

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');

=head2 fetch_by_description

  Arg [1]    : string $description
  Example    : $phenotype = $pheno_adaptor->fetch_all_by_description('diabetes');
  Description: Retrieves a list of Phenotype objects for a phenotype description
               If no phenotype exists undef is returned.
  Returntype : list ref of Bio::EnsEMBL::Variation::Phenotypes
  Exceptions : none
  Caller     : general
=cut
sub fetch_by_description {
    my $self = shift;
    my $desc = shift;
    warn "You need to pass a description as argument for fetch_by_description\n" unless $desc;
    return undef if(!$desc);
    $desc =~ s/'/\\'/g;
    return $self->generic_fetch("p.description = '$desc'");
}

=head2 fetch_by_description_accession_type

  Arg [1]    : string $description
  Arg [2]    : string $mapping_type - default 'is', option 'involves'
  Example    : $phenotype = $pheno_adaptor->fetch_by_description_accession_type('diabetes');
               $phenotype = $pheno_adaptor->fetch_by_description_accession_type('diabetes', 'involves');
  Description: Retrieves a list of Phenotype objects for a phenotype description including
               the corresponding ontology accessions of $mapping_type
               If no phenotype exists undef is returned.
  Returntype : list ref of Bio::EnsEMBL::Variation::Phenotypes
  Exceptions : throw if $description arg is not defined or $mapping_type not supported
  Caller     : general
=cut
sub fetch_by_description_accession_type {
    my $self = shift;
    my $desc = shift;
    my $mapping_type = shift;
    throw("You need to pass a description as argument for fetch_by_description") unless $desc;
    $mapping_type ||= 'is';
    throw("$mapping_type is not a valid mapping type, valid types are: 'is','involves'") unless $mapping_type eq 'is' || $mapping_type eq 'involves';

    $desc =~ s/'/\\'/g;

    my $constraint = qq{ p.description = '$desc'};
    $constraint .= qq{ AND poa.mapping_type = '$mapping_type'} if defined $mapping_type;

    return $self->generic_fetch($constraint);
}

=head2 fetch_all_by_ontology_accession

  Arg [1]    : string ontology accession
  Arg [2]    : string mapping type (is/involves)  optional
  Example    : $phenotype = $pheno_adaptor->fetch_all_by_ontology_accession('EFO:0000712', 'is');
  Description: Retrieves a list of Phenotype objects for an ontology accession
               If no phenotype exists undef is returned.
  Returntype : list ref of Bio::EnsEMBL::Variation::Phenotypes
  Exceptions : throw if dbID arg is not defined
  Caller     : general
  Status     : experimental
=cut
sub fetch_all_by_ontology_accession {
  my $self = shift;
  my $accession = shift;
  my $mapping_type = shift;

  my $constraint = "poa.accession = '$accession'"; 
  $constraint   .= " and mapping_type = '$mapping_type' " if defined $mapping_type;

  return $self->generic_fetch($constraint);
  #return $self->generic_fetch("poa.accession = '$accession'");
}

=head2 fetch_by_OntologyTerm

  Arg [1]    : Bio::EnsEMBL::OntologyTerm
  Arg [2]    : string mapping type (is/involves)  optional
  Example    : $phenotype = $pheno_adaptor->fetch_by_OntologyTerm( $ontologyterm, 'involves');
  Description: Retrieves a Phenotype object via an OntologyTerm
               If no phenotype exists undef is returned.
  Returntype : list ref of Bio::EnsEMBL::Variation::Phenotypes
  Exceptions : throw if dbID arg is not defined
  Caller     : general
  Status     : experimental
=cut
sub fetch_by_OntologyTerm {

  my $self = shift;
  my $term = shift;
  my $mapping_type = shift;

  # Check an OntologyTerm is supplied
  assert_ref($term,'Bio::EnsEMBL::OntologyTerm');

  return $self->fetch_all_by_ontology_accession($term->accession(), $mapping_type);
}



sub _left_join {
  my $self = shift;

  my @lj = ();
  
  push @lj, (
    [ 'phenotype_ontology_accession', 'p.phenotype_id = poa.phenotype_id' ],
    [ 'attrib', 'p.class_attrib_id = a.attrib_id']
  ) ;
  
  return @lj;
}

sub _tables {
    return (['phenotype', 'p'],
            ['phenotype_ontology_accession', 'poa'],
            ['attrib', 'a'] );
}

sub _columns {
    return qw(p.phenotype_id p.stable_id p.name p.description p.class_attrib_id a.value poa.accession poa.mapped_by_attrib poa.mapping_type);
}

sub _objs_from_sth {
    my ($self, $sth) = @_;

    my %row;

    $sth->bind_columns( \( @row{ @{$sth->{NAME_lc} } } ));
   
    ## deal with multiple rows due to multiple phenotype ontology terms 
    while ($sth->fetch) {
        $self->_obj_from_row(\%row);
    }

    # Get the created objects from the temporary hash
    my @objs = values %{ $self->{_temp_objs} };
    delete $self->{_temp_objs};
    
    return \@objs;
}

sub _obj_from_row {

    my ($self, $row) = @_;
    
    # If the object for this phenotype_id hasn't already been created, do that
    my $obj = $self->{_temp_objs}{$row->{phenotype_id}}; 

    unless (defined($obj)) {

      $obj = Bio::EnsEMBL::Variation::Phenotype->new_fast({
              dbID           => $row->{phenotype_id},
              stable_id      => $row->{stable_id},
              name           => $row->{name},
              description    => $row->{description},
              class_attrib_id => $row->{class_attrib_id},
              adaptor        => $self,
            }); 

      $self->{_temp_objs}{$row->{phenotype_id}} = $obj;

    }
    # Add class attrib value if available
    $obj->class_attrib( $row->{value}) if defined $row->{value} ;

    # Add a ontology accession if available
    my $link_source = $self->db->get_AttributeAdaptor->attrib_value_for_id($row->{mapped_by_attrib})
      if $row->{mapped_by_attrib};

    $obj->add_ontology_accession({ accession      => $row->{accession}, 
                                   mapping_source => $link_source,
                                   mapping_type   => $row->{mapping_type} }) if defined $row->{accession} ;
}

=head2 get_all_phenotype_class_types

  Example    : $phenotype_classes = $pheno_adaptor->get_all_phenotype_class_types();
  Description: Retrieves the list of existing phenotype class attribs
               If no phenotype class type exists undef is returned.
  Returntype : list ref of string
  Exceptions : none
  Caller     : general
=cut

sub get_all_phenotype_class_types {
  my $self = shift;

  my $phenos= $self->generic_fetch();
  my %pheno_class_attribs=();
  foreach my $pheno (@$phenos){
    $pheno_class_attribs{$pheno->class_attrib} = 1;
  }

  return values %pheno_class_attribs;
}

sub store{
   my ($self, $pheno) = @_;

   # ensure phenotype class_attrib_id is present if class_attrib is
   if (defined $pheno->{class_attrib} && ! defined $pheno->{class_attrib_id}){
    my $class_attrib_id = $self->db->get_AttributeAdaptor->attrib_id_for_type_value('phenotype_type', $pheno->{class_attrib});
    $pheno->{class_attrib_id} = $class_attrib_id;
   }

    my $dbh = $self->dbc->db_handle;

    my $sth = $dbh->prepare(qq{
        INSERT INTO phenotype (
             name,
             stable_id,
             description,
             class_attrib_id
        ) VALUES (?,?,?,?)
    });

    $sth->execute(        
        $pheno->{name},
        $pheno->{stable_id},
        $pheno->{description},
        $pheno->{class_attrib_id}
    );

    $sth->finish;

    # get dbID
    my $dbID = $dbh->last_insert_id(undef, undef, 'phenotype', 'phenotype_id');
    $pheno->{dbID}    = $dbID;
    $pheno->{adaptor} = $self;

    ## add ontology terms if available
    $self->store_ontology_accessions($pheno) if $pheno->{_ontology_accessions};

    ## maintain previous behaviour
    return $self;
}



sub store_ontology_accessions{

   my ($self, $pheno) = @_;

   my $dbh = $self->dbc->db_handle;

   my $sth = $dbh->prepare(qq{
        INSERT IGNORE INTO phenotype_ontology_accession (
             phenotype_id,
             accession,
             mapped_by_attrib,
             mapping_type
        ) VALUES (?,?,?,?)
    });

  foreach my $link_info (@{$pheno->{_ontology_accessions}} ){

    ## get attrib id for source of link - can this be mandatory?
    my $attrib_id;
    if($link_info->{mapping_source}){
      $attrib_id = $self->db->get_AttributeAdaptor->attrib_id_for_type_value( 'ontology_mapping', $link_info->{mapping_source}); 
      warn "Source type " . $link_info->{mapping_source} . " not supported for linking ontology descriptions to accessions\n" unless $attrib_id;
    }

    $sth->execute(
        $pheno->{dbID},
        $link_info->{accession},
        $attrib_id,
        $link_info->{mapping_type}
    );
   }
    $sth->finish;
}

1;
