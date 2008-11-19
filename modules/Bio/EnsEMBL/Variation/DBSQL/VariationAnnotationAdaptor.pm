
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::VariationAnnotationAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::VariationAnnotationAdaptor

=head1 SYNOPSIS

  $vdb = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(...);
  $db  = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);

  $va = $vdb->get_VariationAdaptor();
  $vaa = $vdb->get_VariationAnnotationAdaptor();

  # Get a VariationAnotation by its internal identifier
  $va = $vaa->fetch_by_dbID(45);

  # fetch all annotation for a particular variation
  $v = $va->fetch_by_name('rs56');

  foreach $va (@{$vaa->fetch_all_by_Variation($v)}) {
    print $va->phenotype_name(), $va->phenotype_description(), $va->source_name(), $va->study_type(), $va->local_stable_id(),"\n";
  }

=head1 DESCRIPTION

This adaptor provides database connectivity between Variation and VariationAnnotation objects.

=head1 AUTHOR - Yuan Chen

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::VariationAnnotationAdaptor;

use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor');

=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL:Variation::Variation $var
  Example    : my @vas = @{$vaa->fetch_all_by_Variation($var)};
  Description: Retrieves all variation annotations for a given variation.
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationAnnotation
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

=cut


sub fetch_all_by_Variation {
  my $self = shift;
  my $var  = shift;

  if(!ref($var) || !$var->isa('Bio::EnsEMBL::Variation::Variation')) {
    throw('Bio::EnsEMBL::Variation::Variation arg expected');
  }

  if(!defined($var->dbID())) {
    throw("Variation arg must have defined dbID");
  }

  return $self->generic_fetch("va.variation_id = ".$var->dbID());
}

=head2 fetch_all_by_phenotype_source_name

  Arg [1]    : string $phenotype_name
  Arg [2]    : string $source_name (optional)
  Example    : $vaa = $va_adaptor->fetch_by_phenotype_source_name('BD','EGA');
  Description: Retrieves a variation annotation object via its phenotype/source name
  Returntype : list of ref of Bio::EnsEMBL::Variation::VariationAnnotation
  Exceptions : throw if phenotype name argument is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_pheonotype_source_name {

  my $self = shift;
  my $phenotype_name  = shift;
  my $source_name = shift;

  throw('phenotype_name argument expected') if(!defined($phenotype_name));

  my $extra_sql = " AND p.phenotype_name = ? ";
  if (defined $source_name ) {
    $extra_sql .= qq( AND s.name = ?);
  }
  
  return $self->generic_fetch("$extra_sql");
  
}

# method used by superclass to construct SQL
sub _tables { return (['variation_annotation', 'va'],
                       [ 'phenotype', 'p'],
                       [ 'source', 's']); 
}

sub _default_where_clause {
  my $self = shift;

  return 'va.phenotype_id = p.phenotype_id AND va.source_id = s.source_id';
}

sub _columns {
  return qw( va.variation_annotation_id va.variation_id 
             p.name p.description 
             s.name va.study_type
             va.local_stable_id
             );
}



sub _objs_from_sth {
  my ($self, $sth) = @_;

  my @features;

  my ($variation_annotation_id,$var_id,$phenotype_name,$phenotype_description,$source_name,$study_type,$local_stable_id);
  $sth->bind_columns(\$variation_annotation_id,\$var_id,\$phenotype_name,\$phenotype_description,\$source_name,\$study_type,\$local_stable_id);

  while($sth->fetch()) {
    push @features, $self->_create_feature_fast('Bio::EnsEMBL::Variation::VariationAnnotation',

    {'dbID' => $variation_annotation_id,
     '_variation_id'         => $var_id,
     'phenotype_name'       => $phenotype_name,
     'phenotype_description'=> $phenotype_description,
     'source_name'          => $source_name,
     'study_type'           => $study_type,
     'local_stable_id'      => $local_stable_id,
     'adaptor'  => $self,
    });
  }

  return \@features;


}


1;
