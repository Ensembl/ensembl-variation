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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::SupportingStructuralVariationAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::SupportingStructuralVariationAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $ssva = $reg->get_adaptor("human","variation","supportingstructuralvariation");

  # fetch a supporting structural variation by its name
  $ssv = $ssva->fetch_by_name('nssv133'); 

  # fetch all supporting evidences for a structural variation
  $sva = $reg->get_adaptor("human","variation","structuralvariation");
  $sv = $sva->fetch_by_dbID(145);
  foreach $ssv (@{$ssva->fetch_all_by_StructuralVariation($sv)}){
    print $ssv->dbID, " - ", $ssv->name ,"\n"; 
  }
  
  # Modify the include_failed_variations flag in DBAdaptor to also return supporting evidences that have been flagged as failed
  $va->db->include_failed_variations(1);

=head1 DESCRIPTION

This adaptor provides database connectivity for SupportingStructuralVariation objects.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::SupportingStructuralVariationAdaptor;

use Bio::EnsEMBL::Variation::DBSQL::BaseStructuralVariationAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::SupportingStructuralVariation;
use DBI qw(:sql_types);

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseStructuralVariationAdaptor');


sub _default_where_clause {
  my $self = shift;
  return 'is_evidence=1';
}


sub _objs_from_sth {
  my ($self, $sth) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches and
  # a fair bit of gymnastics is used.
  #
  my @svs;
  
  my ($struct_variation_id, $variation_name, $validation_status, $source_id, $class_attrib_id,
      $study_id, $is_evidence, $is_somatic, $alias, $clinical_significance, $copy_number);

  $sth->bind_columns(\$struct_variation_id, \$variation_name, \$validation_status, \$source_id, \$class_attrib_id,
                     \$study_id, \$is_evidence, \$is_somatic, \$alias, \$clinical_significance, \$copy_number);
                     
  my $aa  = $self->db->get_AttributeAdaptor;
  
  while($sth->fetch()) {

    my @clin_sig;
    @clin_sig = split(/,/,$clinical_significance) if (defined($clinical_significance));

    push @svs, Bio::EnsEMBL::Variation::SupportingStructuralVariation->new(
       -dbID                  => $struct_variation_id,
       -VARIATION_NAME        => $variation_name,
       -VALIDATION_STATUS     => $validation_status,
       -ADAPTOR               => $self,
       -_SOURCE_ID            => $source_id,
       -CLASS_SO_TERM         => $aa->attrib_value_for_id($class_attrib_id),
       -_STUDY_ID             => $study_id,
       -IS_EVIDENCE           => $is_evidence || 0,
       -IS_SOMATIC            => $is_somatic || 0,
       -ALIAS                 => $alias,
       -CLINICAL_SIGNIFICANCE => \@clin_sig,
       -COPY_NUMBER           => $copy_number
    );
  }
  return \@svs;
}


=head2 fetch_all_by_StructuralVariation

  Arg [1]     : Bio::EnsEMBL::Variation::StructuralVariation $sv
  Example     : my $sv = $sv_adaptor->fetch_by_name('esv9549');
                foreach my $ssv (@{$ssv_adaptor->fetch_all_by_StructuralVariation($sv)}){
                  print $ssv->variation_name,"\n";
                }
  Description : Retrieves all supporting evidences from a specified structural variant
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::SupportingStructuralVariation objects
  Exceptions  : throw if incorrect argument is passed
                warning if provided structural variant does not have a dbID
  Caller      : general
  Status      : Stable

=cut

sub fetch_all_by_StructuralVariation {
  my $self = shift;
  my $sv = shift;

  if(!ref($sv) || !$sv->isa('Bio::EnsEMBL::Variation::StructuralVariation')) {
    throw("Bio::EnsEMBL::Variation::StructuralVariation arg expected");
  }
    
  if(!$sv->dbID()) {
    warning("StructuralVariation does not have dbID, cannot retrieve structural variants");
    return [];
  } 
  
  my $cols = qq{ sa.supporting_structural_variation_id, sv.variation_name, sv.validation_status, sv.source_id, sv.class_attrib_id, 
                 sv.study_id, sv.is_evidence, sv.somatic, sv.alias, sv.clinical_significance, sv.copy_number
               };
  
  my $tables;
  foreach my $t ($self->_tables()) {
    next if ($t->[0] eq 'failed_structural_variation' and !$self->db->include_failed_variations());
    $tables .= ',' if ($tables);
    $tables .= join(' ',@$t);
    # Adds a left join to the failed_structural_variation table
    if ($t->[0] eq 'structural_variation' and !$self->db->include_failed_variations()) {
      $tables .= qq{ LEFT JOIN failed_structural_variation fsv 
                     ON (fsv.structural_variation_id=sv.structural_variation_id)};
    }
  }
  
  # Special case for one human study where SV can be a supporting evidence of an other SV
  my $constraint = '';
  
  # Add the constraint for failed structural variant
  if (!$self->db->include_failed_variations()) {
    $constraint .= ' AND ' . $self->db->_exclude_failed_structural_variations_constraint();
  }
           
  my $sth = $self->prepare(qq{SELECT $cols
                                FROM $tables, structural_variation_association sa
                               WHERE sa.supporting_structural_variation_id=sv.structural_variation_id
                                 AND sa.structural_variation_id = ?$constraint});

  $sth->bind_param(1,$sv->dbID,SQL_INTEGER);
  $sth->execute();

  my $results = $self->_objs_from_sth($sth);

  $sth->finish();

  return $results;
}


=head2 fetch_all_SO_term_by_structural_variation_dbID

  Arg [1]     : int $sv_id
  Example     : my $sv = $sv_adaptor->fetch_by_name('esv9549');
                foreach my $SO_term (@{$ssv_adaptor->fetch_all_SO_term_by_structural_variation_dbID($sv->dbID)}){
                  print $SO_term,"\n";
                }
  Description : Retrieves all supporting evidences classes from a specified structural variant ID
  ReturnType  : reference to list of strings
  Exceptions  : throw if structural variation ID arg is not defined
  Caller      : general
  Status      : Stable

=cut
sub fetch_all_SO_term_by_structural_variation_dbID {
  my $self  = shift;
  my $sv_id = shift;
  
  if(!defined($sv_id)) {
    throw("Structural variation ID arg expected");
  }
  
  my @ssv_SO_list = ();
  
  # No failed SV/SSV in the structural_variation_feature table
  my $sth = $self->prepare(qq{SELECT distinct a.value
                                FROM structural_variation_feature svf, structural_variation_association sa, attrib a
                                WHERE sa.supporting_structural_variation_id=svf.structural_variation_id
                                      AND sa.structural_variation_id = ? 
                                      AND svf.class_attrib_id=a.attrib_id
                              });
 
  $sth->bind_param(1,$sv_id,SQL_INTEGER);
  $sth->execute();                
  
  while (my ($SO_term) = $sth->fetchrow_array) {
    push (@ssv_SO_list, $SO_term);
  }
  return \@ssv_SO_list;
}

1;
