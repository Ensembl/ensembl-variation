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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::SupportingStructuralVariationAdaptor
#
# Copyright (c) 2011 Ensembl
#
# You may distribute this module under the same terms as perl itself
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
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::Variation::SupportingStructuralVariation;
use DBI qw(:sql_types);

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseStructuralVariationAdaptor');


sub _default_where_clause {
  my $self = shift;
  return $self->SUPER::_default_where_clause().' AND is_evidence=1';
}


sub _objs_from_sth {
  my ($self, $sth) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches and
  # a fair bit of gymnastics is used.
  #
	my @svs;
	
  my ($struct_variation_id, $variation_name, $validation_status, $source_name, $source_version, 
		  $source_description, $class_attrib_id, $study_id, $is_evidence, $is_somatic);

  $sth->bind_columns(\$struct_variation_id, \$variation_name, \$validation_status, \$source_name, 
										 \$source_version, \$source_description, \$class_attrib_id, \$study_id, \$is_evidence, 
										 \$is_somatic);
										 
	my $aa  = $self->db->get_AttributeAdaptor;
	my $sta = $self->db->get_StudyAdaptor();
	
  while($sth->fetch()) {
	
		my $study;
		$study = $sta->fetch_by_dbID($study_id) if (defined($study_id));
	
		# Get the validation status
    $validation_status ||= 0;
    my @states = split(/,/,$validation_status);
	
    push @svs, Bio::EnsEMBL::Variation::SupportingStructuralVariation->new(
       -dbID               => $struct_variation_id,
			 -VARIATION_NAME     => $variation_name,
			 -VALIDATION_STATES  => \@states,
       -ADAPTOR            => $self,
       -SOURCE             => $source_name,
       -SOURCE_VERSION     => $source_version,
	     -SOURCE_DESCRIPTION => $source_description,
	     -CLASS_SO_TERM      => $aa->attrib_value_for_id($class_attrib_id),
	     -STUDY              => $study,
			 -IS_EVIDENCE        => $is_evidence || 0,
			 -IS_SOMATIC         => $is_somatic || 0
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
  Status      : At Risk

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
	
	my $cols = qq{ sa.supporting_structural_variation_id, sv.variation_name, sv.validation_status, s.name, 
	               s.version, s.description, sv.class_attrib_id, sv.study_id, sv.is_evidence, sv.somatic };
	
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
	my $constraint = $self->SUPER::_default_where_clause();
	
	# Add the constraint for failed structural variant
	$constraint .= " AND " . $self->db->_exclude_failed_structural_variations_constraint();
						
	my $sth = $self->prepare(qq{SELECT $cols
                                FROM $tables, structural_variation_association sa
                                WHERE $constraint 
                                      AND sa.supporting_structural_variation_id=sv.structural_variation_id
                                      AND sa.structural_variation_id = ?});
	$sth->bind_param(1,$sv->dbID,SQL_INTEGER);
	$sth->execute();

	my $results = $self->_objs_from_sth($sth);

	$sth->finish();

	return $results;
}


1;
