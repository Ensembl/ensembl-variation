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


# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::StructuralVariationAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::StructuralVariationAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $sva = $reg->get_adaptor("human","variation","structuralvariation");
  $sta = $reg->get_adaptor("human","variation","study");
	
  # Get a StructuralVariation by its internal identifier
  $sv = $sva->fetch_by_dbID(145);

	# Get a StructuralVariation by its name
	$sv = $sva->fetch_by_name('esv1285');
	
	# Get all StructuralVariation by a study
	$study = $sta->fetch_by_name('estd1');
	foreach my $sv (@{$sva->fetch_all_by_Study($study)}){
		print $sv->variation_name,"\n";
	}
 
=head1 DESCRIPTION

This adaptor provides database connectivity for StructuralVariation objects.
Genomic locations of structural variations can be obtained from the database using this
adaptor.  See the base class BaseFeatureAdaptor for more information.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::StructuralVariationAdaptor;

use Bio::EnsEMBL::Variation::StructuralVariation;
use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use DBI qw(:sql_types);

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');

# method used by superclass to construct SQL
sub _tables { return (['structural_variation', 'sv'],['source', 's']); }


sub _default_where_clause {
  my $self = shift;
  return 'sv.source_id = s.source_id';
}


sub _columns {
  return qw( sv.structural_variation_id sv.variation_name sv.validation_status s.name s.version 
	           s.description sv.class_attrib_id sv.study_id);
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
		  $source_description, $class_attrib_id, $study_id);

  $sth->bind_columns(\$struct_variation_id, \$variation_name, \$validation_status, \$source_name, 
										 \$source_version, \$source_description, \$class_attrib_id, \$study_id);

	my $aa  = $self->db->get_AttributeAdaptor;
	my $sta = $self->db->get_StudyAdaptor();
	
  while($sth->fetch()) {
	
		my $study;
		$study = $sta->fetch_by_dbID($study_id) if (defined($study_id));
	
		# Get the validation status
    $validation_status ||= 0;
    my @states = split(/,/,$validation_status);
	
    push @svs, Bio::EnsEMBL::Variation::StructuralVariation->new(
       -dbID               => $struct_variation_id,
			 -VARIATION_NAME     => $variation_name,
			 -VALIDATION_STATES  => \@states,
       -ADAPTOR            => $self,
       -SOURCE             => $source_name,
       -SOURCE_VERSION     => $source_version,
	     -SOURCE_DESCRIPTION => $source_description,
	     -CLASS_SO_TERM      => $aa->attrib_value_for_id($class_attrib_id),
	     -STUDY              => $study,
  	);
  }
  return \@svs;


}


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$simple_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all simple features in 
               the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub list_dbIDs {
  my $self = shift;
  return $self->_list_dbIDs('structural_variation');
}


=head2 fetch_by_name

    Args[1]     : string $name
    Example     : my $structural_variation = $sv_adaptor->fetch_by_name('esv263');
    Description : returns the structural variation with the given variation name (or undef if one isn't found).
                  If the name argument is undef this will be converted to NULL in the SQL statement generated.
    ReturnType  : Bio::EnsEMBL::Variation::StructuralVariation
    Exceptions  : thrown if there are multiple objects found with the same variation name
    Caller      : general
    Status      : Stable

=cut

sub fetch_by_name {
    my ($self, $name) = @_;
    
    my $constraint = sprintf('variation_name = %s', $self->dbc->db_handle->quote( $name, SQL_VARCHAR ) );
    my $objs = $self->generic_fetch($constraint);
    throw("Multiple structural variations found with the same name: '$name'") if @$objs > 1;
    return $objs->[0] if @$objs == 1;
}


=head2 fetch_all_by_Study

  Arg [1]     : Bio::EnsEMBL::Variation::Study $study_id
  Example     : my $study = $study_adaptor->fetch_by_name('estd1');
                foreach my $sv (@{$sv_adaptor->fetch_all_by_Study($study)}){
		    		 			print $sv->variation_name,"\n";
                }
  Description : Retrieves all structural variations from a specified study
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::StructuralVariation objects
  Exceptions  : throw if incorrect argument is passed
                warning if provided study does not have a dbID
  Caller      : general
  Status      : At Risk

=cut

sub fetch_all_by_Study {
  my $self = shift;
  my $study = shift;

  if(!ref($study) || !$study->isa('Bio::EnsEMBL::Variation::Study')) {
    throw("Bio::EnsEMBL::Variation::Study arg expected");
  }
    
  if(!$study->dbID()) {
    warning("Study does not have dbID, cannot retrieve structural variants");
    return [];
  } 

  my $cols = join ",", $self->_columns();
	my $tables = '';
	foreach my $tab ($self->_tables()) {
		$tables .= ', ' if ($tables ne '');
		$tables .= $tab->[0].' '.$tab->[1];
	}
	my $clause = $self->_default_where_clause;
	
  my $sth = $self->prepare(qq{SELECT $cols FROM $tables WHERE $clause AND sv.study_id = ?});
  $sth->bind_param(1,$study->dbID,SQL_INTEGER);
  $sth->execute();

  my $results = $self->_objs_from_sth($sth);

  $sth->finish();

  return $results;
}

1;
