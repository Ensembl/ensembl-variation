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
  

=head1 DESCRIPTION

This adaptor provides database connectivity for SupportingStructuralVariation objects.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::SupportingStructuralVariationAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::Variation::SupportingStructuralVariation;

use base qw{Bio::EnsEMBL::DBSQL::BaseAdaptor};


=head2 fetch_by_name

  Arg [1]    : string $name
  Example    : $ssv = $ssv_adaptor->fetch_by_name('nssv133');
  Description: Retrieves a supporting evidence object via its name
  Returntype : Bio::EnsEMBL::Variation::SupportingStructuralVariation
  Exceptions : throw if name argument is not defined
  Caller     : general
  Status     : DEPRECATED: use the fetch_all_by_name method

=cut

sub fetch_by_name {
  my $self = shift;
	deprecate('The supporting structural variation name is no more unique: please use the fetch_all_by_name method instead');
}


=head2 fetch_all_by_name

  Arg [1]    : string $name
  Example    : $ssv = $ssv_adaptor->fetch_all_by_name('nssv133');
  Description: Retrieves a list of supporting evidence objects via its name
  Returntype : listref of Bio::EnsEMBL::Variation::SupportingStructuralVariation
  Exceptions : throw if name argument is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_name {
  my $self = shift;
  my $name = shift;

  throw('name argument expected') if(!defined($name));

	my $cols = join ",", $self->_columns();

  my $sth = $self->prepare(qq{SELECT $cols
                              FROM   supporting_structural_variation 
                              WHERE  name = ?});

  $sth->bind_param(1,$name,SQL_VARCHAR);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);

  $sth->finish();

  return undef if(!@$result);

  return $result;
}

	
=head2 fetch_all_by_dbID_list

  Arg [1]    : listref $list
  Example    : $ssv = $ssv_adaptor->fetch_all_by_dbID_list([907,1132]);
  Description: Retrieves a listref of supporting evidence objects via a list of internal
               dbID identifiers
  Returntype : listref of Bio::EnsEMBL::Variation::SupportingStructuralVariation objects
  Exceptions : throw if list argument is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_dbID_list {
  my $self = shift;
  my $list = shift;

  if(!defined($list) || ref($list) ne 'ARRAY') {
    throw("list reference argument is required");
  }
  
  my $id_str = (@$list > 1)  ? " IN (".join(',',@$list).")"   :   ' = \''.$list->[0].'\'';
	
	my $cols = join ",", $self->_columns();

  my $sth = $self->prepare(qq{SELECT $cols
                              FROM   supporting_structural_variation 
                              WHERE  supporting_structural_variation_id $id_str});
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);

  $sth->finish();

  return undef if(!@$result);

  return $result;
}


=head2 fetch_all_by_StructuralVariation

  Arg [1]     : Bio::EnsEMBL::Variation::StructuralVariation $sv
  Example     : my $sv = $sv_adaptor->fetch_by_name('esv9549');
                foreach my $ssv (@{$ssv_adaptor->fetch_all_by_StructuralVariation($sv)}){
		    		 print $ssv->name,"\n";
                }
  Description : Retrieves all supporting evidences from a specified structural variant
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::SupportingStructuralVariation objects
  Exceptions  : throw if incorrect argument is passed
                warning if provided structural variant does not have a dbID
  Caller      : general
  Status      : At Risk

=cut

sub fetch_all_by_StructuralVariation{
	my $self = shift;
	my $sv = shift;

	if(!ref($sv) || !$sv->isa('Bio::EnsEMBL::Variation::StructuralVariation')) {
		throw("Bio::EnsEMBL::Variation::StructuralVariation arg expected");
	}
    
	if(!$sv->dbID()) {
		warning("StructuralVariation does not have dbID, cannot retrieve structural variants");
		return [];
  } 
	
	my $cols = join ",", $self->_columns();

	my $sth = $self->prepare(qq{SELECT $cols
				                        FROM supporting_structural_variation
				                        WHERE structural_variation_id = ?});
	$sth->bind_param(1,$sv->dbID,SQL_INTEGER);
	$sth->execute();

	my $results = $self->_objs_from_sth($sth);

	$sth->finish();

	return $results;
}

sub _columns {
  return qw( supporting_structural_variation_id name structural_variation_id class_attrib_id );
}

#
# private method, creates supporting evidence objects from an executed statement handle
# ordering of columns must be consistant
#
sub _objs_from_sth {
  my $self = shift;
  my $sth  = shift;

  my @ssvs;

  my ($ssv_id, $name, $structural_variation_id, $class_attrib_id);

	$sth->bind_columns(\$ssv_id, \$name, \$structural_variation_id, \$class_attrib_id);
	
	my $aa = $self->db->get_AttributeAdaptor;
	
  while($sth->fetch()) {
	
    push @ssvs, Bio::EnsEMBL::Variation::SupportingStructuralVariation->new
      (-dbID    => $ssv_id,
       -ADAPTOR => $self,
       -NAME    => $name,
       -STRUCTURAL_VARIATION_ID => $structural_variation_id,
       -CLASS_SO_TERM           => $aa->attrib_value_for_id($class_attrib_id),
			);
  }

  return \@ssvs;
}


1;
