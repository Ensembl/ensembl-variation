=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::GenotypeCodeAdaptor
#
# Copyright (c) 2011 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::GenotypeCodeAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $va = $reg->get_adaptor("human","variation","genotypecode");



=head1 DESCRIPTION

This adaptor provides database connectivity for GenotypeCode objects.
GenotypeCodes may be retrieved from the Ensembl variation database by
several means using this module.

GenotypeCode objects are internal objects utilised by the CompressedGenotype
adaptor - they are not intended for external use.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::GenotypeCodeAdaptor;

use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;

use Bio::EnsEMBL::Variation::GenotypeCode;

use DBI qw(:sql_types);

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');


sub fetch_all_by_dbID_list {
  my ($self, $id_list_ref) = @_;

  if(!defined($id_list_ref) || ref($id_list_ref) ne 'ARRAY') {
    throw("id_list list reference argument is required");
  }

  return [] if(!@$id_list_ref);

  my @out;
  
  # mysql is faster and we ensure that we do not exceed the max query size by
  # splitting large queries into smaller queries of 200 ids
  my $max_size = 200;
  my @id_list = @$id_list_ref;

  while(@id_list) {
    my @ids;
    if(@id_list > $max_size) {
      @ids = splice(@id_list, 0, $max_size);
    } else {
      @ids = splice(@id_list, 0);
    }

    my $id_str;
    if(@ids > 1)  {
      $id_str = " IN (" . join(',', @ids). ")";
    } else {
      $id_str = " = ?";
      $self->bind_param_generic_fetch($ids[0],SQL_INTEGER);
    }

    my $constraint = "gc.genotype_code_id $id_str";

    push @out, @{$self->generic_fetch($constraint)};
  }

  return \@out;
}

# 
sub fetch_all_single_bp {
	my $self = shift;
	
	my $constraint = "length(ac.allele) = 1";
	
	return $self->generic_fetch($constraint);
}

sub _objs_from_sth {
    my $self = shift;
    my $sth = shift;
	
	my $ploidy = $self->ploidy;
	
    my ($gt_code_id, $haplotype_id, $allele);
    my (@result, %gts);
    
    $sth->bind_columns(\$gt_code_id, \$haplotype_id, \$allele);
    
	$gts{$gt_code_id}{$haplotype_id} = $allele while $sth->fetch;
	
	foreach $gt_code_id(keys %gts) {
		my @gt = map {$gts{$gt_code_id}{$_}} sort {$a <=> $b} keys %{$gts{$gt_code_id}};
		
		# splice it down to ploidy size
		@gt = splice @gt, 0, $ploidy;
		
		push @result, Bio::EnsEMBL::Variation::GenotypeCode->new_fast({
			dbID     => $gt_code_id,
			genotype => \@gt,
		});
	}
    
    return \@result;
}

# method used by superclass to construct SQL
sub _tables { 
    return (['genotype_code','gc'],['allele_code','ac']);
}

sub _columns {
  return qw( gc.genotype_code_id gc.haplotype_id ac.allele );
}

sub _default_where_clause {
	return 'gc.allele_code_id = ac.allele_code_id';
}

1;
