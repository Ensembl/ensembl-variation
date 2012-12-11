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


# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor
#
# Copyright (c) 2010 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor

=head1 SYNOPSIS

Abstract class - should not be instantiated.  Implementation of
abstract methods must be performed by subclasses.

=head1 DESCRIPTION

This adaptor provides generic database connectivity for various Variation objects.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use Bio::EnsEMBL::Variation::Utils::Constants qw(%OVERLAP_CONSEQUENCES);

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');

sub AttributeAdaptor {
    my $self = shift;
    
    unless ($self->{_attribute_adaptor}) {
        $self->{_attribute_adaptor} = $self->db->get_AttributeAdaptor if $self->db;
    }
    
    return $self->{_attribute_adaptor};
}

sub _consequence_type_map {
    
    # return a hash mapping between the string terms of a mysql set and
    # the equivalent numerical values

    my ($self, $table, $column) = @_;
    
    my $key = $table.'_'.$column.'_map';

    unless ($self->{$key}) {

        my $map_sth = $self->prepare(qq{SHOW COLUMNS FROM $table LIKE '$column'});
        
        $map_sth->execute;

        my $types = $map_sth->fetchrow_hashref->{Type};

        # Type will look like set('type1','type2'), so tidy it up a bit before splitting
        
        $types =~ s/set\(//;
        $types =~ s/\)//;
        $types =~ s/'//g;

        my $map;

        my $pow = 0;

        # mysql assigns the set values in consecutive powers of 2, so so shall we

        for my $type (split /,/, $types) {
            $map->{$type} = 2**$pow++;
        }

        $self->{$key} = $map;
    }
    
    return $self->{$key};
}

sub _get_consequence_constraint {
    
    my ($self, $table, $query_terms, $without_children, $term_subset) = @_;

    # we build up the numerical value for our query by ORing together all the children of all the terms
    my $query = 0;

    # get a hash mapping consequence terms to numerical values (specifically powers of 2)
    my $cons_map = $self->_consequence_type_map($table, 'consequence_types');

    for my $query_term (@$query_terms) {

        # we allow either an ontology term object, or just a string
        $query_term = UNIVERSAL::can($query_term, 'name') ? $query_term->name : $query_term;
    
        # we store only the most specific consequence term, so we need to get all children of 
        # each query term
        my $terms = $without_children ? [ ($self->_get_term_object($query_term)) ] : $self->_get_child_terms($query_term);

        # and then we OR together all relevant terms

        for my $term (@$terms) {
            next unless $cons_map->{$term->name};
            $query |= $cons_map->{$term->name};
        }
    }

    my $subset_mask;
    if ($term_subset) {
        for my $query_term (@$term_subset) {
    
            # we allow either an ontology term object, or just a string
            $query_term = UNIVERSAL::can($query_term, 'name') ? $query_term->name : $query_term;
        
            my $terms = [ ($self->_get_term_object($query_term)) ];
    
            # and then we OR together all relevant terms
    
            for my $term (@$terms) {
                next unless $cons_map->{$term->name};
                $subset_mask |= $cons_map->{$term->name};
            }
        }
    }

    unless ($self->{_possible_consequences}) {

        # we need a list of the numerical values of all possible 
        # consequence term combinations we have actually observed

        my $sth = $self->dbc->prepare(qq{
            SELECT DISTINCT(consequence_types)
            FROM $table
        });

        $sth->execute;

        my $cons;

        $sth->bind_columns(\$cons);

        my @poss_cons;

        while ($sth->fetch) {
            # construct the numerical value by ORing together each combination
            # (this is much quicker than SELECTing consequence_types+0 above which
            # is what I used to do, but this seems to mean the db can't use the index)
        
            my $bit_val = 0;
            
            for my $term (split /,/, $cons) {
                $bit_val |= $cons_map->{$term};
            }

            push @poss_cons, $bit_val;
        }

        $self->{_possible_consequences} = \@poss_cons;
    }

    # we now find all combinations that include our query by ANDing 
    # the query with all possible combinations and combine these into 
    # our query string

    #my $id_str = join ',', grep { $_ & $query } @{ $self->{_possible_consequences} }; 
    my @cons_vals =  grep { $_ & $query } @{ $self->{_possible_consequences} }; 

    if ($subset_mask) {
        # When only including a subset of types, filter combinations to ones which
        # include at least one of the the subset types.
        @cons_vals =  grep { $_ & $subset_mask } @cons_vals;
    }

    if (!scalar(@cons_vals)) {
      return undef;
    }
   
    my $id_str = join ',', @cons_vals;

    my $constraint = "consequence_types IN ($id_str)"; 

    return $constraint;
}

sub _consequences_for_set_number {
    my ($self, $set_number, $map) = @_;

    my @consequences;

    for my $term (keys %$map) {
        if ($set_number & $map->{$term}) {
            push @consequences, $OVERLAP_CONSEQUENCES{$term};
        }
    }

    return \@consequences;
}

sub _get_term_object {
    my ($self, $term) = @_;

    my $oa = $self->{_ontology_adaptor} ||=
        Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );

    my $terms = $oa->fetch_all_by_name($term, 'SO');

    if (@$terms > 1) {
        warn "Ambiguous term '$term', just using first result";
    }
    elsif (@$terms == 0) {
        warn "Didn't find an ontology term for '$term'";
    }

    return $terms->[0];
}

sub _get_child_terms {
    my ($self, $parent_term) = @_;

    my $parent_obj = $self->_get_term_object($parent_term);

    my $all_terms = $parent_obj->descendants;

    unshift @$all_terms, $parent_obj;

    return $all_terms;
}

sub _get_parent_terms {
    my ($self, $child_term) = @_;

    my $child_obj = $self->_get_term_object($child_term);

    my $all_terms = $child_obj->ancestors;

    unshift @$all_terms, $child_obj;

    return $all_terms;
}

sub _set_number_for_consequences {
    my ($self, $cons_list, $map) = @_;

    my $val = 0;

    for my $cons (@$cons_list) {
        $val |= $map->{$cons->SO_term};
    }

    return $val;
}

sub _transcript_variation_consequences_for_set_number {
    my ($self, $set_number) = @_;
    my $map = $self->_consequence_type_map('transcript_variation', 'consequence_types');
    return $self->_consequences_for_set_number($set_number, $map);
}

sub _variation_feature_consequences_for_set_number {
    my ($self, $set_number) = @_;
    my $map = $self->_consequence_type_map('variation_feature', 'consequence_types');
    return $self->_consequences_for_set_number($set_number, $map);
}

sub _transcript_variation_set_number_for_consequences {
    my ($self, $cons) = @_;
    my $map = $self->_consequence_type_map('transcript_variation', 'consequence_types');
    return $self->_set_number_for_consequences($cons, $map);
}

sub _variation_feature_set_number_for_consequences {
    my ($self, $cons) = @_;
    my $map = $self->_consequence_type_map('variation_feature', 'consequence_types');
    return $self->_set_number_for_consequences($cons, $map);
}

sub _get_all_subsnp_handles_from_variation_ids {
	my $self = shift;
	my $list = shift;
    
    my $in_list = join ",", @$list;
    
	my $sth = $self->dbc->prepare(qq{
		SELECT v.variation_id, sh.handle
		FROM allele a, variation v, subsnp_handle sh
		WHERE a.variation_id = v.variation_id
		AND a.subsnp_id = sh.subsnp_id
		AND v.variation_id IN ($in_list)
		GROUP BY v.variation_id, sh.handle
	});
	
	$sth->execute();
	
	my ($var_id, $handle, %handles);
	$sth->bind_columns(\$var_id, \$handle);
	
	push @{$handles{$var_id}}, $handle while $sth->fetch();
	$sth->finish;
	
	return \%handles;
}

=head2 ploidy

  Arg[1]      : int $ploidy
  Example     : my $ploidy = $adaptor->ploidy();
  Description : Gets/sets the ploidy for this database
  ReturnType  : int
  Exceptions  : None
  Caller      : general
  Status      : At Risk

=cut

sub ploidy {
	my $self = shift;
	my $ploidy = shift;
	
	if(defined($ploidy)) {
		$self->{ploidy} = $ploidy;
	}
	elsif(!defined($self->{ploidy})) {
		my $mc = $self->db->get_MetaContainer;
		throw("Could not retrieve MetaContainer") unless defined($mc);
		
		$self->{ploidy} = $mc->ploidy;
	}
	
	return $self->{ploidy};
}

1;

