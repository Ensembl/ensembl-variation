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


# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor
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
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(raw_freqs_from_gts);

use Scalar::Util qw(weaken);

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
            next unless $term && $cons_map->{$term->name};
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
      return "consequence_types IN ('')";
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
    
    return [] unless $parent_obj;

    my $all_terms = $parent_obj->descendants;

    unshift @$all_terms, $parent_obj;

    return $all_terms;
}

sub _get_parent_terms {
    my ($self, $child_term) = @_;

    my $child_obj = $self->_get_term_object($child_term);
    
    return [] unless $child_obj;

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
  
  return {} unless scalar @$list;
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

sub _generic_fetch_all_by_Variation_from_Genotypes {
  my $self = shift;
  my $variation = shift;
  my $population = shift;
  my $object_type = shift;
  
  # Make sure that we are passed a Variation object
  assert_ref($variation,'Bio::EnsEMBL::Variation::Variation');
  
  # If we got a population argument, make sure that it is a Population object
  assert_ref($population,'Bio::EnsEMBL::Variation::Population') if (defined($population));

  # object type must be "Allele" or "PopulationGenotype"
  throw("Invalid object type \"$object_type\"") unless $object_type eq 'Allele' or $object_type eq 'PopulationGenotype';
  
  # fetch all genotypes
  my $genotypes = $variation->get_all_SampleGenotypes();
  
  return [] unless scalar @$genotypes;
  
  # copy sample ID to save time later
  # also attempt to get hash of missing data, may have been added if retrieving data from VCF
  my $missing = {};
  my %sample_ids = ();
  foreach my $gt(@$genotypes) {
    $gt->{_sample_id} ||= $gt->sample->dbID;
    $sample_ids{$gt->{_sample_id}} = 1;
    $missing = $gt->{_missing_alleles} if defined($gt->{_missing_alleles});
  }

  $sample_ids{$_} = 1 for keys %$missing;
  
  # get populations for samples
  my (@pop_list, $pop_hash);
  
  if(defined($population)) {
    @pop_list = ($population);
    my $pop_id = $population->dbID;
    $pop_hash->{$pop_id}->{$_->dbID} = 1 for @{$population->get_all_Samples};
  }
  else {
    my $pa = $self->db->get_PopulationAdaptor();
    $pop_hash = $pa->_get_sample_population_hash([keys %sample_ids]);
    return [] unless scalar keys %$pop_hash;
  
    @pop_list = @{$pa->fetch_all_by_dbID_list([keys %$pop_hash])};
  }
  
  return [] unless @pop_list and scalar keys %$pop_hash;
  
  my @objs;
  my $class = 'Bio::EnsEMBL::Variation::'.$object_type;

  foreach my $raw(
    @{raw_freqs_from_gts(
      $genotypes,
      [grep {$_->_freqs_from_gts} @pop_list],
      $pop_hash,
      $missing
    )->{$object_type}}
  ) {
    my $obj = $class->new_fast({
      %$raw,
      adaptor   => $self,
      variation => $variation
    });

    weaken($obj->{variation});

    push @objs, $obj;
  }
  
  return \@objs;
}

1;

