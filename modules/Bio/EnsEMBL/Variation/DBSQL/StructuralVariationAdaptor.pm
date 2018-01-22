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


# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::StructuralVariationAdaptor
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
  
  # Modify the include_failed_variations flag in DBAdaptor to also return structural variations that have been flagged as failed
  $va->db->include_failed_variations(1);
 
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
use Bio::EnsEMBL::Variation::DBSQL::BaseStructuralVariationAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use DBI qw(:sql_types);

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseStructuralVariationAdaptor');

my $DEFAULT_ITERATOR_CACHE_SIZE = 10000;

sub _default_where_clause {
  my $self = shift;
  return 'is_evidence=0';
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

    push @svs, Bio::EnsEMBL::Variation::StructuralVariation->new(
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


=head2 fetch_all_by_supporting_evidence

  Arg [1]     : Bio::EnsEMBL::Variation::SupportingStructuralVariation or 
                Bio::EnsEMBL::Variation::StructuralVariation $se
  Example     : my $se = $ssv_adaptor->fetch_by_name('essv2585133');
                foreach my $sv (@{$sv_adaptor->fetch_all_by_supporting_evidence($se)}){
                  print $sv->variation_name,"\n";
                }
  Description : Retrieves all structural variations from a specified supporting evidence
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::StructuralVariation objects
  Exceptions  : throw if incorrect argument is passed
                warning if provided supporting evidence does not have a dbID
  Caller      : general
  Status      : At Risk

=cut

sub fetch_all_by_supporting_evidence {
  my $self = shift;
  my $se = shift;

  if(!ref($se) || (!$se->isa('Bio::EnsEMBL::Variation::SupportingStructuralVariation') &&
                    !$se->isa('Bio::EnsEMBL::Variation::StructuralVariation'))
  ) {
    throw("Bio::EnsEMBL::Variation::SupportingStructuralVariation or Bio::EnsEMBL::Variation::StructuralVariation arg expected");
  }
    
  if(!$se->dbID()) {
    warning("The supporting evidence does not have dbID, cannot retrieve supporting evidence");
    return [];
  } 
  
  my $cols   = join ",", $self->_columns();

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
  
  my $constraint = $self->_default_where_clause();
  
  # Add the constraint for failed structural variant
  $constraint .= " AND " . $self->db->_exclude_failed_structural_variations_constraint();
  
  my $sth = $self->prepare(qq{
                    SELECT $cols 
                    FROM $tables, structural_variation_association sa
                    WHERE $constraint 
                      AND sa.structural_variation_id=sv.structural_variation_id
                      AND sa.supporting_structural_variation_id = ?});
  $sth->bind_param(1,$se->dbID,SQL_INTEGER);
  $sth->execute();

  my $results = $self->_objs_from_sth($sth);

  $sth->finish();

  return $results;
}


sub _generic_fetch_by_VariationSet {
    my $self = shift;
    my $want_iterator = shift;
    my $set = shift;
    
    assert_ref($set,'Bio::EnsEMBL::Variation::VariationSet');

    if(!defined($set->dbID())) {
        warning("Cannot retrieve structural variations for variation set without a dbID");
        return [];
    }
  
    # Get the unique dbIDs for all variations in this set and all of its subsets
    my $dbid_list = $self->fetch_all_dbIDs_by_VariationSet($set);
 
    my $num_vars = @$dbid_list;

    if ($num_vars > 100_000 && !$want_iterator) {
        warn "This set contains a large number ($num_vars) of structural variations, these may not fit".
             "into memory at once, considering using fetch_Iterator_by_VariationSet instead";
    }

    # Use the dbIDs to get all variations and return them
    return $want_iterator ? 
        $self->fetch_Iterator_by_dbID_list($dbid_list) : 
        $self->fetch_all_by_dbID_list($dbid_list);
}


=head2 fetch_all_dbIDs_by_VariationSet

  Arg [1]    : Bio::EnsEMBL::Variation::VariationSet
  Example    : @sv_ids = @{$sva_adaptor->fetch_all_dbIDs_by_VariationSet($vs)};
  Description: Gets an array of internal ids of all structural variations which are present 
               in a specified variation set and its subsets.
  Returntype : listref of integers
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_dbIDs_by_VariationSet {
  my $self = shift;
  my $set  = shift;
  
  # First, get ids for all subsets,
  
  my @var_set_ids = ($set->dbID);
  
  foreach my $var_set (@{$set->adaptor->fetch_all_by_super_VariationSet($set)}) {
    push @var_set_ids, $var_set->dbID;
  }
  
  my $set_str = "(" . join(",",@var_set_ids) .")";

  # Add the constraint for failed structural variations
  my $constraint = $self->_internal_exclude_failed_constraint;
  
  # Then get the dbIDs for all these sets
  my $stmt = qq{
    SELECT DISTINCT
      vssv.structural_variation_id
    FROM
      variation_set_structural_variation vssv LEFT JOIN
      failed_structural_variation fsv ON (
        fsv.structural_variation_id = vssv.structural_variation_id
      )
    WHERE
      vssv.variation_set_id in $set_str
      $constraint
  };

  my $sth = $self->prepare($stmt);
  
  $sth->execute();
  
  my @result;
  my $dbID;
  
  $sth->bind_columns(\$dbID);
  
  while ($sth->fetch()) {
        push @result, $dbID;
  }

  return \@result;
}


=head2 fetch_all_by_VariationSet

  Arg [1]    : Bio::EnsEMBL::Variation::VariationSet
  Example    : @svs = @{$sva_adaptor->fetch_all_by_VariationSet($vs)};
  Description: Retrieves all structural variations which are present in a specified
               variation set and its subsets.
  Returntype : listref of Bio::EnsEMBL::Variation::StructuralVariation
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_VariationSet {
    my $self = shift;
    return $self->_generic_fetch_by_VariationSet(0, @_);
}


=head2 fetch_Iterator_by_VariationSet

  Arg [1]    : Bio::EnsEMBL::Variation::VariationSet
  Example    : $sv_iterator = $sva_adaptor->fetch_Iterator_by_VariationSet($vs);
  Description: Retrieves an iterator for all structural variations which are present 
               in a specified variation set and its subsets.
  Returntype : Bio::EnsEMBL::Utils::Iterator object
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Experimental

=cut

sub fetch_Iterator_by_VariationSet {
    my $self = shift;
    my $set = shift;
    my $cache_size = shift || $DEFAULT_ITERATOR_CACHE_SIZE;
    
    # First, get ids for all subsets,
    my @var_set_ids = ($set->dbID);
    map {push(@var_set_ids,$_->dbID())} @{$set->adaptor->fetch_all_by_super_VariationSet($set)};
    my $var_set_id = join(",",@var_set_ids);
    
    # Prepare a query for getting the span of variation_ids
    my $stmt = qq{
        FROM
            variation_set_structural_variation vssv LEFT JOIN
            failed_structural_variation fsv ON (
              fsv.structural_variation_id = vssv.structural_variation_id
            )
        WHERE
            vssv.variation_set_id IN ($var_set_id)
    };
      
    # Add the constraint for failed structural variations
    my $constraint = $self->_internal_exclude_failed_constraint;
    
    my $sth = $self->prepare(qq{SELECT MIN(vssv.structural_variation_id), MAX(vssv.structural_variation_id) 
                                $stmt $constraint});
    $sth->execute();
    my ($min_sv_id,$max_sv_id);
    $sth->bind_columns(\$min_sv_id,\$max_sv_id);
    $sth->fetch();
    $max_sv_id ||= 0;
    $min_sv_id ||= 1;
    
    # Prepare a statement for getting the ids in a range
    $sth = $self->prepare(qq{SELECT vssv.structural_variation_id $stmt 
                             AND vssv.structural_variation_id BETWEEN ? AND ? $constraint});
    
    # Internally, we keep an Iterator that works on the dbID span we're at
    my $iterator;
        
    return Bio::EnsEMBL::Utils::Iterator->new(sub {

        # If the iterator is empty, get a new chunk of dbIDs, unless we've fetched all dbIDs 
        unless ((defined($iterator) && $iterator->has_next()) || $min_sv_id > $max_sv_id) {
            
            ## check there are ids in the range to return
            my $count_sth = $self->prepare(qq{SELECT count(vssv.structural_variation_id) $stmt 
                                              AND vssv.structural_variation_id BETWEEN ? AND ? $constraint});

            my $done = 0;
            while( $min_sv_id < $max_sv_id && $done == 0  ){
                $count_sth->execute($min_sv_id, $min_sv_id+$cache_size);
                my $count = $count_sth->fetchall_arrayref();
                if ($count->[0]->[0] > 0){
                    $done =1;
                }
                else{
                    $min_sv_id += ($cache_size + 1);
                }
            }

            # Get the next chunk of dbIDs
            $sth->execute($min_sv_id,$min_sv_id+$cache_size);
            $min_sv_id += ($cache_size + 1);
            
            # Use a hash to keep track of the seen dbIDs
            my %seen;
            
            # Loop over the dbIDs and avoid duplicates
            my $dbID;
            my @dbIDs;
            $sth->bind_columns(\$dbID);
            while ($sth->fetch()) {
                push (@dbIDs,$dbID) unless ($seen{$dbID}++);
            }
    
            # Get a new Iterator based on the new dbID span
            $iterator = $self->fetch_Iterator_by_dbID_list(\@dbIDs);
            
        }
        
        return $iterator->next();
    });
}

1;
