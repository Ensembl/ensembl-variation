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


# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::BaseStructuralVariationAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::BaseStructuralVariationAdaptor

=head1 DESCRIPTION

Abstract adaptor class for fetching structural variants. Should not be invoked directly.

By default, the 'fetch_all_by_...'-methods will not return variations
that have been flagged as failed in the Ensembl QC. This behaviour can be modified
by setting the include_failed_variations flag in Bio::EnsEMBL::Variation::DBSQL::DBAdaptor.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::BaseStructuralVariationAdaptor;

use Bio::EnsEMBL::Variation::BaseStructuralVariation;
use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use DBI qw(:sql_types);

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');

my $DEFAULT_ITERATOR_CACHE_SIZE = 10000;

# method used by superclass to construct SQL
sub _tables { 
  my $self = shift;
  my @tables = (['structural_variation', 'sv']);
  
  # If we are excluding failed_structural_variations, add that table
  push(@tables,['failed_structural_variation', 'fsv']) unless ($self->db->include_failed_variations());
  
  return @tables;
}

sub _columns {
  return qw( sv.structural_variation_id sv.variation_name sv.validation_status sv.source_id sv.class_attrib_id
             sv.study_id sv.is_evidence sv.somatic sv.alias sv.clinical_significance sv.copy_number);
}

# Add a left join to the failed_structural_variation table
sub _left_join {
  my $self = shift;
  
  # If we are including failed structural variations, skip the left join
  return () if ($self->db->include_failed_variations());
  return (['failed_structural_variation', 'fsv.structural_variation_id=sv.structural_variation_id']);
}


=head2 fetch_all

  Description: Returns a listref of all germline structural variants
  Returntype : listref of Bio::EnsEMBL::Variation::StructuralVariation or Bio::EnsEMBL::Variation::SupportingStructuralVariation objects
  Status     : Stable

=cut

sub fetch_all {
    my $self = shift;
    my $constraint = 'sv.somatic = 0';
    return $self->generic_fetch($constraint);
}


=head2 fetch_all_somatic

  Description: Returns a listref of all somatic structural variants
  Returntype : listref of Bio::EnsEMBL::Variation::StructuralVariation or Bio::EnsEMBL::Variation::SupportingStructuralVariation objects
  Status     : Stable

=cut

sub fetch_all_somatic {
    my $self = shift;
    my $constraint = 'sv.somatic = 1';
    return $self->generic_fetch($constraint);
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
    ReturnType  : Bio::EnsEMBL::Variation::StructuralVariation or Bio::EnsEMBL::Variation::SupportingStructuralVariation object
    Exceptions  : thrown if there are multiple objects found with the same variation name
    Caller      : general
    Status      : Stable

=cut

sub fetch_by_name {
  my ($self, $name) = @_;

  # hack fix for phencode names with '+' in getting encoded/decoded by the webcode 
  $name =~ s/\s/\+/ if $self->db->species =~ /homo_sapiens/i; 
    
  my $constraint = $self->_internal_exclude_failed_constraint("sv.variation_name='$name'");
  
  my $objs = $self->generic_fetch($constraint);
  throw("Multiple structural variations found with the same name: '$name'") if @$objs > 1;
  return $objs->[0] if @$objs == 1;
}

# alias for fetch_by_name
sub fetch_by_stable_id {
  my $self = shift;
  return $self->fetch_by_name(@_);
}


=head2 fetch_all_by_Study

  Arg [1]     : Bio::EnsEMBL::Variation::Study $study_id
  Example     : my $study = $study_adaptor->fetch_by_name('estd1');
                foreach my $sv (@{$sv_adaptor->fetch_all_by_Study($study)}){
                   print $sv->variation_name,"\n";
                }
  Description : Retrieves all structural variations from a specified study
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::StructuralVariation or 
                Bio::EnsEMBL::Variation::SupportingStructuralVariation objects
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
  
  my $constraint = $self->_internal_exclude_failed_constraint('sv.study_id = '.$study->dbID);
  
  my $result = $self->generic_fetch($constraint);

  return $result;
}


=head2 fetch_all_by_Source

  Arg [1]     : Bio::EnsEMBL::Variation::Source $source_id
  Example     : my $source = $source_adaptor->fetch_by_name('DGVa');
                foreach my $sv (@{$sv_adaptor->fetch_all_by_Source($source)}){
                   print $sv->variation_name,"\n";
                }
  Description : Retrieves all structural variations from a specified source
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::StructuralVariation or 
                Bio::EnsEMBL::Variation::SupportingStructuralVariation objects
  Exceptions  : throw if incorrect argument is passed
                warning if provided source does not have a dbID
  Caller      : general
  Status      : At Risk

=cut

sub fetch_all_by_Source {
  my $self = shift;
  my $source = shift;

  if(!ref($source) || !$source->isa('Bio::EnsEMBL::Variation::Source')) {
    throw("Bio::EnsEMBL::Variation::Source arg expected");
  }
    
  if(!$source->dbID()) {
    warning("Source does not have dbID, cannot retrieve structural variants");
    return [];
  } 
  
  my $constraint = $self->_internal_exclude_failed_constraint('sv.source_id = '.$source->dbID);
  
  my $result = $self->generic_fetch($constraint);

  return $result;
}


=head2 fetch_all_by_dbID_list

  Arg [1]    : listref $list
  Example    : $ssv = $sv_adaptor->fetch_all_by_dbID_list([907,1132]);
  Description: Retrieves a listref of structural variant objects via a list of internal
               dbID identifiers
  Returntype : listref of Bio::EnsEMBL::Variation::StructuralVariation or
               Bio::EnsEMBL::Variation::SupportingStructuralVariation objects
  Exceptions : throw if list argument is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_dbID_list {
  my $self = shift;
  my $list = shift;

  if(!defined($list) || ref($list) ne 'ARRAY') {
    throw("list reference argument is required");
  }
  
  my $id_str = (@$list > 1)  ? " IN (".join(',',@$list).")"   :   ' = \''.$list->[0].'\'';
  
  my $constraint = $self->_internal_exclude_failed_constraint("sv.structural_variation_id $id_str");
  
  my $result = $self->generic_fetch($constraint);

  return $result;
}


=head2 fetch_Iterator_by_dbID_list

  Arg [1]    : reference to list of ints $list
  Example    : $variation_iterator = $va->fetch_Iterator_by_dbID_list([124, 56, 90]);
  Description: Retrieves an iterator over a set of structural variations via their internal identifiers.
  Returntype : Bio::EnsEMBL::Utils::Iterator
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Experimental

=cut

sub fetch_Iterator_by_dbID_list {
    my ($self, $dbid_list, $cache_size) = @_;
    
    unless ((defined $dbid_list) && (ref $dbid_list eq 'ARRAY')) {
        throw("list reference argument is required");
    }

    $cache_size ||= $DEFAULT_ITERATOR_CACHE_SIZE;

    # create an iterator that fetches structural variations in blocks of
    # $cache_size and returns them in turn

    my @object_cache;

    return Bio::EnsEMBL::Utils::Iterator->new(sub {

            if (@object_cache == 0 && @$dbid_list > 0 ) {
                my @dbids = splice @$dbid_list, 0, $cache_size;
                
                # Create a constraint on the dbIDs
                my $id_str = "(" . join(",",@dbids) . ")";
                my $constraint = qq{sv.structural_variation_id IN $id_str};
                
                @object_cache = @{ $self->generic_fetch($constraint) };
            }

            return shift @object_cache;
        }
    );
}


# Exclude the constraint for failed structural variant
sub _internal_exclude_failed_constraint {
  my $self = shift;
  my $constraint = shift;
  $constraint .= " AND " . $self->db->_exclude_failed_structural_variations_constraint();
  
  return $constraint;
}


# API-internal method for getting failed descriptions for an Allele
sub _internal_get_failed_descriptions {
    my $self = shift;
    my $allele = shift;
    my $constraint = shift;
    
    # Assert that the object passed is an Allele
    assert_ref($allele,'Bio::EnsEMBL::Variation::BaseStructuralVariation');
    
    my $stmt = qq{
        SELECT DISTINCT
            fd.description
        FROM
            failed_structural_variation fsv JOIN
            failed_description fd ON (
                fd.failed_description_id = fsv.failed_description_id
            )
        WHERE
            fsv.structural_variation_id = ?
    };
    $stmt .= qq{ AND $constraint } if (defined($constraint));
    
    my $sth = $self->prepare($stmt);
    $sth->execute($allele->dbID());
    my @descriptions;
    my $description;
    $sth->bind_columns(\$description);
    while ($sth->fetch()) {
        push(@descriptions,$description);
    }
    return \@descriptions;
}


=head2 get_all_failed_descriptions

  Arg[1]      : Bio::EnsEMBL::Variation::BaseStructuralVariation $sv
                 The structural variant object to get the failed descriptions for
  Example     : 
                my $failed_descriptions = $adaptor->get_all_failed_descriptions($sv);
                if (scalar(@{$failed_descriptions})) {
                  print "The structural variant'" . $sv->variation_name . "' has been flagged as failed because '" . join("' and '",@{$failed_descriptions}) . "'\n";
                }
    
  Description : Gets the unique descriptions for the reasons why the supplied structural variant has failed.
  ReturnType  : reference to list of strings
  Exceptions  : thrown on incorrect argument
  Caller      : general
  Status      : Stable

=cut

sub get_all_failed_descriptions {
    my $self = shift;
    my $sv = shift;
    
    # Call the internal get method without any constraints
    my $description = $self->_internal_get_failed_descriptions($sv) || [];
    
    return $description;
}


sub store {
    my ($self, $sv) = @_;
    
    my $dbh = $self->dbc->db_handle;
    
    # look up source_id
    if(!defined($sv->{_source_id})) {
      my $sth = $dbh->prepare(q{
           SELECT source_id FROM source WHERE name = ?
      });
      $sth->execute($sv->source_name);
        
      my $source_id;
      $sth->bind_columns(\$source_id);
      $sth->fetch();
      $sth->finish();
      $sv->{_source_id} = $source_id;
    }
    throw("No source ID found for source name ", $sv->source_name) unless defined($sv->{_source_id});
    
    # look up study_id
    if(!defined($sv->{_study_id}) && defined($sv->study)) {
      my $sth = $dbh->prepare(q{
           SELECT study_id FROM study WHERE name = ?
      });
      $sth->execute($sv->study->name);
      
      my $study_id;  
      $sth->bind_columns(\$study_id);
      $sth->fetch();
      $sth->finish();
      $sv->{_study_id} = $study_id;
    }
    
    # look up class_attrib_id
    my $class_attrib_id;
    if(defined($sv->{class_SO_term})) {
      my $sth = $dbh->prepare(q{
           SELECT attrib_id FROM attrib WHERE value = ?
      });
      $sth->execute($sv->{class_SO_term});
        
      $sth->bind_columns(\$class_attrib_id);
      $sth->fetch();
      $sth->finish();
    }
    throw("No class ID found for the class name ", $sv->{class_SO_term}) unless defined($class_attrib_id);
    
    my $sth = $dbh->prepare(q{
        INSERT INTO structural_variation (
            source_id,
            study_id,
            variation_name,
            validation_status,
            class_attrib_id,
            is_evidence,
            somatic,
            alias,
            clinical_significance,
            copy_number
        ) VALUES (?,?,?,?,?,?,?,?,?,?)
    });
    
    $sth->execute(
        $sv->{_source_id},
        $sv->{_study_id} || undef,
        $sv->variation_name,
        $sv->validation_status || undef,
        $class_attrib_id || 0,
        $sv->is_evidence || 0,
        $sv->is_somatic  || 0,
        $sv->alias || undef,
        $sv->get_all_clinical_significance_states ? (join ",", @{$sv->get_all_clinical_significance_states}) : undef,
        $sv->copy_number || undef
    );
    
    $sth->finish;
    
    # get dbID
    my $dbID = $dbh->last_insert_id(undef, undef, 'structural_variation', 'structural_variation_id');
    $sv->{dbID}    = $dbID;
    $sv->{adaptor} = $self;
}

1;
