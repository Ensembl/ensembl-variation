=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::CADDAnnotationAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::DBSQL::CADDAnnotationAdaptor

=head1 SYNOPSIS
  my $reg = 'Bio::EnsEMBL::Registry';

  # set path to configuration file

  ## explicit use

  # get CADD Annotation Adaptor
  my $cadd_annotation_adaptor = $reg->get_adaptor('human', 'variation', 'CADDAnnotation');
  my $va  = $reg->get_adaptor('human', 'variation', 'variation');

  my $v = $va->fetch_by_name('rs699');
  my $vf = $v->get_all_VariationFeatures->[0];

  # iterate over CADD annotation collections
  foreach my $cadd_annotation (@{$cadd_annotation_adaptor->fetch_all})
    # get CADD scores for this VariationFeature
    my $cadd_scores = $cadd_annotation->get_scores_by_VariationFeature($vf);
  }

=head1 DESCRIPTION

This module creates a set of objects that can read CADD scores from tabixed files.

=head1 METHODS

=cut

use strict;
use warnings;


package Bio::EnsEMBL::Variation::DBSQL::CADDAnnotationAdaptor;


use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Variation::CADDAnnotation;
use Bio::EnsEMBL::Variation::Source;
use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;

use Data::Dumper;

use base qw(Bio::EnsEMBL::Variation::DBSQL::BaseAnnotationAdaptor);

our @EXPORT_OK = qw($CONFIG_FILE);

our $CONFIG_FILE;

=head2 new

  Arg [-CONFIG]: string - path to JSON configuration file
  Example    : my $cadd_annotation_adaptor = Bio::EnsEMBL::Variation::CADDAnnotationAdaptor->new(
                 -config => '/path/to/annotation_config.json'
               );

  Description: Constructor.  Instantiates a new CADDAnnotationAdaptor object.
  Returntype : Bio::EnsEMBL::Variation::DBSQL::CADDAnnotationAdaptor
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  $Bio::EnsEMBL::Variation::DBSQL::BaseAnnotationAdaptor::CONFIG_FILE = $CONFIG_FILE if ($CONFIG_FILE); 

  my $self = $class->SUPER::new(@_); 

  my $config = $self->config; 
  my $root_dir = $self->root_dir;
  my $tmpdir = $self->tmpdir;
  foreach my $hash(@{$config->{collections}}) {
    next unless (defined $hash->{annotation_type} && lc $hash->{annotation_type} eq 'cadd');

    # check the species and assembly if we can
    if($self->db) {

      my $species = $hash->{species};
      my $assembly = $hash->{assembly};
      
      # check species
      throw("ERROR: No species defined in collection ".$hash->{id}."\n") unless defined($species);
      next unless lc($self->db->species) =~ lc($species);
      
      # check assembly
      if($assembly) {
        next unless lc($self->db->dnadb->get_CoordSystemAdaptor->fetch_all->[0]->version) eq lc($assembly);
      }
      else {
        warn("WARNING: No assembly defined in collection ".$hash->{id}."\n");
      }
    }

    my $filename_template = $self->_get_filename_template($hash);

    ## create source object if source info available
    my $source = Bio::EnsEMBL::Variation::Source->new_fast({
      name        => $hash->{source_name} || $hash->{id},
      version     => $hash->{source_version}, 
      url         => $hash->{source_url},
      description => $hash->{description} || $hash->{id},
    });

    $self->add_CADDAnnotation(
      Bio::EnsEMBL::Variation::CADDAnnotation->new(
        -id => $hash->{id},
        -description => $hash->{description},
        -type => $hash->{type},
        -use_as_source => $hash->{use_as_source},
        -filename_template => $filename_template,
        -chromosomes => $hash->{chromosomes},
        -assembly  => $hash->{assembly} || undef,
        -source => $source || undef,
        -strict_name_match => $hash->{strict_name_match},
        -use_seq_region_synonyms => $hash->{use_seq_region_synonyms},
        -created =>$hash->{created} || undef,
        -updated =>$hash->{updated} || undef,
        -is_remapped => $hash->{is_remapped} || 0,
        -adaptor => $self,
        -tmpdir => $hash->{tmpdir} || $tmpdir,
    ));
  }
  return $self;
}

=head2 fetch_by_id

  Example    : my $collection = $cadd_annotation_adaptor->fetch_by_id('cadd_version123');
  Description: Fetches CADDAnnotation with given ID
  Returntype : Bio::EnsEMBL::Variation::CADDAnnotation
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_id {
  my $self = shift;
  return $self->SUPER::fetch_by_id(@_);
}


=head2 fetch_all

  Example    : my $collections = $cadd_annotation_adaptor->fetch_all();
  Description: Fetches all configured CADDAnnotations
  Returntype : Arrayref of Bio::EnsEMBL::Variation::CADDAnnotations
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all {
  my $self = shift;
  return $self->SUPER::fetch_all(@_);
}

=head2 add_CADDAnnotation

  Example    : $cadd_annotation_adaptor->add_CADDAnnotation($collection);
  Description: Add a CADDAnnotation to the cached list on this adaptor
  Returntype : Bio::EnsEMBL::Variation::CADDAnnotation
  Exceptions : throw on invalid or missing collection
               throw on attempt to add collection with duplicated ID
  Caller     : general
  Status     : Stable

=cut

sub add_CADDAnnotation {
  my $self = shift;
  my $collection = shift;

  # check class
  assert_ref($collection, 'Bio::EnsEMBL::Variation::CADDAnnotation');

  # check ID
  my $id = $collection->id;
  throw("ERROR: Collection with ID $id already exists\n") if $self->fetch_by_id($id);

  $self->{collections}->{$id} = $collection;
  push @{$self->{order}}, $id;

  return $collection;
}

=head2 remove_CADDAnnotation_by_ID

  Example    : $cadd_annotation_adaptor->remove_CADDAnnotation_by_ID($id);
  Description: Remove a CADDAnnotation from the cached list on this adaptor
  Returntype : bool (0 = ID does not exist; 1 = deleted OK)
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub remove_CADDAnnotation_by_ID {
  my $self = shift;
  my $id = shift;

  return 0 unless $self->fetch_by_id($id);

  delete $self->{collections}->{$id};
  @{$self->{order}} = grep {$_ ne $id} @{$self->{order}};

  return 1;
}

1;
