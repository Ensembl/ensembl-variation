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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::GERPAnnotationAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::DBSQL::GERPAnnotationAdaptor

=head1 SYNOPSIS
  my $reg = 'Bio::EnsEMBL::Registry';

  # set path to configuration file

  ## explicit use

  # get GERP Annotation Adaptor
  my $gerp_annotation_adaptor = $reg->get_adaptor('human', 'variation', 'GERPAnnotation');
  my $va  = $reg->get_adaptor('human', 'variation', 'variation');

  my $v = $va->fetch_by_name('rs699');
  my $vf = $v->get_all_VariationFeatures->[0];

  # iterate over GERP annotation collections
  foreach my $gerp_annotation (@{$gerp_annotation_adaptor->fetch_all})
    # get GERP score for this VariationFeature
    my $gerp_score = $gerp_annotation->get_score_by_VariationFeature($vf);
  }


=head1 DESCRIPTION

This module creates a set of objects that can read GERP scores from BigWig files.

=head1 METHODS

=cut

use strict;
use warnings;


package Bio::EnsEMBL::Variation::DBSQL::GERPAnnotationAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Variation::GERPAnnotation;
use Bio::EnsEMBL::Variation::Source;
use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;

use base qw(Bio::EnsEMBL::Variation::DBSQL::BaseAnnotationAdaptor);

our @EXPORT_OK = qw($CONFIG_FILE);

our $CONFIG_FILE;

=head2 new

  Arg [-CONFIG]: string - path to JSON configuration file
  Example    : my $gaa = Bio::EnsEMBL::Variation::DBSQL::GERPAnnotationAdaptor->new(
                 -config => '/path/to/annotation_config.json'
               );

  Description: Constructor.  Instantiates a new GERPAnnotationAdaptor object.
  Returntype : Bio::EnsEMBL::Variation::DBSQL::GERPAnnotationAdaptor
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
    next unless (defined $hash->{annotation_type} && lc $hash->{annotation_type} eq 'gerp');

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

    $self->add_GERPAnnotation(
      Bio::EnsEMBL::Variation::GERPAnnotation->new(
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

  Example    : my $collection = $gerp_annotation_adaptor->fetch_by_id('gerp_conservation_scores');
  Description: Fetches GERPAnnotation with given ID
  Returntype : Bio::EnsEMBL::Variation::GERPAnnotation
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_id {
  my $self = shift;
  return $self->SUPER::fetch_by_id(@_);
}


=head2 fetch_all

  Example    : my $collections = $gerp_annotation_adaptor->fetch_all();
  Description: Fetches all configured GERPAnnotations
  Returntype : Arrayref of Bio::EnsEMBL::Variation::GERPAnnotations
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all {
  my $self = shift;
  return $self->SUPER::fetch_all(@_);
}


=head2 add_GERPAnnotation

  Example    : $gerp_annotation_adaptor->add_GERPAnnotation($collection);
  Description: Add a GERPAnnotation to the cached list on this adaptor
  Returntype : Bio::EnsEMBL::Variation::GERPAnnotation
  Exceptions : throw on invalid or missing collection
               throw on attempt to add collection with duplicated ID
  Caller     : general
  Status     : Stable

=cut

sub add_GERPAnnotation {
  my $self = shift;
  my $collection = shift;

  # check class
  assert_ref($collection, 'Bio::EnsEMBL::Variation::GERPAnnotation');

  # check ID
  my $id = $collection->id;
  throw("ERROR: Collection with ID $id already exists\n") if $self->fetch_by_id($id);

  $self->{collections}->{$id} = $collection;
  push @{$self->{order}}, $id;

  return $collection;
}

=head2 remove_GERPAnnotation_by_ID

  Example    : $gerp_annotation_adaptor->remove_GERPAnnotation_by_ID($id);
  Description: Remove a GERPAnnotation from the cached list on this adaptor
  Returntype : bool (0 = ID does not exist; 1 = deleted OK)
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub remove_GERPAnnotation_by_ID {
  my $self = shift;
  my $id = shift;

  return 0 unless $self->fetch_by_id($id);

  delete $self->{collections}->{$id};
  @{$self->{order}} = grep {$_ ne $id} @{$self->{order}};

  return 1;
}

=head2 root_dir

  Example    : $root_dir = $gerp_annotation_adaptor->root_dir;
  Description: Sets and returns the root directory.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut


sub root_dir {
  my ($self, $root_dir) = @_;
  if (!$root_dir) {
    if (!$self->{root_dir}) {
      if($self->db && $self->db->gerp_root_dir) {
        $self->{root_dir} = $self->db->gerp_root_dir.'/';
      }
    }
  } else {
    $self->{root_dir} = $root_dir;
  }
  return $self->{root_dir} || '';
}




1;
