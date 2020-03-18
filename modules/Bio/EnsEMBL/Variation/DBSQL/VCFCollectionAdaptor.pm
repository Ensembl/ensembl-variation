=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::DBSQL::VCFCollectionAdaptor

=head1 SYNOPSIS
  my $reg = 'Bio::EnsEMBL::Registry';

  # set path to configuration file

  ## explicit use

  # get VCF Collection Adaptor
  my $vca = $reg->get_adaptor('human', 'variation', 'vcfcollection');
  my $va  = $reg->get_adaptor('human', 'variation', 'variation');

  my $v = $va->fetch_by_name('rs699');
  my $vf = $v->get_all_VariationFeatures->[0];

  # iterate over collections
  foreach my $c(@{$vca->fetch_all})

    # get samples
    my $samples = $c->get_all_Samples;

    # get genotypes for this VariationFeature
    my $gts = $c->get_all_SampleGenotypeFeatures_by_VariationFeature($vf);
  }

  ## implicit use

  # get LD Feature Container Adaptor
  my $ldfca = $reg->get_adaptor('human', 'variation', 'ldfeaturecontainer');
  my $sa    = $reg->get_adaptor('human', 'core', 'slice');

  # tell API to use _only_ VCF to retrieve genotypes (value of 1 uses both DB and VCF)
  $vca->db->use_vcf(2);

  # fetch LD Feature Container
  my $ldfc = $ldfca->fetch_by_Slice($sa->fetch_by_region('chromosome', 1, 230700048, 230720047));

=head1 DESCRIPTION

This module creates a set of objects that can read from tabix-indexed VCF files.

=head1 METHODS

=cut

use strict;
use warnings;


package Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor;

use JSON;
use Cwd;
use Net::FTP;
use URI;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Variation::VCFCollection;
use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;

use base qw(Bio::EnsEMBL::Variation::DBSQL::BaseAnnotationAdaptor);

our @EXPORT_OK = qw($CONFIG_FILE);

our $CONFIG_FILE;

=head2 new

  Arg [-CONFIG]: string - path to JSON configuration file
  Example    : my $vca = Bio::EnsEMBL::Variation::VCFCollectionAdaptor->new(
                 -config => '/path/to/vcf_config.json'
               );

  Description: Constructor.  Instantiates a new VCFCollectionAdaptor object.
  Returntype : Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor
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
#    next unless (defined $hash->{annotation_type} && lc $hash->{annotation_type} eq 'vcf');
    next if (defined $hash->{annotation_type});

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

    ## store populations if available
    my $populations;
    if( defined $hash->{populations}){
      my $prefix = $hash->{population_prefix} || '';
      my $display_group_data = $hash->{population_display_group} || {};

      foreach my $pop_id (keys %{$hash->{populations}}) {

        my $ref = $hash->{populations}->{$pop_id};
        my $hashref;

        if(ref($ref) eq 'HASH') {
          $hashref = $ref;
          
          # we need a name at least
          throw("ERROR: No name given in population from VCF config\n") unless $hashref->{name};

          $hashref->{dbID} ||= $pop_id;
        }
        else {
          $hashref = { name => $ref };
        }

        # add display group data
        $hashref->{$_} = $display_group_data->{$_} for keys %{$display_group_data};

        # add prefix
        if($prefix && !exists($hashref->{_raw_name})) {
          $hashref->{_raw_name} = $hashref->{name};
          $hashref->{name} = $prefix.$hashref->{_raw_name};
        }

        push @{$populations}, Bio::EnsEMBL::Variation::Population->new_fast($hashref);
      }
    }

    $self->add_VCFCollection(Bio::EnsEMBL::Variation::VCFCollection->new(
      -id => $hash->{id},
      -description => $hash->{description},
      -type => $hash->{type},
      -use_as_source => $hash->{use_as_source},
      -filename_template => $filename_template,
      -chromosomes => $hash->{chromosomes},
      -sample_prefix => $hash->{sample_prefix},
      -population_prefix => $hash->{population_prefix},
      -sample_populations => $hash->{sample_populations},
      -populations =>  $populations || undef ,
      -assembly  => $hash->{assembly} || undef,
      -source => $source || undef,
      -strict_name_match => $hash->{strict_name_match},
      -use_seq_region_synonyms => $hash->{use_seq_region_synonyms},
      -created =>$hash->{created} || undef,
      -updated =>$hash->{updated} || undef,
      -is_remapped => $hash->{is_remapped} ||0,
      -adaptor => $self,
      -tmpdir => $hash->{tmpdir} || $tmpdir,
      -ref_freq_index => $hash->{ref_freq_index},
    ));
  }
  return $self;
}

=head2 fetch_by_id

  Example    : my $collection = $vca->fetch_by_id('1000GenomesPhase3');
  Description: Fetches VCFCollection with given ID
  Returntype : Bio::EnsEMBL::Variation::VCFCollection
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_id {
  my $self = shift;
  return $self->SUPER::fetch_by_id(@_);
}


=head2 fetch_all

  Example    : my $collections = $vca->fetch_all();
  Description: Fetches all configured VCFCollections
  Returntype : Arrayref of Bio::EnsEMBL::Variation::VCFCollection
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all {
  my $self = shift;
  return $self->SUPER::fetch_all(@_);
}


=head2 fetch_all_for_web

  Example    : my $collections = $vca->fetch_all_for_web();
  Description: Fetches all configured VCFCollections that can be used
               as sources on the web
  Returntype : Arrayref of Bio::EnsEMBL::Variation::VCFCollection
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_for_web {
  return [grep {$_->use_as_source} @{$_[0]->fetch_all}];
}


=head2 add_VCFCollection

  Example    : $vca->add_VCFCollection($collection);
  Description: Add a VCFCollection to the cached list on this adaptor
  Returntype : Bio::EnsEMBL::Variation::VCFCollection
  Exceptions : throw on invalid or missing collection
               throw on attempt to add collection with duplicated ID
  Caller     : general
  Status     : Stable

=cut

sub add_VCFCollection {
  my $self = shift;
  my $collection = shift;

  # check class
  assert_ref($collection, 'Bio::EnsEMBL::Variation::VCFCollection');

  # check ID
  my $id = $collection->id;
  throw("ERROR: Collection with ID $id already exists\n") if $self->fetch_by_id($id);

  $self->{collections}->{$id} = $collection;
  push @{$self->{order}}, $id;

  return $collection;
}


=head2 remove_VCFCollection_by_ID

  Example    : $vca->remove_VCFCollection_by_ID($id);
  Description: Remove a VCFCollection from the cached list on this adaptor
  Returntype : bool (0 = ID does not exist; 1 = deleted OK)
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub remove_VCFCollection_by_ID {
  my $self = shift;
  my $id = shift;

  return 0 unless $self->fetch_by_id($id);

  delete $self->{collections}->{$id};
  @{$self->{order}} = grep {$_ ne $id} @{$self->{order}};

  return 1;
}

1;
