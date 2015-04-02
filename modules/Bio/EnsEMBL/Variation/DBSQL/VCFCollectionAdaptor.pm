=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
  # optionally it can be set as environment variable $ENSEMBL_VARIATION_VCF_CONFIG_FILE
  $Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor::CONFIG_FILE = '/path/to/vcf_config.json';

  ## explicit use

  # get VCF Collection Adaptor
  my $vca = $reg->get_adaptor('human', 'variation', 'vcfcollection');
  my $va  = $reg->get_adaptor('human', 'variation', 'variation');

  my $v = $va->fetch_by_name('rs699');
  my $vf = $v->get_all_VariationFeatures->[0];

  # iterate over collections
  foreach my $c(@{$vca->fetch_all})

    # get individuals
    my $individuals = $c->get_all_Individuals;

    # get genotypes for this VariationFeature
    my $gts = $c->get_all_IndividualGenotypeFeatures_by_VariationFeature($vf);
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

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Variation::VCFCollection;
use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');

use base qw(Exporter);
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
  
  my $self;
  eval {$self = $class->SUPER::new(shift);};
  $self ||= {};

  my ($config_file) = rearrange([qw(CONFIG_FILE)], @_);
  
  # try and get config file from global variable or ENV
  $config_file ||= $CONFIG_FILE || $ENV{ENSEMBL_VARIATION_VCF_CONFIG_FILE};
  
  # try and find default config file in API dir
  if(!defined($config_file)) {
    my $mod_path  = 'Bio/EnsEMBL/Variation/DBSQL/VCFCollectionAdaptor.pm';
    $config_file  = $INC{$mod_path};
    $config_file =~ s/VCFCollectionAdaptor\.pm/vcf_config\.json/ if $config_file;
  }
  
  throw("ERROR: No config file defined") unless defined($config_file);
  throw("ERROR: Config file $config_file does not exist") unless -e $config_file;
  
  # read config from JSON config file
  open IN, $config_file or throw("ERROR: Could not read from config file $config_file");
  local $/ = undef;
  my $json_string = <IN>;
  close IN;
  
  # parse JSON into hashref $config
  my $config = JSON->new->decode($json_string) or throw("ERROR: Failed to parse config file $config_file");
  
  $self->{config} = $config;
  
  bless($self, $class);
  
  throw("ERROR: No collections defined in config file") unless $config->{collections} && scalar @{$config->{collections}};
  
  $self->{collections} = {};
  $self->{order} = [];
  
  foreach my $hash(@{$config->{collections}}) {
    
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

    ## create source object if source info available
    my $source;
    if( defined $hash->{source_name}){
      $source = Bio::EnsEMBL::Variation::Source->new
          (-name        => $hash->{source_name},
           -version     => $hash->{source_version}, 
           -url         => $hash->{source_url}
        );
    }

    ## store population names if available
    my $populations;
    if( defined $hash->{populations}){
      foreach my $pop_id (keys %{$hash->{populations}}){
 
        my $pop = Bio::EnsEMBL::Variation::Population->new
          (-name        => $hash->{populations}->{$pop_id},
           -dbID        => $pop_id,
        );
        push @{$populations}, $pop;
      }
    }

    
    my $collection = Bio::EnsEMBL::Variation::VCFCollection->new(
      -id => $hash->{id},
      -type => $hash->{type},
      -filename_template => $hash->{filename_template},
      -chromosomes => $hash->{chromosomes},
      -individual_prefix => $hash->{individual_prefix},
      -population_prefix => $hash->{population_prefix},
      -individual_populations => $hash->{individual_populations},
      -populations =>  $populations ||undef ,
      -assembly  => $hash->{assembly} || undef,
      -source => $source || undef,
      -adaptor => $self,
    );
    
    $self->{collections}->{$collection->id} = $collection;
    push @{$self->{order}}, $collection->id;
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
  return $_[0]->{collections}->{$_[1]};
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
  return [values %{$_[0]->{collections} || {}}];
}

1;
