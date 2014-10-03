=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self;
  eval {$self = $class->SUPER::new(shift);};
  $self ||= {};

  my ($config_file) = rearrange([qw(CONFIG_FILE)], @_);
  
  # try and get config file from global variable or ENV
  $config_file ||= $CONFIG_FILE || $ENV{ENSEMBL_VARIATION_VCF_CONFIG_FILE};
  
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
  
  foreach my $hash(@{$config->{collections}}) {
    my $collection = Bio::EnsEMBL::Variation::VCFCollection->new(
      -id => $hash->{id},
      -type => $hash->{type},
      -filename_template => $hash->{filename_template},
      -chromosomes => $hash->{chromosomes},
      -individual_prefix => $hash->{individual_prefix},
      -individual_populations => $hash->{individual_populations},
      -adaptor => $self,
    );
    
    $self->{collections}->{$collection->id} = $collection;
    push @{$self->{order}}, $collection->id;
  }
  
  return $self;
}

sub fetch_by_id {
  return $_[0]->{collections}->{$_[1]};
}

sub fetch_all {
  return [values %{$_[0]->{collections}}];
}

1;
