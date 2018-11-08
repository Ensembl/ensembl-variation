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

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::AnnotationFileAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::DBSQL::AnnotationFileAdaptor

=head1 SYNOPSIS
  my $reg = 'Bio::EnsEMBL::Registry';

  # set path to configuration file
  # optionally it can be set as environment variable $ENSEMBL_VARIATION_ANNOTATION_CONFIG_FILE
  $Bio::EnsEMBL::Variation::DBSQL::AnnotationFileAdaptor::CONFIG_FILE = '/path/to/annotation_config.json';

  ## explicit use

  # get Annotation File Adaptor
  my $annotation_file_adaptor = $reg->get_adaptor('human', 'variation', 'AnnotationFile');
  my $va  = $reg->get_adaptor('human', 'variation', 'variation');
  my $v = $va->fetch_by_name('rs699');
  my $vf = $v->get_all_VariationFeatures->[0];

  my $gerp_file = $annotation_file_adaptor->fetch_by_type('gerp');
  my $cadd_file = $annotation_file_adaptor->fetch_by_type('cadd');
  my $gerp_score = $gerp_file->get_score_by_VariationFeature($vf);
  my $cadd_scores = $gerp_file->get_scores_by_VariationFeature($vf);

  # get scores from variation feature
  my $gerp_score = $vf->get_gerp_score;
  my $cadd_score = $vf->get_cadd_scores;


=head1 DESCRIPTION

This module creates a set of objects that can read annotations from files.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::AnnotationFileAdaptor;

use JSON;
use Cwd;
use Net::FTP;
use URI;
use Data::Dumper;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Variation::CADDFile;
use Bio::EnsEMBL::Variation::GERPFile;
use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');

use base qw(Exporter);
our @EXPORT_OK = qw($CONFIG_FILE);

our $CONFIG_FILE;


=head2 new

  Arg [-CONFIG]: string - path to JSON configuration file
  Example    : my $vca = Bio::EnsEMBL::Variation::AnnotationFileAdaptor->new(
                 -config => '/path/to/annotation_config.json'
               );

  Description: Constructor.  Instantiates a new AnnotationFileAdaptor object.
  Returntype : Bio::EnsEMBL::Variation::DBSQL::AnnotationFileAdaptor
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
  
  bless($self, $class);

  my $config = {};

  # try and get config from DB adaptor
  $config = $self->db->annotation_config if $self->db;


  unless($config && scalar keys %$config) {
    my ($config_file) = rearrange([qw(CONFIG_FILE)], @_);
    
    # try and get config file from global variable or ENV
    $config_file ||= $CONFIG_FILE || ($self->db ? $self->db->annotation_config_file : undef) || $ENV{ENSEMBL_VARIATION_ANNOTATION_CONFIG_FILE};
    
    # try and find default config file in API dir
    if(!defined($config_file)) {
      my $mod_path  = 'Bio/EnsEMBL/Variation/DBSQL/AnnotationFileAdaptor.pm';
      $config_file  = $INC{$mod_path};
      $config_file =~ s/AnnotationFileAdaptor\.pm/annotation_config\.json/ if $config_file;
    }
    throw("ERROR: No config file defined") unless defined($config_file);
    throw("ERROR: Config file $config_file does not exist") unless -e $config_file;
    
    # read config from JSON config file
    open IN, $config_file or throw("ERROR: Could not read from config file $config_file");
    local $/ = undef;
    my $json_string = <IN>;
    close IN;
    
    # parse JSON into hashref $config
    $config = JSON->new->decode($json_string) or throw("ERROR: Failed to parse config file $config_file");
  }
  
  $self->{config} = $config;
  $self->db->annotation_config($config) if $self->db;

  ## set up root dir
  my $root_dir = '';
  if($ENV{ENSEMBL_VARIATION_ANNOTATION_ROOT_DIR}) {
    $root_dir = $ENV{ENSEMBL_VARIATION_ANNOTATION_ROOT_DIR}.'/';
  }
  elsif($self->db && $self->db->annotation_root_dir) {
    $root_dir = $self->db->annotation_root_dir.'/';
  }
  
  $self->{root_dir} = $root_dir;

  my $tmpdir = cwd();
  if($ENV{ENSEMBL_VARIATION_ANNOTATION_TMP_DIR}) {
    $tmpdir = $ENV{ENSEMBL_VARIATION_ANNOTATION_TMP_DIR}.'/';
  }
  elsif($self->db && $self->db->annotation_tmp_dir) {
    $tmpdir = $self->db->annotation_tmp_dir.'/';
  }

  $self->{tmpdir} = $tmpdir;

  throw("ERROR: No annotation data defined in config file") unless $config->{annotation_type} && keys %{$config->{annotation_type}};
 
  if ($self->db) {
    $self->{species} = lc($self->db->species);
    $self->{assembly} = lc($self->db->dnadb->get_CoordSystemAdaptor->fetch_all->[0]->version);
  }

  return $self;
}

=head2 fetch_by_annotation_type
  Example    : my $gerp_file = $annotation_file_adaptor->fetch_by_annotation_type('gerp');
               my $cadd_file = $annotation_file_adaptor->fetch_by_annotation_type('cadd');
  Description: Fetches AnnotationFile by given annotation type
  Returntype : Bio::EnsEMBL::Variation::GERPFile or Bio::EnsEMBL::Variation::CADDFile
  Exceptions : Throws on wrong annotation type
  Caller     : general
  Status     : Stable
=cut

sub fetch_by_annotation_type {
  my $self = shift;
  my $type = shift;

  my $annotation_hash = $self->_get_annotation_hash($type);

  if ($type eq 'gerp') {
    return $self->_get_gerp_annotation_file($annotation_hash);
  } elsif ($type eq 'cadd') {
    return $self->_get_cadd_annotation_file($annotation_hash);
  } else {
    die("Unknown annotation type. Annotation type must be cadd or gerp.");
  }
}

sub _get_annotation_hash {
  my $self = shift;
  my $annotation_type = shift;
  my $config = $self->{config};
  my $species = $self->{species};
  my $assembly = $self->{assembly};
  my @annotations  = grep {lc $_->{species} eq $species && lc $_->{assembly} eq $assembly} @{$config->{annotation_type}->{$annotation_type}};
  if (scalar @annotations == 1) {
    return $annotations[0];
  } elsif (scalar @annotations > 1) {
    die "Found more than one annotation file for annotation type $annotation_type for species $species and assembly $assembly\n";
  } else {
    die "Could not find annotation file for annotation type $annotation_type for species $species and assembly $assembly\n";
  }
}

sub _get_gerp_annotation_file {
  my $self = shift;
  my $hash = shift;
  return Bio::EnsEMBL::Variation::GERPFile->new(
    -id => $hash->{id},
    -type => $hash->{type},
    -filename_template => $self->_get_filename_template($hash),
    -assembly  => $hash->{assembly} || undef,
    -adaptor => $self,
    -tmpdir => $self->{tmpdir},
  );
}

sub _get_cadd_annotation_file {
  my $self = shift;
  my $hash = shift;
  return Bio::EnsEMBL::Variation::CADDFile->new(
    -id => $hash->{id},
    -type => $hash->{type},
    -filename_template => $self->_get_filename_template($hash),
    -assembly  => $hash->{assembly} || undef,
    -adaptor => $self,
    -tmpdir => $self->{tmpdir},
  );
}

sub _get_filename_template {
  my $self = shift;
  my $hash = shift;
  my $root_dir = $self->{root_dir};
  
  my $filename_template = $hash->{filename_template} =~ /(nfs|ftp:)/ ? $hash->{filename_template} : $root_dir.$hash->{filename_template};
  my $file_exists = ($hash->{type} eq 'remote') ? $self->_ftp_file_exists($filename_template) : (-e $filename_template);
  if (!$file_exists) {
    warn("WARNING: Cannot read from file $filename_template for species " . $hash->{species} );
  }
  return $filename_template;
}

# Internal method checking if a remote file exists
sub _ftp_file_exists {
  my $self = shift;
  my $uri = URI->new(shift);

  my $ftp = Net::FTP->new($uri->host) or die "Connection error($uri): $@";
  $ftp->login('anonymous', 'guest') or die "Login error", $ftp->message;
  my $exists = defined $ftp->size($uri->path);
  $ftp->quit;

  return $exists;
}

1;
