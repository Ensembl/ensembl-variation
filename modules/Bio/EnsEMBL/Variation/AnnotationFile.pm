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
# Ensembl module for Bio::EnsEMBL::Variation::AnnotationFile
#
#

=head1 NAME

Bio::EnsEMBL::DBSQL::AnnotationFile

=head1 SYNOPSIS
  my $reg = 'Bio::EnsEMBL::Registry';
  my $annotation_file_adaptor = $reg->get_adaptor('human', 'variation', 'AnnotationFile');
  my $gerp_file = $annotation_file_adaptor->fetch_by_type('gerp');
  my $cadd_file = $annotation_file_adaptor->fetch_by_type('cadd');

=head1 DESCRIPTION
This module represents the base class for CADDFile and GERPFile modules.
It contains methods that are shared by both modules.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::AnnotationFile;

use Cwd;
use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

our %TYPES = (
  'remote' => 1,
  'local'  => 1,
);

=head2 new
  Arg [-ID]:                     string - identifier for this collection
  Arg [-DESCRIPTION]:            string - description for this collection
  Arg [-TYPE]:                   string - "local" or "remote"
  Arg [-FILENAME_TEMPLATE]:      string
  Arg [-ADAPTOR]:                Bio::EnsEMBL::Variation::DBSQL::AnnotationFileAdaptor
  Example    : my $annotation_file = Bio::EnsEMBL::Variation::AnnotationFile->new(
                -id                 => 'test',
                -type               => 'local',
                -filename_template  => '/path/to/annotation_files/annotation_file.gz',
               );
  Description: Constructor.  Instantiates a new AnnotationFile object.
  Returntype : Bio::EnsEMBL::Variation::AnnotationFile
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my (
    $id,
    $description,
    $type,
    $filename_template,
    $assembly,
    $adaptor,
    $tmpdir,
  ) = rearrange(
    [qw(
      ID
      DESCRIPTION
      TYPE
      FILENAME_TEMPLATE
      ASSEMBLY
      ADAPTOR
      TMPDIR
    )],
    @_
  ); 
  
  throw("ERROR: No id defined for collection") unless $id;
  throw("ERROR: Collection type $type invalid") unless $type && defined($TYPES{$type});

  my $self = bless {
    adaptor => $adaptor,
    id => $id,
    type => $type,
    filename_template => $filename_template,
    assembly  => $assembly,
    tmpdir => $tmpdir || cwd(),
  }, $class;

  return $self;
}

=head2 tmpdir
  Arg [1]    : string $tmpdir (optional)
               The new value to set the tmpdir attribute to
  Example    : my $tmpdir = $annotation_file->tmpdir()
  Description: Getter/Setter for the temporary directory path used when
               downloading indexes for remote tabix files.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable
=cut

sub tmpdir {
  my $self = shift;
  $self->{tmpdir} = shift if @_;
  return $self->{tmpdir};
}

=head2 filename_template
  Arg [1]    : string $filename_template (optional)
               The new value to set the filename_template attribute to
  Example    : my $filename_template = $collection->filename_template()
  Description: Getter/Setter for the filename template of this collection.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable
=cut

sub filename_template {
  my $self = shift;
  $self->{filename_template} = shift if @_;
  return $self->{filename_template};
}

sub _seek_by_Slice {
  my $self = shift;
  my $slice = shift;
  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }
  my $parser = $self->_file_parser_obj(); 
  $parser->seek($slice->seq_region_name, $slice->seq_region_start, $slice->seq_region_end);
  return $parser;
}

1;
