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
# Ensembl module for Bio::EnsEMBL::Variation::GERPAnnotation
#

=head1 NAME

Bio::EnsEMBL::DBSQL::GERPAnnotation

=head1 SYNOPSIS
  my $reg = 'Bio::EnsEMBL::Registry';

  my $gerp_annotation_adaptor = $reg->get_adaptor('human', 'variation', 'GERPAnnotation');

  # iterate over all GERPAnnotations
  foreach my $gerp_annotation (@{$gerp_annotation_adaptor->fetch_all})
    my $gerp_score = $gerp_annotation->get_score_by_VariationFeature($vf);
  }

=head1 DESCRIPTION

This module represents a GERP score annotation file or collection of
annotation files, e.g. one for each chromosome. Each file is represented by a
Bio::EnsEMBL::IO::Parser::BigWig object.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::GERPAnnotation;

use Cwd;

use Bio::EnsEMBL::Utils::Exception qw(throw);

use Bio::EnsEMBL::IO::Parser::BigWig;

use base qw(Bio::EnsEMBL::Variation::BaseAnnotation);

my $MAX_OPEN_FILES = 2;

=head2 new

  Arg [-ID]:                     string - identifier for this collection
  Arg [-DESCRIPTION]:            string - description for this collection
  Arg [-TYPE]:                   string - "local" or "remote"
  Arg [-USE_AS_SOURCE]:          boolean
  Arg [-FILENAME_TEMPLATE]:      string
  Arg [-CHROMOSOMES]:            arrayref of strings
  Arg [-STRICT_NAME_MATCH]:      boolean
  Arg [-USE_SEQ_REGION_SYNONYMS]:boolean
  Arg [-ADAPTOR]:                Bio::EnsEMBL::Variation::DBSQL::GERPAnnotationAdaptor

  Example    : my $collection = Bio::EnsEMBL::Variation::GERPAnnotation->new(
                -id                 => 'test',
                -type               => 'local',
                -annotation_type    => 'gerp',
                -filename_template  => '/path/to/gerp.bw'
               );

  Description: Constructor. Instantiates a new GERPAnnotation object.
  Returntype : Bio::EnsEMBL::Variation::DBSQL::GERPAnnotation
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);
  return $self;
}

=head2 get_score_by_VariationFeature
  Arg[1]     : Bio::EnsEMBL::Variation::VariationFeature $vf
  Example    : my $gerp_score = $gerp_file->get_score_by_VariationFeature($vf)
  Description: Return highest GERP score for region overlapped by variant.
  Returntype : float
  Exceptions : throws on wrong argument
  Caller     : Bio::EnsEMBL::Variation::VariationFeature
  Status     : Stable
=cut

sub get_score_by_VariationFeature {
  my $self = shift;
  my $vf = shift;
  if(!ref($vf) || !$vf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
    throw('Bio::EnsEMBL::Variation::VariationFeature arg expected');
  }
  my $parser = $self->_seek_by_VariationFeature($vf);
  return undef if (!defined $parser);
  my $max_score = undef;
  while ($parser->next) {
    my $score = $parser->get_score;
    $max_score = $score if (!defined $max_score);
    $max_score = $score if ($score > $max_score);
  }
  $max_score = sprintf("%.2f", $max_score) if (defined $max_score);
  return $max_score;
}

sub _seek_by_VariationFeature {
  my $self = shift;
  my $vf = shift;
  if(!ref($vf) || !$vf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
    throw('Bio::EnsEMBL::Variation::VariationFeature arg expected');
  }
  my $parser = $self->_file_parser_obj();
  return undef if (!defined $parser);
  my ($start, $end) = ($vf->seq_region_start, $vf->seq_region_end);
  ($start, $end) = ($end, $start) if ($start > $end);
  $parser->seek($vf->seq_region_name, $start - 1, $end);
  return $parser;
}

sub _file_parser_obj {
  my $self = shift;

  # change dir
  my $cwd = cwd();
  chdir($self->tmpdir);
  # open object (remote indexes get downloaded)
  my $obj = Bio::EnsEMBL::IO::Parser::BigWig->open($self->filename_template);
  # change back
  chdir($cwd);

  return $obj;
}


1;
