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
# Ensembl module for Bio::EnsEMBL::Variation::CADDFile
#
#

=head1 NAME

Bio::EnsEMBL::DBSQL::CADDFile

=head1 SYNOPSIS
  my $reg = 'Bio::EnsEMBL::Registry';
  # get Annotation File Adaptor
  my $annotation_file_adaptor = $reg->get_adaptor('human', 'variation', 'AnnotationFile');
  my $va  = $reg->get_adaptor('human', 'variation', 'variation');
  my $v = $va->fetch_by_name('rs699');
  my $vf = $v->get_all_VariationFeatures->[0];

  my $cadd_file = $annotation_file_adaptor->fetch_by_type('cadd');
  my $cadd_scores = $gerp_file->get_scores_by_VariationFeature($vf);
  my $cadd_score = $vf->get_cadd_scores;

=head1 DESCRIPTION
This module represents a CADDFile. It uses the CADDTabix parser for retrieving CADD scores
from a tabixed file.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::CADDFile;

use Cwd;
use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::IO::Parser::CADDTabix;

use base qw(Bio::EnsEMBL::Variation::AnnotationFile);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);  
  return $self;
}

=head2 get_scores_by_VariationFeature
  Arg[1]     : Bio::EnsEMBL::Variation::VariationFeature $vf
  Example    : my $cadd_scores = $cadd_file->get_scores_by_VariationFeature($vf)
  Description: Return CADD phred scores for all possible alternative alleles of
               length 1 for a variant.
  Returntype : hasref of alternative_allele => phred_score
  Exceptions : throws on wrong argument
               warns if length of variation feature is not 1 
  Caller     : Bio::EnsEMBL::Variation::VariationFeature
  Status     : Stable
=cut

sub get_scores_by_VariationFeature {
  my $self = shift;
  my $vf = shift;
  if(!ref($vf) || !$vf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
    throw('Bio::EnsEMBL::Variation::VariationFeature arg expected');
  }

  my $start = $vf->seq_region_start;
  my $end = $vf->seq_region_end;

  if ($end - $start > 0) {
    warn("Can only calculate CADD scores for variants of length 1");
    return {};
  }

  my $parser = $self->_seek_by_VariationFeature($vf);
  my $scores = {};
  while ($parser->next) {
    my $phred_score = $parser->get_phred_score;
    my $alt_allele = $parser->get_alternative;
    $scores->{$alt_allele} = $phred_score;
  } 
  return $scores;
}

sub _seek_by_VariationFeature {
  my $self = shift;
  my $vf = shift;
  if(!ref($vf) || !$vf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
    throw('Bio::EnsEMBL::Variation::VariationFeature arg expected');
  }
  my $parser = $self->_file_parser_obj();
  $parser->seek($vf->seq_region_name, $vf->seq_region_start, $vf->seq_region_end);
  return $parser;
}

sub _file_parser_obj {
  my $self = shift;

  # change dir
  my $cwd = cwd();
  chdir($self->tmpdir);

  # open object (remote indexes get downloaded)
  my $obj = Bio::EnsEMBL::IO::Parser::CADDTabix->open($self->filename_template);

  # change back
  chdir($cwd);

  return $obj;
}

1;
