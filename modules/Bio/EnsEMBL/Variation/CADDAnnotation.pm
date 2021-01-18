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
# Ensembl module for Bio::EnsEMBL::Variation::CADDAnnotation
#
#

=head1 NAME

Bio::EnsEMBL::DBSQL::CADDAnnotation

=head1 SYNOPSIS

  my $reg = 'Bio::EnsEMBL::Registry';

  my $cadd_annotation_adaptor = $reg->get_adaptor('human', 'variation', 'CADDAnnotation');

  # iterate over all CADDAnnotations
  foreach my $cadd_annotation (@{$cadd_annotation_adaptor->fetch_all})
    my $cadd_scores = $cadd_annotation->get_all_scores_by_VariationFeature($vf);
  }
  

=head1 DESCRIPTION

This module represents a CADD score annotation file or collection of
annotation files, e.g. one for each chromosome. Each file is represented by a
Bio::EnsEMBL::IO::Parser::CADDTabix object.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::CADDAnnotation;

use Cwd;

use Bio::EnsEMBL::Utils::Exception qw(throw);

use Bio::EnsEMBL::IO::Parser::CADDTabix;

use base qw(Bio::EnsEMBL::Variation::BaseAnnotation);


=head2 new

  Arg [-ID]:                     string - identifier for this collection
  Arg [-DESCRIPTION]:            string - description for this collection
  Arg [-TYPE]:                   string - "local" or "remote"
  Arg [-USE_AS_SOURCE]:          boolean
  Arg [-FILENAME_TEMPLATE]:      string
  Arg [-CHROMOSOMES]:            arrayref of strings
  Arg [-STRICT_NAME_MATCH]:      boolean
  Arg [-USE_SEQ_REGION_SYNONYMS]:boolean
  Arg [-ADAPTOR]:                Bio::EnsEMBL::Variation::DBSQL::CADDAnnotationAdaptor

  Example    : my $collection = Bio::EnsEMBL::Variation::CADDAnnotation->new(
                -id                 => 'test',
                -type               => 'local',
                -annotation_type    => 'cadd',
                -filename_template  => '/path/to/cadd.tsv.gz',
               );

  Description: Constructor.  Instantiates a new CADDAnnotation object.
  Returntype : Bio::EnsEMBL::Variation::CADDAnnotation
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


=head2 get_all_scores_by_VariationFeature
  Arg[1]     : Bio::EnsEMBL::Variation::VariationFeature $vf
  Example    : my $cadd_scores = $cadd_file->get_all_scores_by_VariationFeature($vf)
  Description: Return CADD phred scores for all possible alternative alleles of
               length 1 for a variant.
  Returntype : hashref of alternative_allele => phred_score
  Exceptions : throws on wrong argument
               warns if length of variation feature is not 1
  Caller     : Bio::EnsEMBL::Variation::VariationFeature
  Status     : Stable
=cut

sub get_all_scores_by_VariationFeature {
  my $self = shift;
  my $vf = shift;
  if(!ref($vf) || !$vf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
    throw('Bio::EnsEMBL::Variation::VariationFeature arg expected');
  }

  my $start = $vf->seq_region_start;
  my $end = $vf->seq_region_end;

  if ($end - $start > 0 || $end < $start) {
    warn("Can only calculate CADD scores for variants of length 1");
    return {};
  }

  my $parser = $self->_seek_by_VariationFeature($vf);
  my $scores = {};
  while ($parser && $parser->next) {
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
  if ($parser && $parser->seek($vf->seq_region_name, $vf->seq_region_start, $vf->seq_region_end)) {
    return $parser;
  }
}

sub _file_parser_obj {
  my $self = shift;

  # change dir
  my $cwd = cwd();
  chdir($self->tmpdir);

  # open object (remote indexes get downloaded) locally, therefore we change the current directory to point to the tmp directory
  my $obj = undef;
  eval { $obj = Bio::EnsEMBL::IO::Parser::CADDTabix->open($self->filename_template); }; warn $@ if $@; 
  # change back
  chdir($cwd);

  return $obj;
}

1;
