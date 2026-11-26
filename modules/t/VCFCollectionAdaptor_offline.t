# Copyright [2016-2026] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;
use File::Temp qw(tempfile);
use JSON qw(encode_json);
use Test::More;

use Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor;

my $selected = {
  id => 'selected.v1_GRCh38',
  species => 'homo_sapiens',
  assembly => 'GRCh38',
  type => 'remote',
  filename_template => 'https://example.org/selected.vcf.gz',
};
my $eva = {
  %$selected,
  id => 'eva_pongo',
  species => 'pongo_abelii',
  filename_template => 'https://example.org/eva.vcf.gz',
};
my $grch37 = {
  %$selected,
  id => 'selected.v1_GRCh37',
  assembly => 'GRCh37',
  filename_template => 'https://example.org/grch37.vcf.gz',
};
my $other_id = {
  %$selected,
  id => 'selected',
  filename_template => 'https://example.org/unrequested.vcf.gz',
};
my @collections = ($selected, $eva, $grch37, $other_id);

my ($config_fh, $config_file) = tempfile(UNLINK => 1);
print {$config_fh} encode_json({ collections => \@collections });
close $config_fh or die "Cannot close $config_file: $!";

my @checked_urls;
{
  no warnings 'redefine';
  local *Bio::EnsEMBL::Variation::DBSQL::BaseAnnotationAdaptor::_get_filename_template = sub {
    my ($adaptor, $collection) = @_;
    push @checked_urls, $collection->{filename_template};
    return $collection->{filename_template};
  };

  my @species = (-species => 'HOMO_SAPIENS');
  my @assembly = (-assembly => 'grch38');
  my @collection_ids = (-collection_ids => ['SELECTED.V1_GRCH38']);
  my @cases = (
    ['all filters', [@species, @assembly, @collection_ids], [$selected]],
    ['species only', \@species, [$selected, $grch37, $other_id]],
    ['assembly only', \@assembly, [$selected, $eva, $other_id]],
    ['collection ID only', \@collection_ids, [$selected]],
    ['no filters', [], \@collections],
  );

  for my $case (@cases) {
    my ($label, $filters, $expected) = @$case;
    @checked_urls = ();
    my $adaptor = Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor->new(
      undef, -config_file => $config_file, @$filters,
    );

    is_deeply([map { $_->id } @{$adaptor->fetch_all}], [map { $_->{id} } @$expected],
      "Expected collections are loaded with $label");
    is_deeply(\@checked_urls, [map { $_->{filename_template} } @$expected],
      "Expected files are checked with $label");
  }
}

done_testing();
