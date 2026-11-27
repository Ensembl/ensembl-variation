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
use Test::Exception;
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

my @selected_args = (undef, -config_file => $config_file, -collection_ids => [$selected->{id}]);

subtest 'Remote checks use the URL scheme' => sub {
  no warnings 'redefine';
  my (@http_urls, @ftp_hosts, @ftp_paths);
  local *HTTP::Tiny::head = sub {
    push @http_urls, $_[1];
    return { success => 1 };
  };
  local *Net::FTP::new = sub {
    my ($class, $host) = @_;
    push @ftp_hosts, $host;
    return bless {}, $class;
  };
  local *Net::FTP::login = sub { 1 };
  local *Net::FTP::size = sub { push @ftp_paths, $_[1]; return 42; };
  local *Net::FTP::quit = sub { 1 };

  my $adaptor = Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor->new(@selected_args);
  my $http_url = 'http://example.org/selected.vcf.gz';
  $adaptor->_get_filename_template({ %$selected, filename_template => $http_url });
  $adaptor->_get_filename_template({ %$selected, filename_template => 'ftp://example.org/selected.vcf.gz' });

  is_deeply(\@http_urls, [$selected->{filename_template}, $http_url], 'HTTPS and HTTP URLs use HEAD');
  is_deeply(\@ftp_hosts, ['example.org'], 'Only the FTP URL opens an FTP connection');
  is_deeply(\@ftp_paths, ['/selected.vcf.gz'], 'The FTP URL is checked with SIZE');
};

subtest 'HTTPS check diagnostics' => sub {
  no warnings 'redefine';
  my $response = { success => 0, status => 599, content => 'Connection timed out' };
  local *HTTP::Tiny::head = sub { return $response };

  throws_ok { Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor->new(@selected_args) }
    qr/\Q$selected->{id}\E.*\Q$selected->{filename_template}\E.*Connection timed out/,
    'A connection error identifies the collection and URL';

  $response = { success => 0, status => 404, reason => 'Not Found' };
  my @warnings;
  local $SIG{__WARN__} = sub { push @warnings, @_ };
  my $adaptor = Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor->new(@selected_args);
  ok($adaptor->fetch_by_id($selected->{id}), 'A missing file allows collection setup to continue');
  like(join('', @warnings), qr/\Q$selected->{filename_template}\E.*\Q$selected->{id}\E/,
    'The missing-file warning identifies the URL and collection');
};

done_testing();
