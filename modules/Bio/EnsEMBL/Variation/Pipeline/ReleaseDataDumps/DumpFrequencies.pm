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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::DumpFrequencies;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

use FileHandle;
use Bio::EnsEMBL::IO::Parser::BedTabix;
use Bio::EnsEMBL::IO::Parser::VCF4Tabix;

sub run {
  my $self = shift;
  my $file_type = $self->param_required('file_type');
  $self->set_variation_col_header;
  $self->set_chr_from_filename;
  $self->set_tabix_parser;
  if ($file_type eq 'gvf') {
    $self->add_frequencies_to_gvf_file();
  } elsif ($file_type eq 'vcf') {
    $self->add_frequencies_to_vcf_file();
  } else {
    die("File type ($file_type) is not supported. File type must be gvf or vcf.");
  }
}

sub set_variation_col_header {
  my $self = shift;
  my $vep_cache_dir = $self->param_required('vep_cache_dir');
  my $release = $self->param_required('ensembl_release');
  my $assembly = $self->param_required('assembly');
  my $info_file = "$vep_cache_dir/homo_sapiens/$release\_$assembly/info.txt";
  if (! -e $info_file) {
    die("File ($info_file) doesn't exist.");
  }
  my $fh = FileHandle->new($info_file, 'r');
  while (<$fh>) {
    next if /^#/;
    chomp;
    my ($param, $value) = split "\t";
    if ($param =~ /variation_col/) {
      my @header = split(',', $value);
      $self->param('variation_col_header', \@header);
    }
  }
  $fh->close;
}

sub set_chr_from_filename {
  my $self = shift;
  my $input_file = $self->param_required('input_file');
  $input_file =~ /.*-chr(.*)\.(gvf|vcf)\.gz/;
  my $chr = $1; 
  if (!$chr) {
    die("Could not get chromosome name from file name ($input_file).");
  }
  $self->param('chr', $chr);
}

sub set_tabix_parser {
  my $self = shift;
  my $vep_cache_dir = $self->param_required('vep_cache_dir');
  my $release = $self->param_required('ensembl_release');
  my $assembly = $self->param_required('assembly');  
  my $chr = $self->param('chr');
  my $tabix_file = "$vep_cache_dir/homo_sapiens/$release\_$assembly/$chr/all_vars.gz";
  if (! -e $tabix_file) {
    die("File ($tabix_file) doesn't exist.");
  }
  my $parser = Bio::EnsEMBL::IO::Parser::BedTabix->open($tabix_file);
  $self->param('parser', $parser);
}

sub add_frequencies_to_gvf_file {
  my $self = shift;
  my $gvf_file = $self->param_required('input_file');
  my $output_file_name = $self->param_required('output_file_name');

  my $homo_sapiens_dump_dir = $self->param_required('homo_sapiens_dump_dir');
  if (! -d $homo_sapiens_dump_dir) {
    die("Directory ($homo_sapiens_dump_dir) doesn't exist");
  }

  my $chr = $self->param('chr');

  my $dir = "$homo_sapiens_dump_dir/gvf/homo_sapiens/";
  if (! -d $dir) {
    die("Directory ($dir) doesn't exist");
  }

  if (! -e "$dir/$gvf_file") {
    die("File ($dir/$gvf_file) doesn't exist.");
  }

  my $gvf_in = FileHandle->new("$dir/$gvf_file", 'r');
  my $gvf_out = FileHandle->new("$dir/$output_file_name-$chr.gvf", 'w');

  my $step_size = $self->param_required('step_size');
  # start and end are used to get data from all_vars tabix file
  my $start = 0;
  my $end = $step_size;
  my $overlap = $self->param_required('overlap');
  # populate a hash $frequencies which is populated from the all_vars tabix file
  # and maps variant to frequencies: {rs372722017 => AMR=0.0588;SAS=0.0038;EUR=0.0167}; 
  # frequencies are only reported for alt alleles
  my $frequencies = $self->get_frequencies($chr, $start, $end);

  while (<$gvf_in>) {
    chomp;
    my @values = split("\t", $_);
    my $chrom = $values[0];
    my $variant_start = $values[3];
    my $variant_end = $values[4];
    if ($variant_end >= $end) {
      $start = $variant_start;
      $end   = $variant_end + $step_size;
      # GVF file is sorted
      # If variant_end becomes bigger than end which we used for getting data from
      # all_vars file we need to update the $frequencies hash with new data
      $frequencies = $self->get_frequencies($chr, $start - $overlap, $end + $overlap);
    }
    my $info = $values[8];
    my @info_values = split(';', $info);
    my $rs = '';
    foreach (@info_values) {
      if (/^Dbxref/) {
        my ($db, $name) = split(':');
        $rs = $name;
      }
    }
    my $id_2_freqs = $frequencies->{$rs};
    next unless ($id_2_freqs);

    $info = $info . ';' . $id_2_freqs;
    print $gvf_out join("\t", $values[0], $values[1], $values[2], $values[3], $values[4], $values[5], $values[6], $values[7], $info), "\n";
  }

  $gvf_in->close;
  $gvf_out->close;
}

sub add_frequencies_to_vcf_file {
  my $self = shift;
  my $vcf_file = $self->param_required('input_file');
  my $vcf_file_with_frequencies = $self->param_required('output_file_name');
  my $homo_sapiens_dump_dir = $self->param_required('homo_sapiens_dump_dir');
  if (! -d $homo_sapiens_dump_dir) {
    die("Directory ($homo_sapiens_dump_dir) doesn't exist");
  }

  my $chr = $self->param('chr');

  my $dir = "$homo_sapiens_dump_dir/vcf/homo_sapiens/";
  if (! -d $dir) {
    die("Directory ($dir) doesn't exist");
  }

  if (! -e "$dir/$vcf_file") {
    die("File ($dir/$vcf_file) doesn't exist.");
  }

  my $vcf_parser = Bio::EnsEMBL::IO::Parser::VCF4Tabix->open("$dir/$vcf_file");

  my $cdba = $self->get_adaptor('homo_sapiens', 'core');
  my $sa = $cdba->get_SliceAdaptor;
  my $slice = $sa->fetch_by_region('chromosome', $chr);
  my $chr_end = $slice->end;
  $vcf_parser->seek($chr, 1, $chr_end);
  my $vcf_out = FileHandle->new("$dir/$vcf_file_with_frequencies-$chr.vcf", 'w');

  my $step_size = $self->param_required('step_size');
  my $start = 0;
  my $end = $step_size;
  my $overlap = $self->param_required('overlap');
  my $frequencies = $self->get_frequencies($chr, $start, $end);

  while ($vcf_parser->next) {
    my @values = @{$vcf_parser->{record}};
    my $variant_start = $values[1];
    my $rs = $values[2];
    my $ref = $values[3];
    my $variant_end = $variant_start + length($ref) - 1;
    if ($variant_end >= $end) {
      $start = $variant_start;
      $end   = $variant_end + $step_size;
      $frequencies = $self->get_frequencies($chr, $start - $overlap, $end + $overlap);
    }
    my $id_2_freqs = $frequencies->{$rs};
    next unless ($id_2_freqs);
    my $info = $values[7];
    $info = $info . ';' . $id_2_freqs;
    print $vcf_out join("\t", $values[0], $values[1], $values[2], $values[3], $values[4], $values[5], $values[6], $info), "\n";
  }
  $vcf_out->close;
  $vcf_parser->close;
  my $parser = $self->param('parser');
  $parser->close;
}

sub get_frequencies {
  my ($self, $chr, $start, $end) = @_;
  my $parser = $self->param('parser');
  my @variation_col_header = @{$self->param('variation_col_header')};
  my @af_keys = @{$self->param('af_keys')};
  $parser->seek($chr, $start, $end);
  my $frequencies = {};
  while ($parser->next) {
    my @record = @{$parser->read_record};
    # create a hash which maps column name to column values for the specific row
    my %data = map {$variation_col_header[$_] => $record[$_]} (0..$#variation_col_header);

    my $id = $data{variation_name};
    # we report frequencies for alt alleles only
    # allele strings are the same for a given variant in database and all_vars file
    my @allele_string = split('/', $data{allele_string});
    my $ref = shift @allele_string;

    my @allele_frequencies = ();
    my $allele_to_frequency = {};
    foreach my $af_key (@af_keys) {
      next if ($data{$af_key} eq '.');
      foreach my $a_to_f (split(',', $data{$af_key})) {
        # frequencies in all_vars are reported per allele: T:0.0038
        my ($a, $f) = split(':', $a_to_f);
        $allele_to_frequency->{$a} = $f;
      }
      my @results = ();
      foreach my $allele (@allele_string) {
      # report frequencies in order of alt alleles
      if (defined $allele_to_frequency->{$allele}) {
          push @results, $allele_to_frequency->{$allele};
        } else {
          push @results, 0;
        }
      }
      push @allele_frequencies, "$af_key=" . join(',', @results);
    }
    if (scalar @allele_frequencies > 0) {
      # for example: AMR=0.0588;SAS=0.0038;EUR=0.0167
      $frequencies->{$id} = join(";", @allele_frequencies);
    }
  }
  return $frequencies;
}

1;
