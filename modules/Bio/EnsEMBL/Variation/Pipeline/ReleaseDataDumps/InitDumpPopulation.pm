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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::InitDumpPopulation;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

use File::Path qw(make_path);
use FileHandle;

sub fetch_input {
  my $self = shift;

  my $populations = {
    'HAPMAP' => {
      'CSHL-HAPMAP:HAPMAP-ASW' => 'ASW',
      'CSHL-HAPMAP:HAPMAP-CHB' => 'CHB',
      'CSHL-HAPMAP:HAPMAP-CHD' => 'CHD',
      'CSHL-HAPMAP:HAPMAP-GIH' => 'GIH',
      'CSHL-HAPMAP:HAPMAP-LWK' => 'LWK',
      'CSHL-HAPMAP:HAPMAP-MEX' => 'MEX',
      'CSHL-HAPMAP:HAPMAP-MKK' => 'MKK',
      'CSHL-HAPMAP:HAPMAP-TSI' => 'TSI',
      'CSHL-HAPMAP:HapMap-CEU' => 'CEU',
      'CSHL-HAPMAP:HapMap-HCB' => 'HCB',
      'CSHL-HAPMAP:HapMap-JPT' => 'JPT',
      'CSHL-HAPMAP:HapMap-YRI' => 'YRI',
    },
    'ESP' => {
      'ESP6500:African_American'  => 'AA',
      'ESP6500:European_American' => 'EA',
    },
  };

  $self->param('populations', $populations);

  my $pipeline_dir = $self->param('pipeline_dir');
  my $file_type = $self->param('file_type');

  if ($file_type eq 'gvf') {    

    my $prefetched_frequencies_dir = $self->param('prefetched_frequencies');
    die "$prefetched_frequencies_dir doesn't exist" unless (-d $prefetched_frequencies_dir);

  }

}

sub run {
  my $self = shift;
  my $file_type = $self->param('file_type');
  if ($file_type eq 'gvf') {
    my $input = $self->get_input_dump_gvf;
    $self->param('input_parameters', $input);
  } elsif ($file_type eq 'vcf') {
    my $input = $self->get_input_dump_vcf;       
    $self->param('input_parameters', $input);
  } else {
    die "File type must be gvf or vcf. File type ($file_type) is not valid.";
  }
}

sub get_input_dump_vcf {
  my $self = shift; 
  my $populations = $self->param('populations');

  my $file_type       = 'vcf';
  my $script_dir      = $self->param('script_dir');
  my $script          = '/misc/release/gvf2vcf.pl';
  my $output_dir      = $self->param('pipeline_dir');
  my $connection_args = '--registry ' . $self->param('registry_file');
  my $species = 'homo_sapiens';
  my @input = ();

  foreach my $group (keys %$populations) { 
    foreach my $population_name (keys %{$populations->{$group}}) {
      my $params = {};
      my $short_name = $populations->{$group}->{$population_name};
      my $file_name = $population_name;
      $file_name =~ s/:/-/;
      my $err = "$output_dir/$file_type/$species/$file_name.err";
      my $out = "$output_dir/$file_type/$species/$file_name.out";
      $params->{'species'}          = $species;
      $params->{'script'}           = "$script_dir/$script";
      $params->{'connection_args'}  = $connection_args;
      $params->{'script_args'}      = "--population $population_name --evidence --ancestral_allele --clinical_significance" ;
      $params->{'err'}              = $err;
      $params->{'out'}              = $out;

      my $gvf_file_name = "$file_name.gvf.gz";
      my $vcf_file_name = "$file_name.vcf";
      my $gvf_file = "$output_dir/gvf/$species/$gvf_file_name";
      die "GVF file $gvf_file not found" unless(-e $gvf_file);
      my $vcf_file = "$output_dir/vcf/$species/$vcf_file_name";
      $params->{'gvf_file'} = "--gvf_file $gvf_file";
      $params->{'vcf_file'} = "--vcf_file $vcf_file";
      push @input, $params;
    }
  }
  return \@input;
}


sub get_gvf_line {
  my $line = shift;
  my $gvf_line = {};
  my @header_names = qw/seq_id source type start end score strand phase/;
  my @header_values = split(/\t/, $line);
  my $attrib = pop @header_values;

  for my $i (0 .. $#header_names) {
    $gvf_line->{$header_names[$i]} = $header_values[$i];
  }

  my @attributes = split(';', $attrib);
  foreach my $attribute (@attributes) {
    my ($key, $value) = split('=', $attribute);
    if ($value) {
      $gvf_line->{attributes}->{$key} = $value;
    }
  }
  return $gvf_line;
}

sub get_input_dump_gvf {
  my $self = shift;
  my @input = ();

  my $populations = $self->param('populations');
  my $pipeline_dir = $self->param('pipeline_dir');
  my $gvf_files_dir = "$pipeline_dir/gvf/homo_sapiens/";
  my $file_type = $self->param('file_type');
  foreach my $group (keys %$populations) {
    foreach my $population_name (keys %{$populations->{$group}}) {
      my $params = {};
      my $short_name = $populations->{$group}->{$population_name};
      my $file_name = $population_name;
      $file_name =~ s/:/-/;
      $params->{population_name} = $population_name;
      $params->{file_name} = $population_name;
      $params->{population_gvf_file} = "$gvf_files_dir/$file_name.gvf";
      $params->{short_name} = $short_name;
      $params->{file_type} = $file_type;
      $params->{group} = $group;
      $params->{gvf_files_dir} = $gvf_files_dir;
      push @input, $params;
    }
  }
  return \@input;
}

sub write_output {
  my $self = shift;
  my $input_parameters = $self->param('input_parameters');
  $self->dataflow_output_id($input_parameters, 2);
  return;
}

1;
