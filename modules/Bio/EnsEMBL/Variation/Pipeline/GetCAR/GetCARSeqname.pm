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

=head1 NAME

Bio::EnsEMBL::Variation::Pipeline::GetCAR::DumpRegionHGVS

=head1 DESCRIPTION

For a given seqname gets all the CAR ids

=cut

package Bio::EnsEMBL::Variation::Pipeline::GetCAR::GetCARSeqname;

use strict;
use warnings;
use File::Basename;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $clingen_url = 'http://reg.genome.network/alleles?file=hgvs&fields=none+@id';

sub fetch_input {
  my $self = shift;
  my $hgvs_dir = $self->param_required('hgvs_dir');
  my $car_lu_dir = $self->param_required('car_lu_dir');
  my $seq_name = $self->param_required('seq_name');

  if (! $car_lu_dir) {
    die("car_lu_dir is not defined");
  }
  if (! $hgvs_dir) {
    die("hgvs_dir is not defined");
  }
  if (! $seq_name) {
    die("seq_name is not defined");
  }
  if ($seq_name !~ /^(\d{1,2}|MT|X|Y)$/) {
    die("seq_name ($seq_name) is invalid");
  }
  $self->warning("Get CAR for ($seq_name)");
}

sub run {
  my $self = shift;

  my $url = $clingen_url;

  my $hgvs_dir = $self->param_required('hgvs_dir');
  my $car_lu_dir = $self->param_required('car_lu_dir');
  my $seq_name = $self->param_required('seq_name');

  # Get a list of files that start hgvs-<seqname>-start-end.tab
  # The files are produced by HGVS dump pipeline. 
  # Column 1 contains the variation_id 
  # Column 6 contains HGVSg
  # Call the lookup file hgvs-<seqname>-car-lu.txt
  # Call the JSON file hgvs-<seqname>-car-lu.json
  # Make the REST API call
  my @hgvs_files;
  opendir(HGVS_DIR, $hgvs_dir) or die("Unable to open $hgvs_dir: $!");
  while (my $hgvs_file = readdir(HGVS_DIR))  {
    next if ($hgvs_file =~ /^\./);
    next if ($hgvs_file !~ /hgvs-${seq_name}-\d{1,}-\d{1,}.tab/);
    if (! -s "$hgvs_dir/$hgvs_file") {
        $self->warning("hgvs_file empty ($hgvs_file)");
        next;
    }
    push @hgvs_files, $hgvs_file;
  }
  close(HGVS_DIR);
  $self->warning(scalar @hgvs_files . ' files to process');
  return if (! @hgvs_files);

  for my $hgvs_file (@hgvs_files) {
    # Set up the output files
    my($filename, $dirs, $suffix) = fileparse($hgvs_file, qr/\.[^.]*/);

    my $lu_file = $filename . '-car-lu' . '.txt';
    my $json_file = $filename . '-car-lu' . '.json';

    # Getting the HGVS value from the full file
    my $cmd = "cut -f6 $hgvs_dir/$hgvs_file > $car_lu_dir/$lu_file";

    my $ret = system($cmd);
    if ($ret) {
        die("Error in getting HGVS from $hgvs_file");
    }
    # Make a REST call
    # Do a sleep between the calls
    # Output file is based on the name of the input_file
    my $input_file = join('/', $car_lu_dir, $lu_file);
    my $out_file = join('/', $car_lu_dir, $json_file);
    my $curl_cmd = ['curl', '-X', 'POST', "\"$url\"",
                    '-H', '"Content-Type: application/json"',
                    '-w',  '"%{http_code}"',
                    '--data-binary', "\@${input_file}",
                    '-o', $out_file, '-s', '-S'];

    $cmd = join(" ",@$curl_cmd,"2>&1");
    my $res = `$cmd`;
    my $exit_code = $?;
    if (($exit_code != 0) || $res !~ /200/) {
      die("Unable to fetch ($hgvs_file) exit_code ($exit_code) res ($res)");
    }
    sleep(30);
  }
}

sub write_output {
  my $self = shift;
}

1;
