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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::FinishDumpPopulation;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

use FileHandle;

sub run {
    my $self = shift;
    my $file_type = $self->param('file_type');
    if ($file_type eq 'gvf') {
        # remove allele_string and variation_id attributes from file    
        # gzip file
        $self->clean_up_gvf_file;
    }
    if ($file_type eq 'vcf') {
        # sort, bgzip, validate
        $self->clean_up_vcf_files;
    }
}

sub clean_up_vcf_files {
  my $self = shift;
  my $pipeline_dir = $self->param('pipeline_dir');
  my $tmp_dir = $self->param('tmp_dir');
  my $human_vcf_dir = "$pipeline_dir/vcf/Homo_sapiens/";
  die "$human_vcf_dir is not a directory" unless (-d $human_vcf_dir);

  opendir(my $dh, "$pipeline_dir/vcf/Homo_sapiens") or die $!;
  my @dir_content = readdir($dh);
  closedir($dh);
  foreach my $file (@dir_content) {
    if ($file =~ m/\.vcf$/) {
      my $vcf_file = "$pipeline_dir/vcf/Homo_sapiens/$file";
      $self->run_cmd("vcf-sort < $vcf_file | bgzip > $vcf_file.gz");
      $self->run_cmd("rm $vcf_file");
    }
    if ($file =~ m/\.out$/ || $file =~ m/\.err$/) {
      $self->run_cmd("rm $pipeline_dir/vcf/Homo_sapiens/$file");
    }
  } 
}

sub clean_up_gvf_file {
  my $self = shift;
  my $pipeline_dir = $self->param('pipeline_dir');
  my $human_gvf_file = "$pipeline_dir/gvf/Homo_sapiens/Homo_sapiens.gvf";
  my $human_gvf_file_before_clean_up = "$pipeline_dir/gvf/Homo_sapiens/Homo_sapiens_before_clean_up.gvf";
  die "Couldn't find $human_gvf_file" unless (-f $human_gvf_file);
  $self->run_cmd("mv $human_gvf_file $human_gvf_file_before_clean_up"); 
  die "Couldn't find $human_gvf_file_before_clean_up" unless (-f $human_gvf_file_before_clean_up);

  my $fh_before = FileHandle->new($human_gvf_file_before_clean_up, 'r');
  my $fh_after = FileHandle->new($human_gvf_file, 'w');

  while (<$fh_before>) {
    chomp;
    if (/^#/) {
      print $fh_after $_, "\n";
    } else {
    my $line = $_;
    my $gvf_line = get_gvf_line($line);
    delete $gvf_line->{attributes}->{variation_id};
    delete $gvf_line->{attributes}->{allele_string};
    $line = join("\t", map {$gvf_line->{$_}} (
      'seq_id',
      'source',
      'type',
      'start',
      'end',
      'score',
      'strand',
      'phase'));
    my $attributes = join(";", map{"$_=$gvf_line->{attributes}->{$_}"} keys %{$gvf_line->{attributes}});
    print $fh_after $line, "\t", $attributes, "\n";
    }
  }

  $fh_before->close();
  $fh_after->close();

  $self->run_cmd("gzip $human_gvf_file");  
  $self->run_cmd("rm $human_gvf_file_before_clean_up");  
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

sub run_cmd {
  my $self = shift;
  my $cmd = shift;
  if (my $return_value = system($cmd)) {
    $return_value >>= 8;
    die "system($cmd) failed: $return_value";
  }
}

1;
