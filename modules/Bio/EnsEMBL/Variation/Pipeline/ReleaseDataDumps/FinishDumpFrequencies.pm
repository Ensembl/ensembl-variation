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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::FinishDumpFrequencies;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

use FileHandle;

sub run {
  my $self = shift;
  my $file_type = $self->param_required('file_type'); 
  if ($file_type eq 'gvf') {
    $self->join_gvf_files();
  } elsif ($file_type eq 'vcf') {
    $self->join_vcf_files();
  } else {
    die("File type ($file_type) is not supported. File type must be gvf or vcf.");
  }
}

sub join_gvf_files {
  my $self = shift;
  my $homo_sapiens_dump_dir = $self->param_required('homo_sapiens_dump_dir');
  my $gvf_dir = "$homo_sapiens_dump_dir/gvf/homo_sapiens/";
  if (! -d $gvf_dir) {
    die("Directory ($gvf_dir) doesn't exist");
  }

  my $gvf_file = "$homo_sapiens_dump_dir/gvf/homo_sapiens/homo_sapiens-chrY.gvf.gz";
  if (! -e $gvf_file) {
    die("File ($gvf_file) doesn't exist.");
  }
  # Create a header file without the ##sequence-region line
  # ##sequence-region lines will be added for all chromosomes
  $self->run_system_command("zgrep ^# $gvf_file > $gvf_dir/gvf_header");
  $self->run_system_command("awk '!/##sequence-region/' $gvf_dir/gvf_header > $gvf_dir/temp");
  $self->run_system_command("mv $gvf_dir/temp $gvf_dir/gvf_header");

  open(my $fh, '>>', "$gvf_dir/gvf_header");
  my $cdba = $self->get_adaptor('homo_sapiens', 'core');
  my $sa = $cdba->get_SliceAdaptor;
  for my $chr (1..22, 'MT', 'X', 'Y') {
    my $slice = $sa->fetch_by_region('chromosome', $chr);  
    print $fh join(' ', '##sequence-region', $slice->seq_region_name, $slice->start, $slice->end), "\n";
  }
  close $fh;

  my $output_file_name = $self->param_required('output_file_name'); 

  my $fh_out = FileHandle->new("$gvf_dir/$output_file_name.gvf", 'w');
  # Start GVF file with header
  my $fh_header = FileHandle->new("$gvf_dir/gvf_header", 'r');
  while (<$fh_header>) {
    chomp;
    print $fh_out "$_\n";
  }
  $fh_header->close;

  my $id = 1;
  for my $chr (1..22, 'MT', 'X', 'Y') {

    my $fh_in = FileHandle->new("$gvf_dir/$output_file_name-$chr.gvf", 'r');

    while (<$fh_in>) {
      chomp;
      next if (/^#/);
      my $line = $_;
      my $gvf_line = get_gvf_line($line);
      $gvf_line->{attributes}->{ID} = $id;
      # GVF files contain a file-wide unique identifier
      # We need to reassign a unique ID value because we are combining several GVF files into one
      $id++;

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
      print $fh_out $line, "\t", $attributes, "\n";
    }
    $fh_in->close;
  }

  $self->run_system_command("rm $gvf_dir/$output_file_name-*.gvf");
  $self->run_system_command("rm $gvf_dir/sorted*");
  $self->run_system_command("rm $gvf_dir/gvf_header");
  $self->run_system_command("gzip $gvf_dir/$output_file_name.gvf");
  $self->compute_checksums_for_directory($gvf_dir);
}

sub join_vcf_files {
  my $self = shift;
  my $homo_sapiens_dump_dir = $self->param_required('homo_sapiens_dump_dir');
  my $vcf_dir = "$homo_sapiens_dump_dir/vcf/homo_sapiens/";
  if (! -d $vcf_dir) {
    die("Directory ($vcf_dir) doesn't exist");
  }

  my $vcf_file = "$homo_sapiens_dump_dir/vcf/homo_sapiens/homo_sapiens-chrY.vcf.gz";
  if (! -e $vcf_file) {
    die("File ($vcf_file) doesn't exist.");
  }

  $self->run_system_command("zgrep ^# $vcf_file > $vcf_dir/vcf_header");
  $self->run_system_command("awk '!/^#CHROM/' $vcf_dir/vcf_header > $vcf_dir/temp");
  $self->run_system_command("mv $vcf_dir/temp $vcf_dir/vcf_header");

  open(my $fh, '>>', "$vcf_dir/vcf_header");
  my $cdba = $self->get_adaptor('homo_sapiens', 'core');
  my $sa = $cdba->get_SliceAdaptor;
  for my $chr (1..22, 'MT', 'X', 'Y') {
    my $slice = $sa->fetch_by_region('chromosome', $chr);  
    my $id = $slice->seq_region_name;
    my $length = $slice->end;
    print $fh "##contig=<ID=$id,length=$length>\n";
  }
  my $af_keys_descriptions = $self->param_required('af_keys_descriptions');
  foreach my $af_key (sort keys %$af_keys_descriptions) {
    my $description = $af_keys_descriptions->{$af_key};
    print $fh "##INFO=<ID=$af_key,Number=A,Type=Float,Description=\"$description\">\n";
  } 
  print $fh join("\t", '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'), "\n";
  close $fh;

  my $output_file_name = $self->param_required('output_file_name'); 
  my @files = ();
  for my $chr (1..22, 'MT', 'X', 'Y') {
    push @files, "$vcf_dir/$output_file_name-$chr.vcf";
  }
  my $files = join(' ', @files);
  $self->run_system_command("cat $vcf_dir/vcf_header $files > $vcf_dir/$output_file_name.vcf");
  $self->run_system_command("rm $vcf_dir/$output_file_name-*.vcf");
  $self->run_system_command("rm $vcf_dir/vcf_header");
  $self->run_system_command("vcf-sort < $vcf_dir/$output_file_name.vcf | bgzip > $vcf_dir/$output_file_name.vcf.gz");
  $self->run_system_command("rm $vcf_dir/$output_file_name.vcf");
  $self->run_system_command("tabix -C -p vcf $vcf_dir/$output_file_name.vcf.gz");   
  $self->compute_checksums_for_directory($vcf_dir);
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

1;
