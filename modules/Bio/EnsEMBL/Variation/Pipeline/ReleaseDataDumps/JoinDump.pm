=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::JoinDump;

use strict;
use warnings;

use FileHandle;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

sub run {
  my $self = shift;

  my $species      = $self->param('species');
  my $pipeline_dir = $self->param('pipeline_dir');
  my $file_type    = 'gvf';	

  my $working_dir = "$pipeline_dir/$file_type/$species/";
  my $files = $self->get_files($working_dir, $file_type);
  if ($self->contains_unjoined_files($files)) {		
    $self->join_gvf($working_dir, $files);
  }
  my @input = ();
  foreach my $file_name (keys %$files) {
    my $params = {};
    $params->{file_name}   = $file_name;
    $params->{working_dir} = $working_dir;
    push @input, $params;
  }
  $self->param('input_for_validation', \@input);	
}

sub get_files {
  my $self = shift;
  my $working_dir = shift;
  my $file_type = shift;	
  opendir(DIR, $working_dir) or die $!; 
  my $files = {};
  my ($range, $file_name);
  while (my $file = readdir(DIR)) {
    next if ($file =~ m/^\./);
    if ($file =~ m/\.$file_type/) {
      $file =~ s/\.$file_type//g;
      my @file_name_components =  split('-', $file);
      if (scalar @file_name_components == 2) {
        $file_name  = shift @file_name_components;
        $range = shift @file_name_components;
        $files->{$file_name}->{$range} = 1;
      } else {
        $file_name = join('-', @file_name_components);
        $files->{$file_name}->{1} = 1;
      }
    } # else .err and .out files
  }
  closedir(DIR);
  return $files;	
}

sub contains_unjoined_files {
  my $self = shift;
  my $files = shift;
  foreach my $type (keys %$files) {
    return ((scalar keys %{$files->{$type}}) > 1);
  }
  return 0;	
} 

sub join_gvf {
  my $self  = shift;
  my $working_dir = shift;
  my $files = shift;
  my $species = $self->param('species');
  my $tmp_dir = $self->param('tmp_dir');	

# check there are the same number of sub_files for each dump_type
  my $first_file_name     = (keys %$files)[0];
  my $sub_file_count = scalar keys %{$files->{$first_file_name}};

  foreach my $file_name (keys %$files) {
    my $file_count = scalar keys %{$files->{$file_name}};
    die "file count differs for type: $file_name" unless ($file_count == $sub_file_count);	
  }

# first assemble header: 
# from first file complete header, then only sequence_region information
  foreach my $file_name (keys %$files) {
    my $fh = FileHandle->new("> $working_dir/$file_name.gvf");
# header ---------------------------------------------------------------------
    my @sub_files = keys %{$files->{$file_name}};
    my $first_file_range = shift @sub_files;
    my $tmp_header_fh = FileHandle->new("< $working_dir/$file_name-$first_file_range.gvf");
    while (<$tmp_header_fh>) {
      chomp;
      my $line = $_;
      if ($line =~ m/^#/) {
        print $fh $line, "\n";
      } else {last;}
    }
    $tmp_header_fh->close();

    foreach my $file_range (@sub_files) {
      $self->warning("<$working_dir/$file_name-$file_range.gvf");		
      $tmp_header_fh = FileHandle->new("< $working_dir/$file_name-$file_range.gvf");
      while (<$tmp_header_fh>) {
        chomp;
        my $line = $_;
        last unless ($line =~ m/^#/);
        if ($line =~ m/^##sequence-region/) {
          print $fh $line, "\n";
        }
      }
      $tmp_header_fh->close();
    }
# header ---------------------------------------------------------------------

# body -----------------------------------------------------------------------
# refresh sub_files
    @sub_files = keys %{$files->{$file_name}};
    my $id_count = 1;
    foreach my $file_range (@sub_files) {
      $tmp_header_fh = FileHandle->new("< $working_dir/$file_name-$file_range.gvf");
      while (<$tmp_header_fh>) {
        chomp;
        my $line = $_;
        next if ($line =~ m/^#/);
# do some corrections here:
        my $gvf_line = get_gvf_line($line);
        $gvf_line->{attributes}->{ID} = $id_count;

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
        print $fh $line, "\t", $attributes, "\n";   
        $id_count++;
      }
      $tmp_header_fh->close();
      system("gzip $working_dir/$file_name-$file_range.gvf");
      system("mv $working_dir/$file_name-$file_range.gvf.gz $tmp_dir");
    }
# body -----------------------------------------------------------------------
    $fh->close();
  }
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

sub write_output {
  my $self = shift;
  $self->dataflow_output_id($self->param('input_for_validation'), 1);
  return;
}

1;
