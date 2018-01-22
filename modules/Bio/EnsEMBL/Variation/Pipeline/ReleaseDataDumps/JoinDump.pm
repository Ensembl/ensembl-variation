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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::JoinDump;

use strict;
use warnings;

use FileHandle;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');


=begin
  Different file joins:
  - join slice pieces to seq_region
  - join seq_regions to complete dump file
  - adjust header lines in dump file:
    - GVF: get correct sequence_regions

=end
=cut

sub run {
  my $self = shift;

  my $species      = $self->param('species');
  my $file_type    = $self->param('file_type');	
  my $config = $self->param('config');
  my $division = $self->param('species_division');

  my $pipeline_dir = $self->data_dir($species);

  my $working_dir = "$pipeline_dir/$file_type/$species/";

  my $mode = $self->param('mode'); 

  if ($mode eq 'join_slice_split') {
    my $files = $self->get_slice_split_files($working_dir, $file_type);  
    $self->join_split_slice_files($working_dir, $files);
    $self->dataflow_output_id({}, 2);
    $self->dataflow_output_id({}, 1);
  } elsif ($mode eq 'final_join') {
    $self->final_join;
  } elsif ($mode eq 'no_join') {
    $self->warning('No final join required');  
  } else {
    die "Unknown mode: $mode in JoinDump";
  }
}

sub final_join {
  my $self = shift;
  my $file_type = $self->param('file_type');
  if ($file_type eq 'gvf') {
    $self->final_join_gvf();
  } else {
    $self->final_join_vcf();
  }
}

sub final_join_gvf {
  my $self = shift;
  my $dir = $self->param('dir');
  my $file_name = $self->param('file_name');
  my $dump_type = $self->param('dump_type');
  my @input_ids = sort @{$self->param('input_ids')};
 
  my $covered_seq_region_ids = $self->get_covered_seq_regions; 
  my $species = $self->param('species');
  my $cdba = $self->get_species_adaptor($species, 'core');
  my $sa = $cdba->get_SliceAdaptor;
  my @sequence_regions = {};
  foreach my $seq_region_id (keys %$covered_seq_region_ids) {
    my $slice = $sa->fetch_by_seq_region_id($seq_region_id);
    push @sequence_regions, {
      'name' => $slice->seq_region_name,
      'start' => $slice->start,
      'end' => $slice->end,
    };
  }

  my $first_file_id = $input_ids[0];
  
  my $fh_join = FileHandle->new("$dir/$file_name.gvf", 'w');

  `gunzip $dir/$dump_type-$first_file_id.gvf.gz`; 

  my $fh = FileHandle->new("$dir/$dump_type-$first_file_id.gvf", 'r');

  # print the header first
  while (<$fh>)  {
    chomp;
    my $line = $_;
    if ($line =~ m/^#/) {
      next if ($line =~ m/^##sequence-region/);
      print $fh_join $line, "\n";
    } 
  }  
  $fh->close();
  `gzip $dir/$dump_type-$first_file_id.gvf`; 

  foreach my $sequence_region (@sequence_regions) {
    print $fh_join join(' ', '##sequence-region', $sequence_region->{name}, $sequence_region->{start}, $sequence_region->{end}), "\n";
  }

  my $id_count = 1; 

  foreach my $file_id (@input_ids) {
    `gunzip $dir/$dump_type-$file_id.gvf.gz`; 
    my $fh = FileHandle->new("$dir/$dump_type-$file_id.gvf", 'r');
    while (<$fh>) {
      chomp;
      my $line = $_;
      next if ($line =~ m/^#/);
      my $gvf_line = get_gvf_line($line, $id_count);
      print $fh_join $gvf_line, "\n";   
      $id_count++;
    }
    $fh->close();
    `rm $dir/$dump_type-$file_id.gvf`;
  }
  $fh_join->close();
}

sub final_join_vcf {
  my $self = shift;
  my $dir = $self->param('dir');
  my $file_name = $self->param('file_name');
  my $dump_type = $self->param('dump_type');
  my @input_ids = sort @{$self->param('input_ids')};
 
  my $first_file_id = shift @input_ids;
  
  my $joined_fh = FileHandle->new("$dir/$file_name.vcf", 'w');

  `gunzip $dir/$dump_type-$first_file_id.vcf.gz`; 
  my $fh = FileHandle->new("$dir/$dump_type-$first_file_id.vcf", 'r');
  while (<$fh>)  {
    chomp;
    print $joined_fh $_, "\n";
  }  
  $fh->close();

  `rm $dir/$dump_type-$first_file_id.vcf`;

  foreach my $file_id (@input_ids) {
    `gunzip $dir/$dump_type-$file_id.vcf.gz`; 
    my $fh = FileHandle->new("$dir/$dump_type-$file_id.vcf", 'r');
    while (<$fh>) {
      chomp;
      my $line = $_;
      next if ($line =~ m/^#/);
      print $joined_fh $line, "\n";
    }
    $fh->close();
    `rm $dir/$dump_type-$file_id.vcf`;
  }
  $joined_fh->close();

  my $vcf_file = "$dir/$file_name.vcf";
  my $cmd = "vcf-sort < $vcf_file | bgzip > $vcf_file.gz";
  $self->run_cmd($cmd); 
  `rm $vcf_file`;
}

sub get_slice_split_files {
  my ($self, $working_dir, $file_type) = @_;
  my $files = {};
  my ($split_slice_range, $file_name);
  opendir(my $dh, $working_dir) or die $!;
  my @dir_content = readdir($dh);
  closedir($dh);
  foreach my $file (@dir_content) {
    next if ($file =~ m/^\./);
    if ($file =~ m/\.$file_type/) {
      $file =~ s/\.$file_type//g;
      my @file_name_components =  split('-', $file);
      if (scalar @file_name_components == 2) {
        $file_name  = shift @file_name_components;
        $split_slice_range = shift @file_name_components;
        my @components = split('_', $split_slice_range);
        if (scalar @components == 3) {
          my ($seq_region_id, $start, $end) = @components; 
          $files->{$file_name}->{$seq_region_id}->{$start} = "$file_name-$seq_region_id\_$start\_$end.$file_type";
        }
      } 
    } # else .err and .out files
  }
  return $files;	
}

sub join_split_slice_files {
  my ($self, $working_dir, $files) = @_;

  my $tmp_dir = $self->param('tmp_dir');

  foreach my $file_type (keys %$files) {
    foreach my $seq_region_id (keys %{$files->{$file_type}}) {
      my $id_count = 1;
      my @start_positions = sort keys %{$files->{$file_type}->{$seq_region_id}};
      my $fh_join = FileHandle->new("$working_dir/$file_type-$seq_region_id.gvf", 'w');
      my $first_start_position = $start_positions[0];
      my $file_name = $files->{$file_type}->{$seq_region_id}->{$first_start_position};
      my $fh = FileHandle->new("$working_dir/$file_name", 'r');
      while (<$fh>)  {
        chomp;
        my $line = $_;
        if ($line =~ m/^#/) {
          next if ($line =~ m/^##sequence-region/);
          print $fh_join $line, "\n";
        }
      }  
      $fh->close();

      foreach my $start_position (@start_positions) {
        my $file_name = $files->{$file_type}->{$seq_region_id}->{$start_position};
        my $fh = FileHandle->new("$working_dir/$file_name", 'r');
        while (<$fh>) {
          chomp;
          my $line = $_;
          next if ($line =~ m/^#/);
          my $gvf_line = get_gvf_line($line, $id_count);
          print $fh_join $gvf_line, "\n";   
          $id_count++;
        }
        $fh->close();
        `gzip $working_dir/$file_name`;
        `mv $working_dir/$file_name.gz $tmp_dir`;
      }
      $fh_join->close();
    }
  }
}

sub get_gvf_line {
  my ($line, $id_count) = @_;
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
  return "$line\t$attributes";
}

sub run_cmd {
  my ($self ,$cmd) = @_;
  if (my $return_value = system($cmd)) {
    $return_value >>= 8;
    die "system($cmd) failed: $return_value";
  }
}

=begin
sub get_covered_seq_regions {
  my $self = shift;
  my $species = $self->param('species');
  my $counts;
  my $vdba = $self->get_species_adaptor($species, 'variation');
  my $dbh = $vdba->dbc->db_handle;
  my $sth = $dbh->prepare(qq{
      SELECT sr.seq_region_id, count(*)
      FROM seq_region sr, variation_feature vf
      WHERE sr.seq_region_id = vf.seq_region_id
      GROUP BY sr.seq_region_id;
      });
  $sth->{'mysql_use_result'} = 1;
  $sth->execute();
  my ($slice_id, $count);
  $sth->bind_columns(\$slice_id, \$count);
  while ($sth->fetch()) {
    if ($count > 0) {
      $counts->{$slice_id} = $count;
    }
  }
  $sth->finish();
  return $counts;
}
=end
=cut

sub get_covered_seq_regions {
  my $self = shift;
  my $species = $self->param('species');
  my $counts;
  my $vdba = $self->get_species_adaptor($species, 'variation');
  my $cdba = $self->get_species_adaptor($species, 'core');
  my $toplevel_seq_region_ids = {};
  my $sa = $cdba->get_SliceAdaptor;
  my $toplevel_slices = $sa->fetch_all('toplevel');
  foreach my $toplevel_slice (@$toplevel_slices) {
    $toplevel_seq_region_ids->{$toplevel_slice->get_seq_region_id} = 1;
  }

  my $dbh = $vdba->dbc->db_handle;
  my $sth = $dbh->prepare(qq{
      SELECT sr.seq_region_id, count(*)
      FROM seq_region sr, variation_feature vf
      WHERE sr.seq_region_id = vf.seq_region_id
      GROUP BY sr.seq_region_id;
      });
  $sth->{'mysql_use_result'} = 1;
  $sth->execute();
  my ($slice_id, $count);
  $sth->bind_columns(\$slice_id, \$count);
  while ($sth->fetch()) {
    if ($count > 0) {
      if ($toplevel_seq_region_ids->{$slice_id}) {
        $counts->{$slice_id} = $count;
      }
    }
  }
  $sth->finish();

  return $counts;
}

sub write_output {
  my $self = shift;
  $self->dataflow_output_id($self->param('input_for_validation'), 1);
  return;
}


1;
