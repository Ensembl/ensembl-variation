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
package Bio::EnsEMBL::Variation::Pipeline::SpliceAI::SplitFiles;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Variation::Pipeline::SpliceAI::BaseSpliceAI');

use FileHandle;
use Bio::EnsEMBL::IO::Parser::BedTabix;
use Bio::EnsEMBL::IO::Parser::VCF4Tabix;

sub run {
  my $self = shift;
  $self->set_chr_from_filename();
  # check = 0 (don't check transcripts); check = 1 (check transcripts)
  my $check = $self->param_required('check_transcripts');

  # check new Mane transcripts
  if($check){
    $self->get_new_transcripts();
    $self->check_split_vcf_file();
  }
  $self->split_vcf_file();
}

sub set_chr_from_filename {
  my $self = shift;
  my $vcf_file = $self->param_required('vcf_file');
  $vcf_file =~ /.*_chr(.*)\.vcf$/;
  my $chr = $1;
  if (!$chr) {
    die("Could not get chromosome name from file name ($vcf_file).");
  }
  $self->param('chr', $chr);
}

sub split_vcf_file {
  my $self = shift;
  my $check = $self->param_required('check_transcripts');
  my $vcf_file = $self->param_required('vcf_file'); # input vcf which is inside the input_directory
  my $split_vcf_no_header_dir = $self->param_required('split_vcf_no_header_dir');
  my $split_vcf_input_dir = $self->param_required('split_vcf_input_dir');
  my $input_dir = $check ? $self->param('input_dir_subset') : $self->param_required('input_directory'); # the input directory is 'input_dir_subset' if check_transcripts = 1 (the input file will only contain the regions for the new transcripts); is the normal 'input_directory' if check_transcripts = 0
  my $split_vcf_output_dir = $self->param_required('split_vcf_output_dir');
  my $step_size = $self->param_required('step_size');

  my $vcf_file_path = $input_dir . '/' . $vcf_file;

  if (! -d $input_dir) {
    die("Directory ($input_dir) doesn't exist");
  }

  if (! -e $vcf_file_path) {
    die("File ($vcf_file_path) doesn't exist.");
  }

  my $chr = $self->param('chr');

  my $split_vcf_no_header_dir_chr = $split_vcf_no_header_dir.'/chr'.$chr; # Create chromosome directory inside 'split_vcf_no_header_dir'
  my $new_file = $split_vcf_no_header_dir_chr.'/all_snps_ensembl_chr'.$chr.'.';

  $self->create_dir($split_vcf_no_header_dir_chr);
  $self->run_system_command("split -l $step_size -d --additional-suffix=.vcf $vcf_file_path $new_file");

  # Files splitted by number of lines (from input)
  # These files contain vcf header
  # Files that are going to be used as input for SpliceAI
  my $split_vcf_input_dir_chr = $split_vcf_input_dir.'/chr'.$chr;

  $self->create_dir($split_vcf_input_dir_chr);

  # Read files from /main_dir/split_vcf (missing header) and write new files to /main_dir/split_vcf_input (ready to be used as input for SpliceAI)
  opendir(my $read_no_header_dir, $split_vcf_no_header_dir_chr) or die $!;

  while(my $split_vcf_no_header_file = readdir($read_no_header_dir)) {
    next if ($split_vcf_no_header_file =~ m/^\./);

    open(my $write, '>', $split_vcf_input_dir_chr . '/' . $split_vcf_no_header_file) or die $!;
    my $header_line = "##fileformat=VCFv4.2\n##contig=<ID=1,length=248956422>\n##contig=<ID=2,length=242193529>\n".
    "##contig=<ID=3,length=198295559>\n##contig=<ID=4,length=190214555>\n##contig=<ID=5,length=181538259>\n".
    "##contig=<ID=6,length=170805979>\n##contig=<ID=7,length=159345973>\n##contig=<ID=8,length=145138636>\n".
    "##contig=<ID=9,length=138394717>\n##contig=<ID=10,length=133797422>\n##contig=<ID=11,length=135086622>\n".
    "##contig=<ID=12,length=133275309>\n##contig=<ID=13,length=114364328>\n##contig=<ID=14,length=107043718>\n".
    "##contig=<ID=15,length=101991189>\n##contig=<ID=16,length=90338345>\n##contig=<ID=17,length=83257441>\n".
    "##contig=<ID=18,length=80373285>\n##contig=<ID=19,length=58617616>\n##contig=<ID=20,length=64444167>\n##contig=<ID=21,length=46709983>\n".
    "##contig=<ID=22,length=50818468>\n##contig=<ID=X,length=156040895>\n##contig=<ID=Y,length=57227415>\n##contig=<ID=MT,length=16569>\n".
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    print $write $header_line;

    open(my $fh, '<:encoding(UTF-8)', $split_vcf_no_header_dir_chr . '/' . $split_vcf_no_header_file)
      or die "Could not open file '$split_vcf_no_header_dir_chr/$split_vcf_no_header_file' $!";

    while (my $row = <$fh>) {
      chomp $row;
      next if($row =~ /^#/);

      my @columns = split /\t/, $row;

      my @alt_splitted = split /,/, $columns[4];
      foreach my $x (@alt_splitted) {
        print $write $columns[0] . "\t" . $columns[1] . "\t" . $columns[2] . "\t" . $columns[3] . "\t" . $x . "\t.\t.\t.\n";
      }
    }
    close($write);
    close($fh);
  }
  close($read_no_header_dir);

  # Create output directory
  my $split_vcf_output_dir_chr = $split_vcf_output_dir.'/chr'.$chr;
  $self->create_dir($split_vcf_output_dir_chr);

  my $output_vcf_files_dir = $split_vcf_output_dir_chr.'/vcf_files';
  $self->create_dir($output_vcf_files_dir);
}

# First checks which variants need to be used and then splits main file using only those variants
sub check_split_vcf_file {
  my $self = shift;
  my $main_dir = $self->param_required('main_dir');
  my $vcf_file = $self->param_required('vcf_file');
  my $input_dir = $self->param_required('input_directory');
  my $step_size = $self->param_required('step_size');
  my $transcripts = $self->param('transcripts');
  my $chr = $self->param('chr');

  my $vcf_file_path = $input_dir . '/' . $vcf_file;

  if (! -d $input_dir) {
    die("Directory ($input_dir) doesn't exist");
  }

  if (! -e $vcf_file_path) {
    die("File ($vcf_file_path) doesn't exist.");
  }

  # prepare input vcf file for tabix
  $self->run_system_command("bgzip $vcf_file_path");
  $self->run_system_command("tabix -p vcf $vcf_file_path.gz");

  # Create directory to store the input file only with the variants of interest (overlapping new transcripts)
  my $vcf_file_path_subset = $main_dir . '/input_vcf_files_subset';
  if (! -e $vcf_file_path_subset) {
    $self->create_dir($vcf_file_path_subset);
  }

  $self->param('input_dir_subset', $vcf_file_path_subset);

  open(my $write, '>', $vcf_file_path_subset . '/' . $vcf_file) or die $!;

  my $positions_of_interest = $transcripts->{$chr};
  foreach my $position (@$positions_of_interest) {
    my @positions = split /-/, $position;
    my $transcript_start = $positions[0];
    my $transcript_end = $positions[1];

    my $pos_string = sprintf("%s:%i-%i", $chr, $transcript_start, $transcript_end);

    (open my $read, "tabix  " . $input_dir . "/" . $vcf_file . ".gz $pos_string |")
      || die "Failed to read from input vcf file " . $input_dir . "/" . $vcf_file . ".gz \n" ;

    while (my $row = <$read>) {
      chomp $row;
      next if($row =~ /^#/);

      my @columns = split /\t/, $row;
      my @alt_splitted = split /,/, $columns[4];
      foreach my $x (@alt_splitted) {
        print $write $columns[0] . "\t" . $columns[1] . "\t" . $columns[2] . "\t" . $columns[3] . "\t" . $x . "\t.\t.\t.\n";
      }
    }
    close($read);
  }
  close($write);

  # Sort new vcf file
  my $vcf_file_subset = $vcf_file_path_subset . '/' . $vcf_file;
  my $vcf_file_subset_sorted = $vcf_file_path_subset . '/sorted_' . $vcf_file;
  $self->run_system_command("sort -t \$'\t' -k1,1 -k2,2n $vcf_file_subset > $vcf_file_subset_sorted");
  $self->run_system_command("mv $vcf_file_subset_sorted $vcf_file_subset");
}

# Check if there are new MANE transcripts since last release
sub get_new_transcripts {
  my $self = shift;

  my %new_transcripts;

  my $registry = 'Bio::EnsEMBL::Registry';
  my $registry_file = $self->param('registry');
  $registry->load_all($registry_file);
  # Only for human
  my $cdba =  $registry->get_DBAdaptor('Homo_sapiens', 'core');
  my $dbh = $cdba->dbc->db_handle;

  # Get new transcripts which need to have the scores re-calculated
  my $sth = $dbh->prepare(qq{ SELECT s.name,t.stable_id,t.version,g.seq_region_start,g.seq_region_end FROM transcript t
                              JOIN transcript_attrib ta ON t.transcript_id = ta.transcript_id
                              JOIN attrib_type atr ON ta.attrib_type_id = atr.attrib_type_id
                              JOIN seq_region s ON t.seq_region_id = s.seq_region_id
                              JOIN gene g ON g.gene_id = t.gene_id
                              WHERE t.stable_id like 'ENST%' and t.biotype = 'protein_coding' and atr.code = 'MANE_Select' and t.modified_date >= DATE_SUB(NOW(), INTERVAL 4 MONTH) });

  $sth->execute();
  while (my $row = $sth->fetchrow_arrayref) {
    my $chr = $row->[0];
    my $transcript_id = $row->[1];
    my $transcript_version = $row->[2];
    my $gene_start = $row->[3];
    my $gene_end = $row->[4];

    if(!$new_transcripts{$chr}) {
      my @positions;
      push @positions, $gene_start.'-'.$gene_end;
      $new_transcripts{$chr} = \@positions;
    }
    else {
      push @{$new_transcripts{$chr}}, $gene_start.'-'.$gene_end;
    }
  }
  $sth->finish();

  $self->param('transcripts', \%new_transcripts);
}

sub count_lines {
  my $self = shift;
  my $filename = shift;

  
}

1;
