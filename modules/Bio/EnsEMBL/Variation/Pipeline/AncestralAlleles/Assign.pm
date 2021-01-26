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

package Bio::EnsEMBL::Variation::Pipeline::AncestralAlleles::Assign;

use strict;
use warnings;

use FileHandle;
use Bio::EnsEMBL::Registry;
use File::Path qw(make_path);
use Data::Dumper;
use Bio::DB::Fasta;
use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub run {
  my $self = shift;
  my $species_name = $self->param('species_name');
  my $batch = $self->param('batch');
  my $id = $batch->{id};
  my $start = $batch->{start};
  my $end = $batch->{end};
  my $ancestral_file = $self->param('ancestral_file');

  $self->setup();
  $self->dump_variation_features();
  $self->assign_ancestral_alleles();
  $self->cleanup();
}

sub setup {
  my $self = shift;
  my $batch = $self->param('batch');
  my $id = $batch->{id};
  my $compara_dir = $self->param('compara_dir'); 
  my $species_dir = $self->param('species_dir');
  my $ancestral_file_targz = $self->param('ancestral_file');
  my $ancestral_file_notargz = $ancestral_file_targz;
  $ancestral_file_notargz =~ s/\.tar\.gz|\.tar\.bz2//;
  # create temp directory for storing and uncompressing ancestral fasta file 
  unless (-d "$species_dir/fasta_$id") {
    my $err;
    make_path("$species_dir/fasta_$id", {error => \$err});
    die "make_path failed: ".Dumper($err) if $err && @$err;
  }
  # cp  ancestral fasta file from compara dir to temp directory
  $self->run_cmd("cp $compara_dir/$ancestral_file_targz $species_dir/fasta_$id");
  # uncompress ancestral fasta file
  if ($ancestral_file_targz =~ /\.tar\.gz/) {
    $self->run_cmd("tar xzf $species_dir/fasta_$id/$ancestral_file_targz -C $species_dir/fasta_$id/");
  } elsif ($ancestral_file_targz =~ /\.tar\.bz2/) {
    $self->run_cmd("tar xjf $species_dir/fasta_$id/$ancestral_file_targz -C $species_dir/fasta_$id/");
  } else {
    die "Could not recognise tar file extension. Supported extensions are tar.gz and tar.bz2";
  }
  # rm bed files
  $self->run_cmd("rm $species_dir/fasta_$id/$ancestral_file_notargz/*.bed");
  my $fasta_files_dir = "$species_dir/fasta_$id/$ancestral_file_notargz"; 
  $self->param('fasta_files_dir', $fasta_files_dir);
  my $fasta_files_dir_root = "$species_dir/fasta_$id/"; 
  $self->param('fasta_files_dir_root', $fasta_files_dir_root);
}

sub dump_variation_features {
  my $self = shift;
  my $species_name = $self->param('species_name');
  my $species_dir = $self->param('species_dir');
  my $batch = $self->param('batch'); 
  my ($start, $end) = ($batch->{start}, $batch->{end});
  my $vf_file = "$species_dir/vf_$start\_$end";
  my $aa_file = "$species_dir/aa_$start\_$end";
  $self->param('vf_file', $vf_file);
  $self->param('aa_file', $aa_file);
  if (! -e "$vf_file") {
    my $fh = FileHandle->new($vf_file, 'w');
    $self->param('species', $species_name);
    my $vdba = $self->get_species_adaptor('variation');
    my $dbc = $vdba->dbc;
    my $sth = $dbc->prepare(qq{
      SELECT vf.variation_feature_id, sr.name, vf.seq_region_start, vf.seq_region_end, vf.seq_region_strand, vf.ancestral_allele
      FROM variation_feature vf, seq_region sr
      WHERE vf.seq_region_id = sr.seq_region_id
      AND vf.variation_feature_id >= ?
      AND vf.variation_feature_id <= ?
      ORDER BY vf.variation_feature_id;}, {mysql_use_result => 1}
    ) or die $dbc->errstr;
    $sth->execute($start, $end) or die $dbc->errstr;
    my ($vf_id, $seq_region, $start, $end, $strand, $aa);
    $sth->bind_columns(\($vf_id, $seq_region, $start, $end, $strand, $aa));
    while ($sth->fetch) {
      $aa ||= 'NULL';
      print $fh join("\t", ($vf_id, $seq_region, $start, $end, $strand, $aa)), "\n";
    }
    $sth->finish();
    $fh->close();
  }
}

sub assign_ancestral_alleles {
  my $self = shift;
  my $fasta_files_dir = $self->param('fasta_files_dir');
  my $db = Bio::DB::Fasta->new($fasta_files_dir, -reindex => 1);
  my @sequence_ids = $db->get_all_ids;
  my %sequence_id_2_chr_number;
  foreach my $sequence_id (@sequence_ids) {
    my @split = split(/:/, $sequence_id);
    $sequence_id_2_chr_number{$split[2]} = $sequence_id;
  }

  my $vf_file = $self->param('vf_file');
  my $aa_file = $self->param('aa_file');
  my $input = FileHandle->new($vf_file, 'r');
  my $fh    = FileHandle->new("$aa_file.out", 'w');
  my $err   = FileHandle->new("$aa_file.err", 'w');
  while (<$input>) {
    chomp;
    my ($vf_id, $chrom, $start, $end, $strand, $old_aa) = split /\t/;
    my $set_aa_to_null = 1;
    if ($start > $end) {
      print $err "Insertion $vf_id\n";
    }
    elsif (($end - $start) >= 50) {
      print $err "Longer than 50bp $vf_id\n";
    } else {
      my $chrom_name = $sequence_id_2_chr_number{$chrom};
      if (!($chrom_name && $start && $end)) {
        print $err "Incomplete coords $vf_id\n";
      } else {
        my $aa = $db->seq("$chrom_name:$start,$end");
        if ($aa && $aa =~ m/^[ACGTacgt]+$/) {
          $aa = uc $aa;
          $set_aa_to_null = 0;
          print $fh "UPDATE variation_feature SET ancestral_allele='$aa' WHERE variation_feature_id=$vf_id;\n" if ($aa ne $old_aa);
        } else {
          print $err "No AA in fasta $vf_id\n";
        }
      }
    }
    if ($set_aa_to_null && $old_aa ne 'NULL') {
      print $fh "UPDATE variation_feature SET ancestral_allele=\\N WHERE variation_feature_id=$vf_id;\n" if ($old_aa ne 'NULL');
    }
  }

  $input->close();
  $fh->close();
  $err->close();
}

sub cleanup {
  my $self = shift;
  my $fasta_files_dir_root = $self->param('fasta_files_dir_root');
  $self->run_cmd("rm -rf $fasta_files_dir_root");
}

1;
