#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

#directly use this script if generate from small genome like zebrafish

use strict;
#use DBH;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Data::Dumper;
use Bio::SeqIO;
use FindBin qw( $Bin );
use Getopt::Long;
use ImportUtils qw(dumpSQL debug create_and_load load);

# try to use named options here and write a sub usage() function
# eg -host -user -pass -port -snp_dbname -core_dbname etc
# optional chromosome name or genomic sequence file
# optional more than one genomic sequence file
# optional a directory or sequence files (for unknown placing)


our ($species,$output_dir,$seq_region_id,$TMP_DIR,$TMP_FILE, $chr_name,$read_file,$generate_input_seq,$read_flank,$registry_file);

GetOptions('species=s'          => \$species,
	   'output_dir=s'       => \$output_dir,
	   'chr_name=s'         => \$chr_name,
	   'seq_region_id=i'    => \$seq_region_id,
	   'tmpdir=s'           => \$ImportUtils::TMP_DIR,
	   'tmpfile=s'          => \$ImportUtils::TMP_FILE,
	   'read_file=s'        => \$read_file,
           'ensembl_registry=s' => \$registry_file,
	   'generate_input_seq' => \$generate_input_seq,
	   'read_flank'         => \$read_flank, ###this is to mapping read_coverage feature
	  );
$registry_file ||= $Bin . "/ensembl.registry";
$output_dir    ||= $ENV{PWD};

-f $registry_file or usage( "Non existent registry file $registry_file" );
-d $output_dir    or usage( "Non existent output dir $output_dir" );

usage('-species argument is required') if(!$species);
usage('-generate_input_seq and/or -read_flank is required' ) 
    if( !$generate_input_seq and !$read_flank );

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');

my $dbCore = $cdb->dbc->db_handle;
my $dbVar = $vdb->dbc->db_handle;

generate_input_seq($cdb, $vdb, $chr_name) if $generate_input_seq;
get_read_flank_seq($cdb, $vdb, $read_file) if $read_flank;

sub generate_input_seq {

  my $cdb = shift;
  my $vdb = shift;
  my $chr_name = shift;

  my $slice_adaptor = $cdb->get_SliceAdaptor();

  my (%variation_ids);

  my $file_size = 5000;
  my $file_count=1;
  my $i = 0;

  if ($chr_name or $seq_region_id) {
    # Only dump sequences corresponding to chr_name or seq_region_id
    if ($chr_name and ! $seq_region_id) {
      # Get the seq_region_id corresponding to the chr_name
      my $seq_region_id_ref = $dbCore->selectall_arrayref(qq{
SELECT seq_region_id 
FROM   seq_region 
WHERE  name = "$chr_name"});
      $seq_region_id = $seq_region_id_ref->[0][0];
    } elsif (! $chr_name and $seq_region_id) {
      # Get the seq_region_name corresponding to the ID
      my $seq_region_name_ref = $dbCore->selectall_arrayref(qq{
SELECT name 
FROM   seq_region
WHERE  seq_region_id = $seq_region_id});
      $chr_name = $seq_region_name_ref->[0][0];
    }
    print "seq_region_name is $chr_name\n";
    print "seq_region_id is $seq_region_id\n";

    if ($chr_name) {
      if ($chr_name =~ /NT|^A|^B|^C|^ZV7/i) {
	$chr_name = "NT";
      }
      if (! -e "$output_dir/$chr_name") {
        # Create per-chromosome output dir
	mkdir "$output_dir/$chr_name" or die "can't make dir $chr_name: $!";
      }
      $output_dir = "$output_dir/$chr_name";
      print "output_dir is $output_dir\n";
      $ImportUtils::TMP_DIR = $output_dir;
      $TMP_DIR = $ImportUtils::TMP_DIR;
      print "tmp_dir is $TMP_DIR\n";
    }
    print "seq_region_id is $seq_region_id\n";
    
    if (! -e "$TMP_DIR/$TMP_FILE" or -z "$TMP_DIR/$TMP_FILE") {
      # Only dump if tmp file does not exist, only dump flanking sequence that belong to $chr
      dumpSQL($dbVar, qq{
SELECT vf.variation_name,
       vf.variation_id,
       vf.seq_region_id,
       vf.seq_region_strand, 
       null, null,
       IF (vf.seq_region_strand =1, 
           vf.seq_region_start-100, 
           vf.seq_region_end+1),
       IF (vf.seq_region_strand =1, 
           vf.seq_region_start-1, 
           vf.seq_region_end+100),
       IF (vf.seq_region_strand =1, 
           vf.seq_region_end+1, 
           vf.seq_region_start-100),
       IF (vf.seq_region_strand =1, 
           vf.seq_region_end+100, 
           vf.seq_region_start-1)
FROM   variation_feature vf,
       yuan_enssnp_mapping_55.failed_mapping_vf_Watson_SARA m
WHERE  vf.seq_region_id = $seq_region_id
AND    vf.variation_id=m.variation_id
#AND    vf.variation_name not like "NT%" ##this is only for mouse
                           }
	     );
    }
  }
  else {
    if (! -e "$TMP_DIR/$TMP_FILE" or -z "$TMP_DIR/$TMP_FILE") {
      # Only dump if tmp file does not exist, or is compressed
      # No chr_name or seq_region_id specified. dump flanking sequence from paralle_processed flanking_sequence table
      #only can dump flanking seq which is mapped
      my $sql = qq(
SELECT vf.variation_name,
       f.variation_id,
       f.seq_region_id,
       f.seq_region_strand,
       f.up_seq,f.down_seq,
       f.up_seq_region_start, 
       f.up_seq_region_end,
       f.down_seq_region_start, 
       f.down_seq_region_end
FROM   variation_feature vf, 
       flanking_sequence f,
       tmp_var m
       #yuan_enssnp_mapping_55.failed_mapping_vf_Celera_SARA m
WHERE  vf.variation_id=f.variation_id
#AND    v.source_id=2
AND    vf.variation_id=m.variation_id
       );
      #dump from paralle_processed/or un-processed flanking sequence which is not mapped
      my $sql2 = qq(
SELECT v.name,
       f.variation_id,
       f.seq_region_id,
       f.seq_region_strand,
       f.up_seq,f.down_seq,
       f.up_seq_region_start, 
       f.up_seq_region_end,
       f.down_seq_region_start, 
       f.down_seq_region_end
FROM flanking_sequence f, variation v
LEFT JOIN variation_feature vf
ON   vf.variation_id=f.variation_id
WHERE v.variation_id = f.variation_id
#AND v.source_id=2
AND vf.variation_id is null);
      #dump all from unprocessed flanking sequences or un-mapped from processed
      my $sql3 = qq(
SELECT v.name,
       f.variation_id,
       f.seq_region_id,
       f.seq_region_strand,
       f.up_seq,f.down_seq,
       f.up_seq_region_start, 
       f.up_seq_region_end,
       f.down_seq_region_start, 
       f.down_seq_region_end
FROM flanking_sequence f, variation v
WHERE v.variation_id = f.variation_id
AND f.seq_region_id is null);

      #dump with which sql???
      dumpSQL($dbVar, $sql);
    }
  }


  open ( FH, "<$TMP_DIR/$TMP_FILE");
  while (<FH>) {
    my ($variation_name,
        $variation_id,
        $seq_region_id,
        $seq_region_strand,
        $up_seq,
        $down_seq,
	$up_seq_region_start, 
        $up_seq_region_end,
	$down_seq_region_start, 
        $down_seq_region_end) = split /\t/, $_;
    #print "$variation_name,$variation_id,$seq_region_id,$up_seq,$down_seq,$up_seq_region_start,$up_seq_region_end,$down_seq_region_start,$down_seq_region_end\n";
    $variation_ids{$variation_id}{'var_name'}=$variation_name;
    $variation_ids{$variation_id}{'up_seq'}=$up_seq;
    $variation_ids{$variation_id}{'down_seq'}=$down_seq;
    $variation_ids{$variation_id}{'seq_region_id'}=$seq_region_id;
    if ($up_seq eq '\N') {
      $variation_ids{$variation_id}{'up_seq_region_start'}=$up_seq_region_start;
      $variation_ids{$variation_id}{'up_seq_region_end'}=$up_seq_region_end;
      $variation_ids{$variation_id}{'seq_region_strand'}=$seq_region_strand;
    }
    if ($down_seq eq '\N') {
      $variation_ids{$variation_id}{'down_seq_region_start'}=$down_seq_region_start;
      $variation_ids{$variation_id}{'down_seq_region_end'}=$down_seq_region_end;
      $variation_ids{$variation_id}{'seq_region_strand'}=$seq_region_strand;
    }
    $i++;
    if( $i == $file_size ){
      print_seqs ($slice_adaptor, $file_count, $output_dir, \%variation_ids);
      %variation_ids = ();
      $file_count++;
      $i = 0;
    }
  }
  close FH;
  ###print the last batch of seq
  print_seqs ($slice_adaptor, $file_count, $output_dir, \%variation_ids);

}

sub print_seqs {

  my ($slice_adaptor, $file_count, $output_dir, $variation_ids) = @_;
  my %variation_ids = %$variation_ids;
  my @ids = keys %variation_ids;

  open OUT, ">$output_dir/$file_count\.query_seq" 
      or die "can't open query_seq file : $!";

  foreach my $var_id (@ids) {
    my ($flanking_sequence,$down_seq,$up_seq);
    if ( $variation_ids{$var_id}{'up_seq'} eq '\N') {
      # Upstream sequence not in query; create from slice
      my $up_tmp_slice = $slice_adaptor->fetch_by_seq_region_id 
          ($variation_ids{$var_id}{'seq_region_id'});
      #make sure the total length of up_flanking_seq is 400 bp
      if ($variation_ids{$var_id}{'up_seq_region_end'} - $variation_ids{$var_id}{'up_seq_region_start'} +1 < 400) {
	my $length_added = 400 - ($variation_ids{$var_id}{'up_seq_region_end'} - $variation_ids{$var_id}{'up_seq_region_start'} +1);
	$variation_ids{$var_id}{'up_seq_region_start'} = $variation_ids{$var_id}{'up_seq_region_start'}-$length_added;
	#sometimes the edge is outside the chromosome length
	if ($variation_ids{$var_id}{'up_seq_region_start'}<1) {
	  $variation_ids{$var_id}{'up_seq_region_start'} = 1;
	}
      }
      my $up_seq_slice = $up_tmp_slice->sub_Slice 
          (
           $variation_ids{$var_id}{'up_seq_region_start'},
           $variation_ids{$var_id}{'up_seq_region_end'},
           $variation_ids{$var_id}{'seq_region_strand'}
           ) if $up_tmp_slice;
      if (! $up_tmp_slice or ! $up_seq_slice) {
	print "variation_id $var_id don't have up_seq_slice ".
              "$variation_ids{$var_id}{'seq_region_id'}, ".
              "$variation_ids{$var_id}{'up_seq_region_start'}, ".
              "$variation_ids{$var_id}{'up_seq_region_end'}, ".
              "$variation_ids{$var_id}{'seq_region_strand'}\n";
      }
      else {
	$up_seq = $up_seq_slice->seq;
      }
    }
    else {
      $up_seq = $variation_ids{$var_id}{'up_seq'};
    }
    if ( $variation_ids{$var_id}{'down_seq'} eq '\N') {
      # Downstream sequence not in query; create from slice
      my $down_tmp_slice = $slice_adaptor->fetch_by_seq_region_id 
          ($variation_ids{$var_id}{'seq_region_id'});
      #make sure the total length of down_flanking_seq is 400 bp
      if ($variation_ids{$var_id}{'down_seq_region_end'} - $variation_ids{$var_id}{'down_seq_region_start'} +1 < 400) {
	my $length_added = 400 - ($variation_ids{$var_id}{'down_seq_region_end'} - $variation_ids{$var_id}{'down_seq_region_start'} +1);
	$variation_ids{$var_id}{'down_seq_region_end'} = $variation_ids{$var_id}{'down_seq_region_end'}+$length_added;
	#edge effect
	if ($variation_ids{$var_id}{'down_seq_region_end'} > $down_tmp_slice->length()) {
	  $variation_ids{$var_id}{'down_seq_region_end'} = $down_tmp_slice->length();
	}
      }
      my $down_seq_slice = $down_tmp_slice->sub_Slice 
          (
           $variation_ids{$var_id}{'down_seq_region_start'},
           $variation_ids{$var_id}{'down_seq_region_end'},
           $variation_ids{$var_id}{'seq_region_strand'}
           ) if $down_tmp_slice;
      if (!$down_seq_slice or !$down_tmp_slice) {
	print "variation_id $var_id don't have down_seq_slice ".
              "$variation_ids{$var_id}{'seq_region_id'}, ".
              "$variation_ids{$var_id}{'down_seq_region_start'}, ".
              "$variation_ids{$var_id}{'down_seq_region_end'}, ".
              "$variation_ids{$var_id}{'seq_region_strand'}\n";
      }
      else {
	$down_seq = $down_seq_slice->seq;
      }
    }
    else {
      $down_seq = $variation_ids{$var_id}{'down_seq'};
    }


#     # Trim flanking sequences to 100 bp
#     if (length($up_seq) >100) {
#       $up_seq = substr($up_seq, -100);
#     }
#     if (length($down_seq) >100) {
#       $down_seq = substr($down_seq,0,100);
#     }
    # Join the flanks
    my $seq = lc($up_seq)."W".lc($down_seq);

    # Print to file
    print OUT ">$variation_ids{$var_id}{'var_name'}\n$seq\n";

    #print OUT "$var_id\t$up_seq\t$down_seq\n";
}

  close OUT;
}

sub get_read_flank_seq {
  my ($cdb, $vdb, $read_file) = @_;
  my $slice_adaptor = $cdb->get_SliceAdaptor();
  my @slices = @{$slice_adaptor->fetch_all('toplevel')};

  my %rec_slice;

  if (@slices) {
    foreach my $slice (@slices) {
      my $seq_region_id=$slice->adaptor->get_seq_region_id($slice);
      $rec_slice{$seq_region_id}=$slice;
    }
  }


  if (! $chr_name and $seq_region_id) {
    my $seq_region_name_ref = $dbCore->selectall_arrayref(qq{
SELECT name 
FROM   seq_region 
WHERE  seq_region_id = $seq_region_id});
    $chr_name = $seq_region_name_ref->[0][0];
  }
  print "chr_name is $chr_name\n";

  if ($chr_name) {
    if ($chr_name =~ /NT|^A|^B|^C/) {
      $chr_name = "NT";
    }
    if (! -e "$output_dir/$chr_name") {
      mkdir "$output_dir/$chr_name" or die "can't make dir $chr_name: $!";
    }
    $output_dir = "$output_dir/$chr_name";
    print "output_dir is $output_dir\n";
    $ImportUtils::TMP_DIR = $output_dir;
    $TMP_DIR = $ImportUtils::TMP_DIR;
    print "tmp_file is $TMP_DIR\n";
  }
  print "seq_region_id is $seq_region_id\n";

  my $file_size = 10000;
  my $file_count=1;
  my $i = 0;

  if (!$read_file) {
    if (! -e "$TMP_DIR/$TMP_FILE" or -z "$TMP_DIR/$TMP_FILE") {
      dumpSQL($vdb->dbc, qq{SELECT * FROM read_coverage WHERE seq_region_id=$seq_region_id});
      open ( FH, "$TMP_DIR/$TMP_FILE" );
    }
  }
  else {
    open FH, "$read_file" or die "can't open $read_file\n";
  }

  open OUT, ">$output_dir/$file_count\_query_seq" or die "can't open output\n";
  while (<FH>) {
    if (/^\d+/) {
      my ($seq_region_id, 
          $seq_region_start, 
          $seq_region_end, 
          $level, 
          $sample_id) = split;

      if ($rec_slice{$seq_region_id}) {
	my $up_seq_start = $rec_slice{$seq_region_id}->sub_Slice 
            ($seq_region_start-100,$seq_region_start-1);
	my $down_seq_start = $rec_slice{$seq_region_id}->sub_Slice 
            ($seq_region_start+1,$seq_region_start+100);
	my $up_seq_end = $rec_slice{$seq_region_id}->sub_Slice 
            ($seq_region_end-100,$seq_region_end-1);
	my $down_seq_end = $rec_slice{$seq_region_id}->sub_Slice 
            ($seq_region_end+1,$seq_region_end+100);
	if ($up_seq_start and $down_seq_start) {
	  print OUT ">$seq_region_id\_$seq_region_start\_$seq_region_end\_$level\_$sample_id\_1\n";
	  print OUT lc($up_seq_start->seq).'W'.lc($down_seq_start->seq),"\n";
	}
	else {
	  print "start not slice: $seq_region_id\_$seq_region_start\_$seq_region_end\_$level\_$sample_id\_1\n";
	}
	if ($up_seq_end and $down_seq_end) {
	  print OUT ">$seq_region_id\_$seq_region_start\_$seq_region_end\_$level\_$sample_id\_2\n";
	  print OUT lc($up_seq_end->seq).'W'.lc($down_seq_end->seq),"\n";
        }
	else {
	  print "end not slice: $seq_region_id\_$seq_region_start\_$seq_region_end\_$level\_$sample_id\_1\n";
	}
	$i=$i+2;
	if ($i == $file_size) {
	  $file_count++;
	  close OUT;
	  open OUT, ">$output_dir/$file_count\_query_seq" or die "can't open output\n";
	  $i=0;
	}
      }
    }
  }
}

sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: generate_input_seq.pl  <options>

options:
    -ensembl_registry <path> 
       Use this ensembl registry file, 
       default: ensembl-variation/scripts/import/ensembl.registry 

    -species  <string>   
       Use DBs for this species in the ensembl registry file

    -chr_name <string>   
       Dump sequences mapping to this chromosome only

    -seq_region_id <int> 
       Dump sequences mapping to this seq_region_id only

    -output_dir <path>
      Dump sequences into this directoty. Def $ENV{PWD}

    -tmpdir <path>       
       Use this dir for tmp files. Def; $ImportUtils::TMP_DIR

    -tmpfile <name>    
       Use this file for tmp file. Def; $ImportUtils::TMP_FILE

    -generate_input_seq
       Whether to generate sequence files

    -read_flank

    -read_file <path>
       Path to file containing flanking sequences

EOF

  die("\n$msg\n\n");
}

###example:
###./generate_input_seq.pl -cdbname homo_sapiens_core_27_35a -vdbname homo_sapiens_variation_27_35a -chr_name 21
