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

#
#./bsub_ssaha2.pl -species human -input_dir /ecs4/scratch2/yuan/hum/mapping_36 -output_dir /ecs4/scratch2/yuan/hum/mapping_36 -target_dir /ecs4/scratch2/yuan/hum/mapping_36/target_dna -start 1 -end 9

use strict;
use SARA::GetSARA;
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

our ($species, $input_dir, $TMP_DIR, $TMP_FILE);

GetOptions('species=s'    => \$species,
	   'input_dir=s'  => \$input_dir,
 	   'tmpdir=s'     => \$ImportUtils::TMP_DIR,
	   'tmpfile=s'    => \$ImportUtils::TMP_FILE,
	  );
my $registry_file ||= $Bin . "/ensembl.registry";

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;


Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdb_var = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $vdb_sara = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'sara_hum');

my $dbCore = $cdb->dbc if $cdb;
my $dbVar = $vdb_var->dbc if $vdb_var;
my $dbSara = $vdb_sara->dbc if $vdb_sara;

#need three files, i.e cigar_file, match_file and ALIGNMENT_file
#&load_cigar_to_db();

my $rat = SARA::GetSARA->new(-dbCore => $dbCore,
			     -dbVar => $dbVar,
			     -dbSara => $dbSara,
			     -tmpdir => $TMP_DIR,
			     -tmpfile => $TMP_FILE,
			     -species => $species,
			    );
$rat->get_sara();

sub load_cigar_to_db {

  debug("Dumping seq_region data");

  #only take toplevel coordinates
  my %seq_region_ids;
  my $sth = $dbCore->prepare(qq{SELECT sr.seq_region_id, sr.name
  		                FROM   seq_region_attrib sra, attrib_type at, seq_region sr
 		                WHERE sra.attrib_type_id=at.attrib_type_id 
	                        AND at.code="toplevel" 
                                AND sr.seq_region_id = sra.seq_region_id 
		               });
  $sth->execute();
  while (my ($seq_region_id,$seq_region_name) = $sth->fetchrow_array()) {
    #print "seq_reion_name is $seq_region_name\n";
    $seq_region_ids{$seq_region_name} = $seq_region_id;
  }

  #my ($base_dir) = $input_dir =~ /^(.*)\/out.*$/;
  my $cigar_file = "$input_dir/cigar_file";
  my $match_file = "$input_dir/match_file";
  my $align_file = "$input_dir/ALIGNMENT_file";
  
  #get_files($input_dir,$cigar_file,"^1");
  #get_files($input_dir,$match_file,"^Matches For Query");

  read_cigar_file($dbSara,$cigar_file,\%seq_region_ids);
  read_match_file($dbSara,$match_file,$align_file);
}

sub get_files {

  my ($input_dir,$file_name,$pattern) = @_;

  debug("Dumping $file_name");

  opendir DIRENTRY, "$input_dir" || die "Failed to open dir : $!";
  open FILE, ">$file_name" or die "can't open file $!";

  my @files = grep /^ssaha_out/, readdir(DIRENTRY);

  foreach my $file (@files) {
    my $call = "grep \"$pattern\" $input_dir/$file >> $file_name";
    system($call);
  }

  close FILE;
  return;
}

sub read_cigar_file {

  my ($dbSara, $cigar_file, $seq_region_ids) = @_;

  debug("Reading cigar_file");

  my %seq_region_ids = %$seq_region_ids;

  open IN, "$cigar_file" or die "can't open cigar_file $!";
  open OUT, ">$TMP_DIR/$TMP_FILE" or die "can't open tmp_file $!";

  while (<IN>) {
    if (/^\S+\s+\d+\s+\d+\s+/) {
    $_ =~ s/^.*\:// if $_ =~/\:/; 
    #gnl|ti|904511325 0 778 + 8-1-129041809 73391177 73391953 + 762 M 754 I 1 M 16 I 1 M 6
    my ($query_name,$query_start,$query_end,$query_strand,$target_name,$target_start,$target_end,$target_strand,$score,@cigars) = split;
    $query_start = $query_start + 1;
    $target_start = $target_start + 1;
    my ($target_name1) = $target_name =~ /^(.*)\-.*\-.*$/;
    my $cigar_string;
    while (@cigars) {
      my $tag = shift @cigars;
      my $base = shift @cigars;
      $cigar_string .= $base.$tag;
    }
    $query_strand = ($query_strand eq "+") ? 1 : -1;
    $target_strand = ($target_strand eq "+") ? 1 : -1;

    print OUT join ("\t", $query_name,$query_start,$query_end,$query_strand,$seq_region_ids{$target_name1},$target_name,$target_start,$target_end,$target_strand,$score,$cigar_string, "SD"),"\n";
  }
  }
  debug("Loading ssahaSNP_feature table");
  load($dbSara, "ssahaSNP_feature", "query_name","query_start","query_end","query_strand","target_seq_region_id","target_name","target_start","target_end","target_strand","score","cigar_string","individual_name");
}

sub read_match_file {

  my ($dbSara, $match_file, $align_file) = @_;
  my %rec_strain_name;
  debug("Reading Match file...");

  open IN, "$match_file" or die "can't open match_file $!";
  open IN1, "$align_file" or die "can't open align_file $!";
  open OUT, ">$TMP_DIR/$TMP_FILE" or die "can't open tmp_file $!";

  while (<IN1>) {
    my @all=split;
    if (@all==13) {
      my $reads_name = $all[2];
      $reads_name =~ s/\.p1k$|\.q1k$//;
      my $strain_name = $all[12];
      $rec_strain_name{$reads_name}=$strain_name;
    }
  }    
  while (<IN>) {
    #match line looks like: Matches For Query 0 (907 bases): gnl|ti|925271746
    $_ =~ /^.*\:Matches For Query\s+\d+\s+\((\d+)\s+bases\)\:\s+(.*)$/;
    my $length = $1;
    my $query_name = $2;
    my $strain_name = ($rec_strain_name{$query_name}) ? $rec_strain_name{$query_name} : '\N';
    print OUT "$query_name\t$length\t$strain_name\n";
  }

  debug("Loading query_match_length table...");
  create_and_load($dbSara,"query_match_length_strain","query_name *","length","strain_name");
}
