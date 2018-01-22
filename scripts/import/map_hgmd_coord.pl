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

# Modified by Laurent the 01/12/2010

use strict;
#use DBH;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Data::Dumper;
use Bio::SeqIO;
use FindBin qw( $Bin );
use Getopt::Long;
use ImportUtils qw(dumpSQL debug create_and_load load);
use DBI qw(:sql_types);
our ($species, $input_file, $source_name, $TMP_DIR, $TMP_FILE, $mapping, $registry_file);

GetOptions('species=s'     => \$species,
             'input_file=s'  => \$input_file,
             'source_name=s' => \$source_name,
             'tmpdir=s'      => \$ImportUtils::TMP_DIR,
            'tmpfile=s'     => \$ImportUtils::TMP_FILE,
           'mapping!'      => \$mapping,
           'registry=s'    => \$registry_file
          );

$species ||= 'human';
die("You need to specify an input file (option -input_file)") if (!defined($input_file));
$source_name ||= 'HGMD-PUBLIC';

my $header = 'HGMD_PUBLIC';
my $ncbi_version = 'NCBI36';
   $registry_file ||= $Bin . "/ensembl.registry";
my $hgmd_url = 'http://www.hgmd.cf.ac.uk/ac/index.php';
my $insertion_label = 'I';

usage('-species argument is required') if(!$species);

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

Bio::EnsEMBL::Registry->load_all( $registry_file );
my $cdb2 = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdb2 = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbCore2 = $cdb2->dbc->db_handle;
my $dbVar2 = $vdb2->dbc->db_handle;

my $sa2 = $cdb2->get_SliceAdaptor();
my $buffer = {};
my $source_id;
my $source_description = 'Variants from HGMD-PUBLIC dataset';

my $month_hash = { 1 => 'March', 2 => 'June', 3 => 'September', 4 => 'December' };

my $new_source_ref = $dbVar2->selectall_arrayref(qq{SELECT source_id FROM source WHERE name = "$source_name"});

$input_file =~ /\/(\d{4})\.(\d+)-hgmd-public\.tsv$/;
my $year      = $1;
my $month_num  = $2;
my $month     = $month_hash->{$month_num} if ($month_num);
my $remapped;
if ($mapping) {
  $remapped = " [remapped from build $ncbi_version]";
}
$source_description .= " $month $year" if ($year and $month);

if ($new_source_ref->[0][0]) {
  $source_id = $new_source_ref->[0][0];
  $dbVar2->do(qq{UPDATE source SET description="$source_description$remapped", version=$year$month_num 
                 WHERE source_id=$source_id });
}
else {
  if ($year and $month) {
    $dbVar2->do(qq{INSERT INTO source(name,description,url,version)
                   values("$source_name","$source_description$remapped","$hgmd_url",$year$month_num)
                  });
  }
  else {
    $dbVar2->do(qq{INSERT INTO source(name,description,url)
                   values("$source_name","$source_description$remapped","$hgmd_url")
                  });
  }
  $source_id = $dbVar2->{'mysql_insertid'};
}


$dbVar2->do(qq{CREATE TABLE IF NOT EXISTS `$header\_variation` (
  `variation_id` int(10) unsigned not null auto_increment, # PK
  `source_id` int(10) unsigned not null,
  `name` varchar(255),
  `type` char(1) NOT NULL,
  
  primary key( `variation_id` ),
  unique ( `name` ),
  key source_idx (`source_id`)
)});


open IN, "$input_file" or die "can't open file $input_file : $!";

while (<IN>) {
  chomp;
  my ($var_name,$gene_symbol,$seq_region_name,$pos,$var_type) = split /\s+/, $_;
  my %data = ('var_name'    => $var_name,
              'gene_symbol' => $gene_symbol,
              'region_name' => $seq_region_name,
              'pos'         => $pos,
              'var_type'    => $var_type,
              'strand'      => 1,
              'status'      => '\N',
              'map_weight'  => 1,
              'allele'      => 'HGMD_MUTATION',
              'flags'       => '\N',
              'consequence' => 'HGMD_MUTATION'
            );
  
  
  if ($mapping) {
    my $slice_newdb_oldasm = $sa2->fetch_by_region('chromosome', $data{region_name}, undef, undef, undef, "$ncbi_version");
    if (!$slice_newdb_oldasm) {
      print_buffered($buffer,"$TMP_DIR/$header\_error",join ("\t",$data{var_name},$data{region_name},$data{pos}) . "\n");
      next;
    }
    
    my $feat = new Bio::EnsEMBL::SimpleFeature(
                  -START  => $data{pos},
                  -END    => $data{pos},
                  -STRAND => $data{strand},
                  -SLICE  => $slice_newdb_oldasm,
                );

    #print $feat->start," ",$feat->strand," ",$feat->slice,"\n";
    my @segments;

    eval {
      @segments = @{$feat->feature_Slice->project('chromosome','GRCh37')};
    };
    if ($@) {
      print "no segments found for ".$data{var_name}."\n" if $@;
      print_buffered($buffer,"$TMP_DIR/$header\_error",join ("\t",$data{var_name},$data{region_name},$data{pos}) . "\n");
      next;
    }
  
    my @slices_newdb_newasm = map { $_->to_Slice }  @segments;
    @slices_newdb_newasm = sort {$a->start<=>$b->start} @slices_newdb_newasm;

    #make new strand is same as old strand 
    my $indel;
    if (@slices_newdb_newasm) {
      $data{start} = $slices_newdb_newasm[0]->start;
      $data{end}   = $slices_newdb_newasm[-1]->end;
      if ($indel) {
        $data{start} ++;
        $data{end} --;
      }

      $data{region_id} = $slices_newdb_newasm[0]->get_seq_region_id;

      print_all_buffers(\%data);
    }
    else {
      print_buffered($buffer,"$TMP_DIR/$header\_error",join ("\t",$var_name,$seq_region_name,$pos) . "\n") if $var_name;
    }
  }
  else {
    $data{start} = $data{pos};
    $data{end}   = $data{pos};
    
    my $slice = $sa2->fetch_by_region('chromosome', $data{region_name}, $data{start},$data{end});
    if (!$slice) {
      print_buffered($buffer,"$TMP_DIR/$header\_error",join ("\t",$data{var_name},$data{region_name},$data{pos}) . "\n");
      next;
    }
    $data{region_id} = $slice->get_seq_region_id;
    
    print_all_buffers(\%data);
  }
}

print_buffered($buffer);

debug("Loading mapping data...");

foreach my $file ("$TMP_DIR/$header\_variation_feature","$TMP_DIR/$header\_variation_annotation","$TMP_DIR/$header\_error") {
  if (-e "$file") {
    system("mv $file  $TMP_DIR/$TMP_FILE") ;
    $file =~ s/$TMP_DIR\///;
    $file =~ s/\-/\_/g;
    my $mapping_ref = $dbVar2->selectall_arrayref(qq{show tables like "$file"});
    if (!$mapping_ref->[0][0]) {
      if ($file =~ /variation_feature/) {
        create_and_load ($dbVar2,"$file","seq_region_id i*","seq_region_start i","seq_region_end i","seq_region_strand i","variation_id i*","allele_string","variation_name","map_weight","flags","source_id i","validation_status","consequence_type");
      }
      elsif ($file =~ /variation_annotation/) {
        create_and_load ($dbVar2,"$file","variation_id i*","variation_name","associated_gene");
      }
      elsif ($file =~ /error/) {
        create_and_load ($dbVar2,"$file","variation_name *","seq_region_name","seq_region_start");
      }
    }
    else {
      if ($file =~ /variation_feature/) {
        load ($dbVar2,"$file","seq_region_id","seq_region_start","seq_region_end","seq_region_strand","variation_id","allele_string","variation_name","map_weight","flags","source_id","validation_status","consequence_type");
      }
      elsif ($file =~ /variation_annotation/) {
        load ($dbVar2,"$file","variation_id","variation_name","associated_gene");
      }
      elsif ($file =~ /error/) {
        load ($dbVar2,"$file","variation_name","seq_region_name","seq_region_start");
      }
    }
  }
}
updated_vf_coordinate_insertion($insertion_label);



# Change the coordinates for the variation where the class is "insertion"
sub updated_vf_coordinate_insertion {
  my $insertion_class = shift;
  my $select_vf_sth = $dbVar2->prepare(qq{
    SELECT distinct v.name,vf.seq_region_end 
    FROM $header\_variation v, $header\_variation_feature vf 
    WHERE v.type='$insertion_class' AND v.name=vf.variation_name
  });

  my $update_vf_sth = $dbVar2->prepare(qq{
    UPDATE $header\_variation_feature  
     SET seq_region_end = ?
     WHERE variation_name = ?;
  });

  $select_vf_sth->execute();
  while(my @res = $select_vf_sth->fetchrow_array){
    my $end = $res[1]-1;
    $update_vf_sth->execute($end, $res[0]) or die "Cannot execute: " . $update_vf_sth->errstr();
  }
}


sub print_buffered {
  my $buffer = shift;
  my $filename = shift;
  my $text = shift;

  local *FH;

  if( ! $filename ) {
    # flush the buffer
    foreach my $file (keys %{$buffer}){
      open( FH, ">>$file" ) or die;
      print FH $buffer->{ $file };
      close FH;
    }
    %{$buffer} = ();

  } 
  else {
    $buffer->{ $filename } .= $text;
    if( length( $buffer->{ $filename } ) > 10_000 ) {
       open( FH, ">>$filename" ) or die;
       print FH $buffer->{ $filename };
       close FH;
       $buffer->{ $filename } = '';
    }
  }
}

sub print_all_buffers {
  my $data = shift;
  my $var_name = $data->{var_name};
  my $var_type = $data->{var_type};
  my $start    = $data->{start};
  
  $dbVar2->do(qq{INSERT INTO $header\_variation(source_id,name,type)values($source_id,"$var_name","$var_type")});
  my $variation_id = $dbVar2->{'mysql_insertid'};

  print_buffered($buffer,"$TMP_DIR/$header\_variation_feature",join ("\t",$data->{region_id},$start,$data->{end},$data->{strand},"$variation_id\t".$data->{allele}."\t$var_name\t".$data->{map_weight}."\t".$data->{flags}."\t$source_id\t".$data->{status}."\t".$data->{consequence}."\n"));
  print_buffered($buffer,"$TMP_DIR/$header\_allele", join ("\t",$variation_id,$data->{allele}) . "\n");
  print_buffered($buffer,"$TMP_DIR/$header\_variation_annotation", join ("\t",$variation_id,$var_name,$data->{gene_symbol}) . "\n");
}
