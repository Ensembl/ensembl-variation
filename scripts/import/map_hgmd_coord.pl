#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2025] EMBL-European Bioinformatics Institute
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
  <http://www.ensembl.org/Help/Contact>.

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
our ($species, $input_file, $source_name, $TMP_DIR, $TMP_FILE, $registry_file);

GetOptions('species=s'     => \$species,
            'input_file=s'  => \$input_file,
            'source_name=s' => \$source_name,
            'tmpdir=s'      => \$ImportUtils::TMP_DIR,
            'tmpfile=s'     => \$ImportUtils::TMP_FILE,
            'registry=s'    => \$registry_file
          );

$species ||= 'human';
die("You need to specify an input file (option -input_file)") if (!defined($input_file));
$source_name ||= 'HGMD-PUBLIC';

my $header = 'HGMD_PUBLIC';
   $registry_file ||= $Bin . "/ensembl.registry";
my $hgmd_url = 'http://www.hgmd.cf.ac.uk/ac/index.php';
my $insertion_label = 'I';

usage('-species argument is required') if(!$species);

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

Bio::EnsEMBL::Registry->load_all( $registry_file );
my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbVar = $vdb->dbc->db_handle;

my $sa = $cdb->get_SliceAdaptor();
my $buffer = {};
my $source_id;
my $source_description = 'Variants from HGMD-PUBLIC dataset';

my $month_hash = { 1 => 'March', 2 => 'June', 3 => 'September', 4 => 'December' };

my $new_source_ref = $dbVar->selectall_arrayref(qq{SELECT source_id FROM source WHERE name = "$source_name"});

$input_file =~ /(\d{4})\.(\d+)-hgmd-public/;
my $year      = $1;
my $month_num = $2;
my $month     = $month_hash->{$month_num} if ($month_num);

$source_description .= " $month $year" if ($year and $month);

if ($new_source_ref->[0][0]) {
  $source_id = $new_source_ref->[0][0];
  $dbVar->do(qq{UPDATE source SET description="$source_description", version=$year$month_num 
                 WHERE source_id=$source_id });
}
else {
  if ($year and $month) {
    $dbVar->do(qq{INSERT INTO source(name,description,url,version)
                   values("$source_name","$source_description","$hgmd_url",$year$month_num)
                  });
  }
  else {
    $dbVar->do(qq{INSERT INTO source(name,description,url)
                   values("$source_name","$source_description","$hgmd_url")
                  });
  }
  $source_id = $dbVar->{'mysql_insertid'};
}


my %seq_regions;

$dbVar->do(qq{CREATE TABLE IF NOT EXISTS `$header\_variation` (
  `variation_id` int(10) unsigned not null auto_increment, # PK
  `source_id` int(10) unsigned not null,
  `name` varchar(255),
  `type` char(1) NOT NULL,
  `class_attrib_id` int(10) unsigned DEFAULT NULL,
  
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
              'start'       => $pos,
              'end'         => $pos,
              'var_type'    => $var_type,
              'strand'      => 1,
              'status'      => '\N',
              'map_weight'  => 1,
              'allele'      => 'HGMD_MUTATION',
              'flags'       => '\N',
              'consequence' => 'HGMD_MUTATION'
            );
    
  my $seq_region_id;

  if(!$seq_regions{$seq_region_name}) {
    my $seq_region_sth = $dbVar->prepare(qq{ SELECT seq_region_id FROM seq_region WHERE name = ? });
    $seq_region_sth->execute($seq_region_name);

    my $result = $seq_region_sth->fetchall_arrayref();
    $seq_region_id = $result->[0]->[0];

    if(!$seq_region_id) {
      print_buffered($buffer,"$TMP_DIR/$header\_error",join ("\t",$data{var_name},$data{region_name},$data{start}) . "\n");
      next;
    }

    $seq_regions{$seq_region_name} = $seq_region_id;
  }
  else {
    $seq_region_id = $seq_regions{$seq_region_name};
  }

  $data{region_id} = $seq_region_id;

  print_all_buffers(\%data);
}

print_buffered($buffer);

debug("Loading mapping data...");

foreach my $file ("$TMP_DIR/$header\_variation_feature","$TMP_DIR/$header\_variation_annotation","$TMP_DIR/$header\_error") {
  if (-e "$file") {
    system("mv $file  $TMP_DIR/$TMP_FILE") ;
    $file =~ s/$TMP_DIR\///;
    $file =~ s/\-/\_/g;
    my $mapping_ref = $dbVar->selectall_arrayref(qq{show tables like "$file"});
    if (!$mapping_ref->[0][0]) {
      if ($file =~ /variation_feature/) {
        create_and_load ($dbVar,"$file","seq_region_id i*","seq_region_start i","seq_region_end i","seq_region_strand i","variation_id i*","allele_string","variation_name","map_weight","flags","source_id i","validation_status","consequence_type", "class_attrib_id");
      }
      elsif ($file =~ /variation_annotation/) {
        create_and_load ($dbVar,"$file","variation_id i*","variation_name","associated_gene");
      }
      elsif ($file =~ /error/) {
        create_and_load ($dbVar,"$file","variation_name *","seq_region_name","seq_region_start");
      }
    }
    else {
      if ($file =~ /variation_feature/) {
        load ($dbVar,"$file","seq_region_id","seq_region_start","seq_region_end","seq_region_strand","variation_id","allele_string","variation_name","map_weight","flags","source_id","validation_status","consequence_type");
      }
      elsif ($file =~ /variation_annotation/) {
        load ($dbVar,"$file","variation_id","variation_name","associated_gene");
      }
      elsif ($file =~ /error/) {
        load ($dbVar,"$file","variation_name","seq_region_name","seq_region_start");
      }
    }
  }
}
updated_vf_coordinate_insertion($insertion_label);
add_class_attrib_id();



# Change the coordinates for the variation where the class is "insertion"
sub updated_vf_coordinate_insertion {
  my $insertion_class = shift;
  my $select_vf_sth = $dbVar->prepare(qq{
    SELECT distinct v.name,vf.seq_region_end 
    FROM $header\_variation v, $header\_variation_feature vf 
    WHERE v.type='$insertion_class' AND v.name=vf.variation_name
  });

  my $update_vf_sth = $dbVar->prepare(qq{
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


# Add class_attrib_id for the tmp variation, variation_feature tables
sub add_class_attrib_id {
  my %attrib = ('M' => 'SNV',
                'D' => 'deletion',
                'I' => 'insertion',
                'X' => 'indel',
                'P' => 'indel',
                'R' => 'sequence_alteration',
                'S' => 'sequence_alteration'
               );

  my %class = ();

  my $select_a_sth = $dbVar->prepare(qq{
    SELECT a.attrib_id FROM attrib a, attrib_type att
    WHERE a.attrib_type_id = att.attrib_type_id AND att.code = 'SO_term' AND a.value = ?;
  });

  my $select_v_sth = $dbVar->prepare(qq{
    SELECT DISTINCT variation_id,type FROM $header\_variation;
  });

  my $update_v_sth = $dbVar->prepare(qq{
    UPDATE $header\_variation SET class_attrib_id = ? WHERE variation_id = ?;
  });

  my $update_vf_sth = $dbVar->prepare(qq{
    UPDATE $header\_variation_feature vf, $header\_variation v SET vf.class_attrib_id = v.class_attrib_id
    WHERE v.variation_id = vf.variation_id;
  });

  while (my ($k,$v) = each (%attrib)) {
    $select_a_sth->execute($v);
    $class{$k} = ($select_a_sth->fetchrow_array)[0];
    print "$k: ".$class{$k}."\n";
  }

  $select_v_sth->execute();
  while (my @res = $select_v_sth->fetchrow_array) {
    my $att = $class{$res[1]};
    if (defined $att) {
      $update_v_sth->execute($att,$res[0]) or die $!;
    }
  }

  $update_vf_sth->execute() or die $!;

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
  
  $dbVar->do(qq{INSERT INTO $header\_variation(source_id,name,type)values($source_id,"$var_name","$var_type")});
  my $variation_id = $dbVar->{'mysql_insertid'};

  print_buffered($buffer,"$TMP_DIR/$header\_variation_feature",join ("\t",$data->{region_id},$start,$data->{end},$data->{strand},"$variation_id\t".$data->{allele}."\t$var_name\t".$data->{map_weight}."\t".$data->{flags}."\t$source_id\t".$data->{status}."\t".$data->{consequence}."\n"));
  print_buffered($buffer,"$TMP_DIR/$header\_allele", join ("\t",$variation_id,$data->{allele}) . "\n");
  print_buffered($buffer,"$TMP_DIR/$header\_variation_annotation", join ("\t",$variation_id,$var_name,$data->{gene_symbol}) . "\n");
}
