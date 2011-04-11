#! /usr/local/bin/perl
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
our ($species, $input_file, $source_name, $TMP_DIR, $TMP_FILE);

GetOptions('species=s'     => \$species,
	   	   	 'input_file=s'  => \$input_file,
	   	   	 'source_name=s' => \$source_name,
	   	   	 'tmpdir=s'      => \$ImportUtils::TMP_DIR,
         	 'tmpfile=s'     => \$ImportUtils::TMP_FILE,
          );

$species ||= 'human';
$input_file ||= '/lustre/scratch103/ensembl/lg10/hgmd/2011.1-hgmd-public.tsv';
$source_name ||= 'HGMD-PUBLIC';

my $header = 'HGMD_PUBLIC';
my $ncbi_version = 'NCBI36';
my $registry_file ||= $Bin . "/ensembl.registry";
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
my $source_description = 'Variants from HGMD-PUBLIC data';

my $month_hash = { 1 => 'January', 2 => 'February', 3 => 'March', 4 => 'April',
									 5 => 'May', 6 => 'June', 7 => 'July', 8 => 'August',
									 9 => 'September', 10 => 'October', 11 => 'November', 12 => 'December'
            		 };

my $new_source_ref = $dbVar2->selectall_arrayref(qq{SELECT source_id FROM source WHERE name = "$source_name"});
if ($new_source_ref->[0][0]) {
  $source_id = $new_source_ref->[0][0];
}
else {
	$input_file =~ /\/(\d{4})\.(\d+)-hgmd-public\.tsv$/;
	my $year      = $1;
	my $month_num	= $2;
	my $month     = $month_hash->{$month_num} if ($month_num);
	if ($year and $month) {
		$source_description .= " set $month $year";
		$month_num = '0'.$month_num if ($month_num<10);
  	$dbVar2->do(qq{INSERT INTO source(name,description,url,version)
									 values("$source_name","$source_description [remapped from build $ncbi_version]","$hgmd_url",$year$month_num)
									});
	}
	else {
		$dbVar2->do(qq{INSERT INTO source(name,description,url)
									 values("$source_name","$source_description [remapped from build $ncbi_version]","$hgmd_url")
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
  my $variation_id;
  
  my $slice_newdb_oldasm = $sa2->fetch_by_region('chromosome', $seq_region_name, undef, undef, undef, "$ncbi_version");
  if (!$slice_newdb_oldasm) {
    print_buffered($buffer,"$TMP_DIR/$header\_error",join ("\t",$var_name,$seq_region_name,$pos) . "\n");
    next;
  }
    
  my $seq_region_strand = 1;
  my $feat = new Bio::EnsEMBL::SimpleFeature(
					      -START => $pos,
					      -END => $pos,
					      -STRAND => $seq_region_strand,
					      -SLICE => $slice_newdb_oldasm,
					    );

  #print $feat->start," ",$feat->strand," ",$feat->slice,"\n";
  my @segments;

  eval {
    @segments = @{$feat->feature_Slice->project('chromosome','GRCh37')};
  };
  if ($@) {
    print "no segments found for $var_name\n" if $@;
    print_buffered($buffer,"$TMP_DIR/$header\_error",join ("\t",$var_name,$seq_region_name,$pos) . "\n");
    next;
  }

  my @slices_newdb_newasm = map { $_->to_Slice }  @segments;
  @slices_newdb_newasm = sort {$a->start<=>$b->start} @slices_newdb_newasm;

  my $validation_status = '\N';
  my $map_weight = 1;
  my $allele_string = 'HGMD_MUTATION';
  my ($allele1,$allele2) = ($allele_string,$allele_string);
  my $flags = '\N';

  #make new strand is same as old strand 
  my $new_strand = $seq_region_strand;
  my ($new_start,$new_end);
  my $consequence = 'HGMD_MUTATION';
  my $indel;
  if (@slices_newdb_newasm) {
    if ($indel) {
      $new_start = $slices_newdb_newasm[0]->start +1;
      $new_end = $slices_newdb_newasm[-1]->end -1;
    }
    else {
      $new_start = $slices_newdb_newasm[0]->start;
      $new_end = $slices_newdb_newasm[-1]->end;
    }

    my $new_seq_region_id = $slices_newdb_newasm[0]->get_seq_region_id;

    # new feature spans start of first projection segment to end of last segment
    #print "old_seq_region_name is ",$old_slice->get_seq_region_id," old_start is ",$feat->start," old_end is ",$feat->end," new_seq_region_name is ",$new_seq_region_id," new_start is ",$new_start," new_end is ",$new_end,"\n" if (@slices_newdb_newasm);
 
    $dbVar2->do(qq{INSERT INTO $header\_variation(source_id,name,type)values($source_id,"$var_name","$var_type")});
    $variation_id = $dbVar2->{'mysql_insertid'};

	 	print_buffered($buffer,"$TMP_DIR/$header\_variation_feature",join ("\t",$new_seq_region_id,$new_start,$new_end,$new_strand,"$variation_id\t$allele_string\t$var_name\t$map_weight\t$flags\t$source_id\t$validation_status\t$consequence\n"));
    print_buffered($buffer,"$TMP_DIR/$header\_flanking_sequence",join ("\t",$variation_id,'\N','\N',$new_start-100,$new_start-1,$new_start+1,$new_start+100,$new_seq_region_id,$new_strand)."\n");
    print_buffered($buffer,"$TMP_DIR/$header\_allele", join ("\t",$variation_id,$allele1) . "\n");
	  print_buffered($buffer,"$TMP_DIR/$header\_variation_annotation", join ("\t",$variation_id,$var_name,$gene_symbol) . "\n");
		
	  #print_buffered($buffer,"$TMP_DIR/$header\_allele", join ("\t",$variation_id,$allele2) . "\n");
    #print OUT join "\t",$slices_newdb_newasm[0]->get_seq_region_id,$new_start,$new_end,$slices_newdb_newasm[0]->strand,"$variation_id\t$allele_string\t$variation_name\t$map_weight\t$flags\t$source_id\t$validation_status\n";
  }
  else {
    print_buffered($buffer,"$TMP_DIR/$header\_error",join ("\t",$var_name,$seq_region_name,$pos) . "\n") if $var_name;
    #print ERR "$variation_id\t$seq_region_id\t$seq_region_start\t$seq_region_end\n";
  }
}


print_buffered($buffer);

debug("Loading mapping data...");

foreach my $file ("$TMP_DIR/$header\_variation_feature","$TMP_DIR/$header\_flanking_sequence","$TMP_DIR/$header\_allele","$TMP_DIR/$header\_variation_annotation","$TMP_DIR/$header\_error") {
  if (-e "$file") {
    system("mv $file  $TMP_DIR/$TMP_FILE") ;
    $file =~ s/$TMP_DIR\///;
    $file =~ s/\-/\_/g;
    my $mapping_ref = $dbVar2->selectall_arrayref(qq{show tables like "$file"});
    if (!$mapping_ref->[0][0]) {
      if ($file =~ /variation_feature/) {
				create_and_load ($dbVar2,"$file","seq_region_id i*","seq_region_start i","seq_region_end i","seq_region_strand i","variation_id i*","allele_string","variation_name","map_weight","flags","source_id i","validation_status","consequence_type");
      }
      elsif ($file =~ /allele/) {
				create_and_load ($dbVar2,"$file","variation_id i*","allele");
      }
      elsif ($file =~ /flanking/) {
				create_and_load ($dbVar2,"$file","variation_id i*","up_seq","down_seq","up_seq_region_start i","up_seq_region_end i","down_seq_region_start i","down_seq_region_end i","seq_region_id i*","seq_region_strand");
      }
			elsif ($file =~ /variation_annotation/) {
				create_and_load ($dbVar2,"$file","variation_id i*","variation_names","associated_gene");
      }
      elsif ($file =~ /error/) {
				create_and_load ($dbVar2,"$file","variation_name *","seq_region_name","seq_region_start");
      }
    }
    else {
      if ($file =~ /variation_feature/) {
				load ($dbVar2,"$file","seq_region_id","seq_region_start","seq_region_end","seq_region_strand","variation_id","allele_string","variation_name","map_weight","flags","source_id","validation_status","consequence_type");
      }
      elsif ($file =~ /allele/) {
				load ($dbVar2,"$file","variation_id","allele");
      }
      elsif ($file =~ /flanking/) {
				load ($dbVar2,"$file","variation_id","up_seq","down_seq","up_seq_region_start","up_seq_region_end","down_seq_region_start","down_seq_region_end","seq_region_id","seq_region_strand");
      }
			elsif ($file =~ /variation_annotation/) {
				load ($dbVar2,"$file","variation_id","variation_names","associated_gene");
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
