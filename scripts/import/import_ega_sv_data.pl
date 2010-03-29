#! /usr/local/bin/perl
#

use strict;
#use DBH;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Data::Dumper;
use FindBin qw( $Bin );
use Getopt::Long;
use ImportUtils qw(dumpSQL debug create_and_load load);
our ($species, $input_file, $source_name, $TMP_DIR, $TMP_FILE, $mapping, $num_gaps, $target_assembly);

GetOptions('species=s'         => \$species,
		   'source_name=s'     => \$source_name,
		   'input_file=s'	   => \$input_file,
           'tmpdir=s'          => \$ImportUtils::TMP_DIR,
           'tmpfile=s'         => \$ImportUtils::TMP_FILE,
		   'mapping'		   => \$mapping,
		   'gaps=i'			   => \$num_gaps,
          );
my $registry_file ||= $Bin . "/ensembl.registry";

$num_gaps = 1 unless defined $num_gaps;
$species ||= 'human';

# set the target assembly
$target_assembly = 'GRCh37';

usage('-species argument is required') if(!$species);
usage('-input_file argument is required') if (!$input_file);

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

# connect to databases
Bio::EnsEMBL::Registry->load_all( $registry_file );
my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbCore = $cdb->dbc->db_handle;
my $dbVar = $vdb->dbc->db_handle;

# run the mapping sub-routine if the data needs mapping
if($mapping) {
  mapping();
}

# otherwise parse the file into the temp file
else {
  open IN, $input_file or die "Could not read from input file $input_file\n";
  open OUT, ">$TMP_DIR/$TMP_FILE" or die "Could not write to output file $TMP_DIR/$TMP_FILE\n";
  
  while(<IN>) {
	next if /STUDY_ID/;		# skip header line
	print OUT $_;
  }
  
  close IN;
  close OUT;
}


read_file();
source();
structural_variation();
meta_coord();

sub read_file{
  debug("Loading file into temporary table");
  
  # drop any table that may be there
  $dbVar->do("DROP TABLE IF EXISTS temp_cnv;");
  
  create_and_load(
	$dbVar, "temp_cnv", "study", "id *", "chr", "outer_start i", "inner_start i",
	"start i", "end i", "inner_end i", "outer_end i", "assembly", "type");
  
  # fix nulls
  foreach my $coord('outer_start', 'inner_start', 'inner_end', 'outer_end') {
	$dbVar->do(qq{UPDATE temp_cnv SET $coord = NULL WHERE $coord = 0;});
  }
}

sub source{
  debug("Inserting into source table");
  $dbVar->do(qq{INSERT INTO source (name) SELECT distinct(concat('DGVa:', study)) FROM temp_cnv;});
}

sub structural_variation{
  
  debug("Inserting into structural_variation table");
  
  # create table if it doesn't exist yet
  $dbVar->do(qq{CREATE TABLE IF NOT EXISTS `structural_variation` (
	`structural_variation_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
	`seq_region_id` int(10) unsigned NOT NULL,
	`seq_region_start` int(11) NOT NULL,
	`seq_region_end` int(11) NOT NULL,
	`seq_region_strand` tinyint(4) NOT NULL,
	`variation_name` varchar(255) DEFAULT NULL,
	`source_id` int(10) unsigned NOT NULL,
	`class` varchar(255) DEFAULT NULL,
	`bound_start` int(11) DEFAULT NULL,
	`bound_end` int(11) DEFAULT NULL,
	PRIMARY KEY (`structural_variation_id`),
	KEY `pos_idx` (`seq_region_id`,`seq_region_start`)
  )});
  
  # now copy the data
  $dbVar->do(qq{
	INSERT INTO structural_variation
	(seq_region_id, seq_region_start, seq_region_end, seq_region_strand,
	variation_name, source_id,
	class, bound_start, bound_end)
	SELECT q.seq_region_id, t.start, t.end, 1,
	t.id, s.source_id, t.type, t.outer_start, t.outer_end
	FROM seq_region q, temp_cnv t, source s
	WHERE q.name = t.chr AND concat('DGVa:', t.study) = s.name;
  });
  
  # cleanup
  $dbVar->do(qq{DROP TABLE temp_cnv;});
}


sub meta_coord{
  
  debug("Adding entry to meta_coord table");
  
  $dbVar->do(qq{INSERT INTO meta_coord(table_name, coord_system_id, max_length) VALUES ('structural_variation', 2, 500);});
}



sub mapping{
  debug("Mapping SV features to $target_assembly");
  
  open IN, $input_file or die "Could not read from input file $input_file\n";
  
  # get a slice adaptor
  my $sa = Bio::EnsEMBL::Registry->get_adaptor('human', 'core', 'slice');
  
  open OUT, ">$TMP_DIR/$TMP_FILE" or die "Could not write to output file $TMP_DIR/$TMP_FILE\n";
  
  my ($num_mapped, $num_not_mapped);
  
  while(<IN>) {
	chomp;
	
	# input file has these columns
	my ($study, $id, $chr, $outer_start, $inner_start, $start, $end, $inner_end, $outer_end, $assembly, $type) = split /\t/, $_;
	
	# skip header line
	next if /STUDY_ID/;
	
	#print "Assembly is $assembly\n";
	
	# try to guess the source assembly
	if($assembly =~ /35/) {
	  $assembly = 'NCBI35';
	}
	elsif($assembly =~ /36/) {
	  $assembly = 'NCBI36';
	}
	elsif($assembly =~ /37/) {
	  # no need to map if it's already on the latest assembly
	  print OUT "$_\n";
	  next;
	}
	else {
	  warn("Could not guess assembly from file - assuming to be NCBI35");
	  $assembly = 'NCBI35';
	}
	
	# get a slice on the old assembly
	my $slice = $sa->fetch_by_region('chromosome', $chr, $start, $end, 1, $assembly);
	
	my $count = 0;
	my ($min_start, $max_end) = (999999999, -999999999);
	my $to_chr;
	
	# project the slice to the latest assembly
	# this may produce several "segments"
	foreach my $segment(@{$slice->project('chromosome', $target_assembly)}) {
	  
	  # get the slice that this segment is on
	  my $to_slice = $segment->to_Slice();
	  
	  # check this segment is on the same chrom
	  if((defined $to_chr) && ($to_chr ne $to_slice->seq_region_name)) {
		$count = 0;
		warn("Segments of $id map to different chromosomes ($to_chr and ", $to_slice->seq_region_name, ")");
		last;
	  }
	  
	  # get the chromosome name
	  $to_chr = $to_slice->seq_region_name;
	  
	  # get the slice start and end
	  my $to_start = $to_slice->start;
	  my $to_end = $to_slice->end;
	  
	  # now calculate the coords of the segment on the slice
	  $to_start = ($to_start + $segment->from_start) - 1;
	  $to_end = ($to_end + $segment->from_start) - 1;
	  
	  # store the min/max
	  $min_start = $to_start if $to_start < $min_start;
	  $max_end = $to_end if $to_end > $max_end;
	  
	  # print this segment
	  #print "$id\t$chr\t$start\t$end\tsize ",($end - $start + 1), "\t\-\>\t", $to_chr, "\t", $to_start, "\t", $to_end, " size ", ($to_end - $to_start + 1), "\n";
	  
	  $count++;
	}
	
	# allow gaps?
	if(($count - 1) <= $num_gaps && $count) {
	  
	  $assembly = $target_assembly;
	  
	  # calculate inner/outer coords
	  $outer_start = $min_start - ($start - $outer_start) if $outer_start >= 1;
	  $inner_start = $min_start + ($inner_start - $start) if $inner_start >= 1;
	  $inner_end = $max_end - ($end - $inner_end) if $inner_end >= 1;
	  $outer_end = $max_end + ($outer_end - $end) if $outer_end >= 1;
	  
	  print OUT (join "\t", ($study, $id, $to_chr, $outer_start, $inner_start, $min_start, $max_end, $inner_end, $outer_end, $assembly, $type));
	  print OUT "\n";
	  
	  $num_mapped++;
	}
	
	else {
	  $num_not_mapped++;
	}
  }
  
  close IN;
  close OUT;
  
  debug("Finished mapping - mapped $num_mapped successfully, $num_not_mapped not mapped");
}


sub usage {
  die shift, qq{

Options -
  -species
  -tmpdir
  -tmpfile
  -input_file	: file containing EGA data dump (required)
  -mapping	: if set, the data will be mapped to $target_assembly using the Ensembl API
  -gaps		: number of gaps allowed in mapping (defaults to 1)
  
  };
}
