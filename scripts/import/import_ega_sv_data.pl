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
our ($species, $input_file, $source_name, $TMP_DIR, $TMP_FILE, $mapping, $num_gaps, $target_assembly, $size_diff);

GetOptions('species=s'         => \$species,
		   'source_name=s'     => \$source_name,
		   'input_file=s'	   => \$input_file,
           'tmpdir=s'          => \$ImportUtils::TMP_DIR,
           'tmpfile=s'         => \$ImportUtils::TMP_FILE,
		   'mapping'		   => \$mapping,
		   'gaps=i'			   => \$num_gaps,
		   'target_assembly=s' => \$target_assembly,
		   'size_diff=i'       => \$size_diff,
          );
my $registry_file ||= $Bin . "/ensembl.registry";

$num_gaps = 1 unless defined $num_gaps;
$species ||= 'human';

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

my $csa = Bio::EnsEMBL::Registry->get_adaptor($species, "core", "coordsystem");
our $default_cs = $csa->fetch_by_name("chromosome");


# set the target assembly
$target_assembly ||= $default_cs->version;

# get the tax_id for the connected species
my $meta_container = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Core', 'MetaContainer' );
our $connected_tax_id = $meta_container->get_taxonomy_id;

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
	$dbVar, "temp_cnv", "study", "id *", "tax_id i", "organism", "chr", "outer_start i", "inner_start i",
	"start i", "end i", "inner_end i", "outer_end i", "assembly", "type");
  
  # fix nulls
  foreach my $coord('outer_start', 'inner_start', 'inner_end', 'outer_end') {
	$dbVar->do(qq{UPDATE temp_cnv SET $coord = NULL WHERE $coord = 0;});
  }
  
  # remove those that are of different species
  $dbVar->do(qq{DELETE FROM temp_cnv WHERE tax_id != $connected_tax_id;});
}

sub source{
  debug("Inserting into source table");
  $dbVar->do(qq{INSERT INTO source (name) SELECT distinct(concat('DGVa:', study)) FROM temp_cnv ORDER BY length(substr(study, 4)), substr(study, 4);});
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
  
  my $max_length_ref = $dbVar->selectall_arrayref(qq{SELECT max(seq_region_end - seq_region_start+1) from structural_variation});
  my $max_length = $max_length_ref->[0][0];
  
  my $cs = $default_cs->dbID;
  
  $dbVar->do(qq{DELETE FROM meta_coord where table_name = 'structural_variation';});
  $dbVar->do(qq{INSERT INTO meta_coord(table_name, coord_system_id, max_length) VALUES ('structural_variation', $cs, $max_length);});
}



sub mapping{
  debug("Mapping SV features to $target_assembly");
  
  open IN, $input_file or die "Could not read from input file $input_file\n";
  
  # get a slice adaptor
  my $sa = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'slice');
  
  open OUT, ">$TMP_DIR/$TMP_FILE" or die "Could not write to output file $TMP_DIR/$TMP_FILE\n";
  
  my (%num_mapped, %num_not_mapped, $no_mapping_needed, $skipped);
  $no_mapping_needed = 0;
  $skipped = 0;
  
  # get the default CS version
  my $cs_version_number = $target_assembly;
  $cs_version_number =~ s/\D//g;
  
  while(<IN>) {
	chomp;
	
	# input file has these columns
	my ($study, $id, $tax_id, $organism, $chr, $outer_start, $inner_start, $start, $end, $inner_end, $outer_end, $assembly, $type) = split /\t/, $_;
	
	next unless $tax_id == $connected_tax_id;
	
	# skip header line
	next if /STUDY_ID/;
	
	#print "Assembly is $assembly\n";
	
	if($species =~ /human|homo/) {
	  # try to guess the source assembly
	  if($assembly =~ /34/) {
		$assembly = 'NCBI34';
	  }
	  elsif($assembly =~ /35/) {
		$assembly = 'NCBI35';
	  }
	  elsif($assembly =~ /36/) {
		$assembly = 'NCBI36';
	  }
	  elsif($assembly =~ /$cs_version_number/) {
		# no need to map if it's already on the latest assembly
		print OUT "$_\n";
		next;
	  }
	  else {
		warn("Could not guess assembly from file - assuming to be NCBI35");
		$assembly = 'NCBI35';
	  }
	}
	
	if($species =~ /mouse|mus.*mus/) {
	  if($assembly =~ /34/) {
		$assembly = 'NCBIM34';
	  }
	  if($assembly =~ /35/) {
		$assembly = 'NCBIM35';
	  }
	  if($assembly =~ /36/) {
		$assembly = 'NCBIM36';
	  }
	  elsif($assembly =~ /$cs_version_number/) {
		# no need to map if it's already on the latest assembly
		print OUT "$_\n";
		$no_mapping_needed++;
		next;
	  }
	}
	
	# get a slice on the old assembly
	my $slice;
	
	eval { $slice = $sa->fetch_by_region('chromosome', $chr, $start, $end, 1, $assembly); };
	
	# check got the slice OK
	if(!defined($slice)) {
	  warn("Unable to map from assembly $assembly, or unable to retrieve slice $chr\:$start\-$end");
	  $skipped++;
	  next;
	}
	
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
	
	my $diff = abs(100 - (100 * (($end - $start + 1) / ($max_end - $min_start + 1))));
	
	#print "Before ",($end - $start + 1), " After ", ($max_end - $min_start + 1), " Diff $diff count $count\n";
	
	# allow gaps?
	if((($count - 1) <= $num_gaps && $count)) {# && $diff < $size_diff) {
	  
	  $num_mapped{$assembly}++;
	  
	  $assembly = $target_assembly;
	  
	  # calculate inner/outer coords
	  $outer_start = $min_start - ($start - $outer_start) if $outer_start >= 1;
	  $inner_start = $min_start + ($inner_start - $start) if $inner_start >= 1;
	  $inner_end = $max_end - ($end - $inner_end) if $inner_end >= 1;
	  $outer_end = $max_end + ($outer_end - $end) if $outer_end >= 1;
	  
	  print OUT (join "\t", ($study, $id, $tax_id, $organism, $to_chr, $outer_start, $inner_start, $min_start, $max_end, $inner_end, $outer_end, $assembly, $type));
	  print OUT "\n";
	}
	
	else {
	  $num_not_mapped{$assembly}++;
	}
  }
  
  close IN;
  close OUT;
  
  debug("Finished mapping\n\tSuccess: ".(join " ", %num_mapped)." not required $no_mapping_needed\n\tFailed: ".(join " ", %num_not_mapped)." skipped $skipped");
}


sub usage {
  die shift, qq{

Options -
  -species
  -target_assembly : assembly version to map to
  -tmpdir
  -tmpfile
  -input_file      : file containing EGA data dump (required)
  -mapping         : if set, the data will be mapped to $target_assembly using the Ensembl API
  -gaps            : number of gaps allowed in mapping (defaults to 1)
  -size_diff       : % difference allowed in size after mapping
  
  };
}
