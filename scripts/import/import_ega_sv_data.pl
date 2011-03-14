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
		   		 'input_file=s'	   	 => \$input_file,
           'tmpdir=s'          => \$ImportUtils::TMP_DIR,
           'tmpfile=s'         => \$ImportUtils::TMP_FILE,
		   		 'mapping'		   		 => \$mapping,
		   		 'gaps=i'			   		 => \$num_gaps,
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

my $study_table = 'study';
my $sv_table = 'structural_variation';
my $ssv_table = 'supporting_structural_variation';

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
my $failed = [];
if($mapping) {
  $failed = mapping();
	print join(',',@{$failed})."\n";
}

# otherwise parse the file into the temp file
else {
  open IN, $input_file or die "Could not read from input file $input_file\n";
  open OUT, ">$TMP_DIR/$TMP_FILE" or die "Could not write to output file $TMP_DIR/$TMP_FILE\n";
  
  while(<IN>) {
	next if /STUDY_ID/;		# skip header line
	chomp;
	print OUT $_ . "\t\n";
  }
  
  close IN;
  close OUT;
}


read_file();
my $source_id = source();
study_table();
structural_variation($failed);
supporting_evidence();
meta_coord();

sub read_file{
  debug("Loading file into temporary table");
  
  # drop any table that may be there
  $dbVar->do("DROP TABLE IF EXISTS temp_cnv;");
  
  create_and_load(
	$dbVar, "temp_cnv", "study", "pmid i", "author", "year", "title l", "id *", "tax_id i", "organism", "chr", "outer_start i", "inner_start i",
	"start i", "end i", "inner_end i", "outer_end i", "assembly", "type", "ssv l", "comment");
  
  # fix nulls
  foreach my $coord('outer_start', 'inner_start', 'inner_end', 'outer_end') {
	$dbVar->do(qq{UPDATE temp_cnv SET $coord = NULL WHERE $coord = 0;});
  }
  
  # remove those that are of different species
  $dbVar->do(qq{DELETE FROM temp_cnv WHERE tax_id != $connected_tax_id;});
  $dbVar->do(qq{DELETE FROM temp_cnv WHERE chr like '%_random';});
}


sub source{
  debug("Inserting into source table");
	
	# Check if the DGVa source already exists, else it create the entry
	if ($dbVar->selectrow_arrayref(qq{SELECT source_id FROM source WHERE name='DGVa';})) {
		$dbVar->do(qq{UPDATE IGNORE source SET description='Database of Genomic Variants Archive',url='http://www.ebi.ac.uk/dgva/',version=201101 where name='DGVa';});
	}
	else {
		$dbVar->do(qq{INSERT INTO source (name,description,url,version) VALUES ('DGVa','Database of Genomic Variants Archive','http://www.ebi.ac.uk/dgva/',201102);});
	}
	my @source_id = @{$dbVar->selectrow_arrayref(qq{SELECT source_id FROM source WHERE name='DGVa';})};
	return $source_id[0];
}

sub study_table{
  debug("Inserting into study table");
  
  # An ugly construct, but let's create a unique key on the study name in order to avoid duplicates but still update the URLs
  # BTW, perhaps we actually should have a unique constraint on the name column??
  my $stmt = qq{
    ALTER TABLE
      $study_table
    ADD CONSTRAINT
      UNIQUE KEY
	name_key (name)
  };
  $dbVar->do($stmt);
  
  # Then insert "new" studys. Update the URL field in case of duplicates
  $stmt = qq{
    INSERT INTO
      $study_table (
	name,
	description,
	url,
	source_id,
	external_reference
      )
    SELECT DISTINCT
      CONCAT(
	study
      ),
      CONCAT(
	author,
	' ',
	year,
	' "',
	title,
	'" PMID:',
	pmid
      ),
      CONCAT(
	'ftp://ftp.ebi.ac.uk/pub/databases/dgva/',
	study,
	'_',
	author,
	'_et_al_',
	year
      ),
			$source_id,
			CONCAT(
			'pubmed/',
			pmid) 
    FROM
      temp_cnv
    ORDER BY
      LENGTH(
	SUBSTR(
	  study,
	  4
	)
      ),
      SUBSTR(
	study,
	4
      )
    ON DUPLICATE KEY UPDATE
      url = CONCAT(
	'http://www.ncbi.nlm.nih.gov/pubmed/',
	pmid
      )
  };
  $dbVar->do($stmt);
  
  # And then drop the name key again
  $stmt = qq{
    ALTER TABLE
      $study_table
    DROP KEY
      name_key
  };
  $dbVar->do($stmt);
  
}


sub structural_variation{
  my $failed = shift;
  
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
	`study_id` int(10) unsigned NOT NULL,
	`class` varchar(255) DEFAULT NULL,
	`bound_start` int(11) DEFAULT NULL,
	`bound_end` int(11) DEFAULT NULL,
	PRIMARY KEY (`structural_variation_id`),
	KEY `pos_idx` (`seq_region_id`,`seq_region_start`),
	KEY `name_idx` (`variation_name`),
	KEY `study_idx` (`study_id`) 
  )});
  
  # The variation name should be unique so for the insert, create a unique key for this column.
  # Should there perhaps in fact be a unique constraint on this column?
	my $stmt = qq{
    ALTER TABLE
      $sv_table
    ADD CONSTRAINT
      UNIQUE KEY
	name_key (variation_name)
  };
  $dbVar->do($stmt);
  
  # now copy the data. If the variation is duplicated, replace the old data
  $stmt = qq{
    REPLACE INTO
      $sv_table (
	seq_region_id,
	seq_region_start,
	seq_region_end,
	seq_region_strand,
	variation_name,
	source_id,
	study_id,
	class,
	bound_start,
	bound_end
      )
    SELECT
      q.seq_region_id,
      t.start,
      t.end,
      1,
      t.id,
      s.source_id,
			st.study_id,
      t.type,
      t.outer_start,
      t.outer_end
    FROM
      seq_region q,
      temp_cnv t,
      source s,
			$study_table st 
    WHERE
      q.name = t.chr AND
      st.source_id=$source_id AND
			st.source_id=s.source_id AND
			st.name=t.study
  };
  $dbVar->do($stmt);
  
  # cleanup
  $stmt = qq{
    ALTER TABLE
      $sv_table
    DROP KEY
      name_key
  };
  $dbVar->do($stmt);

	#	Warn about variations that are on chromosomes for which we don't have a seq region entry. Also, delete them from the structural variation table if they are already present there
  $stmt = qq{
    SELECT
      t.id,
      t.chr,
      t.study,
      t.author,
      t.year
    FROM
      temp_cnv t LEFT JOIN
      seq_region sr ON (
				sr.name = t.chr
      )
    WHERE
      sr.seq_region_id IS NULL
  };
  my $rows = $dbVar->selectall_arrayref($stmt);
	
  while (my $row = shift(@{$rows})) {
    warn ("Structural variant '" . $row->[0] . "' from study '" . $row->[2] . "' (" . $row->[3] . " (" . $row->[4] . ")) is placed on chromosome '" . $row->[1] . "', which could not be found in Ensembl. It will not be imported, if it is already imported with a different (and outdated) location, it will be removed from the database");
    push(@{$failed},qq{'$row->[0]'});
  }
  
  # Remove any variants that have an updated position that could not be determined
	my $condition = join(",",@{$failed});
 	if ($condition ne '') {
		$stmt = qq{
    	DELETE FROM
    	  sv
    	USING
      	$sv_table sv
    	WHERE
      	sv.variation_name IN ( $condition )
  	};
  	$dbVar->do($stmt);
	}
}


sub supporting_evidence {
	
	debug("Inserting into supporting_structural_variation table");
	
	my $stmt = qq{
    SELECT
      sv.structural_variation_id, 
      t.ssv 
    FROM
      temp_cnv t, $sv_table sv
    WHERE
      sv.variation_name=t.id
  };
	$dbVar->do($stmt);
	my %count= ();
  my $rows = $dbVar->selectall_arrayref($stmt);
  while (my $row = shift(@{$rows})) {
		my $sv_id = $row->[0];
		my $ssv_ids = $row->[1];
		my @ssv = split(':',$ssv_ids);
		foreach my $ssv_name (@ssv) {
		# now copy the data. If the variation is duplicated, replace the old data
			$dbVar->do(qq{INSERT INTO $ssv_table (name,structural_variation_id) VALUES ('$ssv_name',$sv_id);});
		}
	}
  
  $dbVar->do($stmt);
  $dbVar->do(qq{DROP TABLE temp_cnv;});
}


sub meta_coord{
  
  debug("Adding entry to meta_coord table");
  
  my $max_length_ref = $dbVar->selectall_arrayref(qq{
    SELECT GREATEST(
      (
	SELECT
	  MAX(seq_region_end - seq_region_start + 1)
	FROM $sv_table
      ),
      (
	SELECT
	  MAX(bound_end - bound_start + 1)
	FROM $sv_table
	     )
    )
  });
  my $max_length = $max_length_ref->[0][0];
	if (!defined($max_length)) { $max_length = 0; }
	
  my $cs = $default_cs->dbID;
  
  $dbVar->do(qq{DELETE FROM meta_coord where table_name = 'structural_variation';});
  print "INSERT INTO meta_coord(table_name, coord_system_id, max_length) VALUES ('structural_variation', $cs, $max_length);\n";
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
  
  # Store and return the variant ids that we could not map. In case these already exist in the database with a different (and thus assumingly outdated) position, we will delete them
  my @failed = ();
  
  while(<IN>) {
	chomp;
	
	# input file has these columns
	my ($study, $pmid, $author, $year, $title, $id, $tax_id, $organism, $chr, $outer_start, $inner_start, $start, $end, $inner_end, $outer_end, $assembly, $type, $ssv) = split /\t/, $_;
	
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
		print OUT "$_\t\n";
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
		print OUT "$_\t\n";
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
	  
	  # calculate inner/outer coords
	  $outer_start = $min_start - ($start - $outer_start) if $outer_start >= 1;
	  $inner_start = $min_start + ($inner_start - $start) if $inner_start >= 1;
	  $inner_end = $max_end - ($end - $inner_end) if $inner_end >= 1;
	  $outer_end = $max_end + ($outer_end - $end) if $outer_end >= 1;
	  
	  print OUT (join "\t", ($study, $pmid, $author, $year, $title, $id, $tax_id, $organism, $to_chr, $outer_start, $inner_start, $min_start, $max_end, $inner_end, $outer_end, $target_assembly, $type, $ssv, "[remapped from build $assembly]"));
	  print OUT "\n";
	}
	
	else {
	  warn ("Structural variant '$id' from study '$study' ($author ($year)) has location '$assembly:$chr:$start\-$end' , which could not be re-mapped to $target_assembly. This variant will not be imported, if it is already present with a different (and outdated) location, it will be removed from the database");
	  push(@failed,qq{'$id'});
	  $num_not_mapped{$assembly}++;
	}
  }
  
  close IN;
  close OUT;
  
  debug("Finished mapping\n\tSuccess: ".(join " ", %num_mapped)." not required $no_mapping_needed\n\tFailed: ".(join " ", %num_not_mapped)." skipped $skipped");
  
  return \@failed;
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
