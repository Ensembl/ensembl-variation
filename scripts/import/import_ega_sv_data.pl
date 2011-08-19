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
use LWP::Simple;

our ($species, $input_file, $source_name, $TMP_DIR, $TMP_FILE, $mapping, $num_gaps, $target_assembly, $size_diff, $version);

GetOptions('species=s'         => \$species,
		  		 'source_name=s'     => \$source_name,
		   		 'input_file=s'	   	 => \$input_file,
           'tmpdir=s'          => \$ImportUtils::TMP_DIR,
           'tmpfile=s'         => \$ImportUtils::TMP_FILE,
		   		 'mapping'		   		 => \$mapping,
		   		 'gaps=i'			   		 => \$num_gaps,
		   		 'target_assembly=s' => \$target_assembly,
		   		 'size_diff=i'       => \$size_diff,
					 'version=i'         => \$version,
          );
my $registry_file ||= $Bin . "/ensembl.registry";

$num_gaps = 1 unless defined $num_gaps;
$species ||= 'human';

usage('-species argument is required')    if(!$species);
usage('-input_file argument is required') if (!$input_file);
usage('-version argument is required')    if (!$version);

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

my $add         = '';
my $study_table = "study$add";
my $sv_table    = "structural_variation$add";
my $ssv_table   = "supporting_structural_variation$add";
my $svf_table   = "structural_variation_feature$add";

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
my %studies = {};
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
verifications();

sub read_file{
  debug("Loading file into temporary table");
  
  # drop any table that may be there
  $dbVar->do("DROP TABLE IF EXISTS temp_cnv;");
  
  create_and_load(
	$dbVar, "temp_cnv", "study", "pmid i", "author", "year", "title l", "organism", "id *", "tax_id i", "type", "chr", "outer_start i", "start i", "inner_start i",
	"inner_end i", "end i", "outer_end i", "assembly", "ssv l", "comment");
	 
	# remove those that are of different species
	$dbVar->do(qq{DELETE FROM temp_cnv WHERE tax_id != $connected_tax_id;});
	if ($species eq 'human') {
		$dbVar->do(qq{DELETE FROM temp_cnv WHERE chr like '%_random';});
	}
	
  # fix nulls
	foreach my $coord('outer_start', 'inner_start', 'inner_end', 'outer_end', 'start', 'end') {
		$dbVar->do(qq{UPDATE temp_cnv SET $coord = NULL WHERE $coord = 0;});
	}
							
	$dbVar->do(qq{UPDATE temp_cnv SET start=outer_start WHERE start is NULL;});
	$dbVar->do(qq{UPDATE temp_cnv SET end=outer_end WHERE end is NULL;});	
	
	$dbVar->do(qq{UPDATE temp_cnv SET start=inner_start WHERE start is NULL;});
	$dbVar->do(qq{UPDATE temp_cnv SET end=inner_end WHERE end is NULL;});
											
	# Case with insertions						
	$dbVar->do(qq{UPDATE temp_cnv SET start=outer_start, end=inner_start, 
	              outer_end=inner_start, inner_start=NULL, inner_end=NULL
								WHERE outer_start=outer_end AND inner_start=inner_end;});																				
}


sub source{
  debug("Inserting into source table");
	
	# Check if the DGVa source already exists, else it create the entry
	if ($dbVar->selectrow_arrayref(qq{SELECT source_id FROM source WHERE name='DGVa';})) {
		$dbVar->do(qq{UPDATE IGNORE source SET description='Database of Genomic Variants Archive',url='http://www.ebi.ac.uk/dgva/',version=$version where name='DGVa';});
	}
	else {
		$dbVar->do(qq{INSERT INTO source (name,description,url,version) VALUES ('DGVa','Database of Genomic Variants Archive','http://www.ebi.ac.uk/dgva/',$version);});
	}
	my @source_id = @{$dbVar->selectrow_arrayref(qq{SELECT source_id FROM source WHERE name='DGVa';})};
	return $source_id[0];
}

sub study_table{
  debug("Inserting into $study_table table");
  
	my $stmt;
	
	# Checks the studies with different assemblies when the option "mapping" is selected.
	if ($mapping) {
		foreach my $st (keys(%studies)) {
			next if (!$studies{$st});
			my $count_assemblies = scalar(@{$studies{$st}});
			next if ($count_assemblies < 2);
			my $assemblies = '';
			my $c_loop = 0;
			foreach my $a_version (@{$studies{$st}}) {
				if ($c_loop == $count_assemblies-1) {
					$assemblies .= ' and ';
				}
				elsif ($assemblies ne '') {
					$assemblies .= ', ';
				}
				$assemblies .= $a_version;
				$c_loop ++;
			}
			$stmt = qq{
				UPDATE temp_cnv SET comment='[remapped from builds $assemblies]' WHERE title='$st'
			};
			$dbVar->do($stmt);
		}
	}
	
	
  # An ugly construct, but let's create a unique key on the study name in order to avoid duplicates but still update the URLs
  # BTW, perhaps we actually should have a unique constraint on the name column??
  
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
	pmid,
	' ',
	comment
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
  
	# Special case for the De Smith url study
	if ($species eq 'human') {
		my @s_data = @{$dbVar->selectrow_arrayref(qq{SELECT study_id,url FROM $study_table WHERE name='estd24';})};
		$s_data[1] =~ s/de\sSmith/De_Smith/;
		$dbVar->do(qq{UPDATE $study_table SET url='$s_data[1]' where study_id=$s_data[0];});
	}
}


sub structural_variation{
  my $failed = shift;
  
  debug("Inserting into $sv_table table");
  
	if ($dbVar->do(qq{show columns from $sv_table like 'tmp_class_name';}) != 1){
		$dbVar->do(qq{ALTER TABLE $sv_table ADD COLUMN tmp_class_name varchar(255);});
	}
  
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
  
	# Structural variations
	$stmt = qq{
    REPLACE INTO
    $sv_table (
	    variation_name,
	    source_id,
	    study_id,
	    tmp_class_name
    )
    SELECT
      t.id,
      s.source_id,
			st.study_id,
      t.type
    FROM
      temp_cnv t,
      source s,
			$study_table st 
    WHERE
			st.source_id=$source_id AND
			st.source_id=s.source_id AND
			st.name=t.study
  };
  $dbVar->do($stmt);
	
	
	# Structural variation features
	debug("Inserting into $svf_table table");
	
	$stmt = qq{
    INSERT INTO $svf_table (
	    seq_region_id,
	    outer_start,
	    seq_region_start,
	    inner_start,
	    inner_end,
	    seq_region_end,
	    outer_end,
	    seq_region_strand,
	    structural_variation_id,
	    variation_name,
	    source_id
    )
    SELECT
      q.seq_region_id,
			t.outer_start,
      t.start,
			t.inner_start,
			t.inner_end,
      t.end,
			t.outer_end,
      1,
			sv.structural_variation_id,
      t.id,
      $source_id
    FROM
      seq_region q,
      temp_cnv t,
			$sv_table sv
    WHERE
      q.name = t.chr AND
			t.id=sv.variation_name AND
			sv.source_id=$source_id
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
      seq_region sr ON ( sr.name = t.chr )
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
	
	debug("Inserting into $ssv_table table");
	
	if ($dbVar->do(qq{show columns from $ssv_table like 'tmp_class_name';}) != 1){
		$dbVar->do(qq{ALTER TABLE $ssv_table ADD COLUMN tmp_class_name varchar(255);});
	}
	
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
		foreach my $ssv_data (@ssv) {
			next if ($ssv_data eq 'None');
			my ($ssv_name, $ssv_class) = split(/\|/,$ssv_data);
			my $exists_ssv = $dbVar->selectall_arrayref(qq{
				SELECT supporting_structural_variation_id FROM $ssv_table WHERE structural_variation_id=$sv_id AND name='$ssv_name' LIMIT 1
  		});
			next if (defined($exists_ssv->[0][0]));
			# now copy the data. If the variation is duplicated, replace the old data
			$dbVar->do(qq{INSERT INTO $ssv_table (name,structural_variation_id,tmp_class_name) VALUES ('$ssv_name',$sv_id,'$ssv_class');});
		}
	}
  
  $dbVar->do($stmt);
  $dbVar->do(qq{DROP TABLE temp_cnv;});
}


sub meta_coord{
  
  debug("Adding entry to meta_coord table");
  
  my $max_length_ref = $dbVar->selectall_arrayref(qq{
	SELECT
	  MAX(seq_region_end - seq_region_start + 1)
	FROM $svf_table
  });
  my $max_length = $max_length_ref->[0][0];
	if (!defined($max_length)) { $max_length = 0; }
	
  my $cs = $default_cs->dbID;
  
  $dbVar->do(qq{DELETE FROM meta_coord where table_name = 'structural_variation_feature';});
  print "INSERT INTO meta_coord(table_name, coord_system_id, max_length) VALUES ('structural_variation_feature', $cs, $max_length);\n";
	$dbVar->do(qq{INSERT INTO meta_coord(table_name, coord_system_id, max_length) VALUES ('structural_variation_feature', $cs, $max_length);});
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
  
	my $sth = $dbVar->prepare(qq{ SELECT seq_region_id FROM seq_region WHERE name=?});
	
	
  while(<IN>) {
	chomp;
	
	# input file has these columns
	my ($study, $pmid, $author, $year, $title, $organism, $id, $tax_id, $type, $chr, $outer_start, $start, $inner_start, $inner_end, $end, $outer_end, $assembly, $ssv) = split /\t/, $_;
	
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
		# Check if the chromosome is stored into the variation database (seq_region table)
		$sth->execute($chr);
		if(!defined(($sth->fetchrow_array)[0])) {
			warn("Unable to find the chromosome '$chr' in the variation database, for the structural variant '$id' from study '$study' ($author ($year))");
	  	$skipped++;
		}
		else {
			# no need to map if it's already on the latest assembly
			print OUT "$_\t\n";
			$no_mapping_needed++;
		}
		next;
	  }
	}
	
	# get a slice on the old assembly
	my $slice;
	my $start_c = $start;
	my $end_c   = $end;
	
	if (!$start) {
		if ($outer_start) {
			$start_c = $outer_start;
		}
		else {
			$start_c = $inner_start;
		}
	}
	if (!$end) {
		if ($outer_end) {
			$end_c = $outer_end;
		}
		else {
			$end_c = $inner_end;
		}
	}
	
	eval { $slice = $sa->fetch_by_region('chromosome', $chr, $start_c, $end_c, 1, $assembly); };
	
	# check got the slice OK
	if(!defined($slice)) {
	  warn("Unable to map from assembly $assembly, or unable to retrieve slice $chr\:$start_c\-$end_c");
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
	
	#my $diff = abs(100 - (100 * (($end - $start + 1) / ($max_end - $min_start + 1))));
	
	#print "Before ",($end - $start + 1), " After ", ($max_end - $min_start + 1), " Diff $diff count $count\n";
	
	# allow gaps?
	if((($count - 1) <= $num_gaps && $count)) {# && $diff < $size_diff) {
	  
	  $num_mapped{$assembly}++;
	  
	  # calculate inner/outer coords
	  $outer_start = $min_start - ($start_c - $outer_start) if $outer_start >= 1;
	  $inner_start = $min_start + ($inner_start - $start_c) if $inner_start >= 1;
	  $inner_end = $max_end - ($end_c - $inner_end) if $inner_end >= 1;
	  $outer_end = $max_end + ($outer_end - $end_c) if $outer_end >= 1;
		
		if ($studies{$title}) {
			push (@{$studies{$title}}, $assembly) if (!grep{$_ eq $assembly} @{$studies{$title}});
		}
		else {
			$studies{$title} = [$assembly];
		}
		
		
	  print OUT (join "\t", ($study, $pmid, $author, $year, $title, $organism, $id, $tax_id, $type, $to_chr, $outer_start, $min_start, $inner_start, $inner_end, $max_end, $outer_end, $target_assembly, $ssv, "[remapped from build $assembly]"));
	  print OUT "\n";
	}
	
	else {
	  warn ("Structural variant '$id' from study '$study' ($author ($year)) has location '$assembly:$chr:$start_c\-$end_c' , which could not be re-mapped to $target_assembly. This variant will not be imported, if it is already present with a different (and outdated) location, it will be removed from the database");
	  push(@failed,qq{'$id'});
	  $num_not_mapped{$assembly}++;
	}
  }
  
  close IN;
  close OUT;
  
  debug("Finished mapping\n\tSuccess: ".(join " ", %num_mapped)." not required $no_mapping_needed\n\tFailed: ".(join " ", %num_not_mapped)." skipped $skipped");
  
  return \@failed;
}


sub verifications {
	debug("Verification of ftp links");
	my $sth = $dbVar->prepare(qq{ SELECT url FROM $study_table WHERE source_id=$source_id});
	$sth->execute();
	my $failed_flag = 0;
	while (my $url = ($sth->fetchrow_array)[0]) {
		if (!head($url)) {
			print STDERR "$url IS NOT VALID\n";
			$failed_flag = 1;
		}
	}
	print "URL verification: OK\n" if ($failed_flag == 0);
}


sub usage {
  die shift, qq{

Options -
  -species         : species name (required)
  -target_assembly : assembly version to map to
  -tmpdir
  -tmpfile
  -input_file      : file containing EGA data dump (required)
  -mapping         : if set, the data will be mapped to $target_assembly using the Ensembl API
  -gaps            : number of gaps allowed in mapping (defaults to 1)
  -size_diff       : % difference allowed in size after mapping
	-version         : version number of the data (required)
  
  };
}
