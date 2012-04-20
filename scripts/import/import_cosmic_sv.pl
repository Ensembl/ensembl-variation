#! /usr/local/bin/perl
#

use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Data::Dumper;
use FindBin qw( $Bin );
use Getopt::Long;
use ImportUtils qw(dumpSQL debug create_and_load load);
use LWP::Simple;

our ($species, $input_file, $input_dir, $source_name, $TMP_DIR, $TMP_FILE, $var_set_id, $mapping, 
     $num_gaps, $target_assembly, $cs_version_number, $size_diff, $version, $registry_file);

GetOptions('species=s'         => \$species,
		  		 'source_name=s'     => \$source_name,
		   		 'input_file=s'	   	 => \$input_file,
					 'input_dir=s'	   	 => \$input_dir,
           'tmpdir=s'          => \$ImportUtils::TMP_DIR,
           'tmpfile=s'         => \$ImportUtils::TMP_FILE,
					 'var_set=i'         => \$var_set_id,
		   		 'mapping'		   		 => \$mapping,
		   		 'gaps=i'			   		 => \$num_gaps,
		   		 'target_assembly=s' => \$target_assembly,
		   		 'size_diff=i'       => \$size_diff,
					 'version=i'         => \$version,
					 'registry=s'        => \$registry_file
          );
$registry_file ||= $Bin . "/ensembl.registry";

$num_gaps = 1 unless defined $num_gaps;
$species ||= 'human';

usage('-species argument is required')    if (!$species);
usage('-input_file or -input_dir argument is required') if (!$input_file and !$input_dir);
usage('-version argument is required')    if (!$version);

my @files;
my $f_count = 1;
my $f_nb    = 1;
if ($input_dir) {
	opendir(DIR, $input_dir) or die $!;
	@files = readdir(DIR); 
	@files = sort(@files);
	close(DIR);
	$f_nb = (scalar(@files)-2);
}
else {
	@files = ($input_file);
}

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

my $clinical_attrib_type = 'dgva_clin_sig';
my $cosmic_url = 'http://www.sanger.ac.uk/genetics/CGP/cosmic/';

my $add          = '';
my $study_table  = "study$add";
my $sv_table     = "structural_variation$add";
my $svf_table    = "structural_variation_feature$add";
my $sva_table    = "structural_variation_association$add";
my $sv_failed    = "failed_structural_variation$add";
my $set_table    = "variation_set_structural_variation$add";
my $svanno_table = "structural_variation_annotation$add";

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

# get the default CS version
$cs_version_number = $target_assembly;
$cs_version_number =~ s/\D//g;

# get the tax_id for the connected species
my $meta_container = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Core', 'MetaContainer' );
our $connected_tax_id = $meta_container->get_taxonomy_id;


# run the mapping sub-routine if the data needs mapping
my (%num_mapped, %num_not_mapped);
my $no_mapping_needed = 0;
my $skipped = 0;
my $failed = [];
my $study_id;

my $source_id = source();


##########
#  Main  #
##########

foreach my $in_file (@files) {
	next if ($in_file !~ /\.gvf$/);
	
	my $msg ="File: $in_file ($f_count/$f_nb)";
	print "$msg\n";
	debug("##########\n$msg\n##########");
	
	my $fname;
	if ($input_dir) {
		$fname = "$input_dir/$in_file";
	}
	else {
		$fname = $in_file;
	}
	
	%num_mapped = ();
	%num_not_mapped = ();
	$no_mapping_needed = 0;
	$skipped = 0;
	$failed = [];
	$study_id;

	parse_gvf($fname);
	read_file();
	structural_variation();
	structural_variation_annotation();
	debug(localtime()." Done!\n");
	
	$f_count ++;
}
# Post processing for mouse annotation (duplicated lines)
if ($species =~ /mouse|mus/i) {
	annotation_post_processing();
}
meta_coord();
verifications();
debug(localtime()." All done!");


#############
#  Methods  #
#############

sub read_file{
  debug(localtime()." Loading file into temporary table");
  
  # drop any table that may be there
  $dbVar->do("DROP TABLE IF EXISTS temp_cnv;");
  
  create_and_load(
	$dbVar, "temp_cnv", "id *", "type", "chr", "outer_start i", "start i", "inner_start i",
	"inner_end i", "end i", "outer_end i", "strand", "parent", "clinic", "phenotype", "sample", "strain",
  "bp_order i", "is_somatic i", "is_ssv i", "is_failed i");		
	
	
  # fix nulls
	foreach my $coord('outer_start', 'inner_start', 'inner_end', 'outer_end', 'start', 'end') {
		$dbVar->do(qq{UPDATE temp_cnv SET $coord = NULL WHERE $coord = 0;});
	}	
	foreach my $col('clinic', 'phenotype', 'sample', 'strain') {
		$dbVar->do(qq{UPDATE temp_cnv SET $col = NULL WHERE $col = '';});
	}
	
	# Case with insertions						
	$dbVar->do(qq{UPDATE temp_cnv SET start=outer_start, end=inner_start, 
	              outer_end=inner_start, inner_start=NULL, inner_end=NULL
								WHERE outer_start=outer_end AND inner_start=inner_end;});	
	# Case with insertions						
	$dbVar->do(qq{UPDATE temp_cnv SET bp_order=null 
								WHERE bp_order=0;});																												
}


sub source{
  debug(localtime()." Inserting into source table");
	
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
	my $data = shift;
	
  debug(localtime()." Inserting into $study_table table");
  
	my $stmt;
	my $study      = $data->{study};
	my $assembly   = $data->{assembly};
	my $study_desc = $data->{desc};
	
	$stmt = qq{ SELECT st.study_id, st.description FROM study st, source s 
		          WHERE s.source_id=st.source_id AND s.name='DGVa' and st.name='$study'};
	my $rows = $dbVar->selectall_arrayref($stmt);		
	
	# UPDATE
	if (scalar (@$rows)) {
		$study_id   = $rows->[0][0];
		$study_desc = $rows->[0][1];
		# Checks the studies with different assemblies when the option "mapping" is selected.
		if ($mapping and $assembly ne $target_assembly) {
			if ($study_desc =~ /^(.+)\[remapped\sfrom\sbuild(s?)\s(.+)\]$/) {
				my $gen_desc = $1;
				my $assembly_desc = $3;
			
				if ($assembly_desc !~ /$assembly/) {
					$assembly_desc .= " and $assembly";
					$stmt = qq{
						UPDATE study SET description='$gen_desc [remapped from builds $assembly_desc]' 
						WHERE study_id=$study_id
					};
					$dbVar->do($stmt);
				}
			}
			else {
				$stmt = qq{
						UPDATE study SET description='$study_desc [remapped from build $assembly]' 
						WHERE study_id=$study_id
				};
				$dbVar->do($stmt);
			}
		}
	}
	# INSERT
	else {
		my $assembly_desc;
		my $pmid         = $data->{pubmed};
		my $first_author = $data->{first_author};
		my $year         = $data->{year};
		my $author       = $data->{author};
		my $author_info  = $data->{author};
		
		$author_info =~ s/_/ /g;
		$author_info =~ /(.+)\set\sal\s(.+)/;
		$first_author = $1 if (!$first_author);
		$year = $2 if (!$year);
		
		my $author_desc = "$first_author $year";
		
		if (length($study_desc)>150) {
			$study_desc = substr($study_desc,0,150);
			$study_desc .= '...';
		}
		$study_desc .= "$author_info $study_desc" if (defined($author_info));
	
		if ($mapping and $assembly ne $target_assembly) {
			$assembly_desc = " [remapped from build $assembly]";
		}
		
		my $pmid_desc = ($pmid) ? "PMID:$pmid" : '';
		my $pubmed    = ($pmid) ? "pubmed/$pmid" : 'NULL';
		
		$stmt = qq{
    	INSERT INTO `$study_table` (
				`name`,
				`description`,
				`url`,
				`source_id`,
				`external_reference`
      )
    	VALUES (
      	'$study',
				CONCAT(
					'$study_desc',
					' $pmid_desc',
					'$assembly_desc'),
				'$cosmic_url',
				$source_id,
				'$pubmed'
			)
		};
  	$dbVar->do($stmt);
		$study_id = $dbVar->{'mysql_insertid'};
	}
}


sub structural_variation{
  
  debug(localtime()." Inserting into $sv_table table");
  
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
  
	# Structural variations & supporting structural variations
	$stmt = qq{
    REPLACE INTO
    $sv_table (
	    variation_name,
	    source_id,
	    study_id,
			class_attrib_id,
	    tmp_class_name,
			is_evidence,
			somatic
    )
    SELECT
      t.id,
      $source_id,
			$study_id,
			a.attrib_id,
      t.type,
			t.is_ssv,
			t.is_somatic
    FROM
      temp_cnv t 
			LEFT JOIN attrib a ON (t.type=a.value)
  };
  $dbVar->do($stmt);
	
	
	# Failed structural variations
	debug(localtime()." Inserting into $sv_failed table for SVs");
	
	$stmt = qq{
    REPLACE INTO
    $sv_failed (
	    structural_variation_id,
	    failed_description_id
    )
    SELECT
      sv.structural_variation_id,
      17
    FROM
			structural_variation sv,
      temp_cnv t
    WHERE
			sv.variation_name=t.id AND
			t.is_ssv=0 AND
      t.is_failed!=0
  };
  $dbVar->do($stmt);
	
	
	# Failed supporting structural variations
	debug(localtime()." Inserting into $sv_failed table for SSVs");
	
	$stmt = qq{
    REPLACE INTO
    $sv_failed (
	    structural_variation_id,
	    failed_description_id
    )
    SELECT
      sv.structural_variation_id,
      18
    FROM
			structural_variation sv,
      temp_cnv t
    WHERE
			sv.variation_name=t.id AND
			t.is_ssv=1 AND
      t.is_failed!=0
  };
  $dbVar->do($stmt);
	
	
	# Structural variation features & supporting structural variation features
	debug(localtime()." Inserting into $svf_table table");
	
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
			class_attrib_id,
	    source_id,
			is_evidence,
			somatic,
			breakpoint_order
    )
    SELECT
			DISTINCT
      q.seq_region_id,
			t.outer_start,
      t.start,
			t.inner_start,
			t.inner_end,
      t.end,
			t.outer_end,
      t.strand,
			sv.structural_variation_id,
      t.id,
			sv.class_attrib_id,
      $source_id,
			t.is_ssv,
			t.is_somatic,
			t.bp_order
    FROM
      seq_region q,
      temp_cnv t,
			$sv_table sv
    WHERE
      q.name = t.chr AND
			t.id=sv.variation_name AND
			t.is_failed!=2 AND
			sv.source_id=$source_id
  };
  $dbVar->do($stmt);
	
	
	# UPDATE breakpoint_order to 1 when more than one SVF and no breakpoint order
	$stmt = qq{
	 SELECT distinct sv.structural_variation_id, count(svf.structural_variation_feature_id) svfc 
	 FROM $sv_table sv, $svf_table svf 
	 WHERE sv.structural_variation_id=svf.structural_variation_id 
	   AND sv.study_id=$study_id
		 AND svf.breakpoint_order is null
	 GROUP BY structural_variation_id having svfc>1
	};
	my $sv_ids = $dbVar->selectall_arrayref($stmt);
	$stmt = qq{
		UPDATE $svf_table SET breakpoint_order=1 WHERE structural_variation_id=?
	};
	foreach my $sv_id (@$sv_ids) {
		my $id = $sv_id->[0];
		$dbVar->do(qq{UPDATE $svf_table SET breakpoint_order=1 WHERE structural_variation_id=$id});
	}
	
	
	# Structural variation associations
	debug(localtime()." Inserting into $sva_table table");
	$stmt = qq{
    REPLACE INTO
    $sva_table (
	    structural_variation_id,
			supporting_structural_variation_id
    )
    SELECT
      sv2.structural_variation_id,
      sv1.structural_variation_id
    FROM
      temp_cnv t1,
			temp_cnv t2,
			structural_variation sv1,
			structural_variation sv2
    WHERE
			t1.id=sv1.variation_name AND
			t1.parent=t2.id AND
			t2.id=sv2.variation_name
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
}


sub structural_variation_annotation {
	my $stmt;
	
	debug(localtime()." Inserting into $svanno_table table");

	my $study_name = ($dbVar->selectrow_arrayref(qq{SELECT name FROM study WHERE study_id=$study_id;}))->[0];

	if ($dbVar->do(qq{show columns from $svanno_table like 'tmp_clinic_name';}) != 1){
		$dbVar->do(qq{ALTER TABLE $svanno_table ADD COLUMN tmp_clinic_name varchar(255);});
	}

	# Create sample entries
	$stmt = qq{ SELECT DISTINCT sample FROM temp_cnv WHERE is_ssv=1 };
  my $rows_samples = $dbVar->selectall_arrayref($stmt);
	foreach my $row (@$rows_samples) {
    my $sample = $row->[0];
		next if ($sample eq	'');
		
		my $srow = $dbVar->selectrow_arrayref(qq{ SELECT sample_id FROM sample WHERE name='$sample' });
		next if (defined($srow));
		$dbVar->do(qq{ INSERT INTO sample (name,description) VALUES ('$sample','Sample from the DGVa study $study_name')});
	}
	
	# Create strain entries
	$stmt = qq{ SELECT DISTINCT strain FROM temp_cnv};
  my $rows_strains = $dbVar->selectall_arrayref($stmt);
	my $strain_type = 1;
	$strain_type = 3 if ($species =~ /human|homo_sapiens/i);
	foreach my $row (@$rows_strains) {
    my $sample = $row->[0];
		next if ($sample eq	'');
		
		my $sample_id;
		my $srow = $dbVar->selectrow_arrayref(qq{SELECT sample_id FROM sample WHERE name='$sample'});
		if (!defined($srow)) {
			$dbVar->do(qq{ INSERT INTO sample (name,description) VALUES ('$sample','Sample from the DGVa study $study_name')});
			$sample_id = $dbVar->{'mysql_insertid'};
		}
		else {
			$sample_id = $srow->[0];
		}
		# Individual entry
		my $indiv_id = $dbVar->selectrow_arrayref(qq{SELECT sample_id FROM individual WHERE sample_id=$sample_id});
		if (!$indiv_id) {
			$dbVar->do(qq{ INSERT INTO individual (sample_id,gender,individual_type_id) VALUES ($sample_id,'Unknown',$strain_type)});
		}	
	}

	
	# Create phenotype entries
	$stmt = qq{ SELECT DISTINCT phenotype FROM temp_cnv};
  my $rows = $dbVar->selectall_arrayref($stmt);
	while (my $row = shift(@$rows)) {
    my $phenotype = $row->[0];
		next if ($phenotype eq '');
		
		$phenotype =~ s/'/\\'/g;
		$stmt = qq{ SELECT phenotype_id FROM phenotype WHERE description='$phenotype' };
		my $srow = $dbVar->selectall_arrayref($stmt);
		if (!scalar(@$srow)) {
			$dbVar->do(qq{ INSERT INTO phenotype (description) VALUES ('$phenotype')});
		}	
	}
	
	$stmt = qq{ SELECT attrib_type_id FROM attrib_type WHERE code='$clinical_attrib_type'};
	my $clinical_attrib_type_id = ($dbVar->selectall_arrayref($stmt))->[0][0];
	my $ext;
	
	# For mouse
	if ($species =~ /mouse|mus/i) {
		$ext = " AND s1.display='UNDISPLAYABLE'";
	}
	
	$stmt = qq{
    REPLACE INTO
    $svanno_table (
	    structural_variation_id,
 			clinical_attrib_id,
			tmp_clinic_name,
	    phenotype_id,
      sample_id,
      strain_id
    )
    SELECT DISTINCT
      sv.structural_variation_id,
			c.attrib_id,
			t.clinic,
			p.phenotype_id,
      s1.sample_id,
			s2.sample_id
    FROM
			structural_variation sv,
      temp_cnv t 
			LEFT JOIN attrib c ON (c.value=t.clinic AND c.attrib_type_id=$clinical_attrib_type_id)
			LEFT JOIN phenotype p ON (p.description=t.phenotype)
			LEFT JOIN sample s1 ON (s1.name=t.sample$ext)
			LEFT JOIN sample s2 ON (s2.name=t.strain)
    WHERE
			sv.variation_name=t.id
  };
  $dbVar->do($stmt);
	
	$stmt = qq{
		DELETE FROM $svanno_table
		WHERE phenotype_id IS NULL 
    	AND sample_id IS NULL 
    	AND strain_id IS NULL
    	AND clinical_attrib_id IS NULL
	};
	
	$dbVar->do($stmt);
  $dbVar->do(qq{DROP TABLE temp_cnv;});
}


sub meta_coord{
  
  debug(localtime()." Adding entries to meta_coord table");
  
	# Structural variation feature
  my $max_length_ref = $dbVar->selectall_arrayref(qq{
		SELECT MAX(seq_region_end - seq_region_start + 1) FROM $svf_table
  });
  my $max_length = $max_length_ref->[0][0];
	if (!defined($max_length)) { $max_length = 0; }
	
  my $cs = $default_cs->dbID;
  
  $dbVar->do(qq{DELETE FROM meta_coord where table_name = '$svf_table';});
	$dbVar->do(qq{INSERT INTO meta_coord(table_name, coord_system_id, max_length) VALUES ('$svf_table', $cs, $max_length);});
}



sub parse_gvf {
	my $infile = shift;
	
	debug(localtime()." Parse GVF file: $infile");
	if ($mapping) {
  	debug(localtime()." Mapping SV features to $target_assembly");
	}
  
  open IN, $infile or die "Could not read from input file $infile\n";
  
  # get a slice adaptor
  my $sa = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'slice');
  
  open OUT, ">$TMP_DIR/$TMP_FILE" or die "Could not write to output file $TMP_DIR/$TMP_FILE\n";
  
	my $sth = $dbVar->prepare(qq{ SELECT seq_region_id FROM seq_region WHERE name=?});
	
	my $header;
	my $assembly;
	
  while(<IN>) {
	
		my $current_line = $_;
		
		# Header (Study)
		if ($current_line =~ /^#/) {
			while ($current_line =~ /^#/) {
				chomp($current_line);
				$header = get_header_info($current_line,$header);
				$current_line = <IN>;
			}
			$assembly = parse_header($header);
		}
	
		chomp ($current_line);
		my $info = parse_line($current_line);
		
		my $is_failed = 0;
		# Flag to know if the entry is a sv or a ssv
		my $is_ssv = ($info->{ID} =~ /ssv/) ? 1 : 0;
		
    ###### Mapping work ######	
		if ($mapping) {
			# get a slice on the old assembly
			my $slice;
			my $start_c = $info->{start};
			my $end_c   = $info->{end};
	
			if (!$info->{start}) {
				if ($info->{outer_start}) {
					$start_c = $info->{outer_start};
				} else {
					$start_c = $info->{inner_start};
				}
			}
			if (!$info->{end}) {
				if ($info->{outer_end}) {
					$end_c = $info->{outer_end};
				} else {
					$end_c = $info->{inner_end};
				}
			}
	
			eval { $slice = $sa->fetch_by_region('chromosome', $info->{chr}, $start_c, $end_c, 1, $assembly); };
	
			# check got the slice OK
			if(!defined($slice)) {
	 			warn("Unable to map from assembly $assembly, or unable to retrieve slice ".$info->{chr}."\:$start_c\-$end_c");
	  		$skipped++;
				$is_failed = 2;
			}
			else {
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
						warn("Segments of ".$info->{ID}." map to different chromosomes ($to_chr and ", $to_slice->seq_region_name, ")");
						$is_failed = 2;
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
				if((($count - 1) <= $num_gaps && $count && !$is_failed)) {# && $diff < $size_diff) {
	  		
	  			$num_mapped{$assembly}++ if (!$is_ssv);
	  		
					$info->{chr}   = $to_chr;
					$info->{start} = $min_start;
					$info->{end}   = $max_end;
				
	  			# calculate inner/outer coords
	  			$info->{outer_start} = $min_start - ($start_c - $info->{outer_start}) if $info->{outer_start} >= 1;
	  			$info->{inner_start} = $min_start + ($info->{inner_start} - $start_c) if $info->{inner_start} >= 1;
	  			$info->{inner_end}   = $max_end - ($end_c - $info->{inner_end}) if $info->{inner_end} >= 1;
	  			$info->{outer_end}   = $max_end + ($info->{outer_end} - $end_c) if $info->{outer_end} >= 1;

				}
				else {
					if (!$is_ssv) {
	  				warn ("Structural variant '".$info->{ID}."' from study '".$header->{study}."' has location '$assembly:".$info->{chr}.":$start_c\-$end_c' , which could not be re-mapped to $target_assembly. This variant will be labelled as failed");
	  				$num_not_mapped{$assembly}++;
					}	
					$is_failed = 1;
				}
			}
		}
		
		
		print OUT (join "\t", ($info->{ID},
													 $info->{SO_term},
		                       $info->{chr}, 
													 $info->{outer_start}, 
													 $info->{start}, 
													 $info->{inner_start}, 
													 $info->{inner_end}, 
													 $info->{end}, 
													 $info->{outer_end},
													 $info->{strand}, 
													 $info->{Parent},
													 $info->{clinical_int},
													 $info->{phenotype},
													 $info->{samples},
													 $info->{strain_name},
													 $info->{bp_order},
													 1,
													 $is_ssv, 
													 $is_failed)
							);
  	print OUT "\n";
  }
  close IN;
  close OUT;
  
  debug(localtime()." Finished SV mapping:\n\t\tSuccess: ".(join " ", %num_mapped)." not required $no_mapping_needed\n\t\tFailed: ".(join " ", %num_not_mapped)." In no existing chromosome: $skipped");
}


sub verifications {
	debug(localtime()." Verification of ftp links");
	my $sth = $dbVar->prepare(qq{ SELECT url FROM $study_table WHERE source_id=$source_id});
	$sth->execute();
	my $failed_flag = 0;
	while (my $url = ($sth->fetchrow_array)[0]) {
		if (!head($url)) {
			print STDERR "\t$url IS NOT VALID\n";
			$failed_flag = 1;
		}
	}
	debug(localtime()." URL verification: OK") if ($failed_flag == 0);
}


sub get_header_info {
	my $line = shift;
	my $h    = shift;
	
	my ($label, $info);
	if ($line =~ /^##/) {
		($label, $info) = split(' ', $line);
	} else {
		($label, $info) = split(':', $line);
	}
	$label =~ s/#//g;
	$label =~ s/^\s//;
	$info =~ s/^\s+//;
	
	$h->{author}       = $info if ($label eq 'Display_name');
	$h->{first_author} = $info if ($label eq 'First_author');
	$h->{assembly}     = $info if ($label eq 'assembly-name');
	$h->{study}        = $info if ($label eq 'Study_accession');
	$h->{pubmed}       = $info if ($label eq 'PMID');
	$h->{year}         = $info if ($label eq 'Publication_year');
	$h->{desc}         = $info if ($label eq 'Paper_title' && $info && $info ne 'None Given');
	$h->{desc}         = 'Catalogue of Somatic Mutations in Cancer (COSMIC)' if ($label eq 'Study_description' && !$h->{desc} && length($info)<256); # Get the Paper title if possible
	
	return $h;
}


sub parse_header {
	my $h = shift;
	
	my $assembly = $h->{assembly};
	
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
	  elsif($assembly !~ /$cs_version_number/) {
			warn("Could not guess assembly from file - assuming to be NCBI35");
			$assembly = 'NCBI35';
	  }
	}
	
	if($species =~ /mouse|mus.*mus/) {
	  if($assembly =~ /34/) {
			$assembly = 'NCBIM34';
	  }
	  elsif($assembly =~ /35/) {
			$assembly = 'NCBIM35';
	  }
	  elsif($assembly =~ /36/) {
			$assembly = 'NCBIM36';
	  }
		elsif($assembly =~ /37/) {
			$assembly = 'NCBIM37';
	  }
	}
	
	$h->{assembly} = $assembly;
	
	study_table($h);
	
	return $assembly;
}

sub parse_line {
	my $line     = shift;
	my @data     = split(/\t/, $line);
	my @last_col = split(';', pop(@data));
	my $info;
	
	foreach my $inf (@last_col) {
	
		if ($inf =~ /Primary site=(\w+)/) {
			$info->{phenotype} = "COSMIC:tumour_site:$1";
		}
		elsif ($inf =~ /sample_name=(\w+)/) {
			$info->{samples} = $1;
		}
		my ($key,$value) = split('=',$inf);
		$info->{$key} = $value;
	}
	
	# Start range
	my @s_range = split(/,/, $info->{Start_range});
	if (scalar @s_range == 2) {
		$info->{outer_start} = $s_range[0] if ($s_range[0] ne '.');
		$info->{inner_start} = $s_range[1] if ($s_range[1] ne '.');
	} 
	# End range
	my @e_range = split(/,/, $info->{End_range});
	if (scalar @e_range == 2) {
		$info->{outer_end} = $e_range[1] if ($e_range[1] ne '.');
		$info->{inner_end} = $e_range[0] if ($e_range[0] ne '.');
	}
	
	
	$data[0] =~ /Chr(.+)$/;
	$data[0] = $1;
	$info->{ID}      = $info->{Name} if ($info->{Name});
	$info->{chr}     = $data[0];
	$info->{SO_term} = $data[2];
	$info->{start}   = $data[3];
	$info->{end}     = $data[4];
  $info->{strand}  = $data[5];

	$info->{bp_order}  = ($info->{submitter_variant_id} =~ /\w_(\d+)$/) ? $1 : undef;
	return $info;
}


sub annotation_post_processing {
	# First round - duplicated combinations sample/strain
	my @sv_list;
	my $sth = $dbVar->prepare(qq{ SELECT structural_variation_id, count(structural_variation_annotation_id) c 
	                              FROM $svanno_table GROUP BY structural_variation_id HAVING c>1
								           });
	$sth->execute();
	while (my @res = ($sth->fetchrow_array)) {
		push (@sv_list,$res[0]);
	}
	my $svs = join(',',@sv_list);
	my $stmt = qq{ DELETE FROM $svanno_table WHERE sample_id=strain_id AND structural_variation_id IN ($svs) };
	print STDRR "QUERY:\n$stmt\n";
	$dbVar->do($stmt);
	
	# Second round - duplicated samples
	@sv_list;
	$sth->execute();
	while (my @res = ($sth->fetchrow_array)) {
		push (@sv_list,$res[0]);
	}
	$svs = join(',',@sv_list);
	$stmt = qq{ DELETE FROM $svanno_table 
	            WHERE sample_id IN (SELECT sample_id FROM individual WHERE individual_type_id=1) 
							AND structural_variation_id IN ($svs) 
						};
	$dbVar->do($stmt);
	
	# Third round - duplicated strains
	@sv_list;
	$sth->execute();
	while (my @res = ($sth->fetchrow_array)) {
		push (@sv_list,$res[0]);
	}
	$svs = join(',',@sv_list);
	$stmt = qq{ DELETE FROM $svanno_table 
	            WHERE strain_id NOT IN (SELECT sample_id FROM individual WHERE individual_type_id=1) 
							AND structural_variation_id IN ($svs) 
						};
	$dbVar->do($stmt);
}


sub usage {
  die shift, qq{

Options -
  -species         : species name (required)
  -target_assembly : assembly version to map to (optional)
  -tmpdir          : (optional)
  -tmpfile         : (optional)
  -input_file      : file containing DGVa data dump (required if no input_dir)
	-input_dir       : directory containing DGVa data dump (required if no input_file)
  -mapping         : if set, the data will be mapped to $target_assembly using the Ensembl API
  -gaps            : number of gaps allowed in mapping (defaults to 1) (optional)
  -size_diff       : % difference allowed in size after mapping (optional)
	-version         : version number of the data (required)
  -registry        : registry file (optional)
  };
}
