#!/usr/bin/env perl

=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

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

our ($species, $input_file, $input_dir, $source_name, $TMP_DIR, $TMP_FILE, $var_set_id, $mapping, $num_gaps,
     $target_assembly, $cs_version_number, $size_diff, $version, $registry_file, $replace_study, $debug);

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
					 'registry=s'        => \$registry_file,
					 'replace!'          => \$replace_study,
					 'debug!'            => \$debug,
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
my %display_name = ( 'COSMIC' => 'http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/' );

my $add          = '';
my $study_table  = "study$add";
my $sv_table     = "structural_variation$add";
my $svf_table    = "structural_variation_feature$add";
my $sva_table    = "structural_variation_association$add";
my $sv_failed    = "failed_structural_variation$add";
my $set_table    = "variation_set_structural_variation$add";
my $svanno_table = "structural_variation_annotation$add";

my $tmp_sv_col  = 'tmp_class_name';
my $tmp_sva_col = 'tmp_clinic_name';
	
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

# variation set
my %var_set = ('pilot1' => 31, 'pilot2' => 32);


# run the mapping sub-routine if the data needs mapping
my (%num_mapped, %num_not_mapped, %samples, %subject_pop, %study_done);
my $no_mapping_needed = 0;
my $skipped = 0;
my $failed = [];
my $study_id;
my $somatic_study;
my $fname;

my $source_id = source();


##########
#  Main  #
##########

foreach my $in_file (@files) {
	next if ($in_file !~ /\.gvf$/);
	
	my $msg ="File: $in_file ($f_count/$f_nb)";
	print "$msg\n";
	debug("##########\n$msg\n##########");
	
	$fname = undef;
	if ($input_dir) {
		$fname = "$input_dir/$in_file";
	}
	else {
		$fname = $in_file;
	}
	
	# Variation set
	if ($species eq 'human' && $fname =~ /pilot(\d)/) { 
		$var_set_id = $var_set{"pilot$1"};
	}
	elsif ($species eq 'human' && $fname =~ /estd199/) {
	  $var_set_id = '1kg';
	}
	else {
		$var_set_id = undef;
	}

	%num_mapped = ();
	%num_not_mapped = ();
	%samples = ();
	$no_mapping_needed = 0;
	$skipped = 0;
	$failed = [];
	$study_id;
	$somatic_study = 0;
	

	# Parsing
	parse_gvf($fname);
	load_file_data();
	
	# Insertion
	structural_variation();
	structural_variation_set() if ($var_set_id);
	structural_variation_annotation();
	debug(localtime()." Done!\n");
	$f_count ++;
}

## Post processings ##

# Post processing for mouse annotation (delete duplicated entries in structural_variation_annotation)
post_processing_annotation() if ($species =~ /mouse|mus/i);
post_processing_feature();
post_processing_sample();

# Finishing methods
meta_coord();
verifications(); # URLs
cleanup();

debug(localtime()." All done!");



#############
#  Methods  #
#############

sub load_file_data{
  debug(localtime()." Loading file into temporary table");
  
  # drop any table that may be there
  $dbVar->do("DROP TABLE IF EXISTS temp_cnv;");
  
  create_and_load(
	$dbVar, "temp_cnv", "id *", "type", "chr", "outer_start i", "start i", "inner_start i",
	"inner_end i", "end i", "outer_end i", "parent", "clinic", "phenotype", "sample", "strain",
  "is_ssv i", "is_failed i", "population", "bp_order i", "is_somatic i", "status");		
	
	
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
	# Case with breakpoints						
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
  
	$data->{desc} = $data->{s_desc} if ($data->{s_desc} and !$data->{desc});
	
	my $stmt;
	my $study      = $data->{study};
	my $assembly   = $data->{assembly};
	my $study_desc = $data->{desc};
	my $study_type = $data->{study_type};
	my $pmid         = $data->{pubmed};
  my $first_author = $data->{first_author};
  my $year         = $data->{year};
	my $author       = $data->{author};
	my $author_info  = $data->{author};
	my $study_type   = ($data->{study_type}) ? $data->{study_type} : 'NULL';
	
	# Author
	$author_info =~ s/_/ /g;
	my $year_desc = $year->[0];
	if ($author_info =~ /(.+)\set\sal\s(.+)/) {
	  $first_author = $1 if (!$first_author);
	  $year_desc = $2 if (defined($2));
	}
				
	my $author_desc;
	$author_desc = "$first_author $year_desc " if ($first_author and defined($year_desc));
	$author = (split('_',$author))[0] if ($author !~ /.+_et_al_.+/);
	
	# Study desc	
	if (length($study_desc)>180) {
		$study_desc = substr($study_desc,0,180);
		$study_desc .= '...';
	}
	
	#	PubMed ID
	my $pmid_desc = '';
	my $pubmed    = 'NULL';
	if (defined($pmid)) {
		my $pmid_len = scalar(@$pmid);
		if ($pmid_len) {
			$pmid_desc = 'PMID:';
      if ($pmid_len > 1 and scalar(@$year) == $pmid_len){
				for (my $i=0; $i<$pmid_len; $i++) {
					$pmid_desc .= ',' if($i>0);
					$pmid_desc .= "$pmid->[$i]($year->[$i])";
				}
			}
			else {		
				$pmid_desc .= join(',',@$pmid);
			}
			$pubmed = "pubmed/".join(",pubmed/", @$pmid);
		}
	}
	my $external_link = ($display_name{$author}) ? $display_name{$author} : $pubmed;
	
	# URL
	$study =~ /(\w+\d+)\.?\d*/;
	my $study_ftp = $1;
	$study_ftp = "ftp://ftp.ebi.ac.uk/pub/databases/dgva/$study_ftp\_$author";
	
	
	my $assembly_desc = " [remapped from build $assembly]" if ($mapping and $assembly ne $target_assembly);
	
	$stmt = qq{ SELECT st.study_id, st.description FROM study st, source s 
		          WHERE s.source_id=st.source_id AND s.name='DGVa' and st.name='$study'};
	my $rows = $dbVar->selectall_arrayref($stmt);		
	my ($assembly_desc,);
	
	
		
	# UPDATE
	if (scalar (@$rows)) {
		$study_id   = $rows->[0][0];
		$study_desc = $rows->[0][1];
		
  	remove_data($study_id) if (defined($replace_study) && defined($study_id) && !defined($study_done{$study}));
		
		# Checks the studies with different assemblies when the option "mapping" is selected.
		if ($mapping and $assembly ne $target_assembly) {
			if ($study_desc =~ /^(.+)\s\[remapped\sfrom\sbuild(s?)\s(.+)\]$/) {
				my $gen_desc = $1;
				$assembly_desc = $3;
			
				if ($assembly_desc !~ /$assembly/) {
					$assembly_desc .= " and $assembly";
					#$stmt = qq{
					#	UPDATE study SET description='$gen_desc [remapped from builds $assembly_desc]' 
					#	WHERE study_id=$study_id
					#};
					#$dbVar->do($stmt);
					$study_desc = "$gen_desc [remapped from builds $assembly_desc]";
				}
			}
			else {
				#$stmt = qq{
				#		UPDATE study SET description='$study_desc [remapped from build $assembly]' 
				#		WHERE study_id=$study_id
				#};
				#$dbVar->do($stmt);
				$study_desc = "$study_desc [remapped from build $assembly]";
			}
		}
		elsif ($assembly eq $target_assembly && $study_desc =~ /^(.+)\s\[remapped/) {
			$study_desc = $1;
		}
		
		$stmt = qq{ UPDATE $study_table SET 
		              description='$study_desc',
								  external_reference='$external_link',
								  study_type='$study_type',
								  url='$study_ftp'
								WHERE study_id=$study_id
							};
    $dbVar->do($stmt);
	}
	# INSERT
	else {
	
		$stmt = qq{
    	INSERT INTO `$study_table` (
				`name`,
				`description`,
				`url`,
				`source_id`,
				`external_reference`,
				`study_type`
      )
    	VALUES (
      	'$study',
				CONCAT(
					'$author_desc',
					'"',
					'$study_desc',
					'" $pmid_desc',
					'$assembly_desc'),
				'$study_ftp',
				$source_id,
				'$external_link',
				'$study_type'
			)
		};
    $dbVar->do($stmt);
		$study_id = $dbVar->{'mysql_insertid'};
	}
	
	# Avoid to replace a study twice, when a second file with the patched assembly is also imported.
	$study_done{$study} = 1;
}


sub structural_variation {
  
  debug(localtime()." Inserting into $sv_table table");
  
	if ($dbVar->do(qq{show columns from $sv_table like '$tmp_sv_col';}) != 1){
		$dbVar->do(qq{ALTER TABLE $sv_table ADD COLUMN $tmp_sv_col varchar(255);});
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
    INSERT IGNORE INTO
    $sv_table (
	    variation_name,
	    source_id,
	    study_id,
			class_attrib_id,
	    $tmp_sv_col,
			is_evidence,
			somatic,
			validation_status
    )
    SELECT
      t.id,
      $source_id,
			$study_id,
			a.attrib_id,
      t.type,
			t.is_ssv,
			t.is_somatic,
			t.status
    FROM
      temp_cnv t 
			LEFT JOIN attrib a ON (t.type=a.value)
		WHERE
			t.chr IN (SELECT name FROM seq_region)
  };
  $dbVar->do($stmt);
	
	
	# Failed structural variations
	debug(localtime()." Inserting into $sv_failed table for SVs");
	
	$stmt = qq{
    INSERT IGNORE INTO
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
    INSERT IGNORE INTO
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
    INSERT IGNORE INTO $svf_table (
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
      1,
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
			t.is_failed=0 AND
			sv.source_id=$source_id
  };
  $dbVar->do($stmt);
	
	
	# UPDATE breakpoint_order to 1 when more than one SVF and no breakpoint order
	if ($somatic_study) {
	  $stmt = qq{
	    SELECT distinct sv.structural_variation_id, count(svf.structural_variation_feature_id) svfc 
	    FROM $sv_table sv, $svf_table svf 
	    WHERE sv.structural_variation_id=svf.structural_variation_id 
	      AND sv.study_id=$study_id
				AND sv.somatic=1
		    AND svf.breakpoint_order is null
	      GROUP BY structural_variation_id having svfc>1
	  };
	  my $sv_ids = $dbVar->selectall_arrayref($stmt);
	  foreach my $sv_id (@$sv_ids) {
		  my $id = $sv_id->[0];
		  $dbVar->do(qq{UPDATE $svf_table SET breakpoint_order=1 WHERE structural_variation_id=$id AND somatic=1});
	  }
	}
	
	# Structural variation associations
	debug(localtime()." Inserting into $sva_table table");
	$stmt = qq{
    INSERT IGNORE INTO
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


sub structural_variation_set {
	debug(localtime()." Inserting into $set_table table");
	
	my $stmt;
	my @var_set_list;
	
	if ($var_set_id eq '1kg') {
	
    my %pop_1kg_p1 = (
      'ASW' => 'AFR', 'LWK' => 'AFR', 'YRI' => 'AFR',                                 # AFR
      'CLM'	=> 'AMR', 'PEL'	=> 'AMR', 'PUR'	=> 'AMR', 'MXL'	=> 'AMR',                 # AMR
      'CEU'	=> 'EUR', 'TSI'	=> 'EUR', 'FIN'	=> 'EUR', 'GBR'	=> 'EUR', 'IBS'	=> 'EUR', # EUR
      'CHB'	=> 'ASN', 'JPT'	=> 'ASN', 'CHS'	=> 'ASN', 'CDX'	=> 'ASN', 'KHV'	=> 'ASN', # ASN
		);
		
	  while (my($sub_pop, $sup_pop) = each (%pop_1kg_p1)) {
			$dbVar->do(qq{UPDATE temp_cnv SET population='$sup_pop' WHERE population='$sub_pop'});
		}
		
	  $stmt = qq{
   	  INSERT IGNORE INTO
   	  $set_table (
	      structural_variation_id,
			  variation_set_id
      )
      SELECT DISTINCT
     	  sva.structural_variation_id,
     	  vs.variation_set_id
      FROM
     	  temp_cnv t,
			  structural_variation sv,
				structural_variation_association sva,
			  variation_set vs
      WHERE
			  t.id=sv.variation_name AND
				sva.supporting_structural_variation_id=sv.structural_variation_id AND
				t.population is not null AND
			  vs.name=CONCAT('1000 Genomes - ',t.population)
    };
	  $dbVar->do($stmt);
		
		# Variation set "All"
		$stmt = qq{
   	  INSERT IGNORE INTO
   	  $set_table (
	      structural_variation_id,
			  variation_set_id
      )
      SELECT DISTINCT
     	  sva.structural_variation_id,
     	  vs.variation_set_id
				
		  FROM
     	  temp_cnv t,
			  structural_variation sv,
				structural_variation_association sva,
			  variation_set vs
      WHERE
			  t.id=sv.variation_name AND
				sva.supporting_structural_variation_id=sv.structural_variation_id AND
				t.population is not null AND
			  vs.name='1000 Genomes - All'
    };
		$dbVar->do($stmt);
		
		my $sets_1kg = $dbVar->selectrow_arrayref(qq{
				 SELECT DISTINCT vssv.variation_set_id 
				 FROM variation_set_structural_variation vssv, structural_variation sv, temp_cnv t 
				 WHERE vssv.structural_variation_id=sv.structural_variation_id AND sv.variation_name=t.id
			 });
		foreach my $s_id (@$sets_1kg) {
		  push (@var_set_list,$s_id);
		}	 
		
	}
	else {
	
	  $stmt = qq{
   	  INSERT IGNORE INTO
   	  $set_table (
	      structural_variation_id,
			  variation_set_id
      )
      SELECT
     	  sv.structural_variation_id,
     	  $var_set_id
      FROM
     	  temp_cnv t,
			  structural_variation sv
      WHERE
			  t.is_ssv=0 AND
			  t.id=sv.variation_name
    };
		$dbVar->do($stmt);
		
		push(@var_set_list,$var_set_id);
  }  
	
	my $set_id_list = join(',',@var_set_list);
	
  # Populate the meta table
	my $sets = $dbVar->selectall_arrayref(qq{SELECT s.name,a.value FROM variation_set s, attrib a 
	                                        WHERE s.variation_set_id IN ($set_id_list) 
																					AND a.attrib_id=s.short_name_attrib_id;}
																			 );
	foreach my $set (@$sets) {								
	  my $meta_value = "sv_set#".$set->[0]."#".$set->[1]; 
	  # Check if the meta entry already exists, else it create the entry
	  if (!$dbVar->selectrow_arrayref(qq{SELECT meta_id FROM meta WHERE meta_key='web_config' 
	      															 AND meta_value='$meta_value';})) {
		  $dbVar->do(qq{INSERT INTO meta (meta_key,meta_value) VALUES ('web_config','$meta_value');});
	  }
	}
}


sub structural_variation_annotation {
	my $stmt;
	
	debug(localtime()." Inserting into $svanno_table table");

	my $study_name = ($dbVar->selectrow_arrayref(qq{SELECT name FROM study WHERE study_id=$study_id;}))->[0];

	if ($dbVar->do(qq{show columns from $svanno_table like '$tmp_sva_col';}) != 1){
		$dbVar->do(qq{ALTER TABLE $svanno_table ADD COLUMN $tmp_sva_col varchar(255);});
	}

	# Create sample entries
	$stmt = qq{ SELECT DISTINCT sample FROM temp_cnv WHERE is_ssv=1 
	            AND sample NOT IN (SELECT DISTINCT name from sample)
						};
  my $rows_samples = $dbVar->selectall_arrayref($stmt);
	foreach my $row (@$rows_samples) {
    my $sample = $row->[0];
		next if ($sample eq	'');
		$dbVar->do(qq{ INSERT INTO sample (name,description) VALUES ('$sample','Subject from the DGVa study $study_name')});
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
		my $srow = $dbVar->selectrow_arrayref(qq{SELECT s.sample_id FROM sample s WHERE s.name='$sample' LIMIT 1});
		if (!defined($srow)) {
			# Strain entry
			$dbVar->do(qq{ INSERT INTO sample (name,description) VALUES ('$sample','Strain from the DGVa study $study_name')});
			$sample_id = $dbVar->{'mysql_insertid'};
		}
		else {
			$sample_id = $srow->[0];
		}
	}

	
	# Create phenotype entries
	$stmt = qq{ SELECT DISTINCT phenotype FROM temp_cnv 
	            WHERE phenotype NOT IN (SELECT description FROM phenotype)
						};
  my $rows = $dbVar->selectall_arrayref($stmt);
	while (my $row = shift(@$rows)) {
    my $phenotype = $row->[0];
		next if ($phenotype eq '');
		
		$phenotype =~ s/'/\\'/g;
		$dbVar->do(qq{ INSERT INTO phenotype (description) VALUES ('$phenotype')});	
	}
	
	$stmt = qq{ SELECT attrib_type_id FROM attrib_type WHERE code='$clinical_attrib_type'};
	my $clinical_attrib_type_id = ($dbVar->selectall_arrayref($stmt))->[0][0];
	my $ext;
	
	# For mouse
	if ($species =~ /mouse|mus/i) {
		$ext = " AND s1.display='UNDISPLAYABLE'";
	}
	
	$stmt = qq{
    INSERT IGNORE INTO
    $svanno_table (
	    structural_variation_id,
 			clinical_attrib_id,
			$tmp_sva_col,
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
  $dbVar->do(qq{DROP TABLE temp_cnv;}) if (!defined($debug));
}




#### Parsing methods ####


sub parse_gvf {
	my $infile = shift;
	
	debug(localtime()." Parse GVF file: $infile");
  debug(localtime()." Start mapping SV features to $target_assembly if needed") if $mapping;
  
  open IN, $infile or die "Could not read from input file $infile\n";
  
  # get a slice adaptor
  my $sa = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'slice');
  
  open OUT, ">$TMP_DIR/$TMP_FILE" or die "Could not write to output file $TMP_DIR/$TMP_FILE\n";
	
	my $header;
	my $assembly;
	
  while(<IN>) {
	
		my $current_line = $_;
		
		###### Header (Study) ######
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
		
    ###### Mapping work ######	
		if ($mapping && $assembly !~ /$cs_version_number/) {
			
			
			
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
	  		
	  			$num_mapped{$assembly}++ if (!$info->{is_ssv});
	  		
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
					if (!$info->{is_ssv}) {
	  				warn ("Structural variant '".$info->{ID}."' from study '".$header->{study}."' has location '$assembly:".$info->{chr}.":$start_c\-$end_c' , which could not be re-mapped to $target_assembly. This variant will be labelled as failed");
	  				$num_not_mapped{$assembly}++;
					}	
					$is_failed = 1;
				}
			}
		}
		
		$info->{is_failed} = $is_failed;
							 
		my $data = generate_data_row($info);					 
							 
		print OUT (join "\t", @{$data})."\n";
  }
  close IN;
  close OUT;
	
	if ($mapping && $assembly !~ /$cs_version_number/) {
		debug(localtime()." Finished SV mapping:\n\t\tSuccess: ".(join " ", %num_mapped)." not required $no_mapping_needed\n\t\tFailed: ".(join " ", %num_not_mapped)." In no existing chromosome: $skipped");
	}
}


sub get_header_info {
	my $line = shift;
	my $h    = shift;
	
	chomp($line);
	
	my ($label, $info);
	if ($line =~ /^##/) {
		($label, $info) = split(' ', $line);
	} else {
		($label, $info) = split(':', $line);
	}
	$label =~ s/#//g;
	$label =~ s/^\s//;
	$info  =~ s/^\s+//;
	
	$h->{author}       = $info if ($label =~ /Display.+name/i);
	$h->{first_author} = $info if ($label =~ /First.+author/i);
	$h->{assembly}     = $info if ($label =~ /Assembly.+name/i);
	$h->{study}        = $info if ($label =~ /Study.+accession/i);
	$h->{study_type}   = $info if ($label =~ /Study.+type/i);
	
	# Publication information
	if ($label =~ /Publication/i && $info !~ /Not.+applicable/i) {
		foreach my $pub (split(';',$info)) {
			my ($p_label,$p_info) = split('=',$pub);
			 	
			push(@{$h->{pubmed}}, $p_info) if ($p_label =~ /PMID/i);
			push(@{$h->{year}}, $p_info) if ($p_label =~ /Publication.+year/i);
			$h->{desc} = $p_info if ($p_label =~ /Paper.+title/i && $p_info && $p_info ne 'None Given');
		}
	}
	
	# Study information
	elsif ($label eq 'Study') {
		foreach my $st (split(';',$info)) {
			my ($s_label,$s_info) = split('=',$st);
			if ($s_label =~ /First.+author/i) {
				$s_info =~ /(\S+)\s*(\w+)/;
				$h->{first_author} = ($2) ? $2 : $s_info;
			}
			$h->{s_desc} = $s_info if ($s_label =~ /Description/i);
			$somatic_study = 1 if ($s_label =~ /Contains/i && $s_info =~ /somatic/i);
		}
	}
	
	# Sample information
	elsif ($label eq 'sample') {
		my ($sample,$subject,$tissue);
		foreach my $s_info (split(';',$info)) {
			my ($key,$value) = split('=',$s_info);
			$sample  = $value if ($key eq 'sample_name');
			$subject = $value if ($key eq 'subject_name');
			
			# Tissue (human)
			if ($key eq 'sample_cell_type'){
        $s_info =~ /Primary site=(.+),/;
				$tissue = $1;
			}
		}
		$tissue = ($tissue ne 'NS') ? $tissue : undef;
		
		# COSMIC phenotype
		if (defined($tissue) && ($h->{desc} =~ /cosmic/i || $h->{s_desc} =~ /cosmic/i)){
		  $tissue = "COSMIC:tumour_site:$tissue";
		}
		
		if (defined($sample) and defined($subject)){
			$samples{$sample}{name} = $subject;
			$samples{$sample}{tissue} = $tissue;
		}
	}
	
	# Subject information
	elsif ($label eq 'subject') {
		my ($subject,$pop,$phenotype);
		foreach my $s_info (split(';',$info)) {
			my ($key,$value) = split('=',$s_info);
			$pop = $value if ($key eq 'subject_population');
			$subject = $value if ($key eq 'subject_name');
			$phenotype = $value if ($key eq 'phenotype_description');
		}
		if (defined($pop) and defined($subject)){
			$samples{$subject}{population} = $pop;
			$samples{$subject}{phenotype}  = $phenotype;
		}
	}
	
	$h->{author} =~ s/\s/_/g;
	
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
			warn("Could not guess assembly from file - assuming to be GRCh37");
			$assembly = 'GRCh37';
	  }
		# For patches
		elsif($assembly =~ /(\w+\d+)\.p\d+/) {
			$assembly = $1;
	  }
	}
	elsif($species =~ /mouse|mus.*mus/) {
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
		elsif($assembly =~ /38/) {
			$assembly = 'GRCm38';
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
	
	$data[0] =~ /Chr(.+)$/;
	$info->{chr} = defined($1) ? $1 : $data[0];
	$info->{SO_term} = $data[2];
	$info->{start}   = $data[3];
	$info->{end}     = $data[4];

  $info = parse_9th_col($info,\@last_col);

	return $info;
}


sub parse_9th_col {
	my $info     = shift;
	my $last_col = shift;

	foreach my $inf (@$last_col) {
		my ($key,$value) = split('=',$inf);
		$info->{$key} = $value; # Default
		
		$info->{ID} = $value if ($key eq 'Name');
		
		# Start range
		if ($key eq 'Start_range') {
			my @s_range = split(/,/, $value);
			if (scalar @s_range == 2) {
				$info->{outer_start} = $s_range[0] if ($s_range[0] ne '.');
				$info->{inner_start} = $s_range[1] if ($s_range[1] ne '.');
			} 
		}
		# End range
		if ($key eq 'End_range') {
			my @e_range = split(/,/, $value);
			if (scalar @e_range == 2) {
				$info->{outer_end} = $e_range[1] if ($e_range[1] ne '.');
				$info->{inner_end} = $e_range[0] if ($e_range[0] ne '.');
			}
		}
			
		$info->{clinical}  = $value if ($key eq 'clinical_significance');
		
		if ($value !~ /Unknown/i) {
		  if ($key eq 'sample_name' || ($key eq 'subject_name' && !$info->{sample})) {
			  $info->{sample}     = ($samples{$value}{name}) ? $samples{$value}{name} : $value;
				$info->{population} = ($samples{$value}{population}) ? $samples{$value}{population} : undef; # 1000 Genomes study
		    $info->{phenotype}  = ($samples{$value}{phenotype}) ? decode_text($samples{$value}{phenotype}) : $samples{$value}{tissue};
			} 
		}
		
		$info->{strain_name} = $info->{population} if ($species =~ /mouse|mus_musculus/i);
		
		$info->{parent}     = $value if ($key eq 'Parent'); # Check how the 'parent' key is spelled
		$info->{is_somatic} = 1 if ($key eq 'var_origin' && $value =~ /somatic/i);
		$info->{bp_order}   = ($info->{submitter_variant_id} =~ /\w_(\d+)$/) ? $1 : undef;
		
		$info->{status}     = 'High quality' if ($key eq 'variant_region_description' && $value =~ /high.quality/i);
		
	}
	
	# Flag to know if the entry is a sv or a ssv
	$info->{is_ssv} = ($info->{ID} =~ /ssv/) ? 1 : 0;
	
	# Check if a SV has its SSV(s) labelled as 'somatic'
	if ($somatic_study && $info->{is_ssv} == 0) {
	  my $cmd  = "grep 'parent=".$info->{ID}."' $fname";
	  my $text = `$cmd`;
		$info->{is_somatic} = 1 if ($text =~ /var_origin=somatic/i || !$text);
	}
	
	return $info;
}

sub decode_text {
	my $text = shift;
	
	$text  =~ s/%3B/;/g;		
	$text  =~ s/%3D/=/g;
	$text  =~ s/%25/%/g;
	$text  =~ s/%26/&/g;
	$text  =~ s/%2C/,/g;
	
	return $text;
}




#### Post processing ####

sub post_processing_annotation {
	debug(localtime()." Post processing of the table $svanno_table");
	
	# First round - duplicated combinations sample/strain
	my @sv_list;
	my $sth = $dbVar->prepare(qq{ SELECT structural_variation_id, count(structural_variation_annotation_id) c 
	                              FROM $svanno_table GROUP BY structural_variation_id HAVING c>1
								           });
	$sth->execute();
	while (my @res = ($sth->fetchrow_array)) {
		push (@sv_list,$res[0]);
	}
	if (scalar(@sv_list)>0) {
		my $svs  = join(',',@sv_list); 
		my $stmt = qq{ DELETE FROM $svanno_table WHERE sample_id=strain_id AND structural_variation_id IN ($svs) };
		$dbVar->do($stmt);
	}
	
	# Second round - duplicated samples
	@sv_list;
	$sth->execute();
	while (my @res = ($sth->fetchrow_array)) {
		push (@sv_list,$res[0]);
	}
	if (scalar(@sv_list)>0) {
		my $svs  = join(',',@sv_list);
		my $stmt = qq{ DELETE FROM $svanno_table 
	                 WHERE sample_id IN (SELECT sample_id FROM individual WHERE individual_type_id=1) 
							       AND structural_variation_id IN ($svs) 
						     };
		$dbVar->do($stmt);
	}
	
	# Third round - duplicated strains
	@sv_list;
	$sth->execute();
	while (my @res = ($sth->fetchrow_array)) {
		push (@sv_list,$res[0]);
	}
	if (scalar(@sv_list)>0) {
		my $svs  = join(',',@sv_list);
		my $stmt = qq{ DELETE FROM $svanno_table 
	                 WHERE strain_id NOT IN (SELECT sample_id FROM individual WHERE individual_type_id=1) 
							       AND structural_variation_id IN ($svs) 
						     };
		$dbVar->do($stmt);
	}
}


sub post_processing_feature {
	debug(localtime()." Post processing of the table $svf_table");
	
	my $stmt_s = qq{ UPDATE $svf_table SET outer_start=NULL, inner_start=NULL 
	                 WHERE outer_start=seq_region_start AND inner_start=seq_region_start
								 };
	$dbVar->do($stmt_s);												 

	my $stmt_e = qq{ UPDATE $svf_table SET outer_end=NULL, inner_end=NULL 
	                 WHERE outer_end=seq_region_end AND inner_end=seq_region_end
								 };
	$dbVar->do($stmt_e);
}


# Change the sample description if a subject is associated to several studies
sub post_processing_sample {
	debug(localtime()." Post processing of the table sample");
	
	my $sth = $dbVar->prepare(qq{ SELECT s.sample_id, s.description, count(DISTINCT sv.study_id) c 
	                              FROM $sv_table sv, sample s, $svanno_table sva
	                              WHERE sv.structural_variation_id=sva.structural_variation_id
																  AND s.sample_id=sva.sample_id
																GROUP BY sva.sample_id HAVING c>1
								           });
	
	$sth->execute();
	
	while (my @res = ($sth->fetchrow_array)) {
		if ($res[1] =~ /^(Sample|Subject) from the DGVa study/) {
			my $stmt_s = qq{ UPDATE sample SET description='Subject from several DGVa studies'
	                     WHERE sample_id=$res[0]
								 		 };
			$dbVar->do($stmt_s);
		}
	}
}




#### Finishing methods ####

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


sub cleanup {
  debug(localtime()." Cleanup temporary columns");

  # Column tmp_class_name" in structural_variation
	my $sth1 = $dbVar->prepare(qq{ SELECT count(*) FROM $sv_table WHERE source_id=$source_id AND class_attrib_id=0});
  $sth1->execute();
	my $sv_count = ($sth1->fetchrow_array)[0];
	$sth1->finish;
  if ($sv_count != 0) {
	  print STDERR "\tThe table $sv_table has $sv_count variants with no class_attrib_id defined!\nPlease, look at the column '$tmp_sv_col'\n";
	} 
	else {
	  $dbVar->do(qq{ALTER TABLE $sv_table DROP COLUMN $tmp_sv_col});
		debug(localtime()." Table $sv_table: cleaned");
	}
	
	# Column tmp_clinic_name" in structural_variation_annotation
	my $sth2 = $dbVar->prepare(qq{ SELECT count(*) FROM $svanno_table WHERE clinical_attrib_id is NULL AND $tmp_sva_col is not NULL});
  $sth2->execute();
	my $sva_count = ($sth2->fetchrow_array)[0];
	$sth2->finish;
  if ($sva_count != 0) {
	  print STDERR "\tThe table $svanno_table has $sva_count variants with no clinical_attrib_id value where a clinical significance is defined!\nPlease, look at the column '$tmp_sva_col'\n";
	} 
	else {
	  $dbVar->do(qq{ALTER TABLE $svanno_table DROP COLUMN $tmp_sva_col});
		debug(localtime()." Table $sva_table: cleaned");
	}
}




#### Other methods ####

sub remove_data {
  my $study_id = shift;
	
  debug(localtime()." Remove study data");

  # structural_variation_annotation
  $dbVar->do(qq{ DELETE from $svanno_table WHERE structural_variation_id IN 
	               (SELECT structural_variation_id FROM $sv_table WHERE study_id=$study_id);
						});
	# structural_variation_association					
  $dbVar->do(qq{ DELETE from $sva_table WHERE structural_variation_id IN 
	               (SELECT structural_variation_id FROM $sv_table WHERE study_id=$study_id)
						});	
	# structural_variation_feature									
  $dbVar->do(qq{ DELETE from $svf_table WHERE structural_variation_id IN 
	               (SELECT structural_variation_id FROM $sv_table WHERE study_id=$study_id)
						});
	# failed_structural_variation
	$dbVar->do(qq{ DELETE from $sv_failed WHERE structural_variation_id IN 
	               (SELECT structural_variation_id FROM $sv_table WHERE study_id=$study_id)
						});
	# structural_variation
	$dbVar->do(qq{ DELETE from $sv_table WHERE study_id=$study_id });									
}


sub generate_data_row {
  my $info = shift;
	my $somatic = shift;
 
  my $bp_order = $info->{bp_order};
  if ($info->{is_somatic} == 1 || $somatic) {
    $bp_order = 1 if(!defined($bp_order));
	} else {
	  $bp_order = undef;
	}
 
  my @row = ($info->{ID},
						 $info->{SO_term},
		         $info->{chr}, 
						 $info->{outer_start}, 
						 $info->{start}, 
						 $info->{inner_start}, 
						 $info->{inner_end}, 
						 $info->{end}, 
						 $info->{outer_end}, 
						 $info->{parent},
						 $info->{clinical},
						 $info->{phenotype},
						 $info->{sample},
						 $info->{strain_name},
						 $info->{is_ssv}, 
						 $info->{is_failed},
						 $info->{population},
						 $bp_order,
						 $info->{is_somatic},
						 $info->{status}
			  		);
  return \@row;
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
  -replace         : flag to remove the existing study data from the database before import them (optional)
  -debug           : flag to keep the temp_cnv table (optional)
  };
}
