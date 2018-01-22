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
     $target_assembly, $cs_version_number, $size_diff, $version, $registry_file, $medgen_file, $hpo_file, $replace_study, $debug);

GetOptions('species=s'         => \$species,
           'source_name=s'     => \$source_name,
           'input_file=s'      => \$input_file,
           'input_dir=s'       => \$input_dir,
           'tmpdir=s'          => \$ImportUtils::TMP_DIR,
           'tmpfile=s'         => \$ImportUtils::TMP_FILE,
           'var_set=i'         => \$var_set_id,
           'mapping!'          => \$mapping,
           'gaps=i'            => \$num_gaps,
           'target_assembly=s' => \$target_assembly,
           'size_diff=i'       => \$size_diff,
           'version=i'         => \$version,
           'registry=s'        => \$registry_file,
           'medgen_file=s'     => \$medgen_file,
           'hpo_file=s'        => \$hpo_file,
           'replace!'          => \$replace_study,
           'debug!'            => \$debug,
          );
$registry_file ||= $Bin . "/ensembl.registry";

$num_gaps = 1 unless defined $num_gaps;

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

my %study_to_skip =  ( 'nstd8'  => 1, # Human & chimp
                       'nstd9'  => 1, # Chimp
                       'nstd82' => 1, # Several monkey species
                       #'nstd95' => 1, # Mus musculus domesticus
                     );

my @attribs = ('ID','SO_term','chr','outer_start','start','inner_start','inner_end','end',
               'outer_end','strand','parent','clinical','subject','sample','gender','is_ssv',
               'is_failed','population','bp_order','is_somatic','status','alias','length','copy_number','zygosity');

my %attribs_col = ('ID'          => 'id *',
                   'SO_term'     => 'type',
                   'chr'         => 'chr *',
                   'outer_start' => 'outer_start i',
                   'start'       => 'start i',
                   'inner_start' => 'inner_start i',
                   'inner_end'   => 'inner_end i',
                   'end'         => 'end i',
                   'outer_end'   => 'outer_end i',
                   'strand'      => 'strand i',
                   'parent'      => 'parent *',
                   'clinical'    => 'clinic',
                   'subject'     => 'subject *', # or strain
                   'sample'      => 'sample *',
                   'gender'      => 'gender',
                   'is_ssv'      => 'is_ssv i',
                   'is_failed'   => 'is_failed i',
                   'population'  => 'population',
                   'bp_order'    => 'bp_order i',
                   'is_somatic'  => 'is_somatic i',
                   'status'      => 'status',
                   'alias'       => 'alias',
                   'length'      => 'length i',
                   'copy_number' => 'copy_number', # Not 'i' because we need to differenciate between "0 copy" and the "default value (0)"
                   'zygosity'    => 'zygosity i'
                  );


my %display_name = ( 'COSMIC' => 'http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/' );

my $add             = '';
my $study_table     = "study$add";
my $sv_table        = "structural_variation$add";
my $svf_table       = "structural_variation_feature$add";
my $sva_table       = "structural_variation_association$add";
my $sv_failed       = "failed_structural_variation$add";
my $set_table       = "variation_set_structural_variation$add";
my $svs_table       = "structural_variation_sample$add";
my $pf_table        = "phenotype_feature$add";
my $temp_table      = "temp_sv$add";
my $temp_phen_table = "temp_sv_phenotype$add";
my $temp_clin_table = "temp_sv_clin_sign$add";

my $tmp_sv_col      = 'tmp_class_name';
my $tmp_sv_clin_col = 'tmp_clinic_name';
  
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
my (%num_mapped, %num_not_mapped, %samples, %subjects, %study_done);
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

pre_processing();

foreach my $in_file (@files) {
  next if ($in_file !~ /\.gvf$/);
  
  $in_file =~ /^(\w{4}\d+)_/;
  if ($study_to_skip{$1}) {
    debug("Study $1 (".$in_file.") skipped because it contains non ".$species." data");
    next;
  }

  my $msg ="File: $in_file ($f_count/$f_nb)";
  print "$msg\n";
  debug("\n##############################\n$msg\n##############################");
  
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
  failed_structural_variation();
  structural_variation_feature();
  structural_variation_association();
  somatic_study_processing() if ($somatic_study == 1);
  structural_variation_set() if ($var_set_id);
  structural_variation_sample();
  phenotype_feature();
  drop_tmp_table() if (!defined($debug));
  debug(localtime()." Done!\n");
  $f_count ++;
}

## Post processings ##

# Post processing for mouse annotation (delete duplicated entries in structural_variation_sample)
post_processing_annotation() if ($species =~ /mouse|mus/i);
post_processing_feature();
post_processing_sample();
post_processing_phenotype();
post_processing_clinical_significance();
post_processing_failed_variants();

# Finishing methods
meta_coord();
verifications(); # URLs
cleanup() if (!defined($debug));

debug(localtime()." All done!");



#############
#  Methods  #
#############

sub load_file_data{
  debug(localtime()." Loading file into temporary table");
  
  # drop any table that may be there
  $dbVar->do("DROP TABLE IF EXISTS $temp_table;");
  $dbVar->do("DROP TABLE IF EXISTS $temp_phen_table;");
  
  my @cols = map { $attribs_col{$_}} @attribs;
  create_and_load($dbVar, "$temp_table", @cols);
  
  # Move the 2nd file in $TMP_DIR/$TMP_FILE before loading the 2nd table
  `mv $TMP_DIR/$TMP_FILE\2 $TMP_DIR/$TMP_FILE`;
  create_and_load( $dbVar, "$temp_phen_table", "id *", "phenotype *");
  
  # fix nulls
  if (($dbVar->selectrow_arrayref(qq{SELECT count(*) FROM $temp_table WHERE (outer_start=0 AND start=0 AND inner_start=0) OR (inner_end=0 AND end=0 AND outer_end=0);}))->[0] != 0) {
    warn "Structural variants with start and/or end coordinates equal to 0 found in the file!";
    $dbVar->do(qq{UPDATE $temp_table SET outer_start=outer_end, start=end, inner_start=inner_end WHERE outer_start=0 AND start=0 AND inner_start=0});
    $dbVar->do(qq{UPDATE $temp_table SET outer_end=outer_start, end=start, inner_end=inner_start WHERE outer_end=0 AND end=0 AND inner_end=0});
  }
  
  foreach my $coord('outer_start', 'inner_start', 'inner_end', 'outer_end', 'start', 'end', 'bp_order') {
    $dbVar->do(qq{UPDATE $temp_table SET $coord = NULL WHERE $coord = 0;});
  }
  
  
  foreach my $col('clinic', 'sample', 'subject', 'status', 'alias', 'length', 'copy_number', 'zygosity') {
    $dbVar->do(qq{UPDATE $temp_table SET $col = NULL WHERE $col = '';});
  }
  
  # Case with insertions            
  $dbVar->do(qq{UPDATE $temp_table SET start=outer_start, end=inner_start, 
                outer_end=inner_start, inner_start=NULL, inner_end=NULL
                WHERE outer_start=outer_end AND inner_start=inner_end;});
}                

sub source {
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
  my $study        = $data->{study};
  my $assembly     = $data->{assembly};
  my $study_desc   = $data->{desc};
  my $study_type   = $data->{study_type};
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
  
  #  PubMed ID
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
      $pubmed = "PMID:".join(",PMID:", @$pmid);
    }
  }
  my $external_link = ($display_name{$author}) ? $display_name{$author} : $pubmed;
  
  # URL
  $study =~ /(\w+\d+)\.?\d*/;
  my $study_ftp = $1;
  $study_ftp = "ftp://ftp.ebi.ac.uk/pub/databases/dgva/$study_ftp\_$author";
  
  
  my $assembly_desc = " [remapped from build $assembly]" if ($mapping and $assembly ne $target_assembly);
  
  $stmt = qq{ SELECT st.study_id, st.description, st.external_reference FROM study st, source s 
              WHERE s.source_id=st.source_id AND s.name='DGVa' and st.name='$study'};
  my $rows = $dbVar->selectall_arrayref($stmt);    
  
  my $assembly_desc;
  
    
  # UPDATE
  if (scalar (@$rows)) {
    $study_id   = $rows->[0][0];
    $study_desc = $rows->[0][1];
    my $study_xref = $rows->[0][2];
    
    $study_desc =~ s/'/\\'/g;

    remove_data($study_id) if (defined($replace_study) && defined($study_id) && !defined($study_done{$study}));
    
    # Checks the studies with different assemblies when the option "mapping" is selected.
    if ($mapping and $assembly ne $target_assembly) {
      if ($study_desc =~ /^(.+)\s\[remapped\sfrom\sbuild(s?)\s(.+)\]$/) {
        my $gen_desc = $1;
        $assembly_desc = $3;
      
        if ($assembly_desc !~ /$assembly/) {
          $assembly_desc .= " and $assembly";
          $study_desc = "$gen_desc [remapped from builds $assembly_desc]";
        }
      }
      else {
        $study_desc = "$study_desc [remapped from build $assembly]";
      }
    }
    elsif ($assembly eq $target_assembly && $study_desc =~ /^(.+)\s\[remapped/) {
      $study_desc = $1;
    }
    
    my $external_link_sql;
    if ($external_link =~ /NULL/i) {
      if ($study_xref && $study_xref !~ /NULL/i && $study_xref ne '') {
        $external_link_sql = '';
      }
      else {
        $external_link_sql = "external_reference='$external_link',";
      }
    }
    else {
      $external_link_sql = "external_reference='$external_link',";
    }

    $stmt = qq{ UPDATE $study_table SET 
                  description='$study_desc',
                  $external_link_sql
                  study_type="$study_type",
                  url="$study_ftp"
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
          "$study_desc",
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


# Structural variations & supporting structural variations
sub structural_variation {
  
  debug(localtime()." Inserting into $sv_table table");
  
  my $stmt = qq{
    INSERT IGNORE INTO
    $sv_table (
      variation_name,
      source_id,
      study_id,
      class_attrib_id,
      $tmp_sv_col,
      clinical_significance,
      $tmp_sv_clin_col,
      is_evidence,
      somatic,
      validation_status,
      alias,
      copy_number
    )
    SELECT
      DISTINCT
      t.id,
      $source_id,
      $study_id,
      a1.attrib_id,
      t.type,
      t.clinic,
      t.clinic,
      t.is_ssv,
      t.is_somatic,
      t.status,
      t.alias,
      t.copy_number
    FROM
      $temp_table t 
      LEFT JOIN attrib a1 ON (t.type=a1.value)
  };
  $dbVar->do($stmt);
}

  
# Failed structural variations
sub failed_structural_variation   {  
  debug(localtime()." Inserting into $sv_failed table for SVs");
  
  my $stmt = qq{
    INSERT IGNORE INTO
    $sv_failed (
      structural_variation_id,
      failed_description_id
    )
    SELECT
      sv.structural_variation_id,
      17
    FROM
      $sv_table sv,
      $temp_table t
    WHERE
      sv.variation_name=t.id AND
      t.is_ssv=0 AND 
      (t.is_failed!=0 OR t.chr NOT IN (SELECT name FROM seq_region))
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
      $sv_table sv,
      $temp_table t
    WHERE
      sv.variation_name=t.id AND
      t.is_ssv=1 AND 
      (t.is_failed!=0 OR t.chr NOT IN (SELECT name FROM seq_region))
  };
  $dbVar->do($stmt);
}
  
  
# Structural variation features & supporting structural variation features
sub structural_variation_feature {
  debug(localtime()." Inserting into $svf_table table");
  
  my $stmt = qq{
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
      study_id,
      is_evidence,
      somatic,
      breakpoint_order,
      length
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
      $study_id,
      t.is_ssv,
      t.is_somatic,
      t.bp_order,
      t.length
    FROM
      seq_region q,
      $temp_table t,
      $sv_table sv
    WHERE
      q.name = t.chr AND
      t.id=sv.variation_name AND
      t.is_failed=0 AND
      sv.source_id=$source_id
  };
  $dbVar->do($stmt);
}

  
# Structural variation association
sub structural_variation_association {
  debug(localtime()." Inserting into $sva_table table");
  my $stmt = qq{
    INSERT IGNORE INTO
    $sva_table (
      structural_variation_id,
      supporting_structural_variation_id
    )
    SELECT
      sv2.structural_variation_id,
      sv1.structural_variation_id
    FROM
      $temp_table t1,
      $temp_table t2,
      $sv_table sv1,
      $sv_table sv2
    WHERE
      t1.id=sv1.variation_name AND
      t1.parent=t2.id AND
      t2.id=sv2.variation_name
  };
  $dbVar->do($stmt);
}


# Somatic study processing  
sub somatic_study_processing {  
  # UPDATE breakpoint_order to 1 when more than one SVF and no breakpoint order
  my $stmt = qq{
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
    
  # Move alias from SSV to SV
  $stmt = qq{ 
      UPDATE $sv_table sv, $sv_table ssv, $sva_table sva SET sv.alias=ssv.alias
      WHERE sv.structural_variation_id=sva.structural_variation_id
        AND ssv.structural_variation_id=sva.supporting_structural_variation_id
        AND sv.study_id=$study_id
        AND sv.is_evidence=0
        AND ssv.alias is not null
  };
  $dbVar->do($stmt);
  $stmt = qq{
    UPDATE $sv_table SET alias=NULL WHERE study_id=$study_id AND is_evidence=1;
  };
  $dbVar->do($stmt);
}


sub structural_variation_set {
  debug(localtime()." Inserting into $set_table table");
  
  my $stmt;
  my @var_set_list;
  
  if ($var_set_id eq '1kg') {
  
    my %pop_1kg_p1 = (
      'ASW' => 'AFR', 'LWK' => 'AFR', 'YRI' => 'AFR',                                      # AFR
      'CLM'  => 'AMR', 'PEL'  => 'AMR', 'PUR'  => 'AMR', 'MXL'  => 'AMR',                  # AMR
      'CEU'  => 'EUR', 'TSI'  => 'EUR', 'FIN'  => 'EUR', 'GBR'  => 'EUR', 'IBS'  => 'EUR', # EUR
      'CHB'  => 'ASN', 'JPT'  => 'ASN', 'CHS'  => 'ASN', 'CDX'  => 'ASN', 'KHV'  => 'ASN', # ASN
    );
    
    while (my($sub_pop, $sup_pop) = each (%pop_1kg_p1)) {
      $dbVar->do(qq{UPDATE $temp_table SET population='$sup_pop' WHERE population='$sub_pop'});
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
         $temp_table t,
        $sv_table sv,
        $sva_table sva,
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
         $temp_table t,
        $sv_table sv,
        $sva_table sva,
        variation_set vs
      WHERE
        t.id=sv.variation_name AND
        sva.supporting_structural_variation_id=sv.structural_variation_id AND
        t.population is not null AND
        vs.name='1000 Genomes - All'
    };
    $dbVar->do($stmt);
    
    # Variation set "High quality"
    $stmt = qq{
       INSERT IGNORE INTO
       $set_table (
        structural_variation_id,
        variation_set_id
      )
      SELECT DISTINCT
         sv.structural_variation_id,
         vs.variation_set_id
      FROM
         $temp_table t,
        $sv_table sv,
        variation_set vs
      WHERE
        t.id=sv.variation_name AND
        sv.is_evidence=0 AND
        t.population is not null AND
        t.status = 'High quality' AND
        vs.name='1000 Genomes - High quality'
    };
    $dbVar->do($stmt);
    
    
    my $sets_1kg = $dbVar->selectrow_arrayref(qq{
         SELECT DISTINCT vssv.variation_set_id 
         FROM variation_set_structural_variation vssv, structural_variation sv, $temp_table t 
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
         $temp_table t,
        $sv_table sv
      WHERE
        t.is_ssv=0 AND
        t.id=sv.variation_name
    };
    $dbVar->do($stmt);
    
    push(@var_set_list,$var_set_id);
  }  
  
  my $set_id_list = join(',',@var_set_list);
  
  # Populate the meta table
  my $sets = $dbVar->selectall_arrayref(qq{SELECT vs.name,a.value FROM variation_set vs, attrib a 
                                          WHERE vs.variation_set_id IN ($set_id_list) 
                                          AND a.attrib_id=vs.short_name_attrib_id;}
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


sub structural_variation_sample {
  my $stmt;
  
  debug(localtime()." Inserting into $svs_table table");

  my $study_name = ($dbVar->selectrow_arrayref(qq{SELECT name FROM study WHERE study_id=$study_id;}))->[0];

  # Create strain entries
  if ($species =~ /mouse|mus/i) {
    $stmt = qq{ SELECT DISTINCT subject FROM $temp_table 
                WHERE subject NOT IN (SELECT DISTINCT name from individual WHERE individual_type_id=1)
              };
    my $rows_strains = $dbVar->selectall_arrayref($stmt);
    foreach my $row (@$rows_strains) {
      my $strain = $row->[0];
      next if ($strain eq  '');
    
      if (!$dbVar->selectrow_arrayref(qq{SELECT individual_id FROM individual WHERE name='$strain' LIMIT 1})) {
        $dbVar->do(qq{ INSERT IGNORE INTO individual (name,description,gender,individual_type_id) VALUES ('$strain','Strain from the DGVa study $study_name','Unknown',1)});
      }
      else{
        if ($dbVar->selectrow_arrayref(qq{SELECT individual_id FROM individual WHERE name='$strain' AND individual_type_id!=1})) {
          $dbVar->do(qq{UPDATE individual SET individual_type_id=1 WHERE name='$strain' AND individual_type_id!=1});
        }
      }
    }
    
    # Create sample entries
    $stmt = qq{ SELECT DISTINCT sample, subject FROM $temp_table WHERE is_ssv=1 AND sample NOT IN (SELECT DISTINCT name from sample WHERE study_id=$study_id AND display!="UNDISPLAYABLE")};
    my $rows_samples = $dbVar->selectall_arrayref($stmt);
    foreach my $row (@$rows_samples) {
      my $sample  = $row->[0];
      my $subject = $row->[1];
      next if ($sample eq  '' || $subject eq '');

      $dbVar->do(qq{ INSERT IGNORE INTO sample (name,description,study_id,display,individual_id) SELECT '$sample','Sample from the DGVa study $study_name', $study_id,"MARTDISPLAYBLE",min(individual_id) FROM individual WHERE name='$subject' LIMIT 1});
    }
  }
  # Create individual entries (not for mouse)
  else { 
    $stmt = qq{ SELECT DISTINCT subject, gender FROM $temp_table WHERE is_ssv=1 AND subject NOT IN (SELECT DISTINCT name from individual)};
    my $rows_subjects = $dbVar->selectall_arrayref($stmt);
    foreach my $row (@$rows_subjects) {
      my $subject = $row->[0];
      my $gender = ($row->[1] =~ /\w+/) ? $row->[1] : 'Unknown';
      next if ($subject eq  '');

      my $itype_val = ($species =~ /human|homo/i) ? 3 : 2;
    
      $dbVar->do(qq{ INSERT IGNORE INTO individual (name,description,gender,individual_type_id) VALUES ('$subject','Subject from the DGVa study $study_name','$gender',$itype_val)});
    }
  
    # Create sample entries
    $stmt = qq{ SELECT DISTINCT sample, subject FROM $temp_table WHERE is_ssv=1 AND sample NOT IN (SELECT DISTINCT name from sample)};
    my $rows_samples = $dbVar->selectall_arrayref($stmt);
    foreach my $row (@$rows_samples) {
      my $sample  = $row->[0];
      my $subject = $row->[1];
      next if ($sample eq  '' || $subject eq '');

      #$dbVar->do(qq{ INSERT IGNORE INTO sample (name,description,study_id,individual_id) SELECT '$sample','Sample from the DGVa study $study_name', $study_id, min(individual_id) FROM individual WHERE name='$subject'});
      $dbVar->do(qq{ INSERT IGNORE INTO sample (name,description,individual_id) SELECT '$sample','Sample from the DGVa study $study_name', min(individual_id) FROM individual WHERE name='$subject'});
    }
  }

  $dbVar->do(q{OPTIMIZE TABLE sample});
  $dbVar->do(q{OPTIMIZE TABLE individual});

  my $ext1;
  my $ext2;
  
  # For mouse
  if ($species =~ /mouse|mus/i) {
    $stmt = qq{
      INSERT IGNORE INTO $svs_table (
        structural_variation_id,
        sample_id,
        zygosity
      )
      SELECT DISTINCT
        sv.structural_variation_id,
        s.sample_id,
        t.zygosity
      FROM
        $sv_table sv,
        $temp_table t
        LEFT JOIN sample s ON (s.name=t.sample AND s.sample_id IN (select min(sample_id) from sample s2 where s.name=s2.name) AND s.display!='UNDISPLAYABLE')
      WHERE
        sv.variation_name=t.id AND
        (t.sample is not null OR t.subject is not null)
      };
  }
  else {
    $stmt = qq{
      INSERT IGNORE INTO $svs_table (
        structural_variation_id,
        sample_id,
        zygosity
      )
      SELECT DISTINCT
        sv.structural_variation_id,
        s.sample_id,
        t.zygosity
      FROM
        $sv_table sv,
        $temp_table t
        LEFT JOIN sample s ON (s.name=t.sample AND s.sample_id IN (select min(sample_id) from sample s2 where s.name=s2.name))
      WHERE
        sv.variation_name=t.id AND
        t.sample is not null
      };
  }
  $dbVar->do($stmt);
}


sub phenotype_feature {
  my $stmt;
  
  # Check if there is some phenotype entries
  $stmt = qq{ SELECT count(phenotype) FROM $temp_phen_table WHERE phenotype is not null };
  if (($dbVar->selectall_arrayref($stmt))->[0][0] == 0) {
    debug(localtime()." Inserting into $pf_table table skipped: no phenotype data for this study");
    return;
  }
  
  debug(localtime()." Inserting into $pf_table table");
  
  # Create phenotype entries
  $stmt = qq{ SELECT DISTINCT phenotype FROM $temp_phen_table 
              WHERE phenotype NOT IN (SELECT description FROM phenotype)
            };
  my $rows = $dbVar->selectall_arrayref($stmt);
  while (my $row = shift(@$rows)) {
    my $phenotype = $row->[0];
    next if ($phenotype eq '');
    
    $phenotype =~ s/'/\\'/g;
    $dbVar->do(qq{ INSERT INTO phenotype (description) VALUES ('$phenotype')});  
  }
    
  $stmt = qq{
    INSERT IGNORE INTO 
    $pf_table (
      phenotype_id,
      object_id,
      source_id,
      study_id,
      type,
      seq_region_id,
      seq_region_start,
      seq_region_end,
      seq_region_strand
    )
    SELECT DISTINCT
      p.phenotype_id,
      t.id,
      $source_id,
      $study_id,
      IF(t.is_ssv=1,'SupportingStructuralVariation','StructuralVariation'),
      svf.seq_region_id,
      svf.seq_region_start,
      svf.seq_region_end,
      svf.seq_region_strand
    FROM
      $temp_table t, 
      $temp_phen_table tp,
      $svf_table svf,
      phenotype p
    WHERE 
      svf.variation_name=t.id AND
      t.id=tp.id AND
      p.description=tp.phenotype
  };
  $dbVar->do($stmt);
}  


sub drop_tmp_table {
  $dbVar->do(qq{DROP TABLE $temp_table;}) unless ($debug);
  $dbVar->do(qq{DROP TABLE $temp_phen_table;}) unless ($debug);
}




#### Parsing methods ####


sub parse_gvf {
  my $infile = shift;
  
  debug(localtime()." Parse GVF file: $infile");
  debug(localtime()." Start mapping SV features to $target_assembly if needed") if $mapping;
  
  open IN, $infile or die "Could not read from input file $infile\n";
  
  # get a slice adaptor
  my $sa = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'slice');
  
  open OUT,  ">$TMP_DIR/$TMP_FILE" or die "Could not write to output file $TMP_DIR/$TMP_FILE\n";
  open OUT2, ">$TMP_DIR/$TMP_FILE\2" or die "Could not write to output file $TMP_DIR/$TMP_FILE\2\n";
  
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
      study_table($header);
    }
  
    chomp ($current_line);
    my $line_info = parse_line($current_line);
    
    my $infos = check_breakpoint_coordinates($line_info);
    
    my $is_failed = 0;
    
    
    foreach my $info (@$infos) {
    
      ###### Mapping work ######  
      if ($mapping && $assembly !~ /$cs_version_number/) {
        
        # get a slice on the old assembly
        my $slice;
        my $start_c = $info->{start};
        my $end_c   = $info->{end};
  
        if (!$info->{start}) {
          $start_c = ($info->{outer_start}) ? $info->{outer_start} : $info->{inner_start};
        }
        if (!$info->{end}) {
          $end_c = ($info->{outer_end}) ? $info->{outer_end} : $info->{inner_end};
        }
  
        eval { $slice = $sa->fetch_by_region('chromosome', $info->{chr}, $start_c, $end_c, 1, $assembly); };
  
        # check got the slice OK
        if(!defined($slice)) {
           warn("Structural variant '".$info->{ID}."' (study ".$header->{study}."): Unable to map from assembly $assembly or unable to retrieve slice ".$info->{chr}."\:$start_c\-$end_c");
          $skipped++;
          $num_not_mapped{$assembly}++;
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
              warn ("Structural variant '".$info->{ID}."' (study ".$header->{study}.") has location '$assembly:".$info->{chr}.":$start_c\-$end_c' , which could not be re-mapped to $target_assembly. This variant will be labelled as failed");
              $num_not_mapped{$assembly}++;
            }  
            $is_failed = 1;
          }
        }
      }
    
      $info->{is_failed} = $is_failed;
               
      my $data = generate_data_row($info);           
      print OUT (join "\t", @{$data})."\n";
      
      if ($info->{phenotype}) {
        my @phenotypes = map { decode_text($_) } keys(%{$info->{phenotype}});
        foreach my $phenotype (@phenotypes) {
          my $phen_data = generate_phenotype_data_row($info,$phenotype);
          print OUT2 (join "\t", @{$phen_data})."\n";
        }
      }
    }
  }
  close IN;
  close OUT;
  close OUT2;
  
  if ($mapping && $assembly !~ /$cs_version_number/) {
    debug(localtime()." Finished SV mapping:\n\t\tSuccess: ".(join " ", %num_mapped)." not required $no_mapping_needed\n\t\tFailed: ".(join " ", %num_not_mapped)." (In no existing chromosome: $skipped)");
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
    $line =~ /^(.+)\:\s+(.+)$/;
    $label = $1;
    $info  = $2;
  }
  
  $label =~ s/#//g;
  $label =~ s/^\s//;
  $info  =~ s/^\s+//;
  
  $h->{author}       = $info if ($label =~ /Display.+name/i);
  $h->{first_author} = $info if ($label =~ /First.+author/i);
  $h->{assembly}     = $info if ($label =~ /Assembly.+name/i);
  $h->{study}        = $info if ($label =~ /Study.+accession/i);
  $h->{study_type}   = $info if ($label =~ /Study.+type/i);
  
  $somatic_study = 1 if ($h->{study_type} =~ /(somatic)|(tumor)/i);
  
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
    my ($sample,$subject,$tissue, $pop);
    foreach my $s_info (split(';',$info)) {
      my ($key,$value) = split('=',$s_info);
      $sample  = $value if ($key eq 'sample_name');
      $subject = $value if ($key eq 'subject_name');
      
      # Population (mouse)
      if ($key eq 'links' && $value =~ /^JAX:(.+)$/i) {
        $pop = $1;
      } 
      
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
    
    
    if (defined($sample)) {
      $subject = defined($subject) ? $subject : $sample;
      
      $samples{$sample}{subject}    = $subject;
      $samples{$sample}{tissue}     = $tissue;
      $samples{$sample}{gender}     = $subjects{$subject}->{subject_sex};
      $samples{$sample}{phenotype}  = $subjects{$subject}->{phenotype_description};
      $samples{$sample}{population} = $subjects{$subject}->{subject_population};
      
      if (!defined($subjects{$subject}->{subject_population}) && $species =~ /mouse|mus.*mus/ && defined($pop)){
        $samples{$sample}{population} = $pop;
      }
    }
  }
  
  # Subject information
  elsif ($label eq 'subject') {
    my %subject_data;
    
    foreach my $s_info (split(';',$info)) {
      my ($key,$value) = split('=',$s_info);
      $subject_data{$key} = $value;
    }
    
    $subject_data{subject_sex} = $subject_data{sex} if ($subject_data{sex});
    $subject_data{subject_population} = $subject_data{population} if ($subject_data{population});
    
    if ($subject_data{subject_sex}) {
      $subject_data{subject_sex} = ($subject_data{subject_sex} =~ /^m/i) ? 'Male' : 'Female';
    }
    
    my $subject = $subject_data{subject_name};
    if (defined($subject)) {
      $subjects{$subject} = \%subject_data;
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
    elsif($assembly =~ /37/) {
      $assembly = 'GRCh37';
    }
    elsif($assembly !~ /$cs_version_number/) {
      warn("Could not guess assembly from file - assuming to be GRCh38");
      $assembly = 'GRCh38';
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
  elsif($species =~ /dog|can.*fam/) {
    $assembly = 'BROADD2' if ($assembly =~ /2\.0/);
  }
  $h->{assembly} = $assembly;
  
  return $assembly;
}

sub parse_line {
  my $line     = shift;
  my @data     = split(/\t/, $line);
  my @last_col = split(';', pop(@data));
  my $info;
  
  $data[0] =~ /Chr(.+)$/;
  $info->{chr} = defined($1) ? $1 : $data[0];
  $info->{SO_term} = ($data[2] eq 'alu_insertion') ? ucfirst($data[2]) : $data[2];
  $info->{start}   = $data[3];
  $info->{end}     = $data[4];
  $info->{strand}  = ($data[6] eq '.') ? 1 : ($data[6] eq '-') ? -1 : 1;

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
    
    # Sample information
    if ($value !~ /Unknown/i) {
      if ($key eq 'sample_name' || ($key eq 'subject_name' && !$info->{sample})) {
        if ($species =~ /dog|can.*fam/) {
          $info->{sample} = ($samples{$value}{population}) ? $samples{$value}{population} : $value;
        } else {
          $info->{sample} = ($samples{$value}{subject}) ? $samples{$value}{subject} : $value;
        }
        if (scalar(keys(%subjects)) == 1 ) {
          my $subject_name = (keys(%subjects))[0];
          $info->{subject} = $subject_name;
        } else {
          $info->{subject} = ($samples{$value}{subject}) ? $samples{$value}{subject} : $info->{sample};
        }
        $info->{population}  = ($samples{$value}{population}) ? $samples{$value}{population} : undef; # 1000 Genomes study
        $info->{gender}      = ($samples{$value}{gender}) ? $samples{$value}{gender} : undef;
        
        my $s_phenotype = ($samples{$value}{phenotype}) ? $samples{$value}{phenotype} : $samples{$value}{tissue};
        $info->{phenotype}{$s_phenotype} = 1;
      }
    }
    
    $info->{clinical}    = lc($value) if ($key eq 'clinical_significance');
    $info->{parent}      = $value if ($key eq 'Parent'); # Check how the 'parent' key is spelled
    $info->{is_somatic}  = 1 if ($key eq 'var_origin' && $value =~ /somatic/i);
    $info->{bp_order}    = ($info->{submitter_variant_id} =~ /\w_(\d+)$/) ? $1 : undef;
    $info->{bp_order}    = 1 if ($info->{SO_term} =~ /translocation/i);
    $info->{status}      = 'High quality' if ($key eq 'variant_region_description' && $value =~ /high.quality/i);
    $info->{alias}       = $value if ($key eq 'Alias' && $value !~ /^\d+$/);
    $info->{length}      = $value if ($key eq 'variant_call_length');
    
    # Zygosity
    if ($key eq 'Zygosity') {
      if ($value eq 'heterozygous') {
        $info->{zygosity} = 1;
      }
      elsif ($value eq 'homozygous') {
        $info->{zygosity} = 2;
      }
    }
    
    # Copy number
    if ($key eq 'copy_number') {
      my $copy = (split(/[-,\.\~]/,$value))[0];
         $copy =~ s/\D//g; # Remove non numeric character
      $info->{copy_number} = $copy if ($copy ne '');
    }

    # Breakpoint definition
    if ($species =~ /^(homo|human)/i && $info->{is_somatic}) {
      $info->{_bp_detail} = $value if ($key =~ /Breakpoint_detail/i);
      $info->{_bp_range}  = $value if ($key =~ /Breakpoint_range/i);
    }
    
    # Phenotype
    $info->{phenotype}{$value} = 1 if ($key eq 'phenotype_description');
    $info->{phenotype_link} = $value if ($key eq 'phenotype_link');
  }
  
  
  ## Phenotype(s) ##
  
  if ($info->{phenotype_link}) {
    
    my @phenotype_links = split(',',$info->{phenotype_link});

    foreach my $p_link (@phenotype_links) {
      # Look at the MedGen data
      if ($medgen_file && $p_link =~ /^MedGen:\w+$/i) {
        $p_link =~ /^MedGen:(\w+)$/i;
        my $phenotype_value = get_medgen_phenotype($1);
        $info->{phenotype}{$phenotype_value} = 1 if (defined($phenotype_value) && $phenotype_value ne '');
      }
    
      # Look at the HP ontology data
      if ($hpo_file && $p_link =~ /^HP:\d+$/) {
        my $phenotype_value = get_hpo_phenotype($p_link);
        $info->{phenotype}{$phenotype_value} = 1 if (defined($phenotype_value) && $phenotype_value ne '');
      }
    }
  }
  
  # Flag to know if the entry is a sv or a ssv
  $info->{is_ssv} = ($info->{ID} =~ /ssv/) ? 1 : 0;
  
  # Somatic alias
  if ($info->{is_somatic} == 1 && $info->{alias} =~ /^(.+)_\d+$/) {
    $info->{alias} = $1;
  }
  
  # Check if a SV has its SSV(s) labelled as 'somatic'
  if ($somatic_study && $info->{is_ssv} == 0) {
    my $cmd  = "grep 'parent=".$info->{ID}."' $fname";
    my $text = `$cmd`;
    $info->{is_somatic} = 1 if ($text =~ /var_origin=Somatic/i || !$text);
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


#### Pre processing ####
sub pre_processing {

  debug(localtime()." Prepare the structural variation tables by adding temporary keys and columns");
  debug(localtime()."\t - Add temporary unique keys in $svf_table and $svs_table tables");
  
  
  # Prepare the structural_variation_feature table
  if ($dbVar->do(qq{SHOW KEYS FROM $svf_table WHERE Key_name='name_coord_key';}) < 1){
    $dbVar->do(qq{ALTER TABLE $svf_table ADD CONSTRAINT  UNIQUE KEY
                  `name_coord_key` (`variation_name`,`seq_region_id`,`seq_region_start`,`seq_region_end`)
                 });
  }
  if ($dbVar->do(qq{SHOW KEYS FROM $svf_table WHERE Key_name='name_key';}) < 1){
    $dbVar->do(qq{ALTER TABLE $svf_table ADD KEY `name_key` (`variation_name`)});
  }
  
  # Prepare the structural_variation_sample table
  if ($dbVar->do(qq{SHOW KEYS FROM $svs_table WHERE Key_name='sv_id_key';}) < 1){
    $dbVar->do(qq{ALTER TABLE $svs_table ADD CONSTRAINT UNIQUE KEY `sv_id_key` (`structural_variation_id`,`sample_id`)});
  }
  
  # Prepare the individual table
  debug(localtime()."\t - Add temporary unique keys in the individual table");
  if ($dbVar->do(qq{SHOW KEYS FROM individual WHERE Key_name='name_key';}) < 1){
    $dbVar->do(qq{ALTER TABLE individual ADD KEY `name_key` (`name`)});
  }
  
  # Prepare the sample table
  debug(localtime()."\t - Add temporary unique keys in the sample table");
  if ($dbVar->do(qq{SHOW KEYS FROM sample WHERE Key_name='name_key';}) < 1){
    $dbVar->do(qq{ALTER TABLE sample ADD KEY `name_key` (`name`)});
  }

  debug(localtime()."\t - Add temporary columns in the $sv_table table");
  
  if ($dbVar->do(qq{show columns from $sv_table like '$tmp_sv_col';}) != 1){
    $dbVar->do(qq{ALTER TABLE $sv_table ADD COLUMN $tmp_sv_col varchar(255);});
  }
  
  if ($dbVar->do(qq{show columns from $sv_table like '$tmp_sv_clin_col';}) != 1){
    $dbVar->do(qq{ALTER TABLE $sv_table ADD COLUMN $tmp_sv_clin_col varchar(255);});
  }
  
}


#### Post processing ####

sub post_processing_annotation {
  debug(localtime()." Post processing of the table $svs_table");
  
  my $stmt = qq{ DELETE FROM $svs_table WHERE sample_id IS NULL };
  $dbVar->do($stmt);
  
  # Duplicate strains at the SV level (only for mouse)
  if ($species =~ /mouse|mus/i) {
    $stmt = qq{
      INSERT IGNORE INTO $svs_table (structural_variation_id,sample_id) 
      SELECT distinct sva.structural_variation_id,svs.sample_id FROM $svs_table svs, $sva_table sva, sample s 
      WHERE sva.supporting_structural_variation_id=svs.structural_variation_id AND s.sample_id=svs.sample_id
        AND s.display!='UNDISPLAYABLE'
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
  debug(localtime()." Post processing of the table individual");
  my $stmt;
  
  # Sample
  my $sth = $dbVar->prepare(qq{ SELECT s.sample_id, s.description, count(DISTINCT sv.study_id) c 
                                FROM $sv_table sv, $svs_table svs, sample s
                                WHERE sv.structural_variation_id=svs.structural_variation_id AND s.sample_id=svs.sample_id
                                GROUP BY svs.sample_id HAVING c>1
                              });
  $sth->execute();
  while (my @res = ($sth->fetchrow_array)) {
    if ($res[1] =~ /^(Sample|Subject) from the DGVa study/) {
      $stmt = qq{UPDATE sample SET description='Sample from several DGVa studies', study_id=null WHERE sample_id=$res[0]};
      $dbVar->do($stmt);
    }
  }
  $sth->finish;

  # Individual
  my $type = ($species =~ /mouse|mus/i) ? 'Strain' : 'Subject';
  $sth = $dbVar->prepare(qq{ SELECT i.individual_id, i.description, count(DISTINCT sv.study_id) c 
                                FROM $sv_table sv, individual i, $svs_table svs, sample s
                                WHERE sv.structural_variation_id=svs.structural_variation_id AND s.sample_id=svs.sample_id
                                  AND s.individual_id=i.individual_id
                                GROUP BY i.individual_id HAVING c>1
                              });
  $sth->execute();
  while (my @res = ($sth->fetchrow_array)) {
    if ($res[1] =~ /^(Sample|Subject) from the DGVa study/) {
      $stmt = qq{UPDATE individual SET description='$type from several DGVa studies' WHERE individual_id=$res[0]};
      $dbVar->do($stmt);
    }
  }
  $sth->finish;
}


# Create entries in phenotype_feature for the structural variation (most of the time, the phenotypes are at the SSV level)
sub post_processing_phenotype {
  debug(localtime()." Post processing of the table $pf_table");

  $dbVar->do(q{OPTIMIZE TABLE phenotype});

  my $stmt;
  
  $stmt = qq{
    INSERT IGNORE INTO 
    $pf_table (
      phenotype_id,
      object_id,
      source_id,
      study_id,
      type,
      seq_region_id,
      seq_region_start,
      seq_region_end,
      seq_region_strand
    )
    SELECT DISTINCT
      pf1.phenotype_id,
      svf.variation_name,
      $source_id,
      svf.study_id,
      'StructuralVariation',
      svf.seq_region_id,
      svf.seq_region_start,
      svf.seq_region_end,
      svf.seq_region_strand
    FROM
      $sv_table sv,
      $svf_table svf,
      $sva_table sva,
      $pf_table pf1
    WHERE 
      svf.structural_variation_id=sva.structural_variation_id AND
      svf.is_evidence=0 AND
      sva.supporting_structural_variation_id=sv.structural_variation_id AND
      sv.is_evidence=1 AND
      sv.variation_name=pf1.object_id AND
      pf1.type='SupportingStructuralVariation' AND
      NOT EXISTS (SELECT pf2.phenotype_feature_id FROM $pf_table pf2 WHERE 
                    pf2.object_id=svf.variation_name AND
                    pf2.phenotype_id=pf1.phenotype_id AND
                    pf2.source_id=$source_id AND
                    pf2.study_id=svf.study_id AND
                    pf2.seq_region_id=svf.seq_region_id AND 
                    pf2.seq_region_start=svf.seq_region_start AND
                    pf2.seq_region_end=svf.seq_region_end AND
                    pf2.seq_region_strand=svf.seq_region_strand
                 )
  };
  $dbVar->do($stmt);
}  
  

# Update clinical significances at the SV level (the clinical significance is defined at the SSV level)
sub post_processing_clinical_significance {
  debug(localtime()." Post processing of the table $sv_table for clinical significance");

  # Create temporary table
  $dbVar->do(qq{DROP TABLE IF EXISTS $temp_clin_table});
  $dbVar->do(qq{CREATE TABLE $temp_clin_table (sv_id int(10), clin_sign varchar(255), primary key (sv_id))});

  my $stmt = qq{
    INSERT INTO
    $temp_clin_table (sv_id, clin_sign)
    SELECT
      sva.structural_variation_id,
      GROUP_CONCAT(distinct sv.clinical_significance ORDER BY sv.clinical_significance DESC SEPARATOR ',')
    FROM
      $sv_table sv,
      $sva_table sva
    WHERE
      sv.structural_variation_id=sva.supporting_structural_variation_id AND
      sv.clinical_significance is not null
    GROUP BY sva.structural_variation_id;
  };
  $dbVar->do($stmt);

  my $stmt2 = qq{
    UPDATE $sv_table sv, $temp_clin_table t SET sv.clinical_significance=t.clin_sign WHERE sv.structural_variation_id=t.sv_id;
  };
  $dbVar->do($stmt2);

  $dbVar->do(qq{DROP TABLE $temp_clin_table});
}


# Unfail the variants which have at least one good mapping (e.g. Mapping to chromosome 14 'OK' and mapping to sequence NT_187600.1 'Not OK').
sub post_processing_failed_variants {
  debug(localtime()." Post processing of the table $sv_failed: unflag the variants which have at least one entry in the $svf_table table");
  
  my $stmt = qq{DELETE FROM $sv_failed WHERE structural_variation_id IN (SELECT structural_variation_id FROM $svf_table) and failed_description_id IN (17,18)};
  
  $dbVar->do($stmt);
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
  debug(localtime()." Cleanup temporary columns and keys");

  # Drop constraints in structural_variation_feature;
  $dbVar->do(qq{ALTER TABLE $svf_table DROP KEY name_coord_key});
  $dbVar->do(qq{ALTER TABLE $svf_table DROP KEY name_key});
  debug(localtime()."\t - Table $svf_table: cleaned");

  # Drop a unique constraint in structural_variation_sample;
  $dbVar->do(qq{ALTER TABLE $svs_table DROP KEY sv_id_key});
  debug(localtime()."\t - Table $svs_table: cleaned");
  
  # Drop a unique constraint in individual;
  $dbVar->do(qq{ALTER TABLE individual DROP KEY name_key});
  debug(localtime()."\t - Table individual: cleaned");

  # Drop a unique constraint in sample;
  $dbVar->do(qq{ALTER TABLE sample DROP KEY name_key});
  debug(localtime()."\t - Table sample: cleaned");

  # structural_variation table
  
  my $sv_flag = 0;
  # Column "tmp_class_name" in structural_variation
  my $sth1 = $dbVar->prepare(qq{ SELECT count(*) FROM $sv_table WHERE source_id=$source_id AND class_attrib_id=0});
  $sth1->execute();
  my $sv_count = ($sth1->fetchrow_array)[0];
  $sth1->finish;
  if ($sv_count != 0) {
    print STDERR "\tThe table $sv_table has $sv_count variants with no class_attrib_id defined!\nPlease, look at the column '$tmp_sv_col'\n";
    $sv_flag = 1;
  } 
  else {
    $dbVar->do(qq{ALTER TABLE $sv_table DROP COLUMN $tmp_sv_col});
  }
  
  # Column "tmp_clinic_name" in structural_variation
  my $sth2 = $dbVar->prepare(qq{ SELECT count(*) FROM $sv_table WHERE clinical_significance is NULL AND $tmp_sv_clin_col is not NULL});
  $sth2->execute();
  my $sv_clin_count = ($sth2->fetchrow_array)[0];
  $sth2->finish;
  if ($sv_clin_count != 0) {
    print STDERR "\tThe table $sv_table has $sv_clin_count variants with clinical_significance values which don't match the authorized values!\nPlease, look at the column '$tmp_sv_clin_col'\n";
    $sv_flag = 1;
  } 
  else {
    $dbVar->do(qq{ALTER TABLE $sv_table DROP COLUMN $tmp_sv_clin_col});
  }
  
  debug(localtime()."\t - Table $sv_table: cleaned") if ($sv_flag == 0);

}




#### Other methods ####

sub remove_data {
  my $study_id = shift;
  
  debug(localtime()." Remove the structural variation data of this study");

  # structural_variation_sample
  $dbVar->do(qq{ DELETE from $svs_table WHERE structural_variation_id IN 
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
  # variation_set_structural_variation          
  $dbVar->do(qq{ DELETE from $set_table WHERE structural_variation_id IN 
                 (SELECT structural_variation_id FROM $sv_table WHERE study_id=$study_id)
            });          
  # phenotype_feature                  
  $dbVar->do(qq{ DELETE from $pf_table WHERE study_id=$study_id });
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
 
  if(!defined($info->{bp_order})) {
    if ($info->{is_somatic} == 1 || $somatic) {
      $info->{bp_order} = 1;
    } else  {
      $info->{bp_order} = undef;
    }
  }
  $info->{phenotype} = decode_text($info->{phenotype}); 
 
  my @row = map { $info->{$_} } @attribs;

  return \@row;
}

sub generate_phenotype_data_row {
  my $info      = shift;
  my $phenotype = shift;
  
  my @row = ($info->{ID}, $phenotype);
  
  return \@row;
}


sub check_breakpoint_coordinates {
  my $info = shift;
  
  my @data = ($info);

  if (defined($info->{_bp_detail})) {
  
    my %bp_info = %$info;
  
    my ($chr,$pos,$strand) = split(':',$info->{_bp_detail});
    my ($start,$end) = split(',',$info->{_bp_range});
    
    $bp_info{chr} = $chr;
    $bp_info{start} = (defined($start) && $start!=$pos) ? $start : $pos;
    $bp_info{end} = (defined($end) && $end!=$pos) ? $end : $pos;
    foreach my $coord ('outer_start','inner_start','inner_end','outer_end') {
      $bp_info{$coord} = 0;
    }
    if (defined($strand)) {
      $bp_info{strand} = ($strand eq '-') ? -1 : 1;
    }
    push @data, \%bp_info;
  }
  
  return \@data;
}


## Get phenotype descriptions ##

sub get_medgen_phenotype {
  my $id = shift;
  
  return undef unless (-e $medgen_file && defined($id));
  
  my $res = `grep $id $medgen_file`;
  
  if ($res) {
    my $phen_desc;
    if ($res =~ /"(.+)"/) {
      $phen_desc = $1;
    }
    else {
      $phen_desc = (split(',',$res))[1];
      $phen_desc =~ s/"//g;
    }
    return $phen_desc;
  }

  return undef;
}

sub get_hpo_phenotype {
  my $id = shift;
  
  return undef unless (-e $hpo_file && defined($id));
  
  my $res = `grep -A1 -w "id: $id" $hpo_file`;
  
  if ($res) {
    $res =~ /name:\s(.+)$/;
    return $1 if ($1);
  }
  
  return undef;
}


sub usage {
  die shift, qq{

Options:
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
  -medgen_file     : Path to the unzipped MedGen file (see on ftp://ftp.ncbi.nlm.nih.gov/pub/medgen/csv/NAMES.csv.gz)
  -hpo_file        : Path to the Human Phenotype Ontology (HPO) file (see on http://compbio.charite.de/hudson/job/hpo/lastSuccessfulBuild/artifact/hp/hp.obo)
  -replace         : flag to remove the existing study data from the database before import them (optional)
  -debug           : flag to keep the $temp_table table (optional)
  };
}
