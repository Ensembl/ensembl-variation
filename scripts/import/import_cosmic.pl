#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use ImportUtils qw(create);
use DBI qw(:sql_types);
use Text::CSV;

my ( $infile, $registry_file, $version, $help );

GetOptions(
  "import|i=s"   => \$infile,
  "registry|r=s" => \$registry_file,
  "version=s"    => \$version,
  "help|h"       => \$help,
);

unless (defined($registry_file) && defined($infile) && defined($version)) {
    print "Must supply an import file, a registry file and a version ...\n" unless $help;
    $help = 1;
}
if ($help) {
    print "Usage: $0 --import <input_file> --registry <reg_file> --version <cosmic_version>\n";
    exit(0);
}

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_all($registry_file);

my $dbh = $registry->get_adaptor(
    'human', 'variation', 'variation'
)->dbc;
my $dbVar = $dbh->db_handle;


my $source_name = 'COSMIC';
my $source_id = get_source_id(); # COSMIC source_id
my $variation_set_cosmic = get_variation_set_id($source_name); # COSMIC variation set
my $variation_set_pheno = get_variation_set_id("All phenotype/disease-associated variants"); #All phenotype/disease variants
my $pheno_evidence_id = get_attrib_id('evidence','Phenotype_or_Disease');
my $pheno_class_attrib_id = get_attrib_id('phenotype_type', 'tumour');

my $temp_table      = 'MTMP_tmp_cosmic';
my $temp_phen_table = 'MTMP_tmp_cosmic_phenotype';
my $temp_varSyn_table = 'MTMP_tmp_cosmic_synonym';

my $default_class = 'sequence_alteration'; 
my %class_mapping = ( 'Substitution' => 'SNV',
                      'Indel'        => 'indel',
                      'Insertion'    => 'insertion',
                      'Deletion'     => 'deletion',
                    );

my $default_strand = 1;
my $somatic = 1;
my $allele  = 'COSMIC_MUTATION';
my $phe_suffix = 'tumour';
  
$dbVar->do("DROP TABLE IF EXISTS $temp_table;");
$dbVar->do("DROP TABLE IF EXISTS $temp_phen_table;");
$dbVar->do("DROP TABLE IF EXISTS $temp_varSyn_table;");
  
my @cols = ('name *', 'seq_region_id i*', 'seq_region_start i', 'seq_region_end i', 'class i');
create($dbVar, "$temp_table", @cols);
$dbVar->do("ALTER TABLE $temp_table ADD PRIMARY KEY (name, seq_region_id, seq_region_start, seq_region_end);");

my @cols_phen = ('name *', 'phenotype_id i*');
create($dbVar, "$temp_phen_table", @cols_phen);
$dbVar->do("ALTER TABLE $temp_phen_table ADD PRIMARY KEY (name, phenotype_id);");

my @cols_syn = ('name *', 'old_name *');
create($dbVar, "$temp_varSyn_table", @cols_syn);
$dbVar->do("ALTER TABLE $temp_varSyn_table ADD PRIMARY KEY (name, old_name);");

my $cosmic_ins_stmt = qq{
    INSERT IGNORE INTO
      $temp_table (
        name,
        seq_region_id,
        seq_region_start,
        seq_region_end,
        class
      )
      VALUES (
        ?,
        ?,
        ?,
        ?,
        ?
      )
};
my $cosmic_ins_sth = $dbh->prepare($cosmic_ins_stmt);

my $cosmic_phe_ins_stmt = qq{
    INSERT IGNORE INTO
      $temp_phen_table (
        name,
        phenotype_id
      )
      VALUES (
        ?,
        ?
      )
};
my $cosmic_phe_ins_sth = $dbh->prepare($cosmic_phe_ins_stmt);

my $cosmic_syn_ins_stmt = qq{
    INSERT IGNORE INTO
      $temp_varSyn_table (
        name,
        old_name
      )
      VALUES (
        ?,
        ?
      )
};
my $cosmic_syn_ins_sth = $dbh->prepare($cosmic_syn_ins_stmt);

my $class_attrib_ids = get_class_attrib_ids();
my $seq_region_ids   = get_seq_region_ids();
my $phenotype_ids    = get_phenotype_ids();

my %chr_names = ( '23' => 'X',
                  '24' => 'Y',
                  '25' => 'MT');


if ($infile =~ /gz$/) {
  open IN, "zcat $infile |" or die ("Could not open $infile for reading");
}
else {
  open(IN,'<',$infile) or die ("Could not open $infile for reading");
}

my $csvP = Text::CSV->new({ sep_char => ',' });

# Read through the file and parse out the desired fields
while (<IN>) {
  chomp;
  if (!$csvP->parse($_)){
    print STDERR "WARNING: could not parse line: $_\n";
    next;
  }
  my @line = $csvP->fields();

  my $chr = shift(@line);
     $chr = $chr_names{$chr} if ($chr_names{$chr});
  my $start         = shift(@line);
  my $end           = shift(@line);
  my $cosv_id       = shift(@line);
  my $cosmic_id     = shift(@line);
  shift(@line); # Skip this column
  my @phenos        = split(',', shift(@line));
  my $cosmic_class  = pop(@line);
  
  my $class = get_equivalent_class($cosmic_class,$start,$end);
  
  my $seq_region_id = $seq_region_ids->{$chr};

  if (!$seq_region_id) {
    print STDERR "COSMIC $cosmic_id: chromosome '$chr' not found in ensembl. Entry skipped.\n";
    next;
  }

  my $class_attrib_id = $class_attrib_ids->{$class};
  
  $cosmic_ins_sth->bind_param(1,$cosv_id,SQL_VARCHAR);
  $cosmic_ins_sth->bind_param(2,$seq_region_id,SQL_INTEGER);
  if ($class eq 'insertion' ){
    $cosmic_ins_sth->bind_param(3,$end+1,SQL_INTEGER);
    $cosmic_ins_sth->bind_param(4,$end,SQL_INTEGER);
  } else {
    $cosmic_ins_sth->bind_param(3,$start,SQL_INTEGER);
    $cosmic_ins_sth->bind_param(4,$end,SQL_INTEGER);
  }
  $cosmic_ins_sth->bind_param(5,$class_attrib_id,SQL_INTEGER);
  $cosmic_ins_sth->execute();

  if ( $cosmic_id =~ /COSM/ ){
    $cosmic_syn_ins_sth->bind_param(1,$cosv_id,SQL_VARCHAR);
    $cosmic_syn_ins_sth->bind_param(2,$cosmic_id,SQL_VARCHAR);
    $cosmic_syn_ins_sth->execute();
  }
  
  foreach my $phenotype (@phenos) {
    $phenotype =~ s/_/ /g;
    $phenotype = ucfirst($phenotype)." $phe_suffix";

    my $phenotype_id = $phenotype_ids->{$phenotype};
    
    if (!$phenotype_id) {
      $phenotype_id = add_phenotype($phenotype);
      print STDERR "COSMIC $cosmic_id: phenotype '$phenotype' not found in ensembl. Phenotype added.\n";
    }
    
    $cosmic_phe_ins_sth->bind_param(1,$cosv_id,SQL_VARCHAR);
    $cosmic_phe_ins_sth->bind_param(2,$phenotype_id,SQL_INTEGER);
    $cosmic_phe_ins_sth->execute();
  }
}
close(IN);
$cosmic_ins_sth->finish();
$cosmic_phe_ins_sth->finish();
$cosmic_syn_ins_sth->finish();

# Insert COSMIC in the latest release which are not in COSMIC 71
insert_cosmic_entries();


sub get_equivalent_class {
  my $type  = shift;
  my $start = shift;
  my $end   = shift;

  my @type_parts = split(',|\s',$type);
  $type = $type_parts[0];

  my $class = $default_class;
  # map the COSMIC class type into the predefined set of class types in %class_mapping
  if ($type eq 'Substitution') {
    $class = ($start == $end) ? $class_mapping{$type} : $class_mapping{'Indel'};
  }
  elsif ($class_mapping{$type}) {
    $class = $class_mapping{$type};
  }
  return $class;
}

sub get_class_attrib_ids {
  my %class_attrib_ids;
  my $get_class_attrib_ids_sth = $dbh->prepare(
  qq{
    SELECT a.value, a.attrib_id
    FROM   attrib a, attrib_type at
    WHERE  a.attrib_type_id = at.attrib_type_id
    AND    at.code = 'SO_term'
  });
  $get_class_attrib_ids_sth->execute;
  while (my ($value, $attrib_id) = $get_class_attrib_ids_sth->fetchrow_array) {
    $class_attrib_ids{$value} = $attrib_id;
  }
  $get_class_attrib_ids_sth->finish();
  
  return \%class_attrib_ids;
}

sub get_seq_region_ids {
  
  my $sth = $dbh->prepare(
  qq{
    SELECT seq_region_id, name
    FROM seq_region
  });
  $sth->execute;
  
  my (%seq_region_ids, $id, $name);
  $sth->bind_columns(\$id, \$name);
  $seq_region_ids{$name} = $id while $sth->fetch();
  $sth->finish;
  
  return \%seq_region_ids;
}

sub get_phenotype_ids {
  
  my $sth = $dbh->prepare(
  qq{
    SELECT phenotype_id, description
    FROM phenotype
  });
  $sth->execute;
  
  my (%phenotype_ids, $id, $desc);
  $sth->bind_columns(\$id, \$desc);
  $phenotype_ids{$desc} = $id while $sth->fetch();
  $sth->finish;
  
  return \%phenotype_ids;
}


sub get_source_id {
  
  # Check if the COSMIC source already exists, else it create the entry
  if ($dbVar->selectrow_arrayref(qq{SELECT source_id FROM source WHERE name="$source_name";})) {
    $dbVar->do(qq{UPDATE IGNORE source SET version=$version where name="$source_name";});
  }
  else {
    $dbVar->do(qq{INSERT INTO source (name,description,url,version,somatic_status,data_types) VALUES ("$source_name",'Somatic mutations found in human cancers from the COSMIC project - Public version','https://cancer.sanger.ac.uk/cosmic/',$version,'somatic','variation,variation_synonym,phenotype_feature');});
  }
  my @source_id = @{$dbVar->selectrow_arrayref(qq{SELECT source_id FROM source WHERE name="$source_name";})};
  return $source_id[0];
}

# get attributes based on attrib and attrib_type tables
sub get_attrib_id {
  my ($type, $value) = @_;

  #GET attrib of specific type
  my $aid = $dbVar->selectrow_arrayref(qq{
    SELECT a.attrib_id
    FROM attrib a JOIN attrib_type att
      ON a.attrib_type_id = att.attrib_type_id
    WHERE att.code="$type" and a.value = "$value";});

  if (!$aid->[0]){
    die("Couldn't find the $value attrib of $type type\n");
  } else {
    return $aid->[0];
  }
}

sub add_phenotype {
  my $phenotype = shift;
  $dbVar->do(qq{INSERT IGNORE INTO phenotype (description, class_attrib_id)
                VALUES ("$phenotype", $pheno_class_attrib_id )});
  my $phenotype_id = $dbVar->selectrow_arrayref(qq{SELECT phenotype_id FROM phenotype WHERE description="$phenotype"});

  # Update the list of phenotypes
  $phenotype_ids->{$phenotype} = $phenotype_id->[0];

  return $phenotype_id->[0];
}

sub get_variation_set_id {
  my $type = shift;

  my $variation_set_ids = $dbVar->selectrow_arrayref(qq{SELECT variation_set_id FROM variation_set WHERE name LIKE '$type%'});

  if (!$variation_set_ids) {
    die("Couldn't find the '$type' variation set");
  }
  else {
    return $variation_set_ids->[0];
  }
}


sub insert_cosmic_entries {
  
  # Insert Var
  my $stmt_var = qq{INSERT IGNORE INTO variation
                    (name, source_id, class_attrib_id, somatic,
                     evidence_attribs, display)
                    SELECT name, ?, class, ?, '$pheno_evidence_id', 1
                    FROM $temp_table};
  my $sth_var  = $dbh->prepare($stmt_var);
  $sth_var->execute($source_id, $somatic);

  # Insert VF
  my $stmt_vf = qq{INSERT IGNORE INTO variation_feature 
                   (variation_id, variation_name, source_id, variation_set_id, class_attrib_id, somatic, allele_string, seq_region_id, seq_region_start, seq_region_end, seq_region_strand,
                   evidence_attribs, display)
                   SELECT v.variation_id, v.name, v.source_id,'$variation_set_pheno,$variation_set_cosmic', v.class_attrib_id, v.somatic, ?, c.seq_region_id, c.seq_region_start, c.seq_region_end, ?,
                   '$pheno_evidence_id', 1
                   FROM variation v, $temp_table c WHERE v.name=c.name};
  my $sth_vf  = $dbh->prepare($stmt_vf);
  $sth_vf->execute($allele, $default_strand);
  
  # Insert variation synonym
  my $stmt_vs = qq{INSERT IGNORE INTO variation_synonym
                   (variation_id, source_id, name)
                   SELECT v.variation_id, v.source_id, c.old_name
                   FROM variation v, $temp_varSyn_table c WHERE v.name=c.name};
  my $sth_vs  = $dbh->prepare($stmt_vs);
  $sth_vs->execute();

  # Insert PF
  my $stmt_pf = qq{INSERT IGNORE INTO phenotype_feature 
                   (object_id, type, source_id, phenotype_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand)
                   SELECT v.name, "Variation", v.source_id, pc.phenotype_id, c.seq_region_id, c.seq_region_start, c.seq_region_end, ? 
                   FROM variation v, $temp_table c, $temp_phen_table pc WHERE v.name=c.name AND c.name=pc.name};
  my $sth_pf  = $dbh->prepare($stmt_pf);
  $sth_pf->execute($default_strand);
  
  # Insert Set
  my $stmt_set = qq{INSERT IGNORE INTO variation_set_variation (variation_id, variation_set_id)
                    SELECT variation_id, ? FROM variation WHERE source_id=?};
  my $sth_set  = $dbh->prepare($stmt_set);
  $sth_set->execute($variation_set_cosmic, $source_id);
}

