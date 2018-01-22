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
use Getopt::Long;
use DBI qw(:sql_types);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use ImportUtils;

my $infile;
my $host;
my $dbname;
my $user;
my $pass;
my $port;
our $verbose;
my $skip_synonyms;
my $skip_phenotypes;
my $skip_sets;
my $source;
my $source_version;
my $set;
my $help;

my $UNIPROT_SOURCE_NAME = "Uniprot";
my $UNIPROT_SOURCE_DESCRIPTION = "Variants with protein annotation imported from Uniprot";
my $UNIPROT_SOURCE_URL = "http://www.uniprot.org/";
my $UNIPROT_SET_NAME = "ph_uniprot";

my $OMIM_SOURCE_NAME = "OMIM";
my $OMIM_SOURCE_DESCRIPTION = "Variations linked to entries in the Online Mendelian Inheritance in Man (OMIM) database";
my $OMIM_SOURCE_URL = "http://www.omim.org/";
my $OMIM_SET_NAME = "ph_omim";

my $NHGRI_SOURCE_NAME = "NHGRI_GWAS_catalog";
my $NHGRI_SOURCE_DESCRIPTION = "Variants associated with phenotype data from the NHGRI GWAS catalog";
my $NHGRI_SOURCE_URL = "http://www.genome.gov/gwastudies/";
my $NHGRI_SET_NAME = "ph_nhgri";

my $EGA_SOURCE_NAME = "EGA";
my $EGA_SOURCE_DESCRIPTION = "Variants imported from the European Genome-phenome Archive with phenotype association";
my $EGA_SOURCE_URL = "http://www.ebi.ac.uk/ega/";

usage() if (!scalar(@ARGV));
 
GetOptions(
    'infile=s' => \$infile,
    'host=s' => \$host,
    'dbname=s' => \$dbname,
    'user=s' => \$user,
    'pass=s' => \$pass,
    'port=s' => \$port,
    'source=s' => \$source,
    'version=i' => \$source_version,
    'verbose!' => \$verbose,
    'skip_synonyms!' => \$skip_synonyms,
    'skip_phenotypes!' => \$skip_phenotypes,
    'skip_sets!' => \$skip_sets,
    'help!' => \$help
);

usage() if ($help);

die ("An input file is required") unless (defined($infile));
die ("Database credentials are required") unless (defined($host) && defined($dbname) && defined($user) && defined($pass));
die ("A source version is required") unless (defined($source_version));

$port ||= 3306;

my $result;
my $source_name;
my $source_description;
my $source_url;

=head

    The parser subroutines parse the input file into a common data structure.
    Currently what is returned should be a reference to a hash. The hash should
    contain the key 'phenotypes' with the value being a reference to an array of
    phenotype data objects. The phenotype data object is a reference to a hash
    where the keys correspond to the column names in the variation_annotation table.
    In addition, there are the keys 'rsid' which holds the rs-id that the phenotype
    annotates and 'description' and 'name' which correspond to the columns in the
    phenotype table.
    
    In addition, the hash can also contain the key 'synonyms' and the value is a
    reference to a hash where the keys are rs-ids and each respective value is a
    reference to an array of synonyms for the rs-id.

=cut

# Make sure that the input file is XML compliant
ImportUtils::make_xml_compliant($infile);

# Remove carriage return in the input file
remove_carriage_return($infile);


# Connect to the variation database
print STDOUT localtime() . "\tConnecting to database $dbname\n" if ($verbose);
my $db_adaptor = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host => $host,
  -user => $user,
  -pass => $pass,
  -port => $port,
  -dbname => $dbname
) or die("Could not get a database adaptor for $dbname on $host:$port");
print STDOUT localtime() . "\tConnected to $dbname on $host:$port\n" if ($verbose);

my @rsids;

# Parse the input files into a hash
if ($source =~ m/uniprot/i) {
  $result = parse_uniprot($infile);
	$source_name = $UNIPROT_SOURCE_NAME;
	$source_description = $UNIPROT_SOURCE_DESCRIPTION;
	$source_url = $UNIPROT_SOURCE_URL;
	$set = $UNIPROT_SET_NAME;
}
elsif ($source =~ m/nhgri/i) {
	$result = parse_nhgri($infile);
	$source_name = $NHGRI_SOURCE_NAME;
	$source_description = $NHGRI_SOURCE_DESCRIPTION;
	$source_url = $NHGRI_SOURCE_URL;
	$set = $NHGRI_SET_NAME;
}
elsif ($source =~ m/omim/i) {
	$result = parse_omim($infile);
	$source_name = $OMIM_SOURCE_NAME;
	$source_description = $OMIM_SOURCE_DESCRIPTION;
	$source_url = $OMIM_SOURCE_URL;
	$set = $OMIM_SET_NAME;
}
elsif ($source =~ m/dbsnp/i) {
	$result = parse_dbsnp_omim($infile);
	$source_name = $OMIM_SOURCE_NAME;
	$source_description = $OMIM_SOURCE_DESCRIPTION;
	$source_url = $OMIM_SOURCE_URL;
	$set = $OMIM_SET_NAME;
}
elsif ($source =~ m/ega/i) {
	$source_name = $EGA_SOURCE_NAME;
	$source_description = $EGA_SOURCE_DESCRIPTION;
	$source_url = $EGA_SOURCE_URL;
	my $source_id = get_or_add_source($source_name,$source_description,$source_url,$db_adaptor);
	print STDOUT "$source source_id is $source_id\n" if ($verbose);
	parse_ega($infile,$source_id);
	exit(0);
}
else {
	die("Source $source is not recognized");
}

my %synonym;
my @phenotypes;
if (exists($result->{'synonyms'})) {
	%synonym = %{$result->{'synonyms'}};
  # To get all the rsids of the source (Uniprot)
  @rsids = keys(%synonym);
}
if (exists($result->{'phenotypes'})) {
  @phenotypes = @{$result->{'phenotypes'}};
}

# Get internal variation ids for the rsIds
if (scalar @rsids == 0) {
  @rsids = map {$_->{'rsid'}} @phenotypes;
}
my $variation_ids = get_dbIDs(\@rsids,$db_adaptor);

# Get or add a source
my $source_id = get_or_add_source($source_name,$source_description,$source_url,$db_adaptor);
print STDOUT "$source source_id is $source_id\n" if ($verbose);

# Add the synonyms if required
add_synonyms(\%synonym,$variation_ids,$source_id,$db_adaptor) unless ($skip_synonyms);

# Now, insert phenotypes
add_phenotypes(\@phenotypes,$variation_ids,$source_id,$db_adaptor) unless ($skip_phenotypes);

# Add the variation sets if required
add_set($set,$source_id,$db_adaptor) unless ($skip_sets);


# Loop over the remaining rsids (the ones that could not be find in the db) and print them out
while (my ($rs_id,$var_id) = each(%{$variation_ids})) {
    next if (defined($var_id->[0]));
    print STDOUT "$rs_id could not be found in $dbname";
    if (defined($synonym{$rs_id})) {
        print STDOUT " (Synonyms: " . join(", ",@{$synonym{$rs_id}}) . ")";
    }
    print STDOUT "\n";
}



###########
# METHODS #
###########

sub parse_uniprot {
  my $infile = shift;
  
  my %synonym;
  my @phenotypes;

  # Open the input file for reading
  open(IN,'<',$infile) or die ("Could not open $infile for reading");
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;
    
    # A regexp to catch the meta information in the header. Just echo this to the stdout for logging purposes
    if ($_ =~ m/^\s*(Description|Name|Release)\:/) {
      print STDOUT $_ . "\n";
    }
    
    # Main regexp to extract relevant variation information
    if ($_ =~ m/^(\S+)\s+\S+\s+(VAR\_\d+)\s+\w\.\S+\s+(Disease|Polymorphism|Unclassified)\s+(\-|rs\d*)\s*(.*)$/) {
      
      # Get the data that was caught by the regexp
      my $gene = $1;
      my $uniprot_id = $2;
      my $rs_id = $4;
      my $phenotype = $5;
      
      # If no rsId was given, will attempt to get one by looking up the Uniprot id in the synonym table
      if ($rs_id ne '-') {
        push(@{$synonym{$rs_id}},$uniprot_id);
      }
      else {
        $rs_id = $uniprot_id;
      }
      
      $phenotype ||= '-';
      
      # Try to further split the phenotype into a short name, description and possibly MIM id
      if ($phenotype ne '-') {
        my $description;
        my $name;
        my $mim_id;
        
        ($description,$mim_id) = $phenotype =~ m/^([^\[]+)(?:\[(MIM\:.+?)\])?$/;
        ($description,$name) = $description =~ m/^(.+?)\s*(?:\((.+?)\))?\s*$/;
        
        $mim_id &&= join(",MIM:",split(",",$mim_id));
        $mim_id =~ s/\s+//g if (defined($mim_id));
        
        push(
          @phenotypes,
          {
            "rsid"      => $rs_id,
            "associated_gene" => $gene,
            "description"   => $description,
            "name"      => $name,
            "variation_names" => $rs_id,
            "study"       => $mim_id,
            "study_type"    => 'GWAS'
          }
        );
      }
    }
  }
  close(IN);
  
  print STDOUT "Parsed " . scalar(keys(%synonym)) . " rs-ids with Uniprot synonyms\n" if ($verbose);
  print STDOUT scalar(@phenotypes) . " phenotypes were found linked to rs-ids\n" if ($verbose);
  
  my %result = ('synonyms' => \%synonym, 'phenotypes' => \@phenotypes);
  return \%result;
}

sub parse_nhgri {
  my $infile = shift;
  
  my @phenotypes;
  
  # Open the input file for reading
  open(IN,'<',$infile) or die ("Could not open $infile for reading");
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;
    
  my @content = split(/\t/,$_);
    
    my $pubmed_id    = $content[1];
    my $study      = $content[6];
    my $phenotype    = $content[7];
    my $gene       = $content[13];
    my $rs_risk_allele = $content[20];
    my $rs_id      = $content[21];
    my $risk_frequency = $content[26];
    my $pvalue     = $content[27];
    
    my %data = (
      'study_type' => 'GWAS',
      'description' => $phenotype,
      'associated_gene' => $gene,
      'associated_variant_risk_allele' => $rs_risk_allele,
      'risk_allele_freq_in_controls' => $risk_frequency,
      'p_value' => $pvalue,
      'study_description' => $study,
    );
    
    # Parse the rsids
    my @rsids;
    $rs_id ||= "";
    while ($rs_id =~ m/(rs[0-9]+)/g) {
      push(@rsids,$1);
    }
    $data{'variation_names'} = join(',',@rsids);
    $data{'study'} = 'pubmed/' . $pubmed_id if (defined($pubmed_id));
    
    # If we didn't get any rsIds, skip this row (this will also get rid of the header)
    warn("Could not parse any rsIds from string '$rs_id'") if (!scalar(@rsids));
    next if (!scalar(@rsids));
    
    map {
      my %t_data = %{\%data};
      $t_data{'rsid'} = $_;
      push(@phenotypes,\%t_data)
    } @rsids;
  }
  close(IN);

  my %result = ('phenotypes' => \@phenotypes);
  
  return \%result;
}

sub parse_omim {
  my $infile = shift;
  
  my @phenotypes;
  
  # Open the input file for reading
  open(IN,'<',$infile) or die ("Could not open $infile for reading");
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;
    
    my @attributes = split(/\t/);
    
    next if (!$attributes[4] or $attributes[0] =~ /^\D/);
    
    my ($study,$allele) = split(/\./,$attributes[0]);
    my $rsids = $attributes[4];
    
    my @rs_list = split(',',$rsids);
    
    #  Get one line for each variation_names
    foreach my $rs (@rs_list) {
    
       # Skip the risk allele if the variant is "0000"
      my $data = {
        'rsid'      => $rs,
        'study'       => 'MIM:'.$study,
        'associated_gene' => $attributes[2],
        'variation_names' => $rsids,
        'description'   => $attributes[1],
        'associated_variant_risk_allele' => ($allele !~ m/^\s*0+\s*$/ ? $allele : undef),
      };
    
      push(@phenotypes,$data);
    }
  }
  
  close(IN);
  
  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}


sub parse_dbsnp_omim {
  my $infile = shift;
  
  my @phenotypes;
  my @attribute_keys = (
    'ID',
    'Phenotype_study',
    'Phenotype_associated_variant_risk_seq',
    'Omim_title',
    'Allele_title',
    'Gene_names'
  );
  
  # Open the input file for reading
  open(IN,'<',$infile) or die ("Could not open $infile for reading");
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;
    
    my @attributes = split(/\t/);
    
    # Skip the risk allele if the variant is "0000"
    my $data = {
      'rsid'      => 'rs' . $attributes[0],
      'study'       => 'MIM:' . $attributes[1],
      'associated_gene' => $attributes[5],
      'variation_names' => 'rs' . $attributes[0],
      'associated_variant_risk_allele' => ($attributes[2] !~ m/^\s*0+\s*$/ ? $attributes[2] : undef),

    };
    
    # If available, use the variant title, else use the omim record title
    if (defined($data->{'associated_variant_risk_allele'}) and $attributes[4] ne '') {
      $data->{'description'} = $attributes[4];
    } else {
      $data->{'description'} = $attributes[3];
    }
    
    # If possible, try to extract the last comma-separated word as this should be the short name for the phenotype
    @attributes = split(/;/,$data->{'description'});
    if (scalar(@attributes) > 1) {
      ($data->{'name'}) = pop(@attributes) =~ m/(\S+)/;
      $data->{'description'} = join(';',@attributes);
    }
    
    push(@phenotypes,$data);
  }
  
  close(IN);
  
  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}


sub parse_ega {
  my $infile = shift;
  my $source_id = shift;
  
  my $study_check_stmt = qq{
    SELECT
      study_id
    FROM
      study
    WHERE
      name=? AND source_id=$source_id
    LIMIT 1
  };
  my $nhgri_check_stmt = qq{
    SELECT
      st.study_id,st.study_type
    FROM
      study st, source s
    WHERE
      external_reference=? AND s.name like '%nhgri%'
      AND s.source_id=st.source_id
    LIMIT 1
  };
  my $study_ins_stmt = qq{
    INSERT INTO
      study (
      name,
      source_id,
      external_reference,
      url,
      study_type
      )
    VALUES (
      ?,
      $source_id,
      ?,
      ?,
      ?
    )
  };
  # NHGRI and EGA associated studies
  my $asso_study_check_stmt = qq{
    SELECT
      study1_id
    FROM
      associate_study
    WHERE
      study1_id = ? AND study2_id = ?
    LIMIT 1
  };
  my $asso_study_ins_stmt = qq{
    INSERT INTO
      associate_study (study1_id,study2_id)
    VALUES (?,?)
  };
  
  my $nhgri_check_sth = $db_adaptor->dbc->prepare($nhgri_check_stmt);
  my $study_check_sth = $db_adaptor->dbc->prepare($study_check_stmt);
  my $study_ins_sth   = $db_adaptor->dbc->prepare($study_ins_stmt);
  my $asso_study_check_sth = $db_adaptor->dbc->prepare($asso_study_check_stmt);
  my $asso_study_ins_sth   = $db_adaptor->dbc->prepare($asso_study_ins_stmt);
  
  # Open the input file for reading
  open(IN,'<',$infile) or die ("Could not open $infile for reading");
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp $_;
    my @attributes = split(",",$_);
    next if ($attributes[1] eq '');
    my $name = $attributes[0];
    my $pubmed = 'pubmed/'.$attributes[1];
    my $url = $attributes[2];
    
    # NHGRI study
    my $nhgri_study_id;
    my $study_type;
    $nhgri_check_sth->bind_param(1,$pubmed,SQL_VARCHAR);
    $nhgri_check_sth->execute();
    $nhgri_check_sth->bind_columns(\$nhgri_study_id,\$study_type);
    $nhgri_check_sth->fetch();
    
    if (!defined($nhgri_study_id)) {
      print "No NHGRI study found for the EGA $name | $pubmed !\n";
      next;
    }
    
    # EGA study
    my $study_id;
    $study_check_sth->bind_param(1,$name,SQL_VARCHAR);
    $study_check_sth->execute();
    $study_check_sth->bind_columns(\$study_id);
    $study_check_sth->fetch();
    if (!defined($study_id)) {
      $study_ins_sth->bind_param(1,$name,SQL_VARCHAR);
      $study_ins_sth->bind_param(2,$pubmed,SQL_VARCHAR);
      $study_ins_sth->bind_param(3,$url,SQL_VARCHAR);
      $study_ins_sth->bind_param(4,$study_type,SQL_VARCHAR);
      $study_ins_sth->execute();
      $study_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
    }
    
    my $is_associated;
    $asso_study_check_sth->bind_param(1,$nhgri_study_id,SQL_INTEGER);
    $asso_study_check_sth->bind_param(2,$study_id,SQL_INTEGER);
    $asso_study_check_sth->execute();
    $asso_study_check_sth->bind_columns(\$is_associated);
    $asso_study_check_sth->fetch();
    
    if (!defined($is_associated)) {
      $asso_study_ins_sth->bind_param(1,$nhgri_study_id,SQL_INTEGER);
      $asso_study_ins_sth->bind_param(2,$study_id,SQL_INTEGER);
      $asso_study_ins_sth->execute();
    }
  }
  close(IN);
}

sub get_dbIDs {
  my $rs_ids = shift;
  my $db_adaptor = shift;
  
  my $id_stmt = qq{
    SELECT
      DISTINCT
      v.variation_id,
      v.name
    FROM
      variation v,
      variation_feature vf
      
    WHERE
      v.name = ? AND
      v.variation_id=vf.variation_id
    LIMIT 1
  };
  my $syn_stmt = qq{
    SELECT
      DISTINCT
      v.variation_id,
      v.name
    FROM
      variation_feature vf,
      variation_synonym vs JOIN
      variation v ON vs.variation_id = v.variation_id
    WHERE
      vs.name = ? AND
      v.variation_id=vf.variation_id
    LIMIT 1
  };
  my $id_sth = $db_adaptor->dbc->prepare($id_stmt);
  my $syn_sth = $db_adaptor->dbc->prepare($syn_stmt);
  
  my %mapping;
  
  foreach my $rs_id (@{$rs_ids}) {
     $id_sth->bind_param(1,$rs_id,SQL_VARCHAR);
    $id_sth->execute();
    my ($var_id,$var_name);
    $id_sth->bind_columns(\$var_id,\$var_name);
    $id_sth->fetch();
    
    # If we couldn't find the rs_id, look in the synonym table
    if (!defined($var_id)) {
      $syn_sth->bind_param(1,$rs_id,SQL_VARCHAR);
      $syn_sth->execute();
      $syn_sth->bind_columns(\$var_id,\$var_name);
      $syn_sth->fetch();
    }
    $mapping{$rs_id} = [$var_id,$var_name];
  }
  
  return \%mapping;
}

sub get_or_add_source {
  my $source_name = shift;
  my $source_description = shift;
  my $source_url = shift;
  my $db_adaptor = shift;
  
  my $stmt = qq{
    SELECT
      source_id
    FROM
      source
    WHERE
      name = '$source_name'
    LIMIT 1
  };
  my $sth = $db_adaptor->dbc->prepare($stmt);
  $sth->execute();
  my $source_id;
  $sth->bind_columns(\$source_id);
  $sth->fetch();
  
  if (!defined($source_id)) {
    $stmt = qq{
      INSERT INTO
        source (
          name,
          description,
          url
        )
      VALUES (
        '$source_name',
        '$source_description',
        '$source_url'
      )
    };
    $db_adaptor->dbc->do($stmt);
    $source_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
    
    print STDOUT "Added source for $source_name (source_id = $source_id)\n" if ($verbose);
  }
  else {
    $stmt = qq{
      UPDATE
        source 
      SET name=?,
        description=?,
        url=?,
        version=?
      WHERE
        source_id=?
    };
    my $update_source_sth =$db_adaptor->dbc->prepare($stmt);
    $update_source_sth->execute($source_name,$source_description,$source_url,$source_version,$source_id);
  }

  return $source_id;
}

sub add_phenotypes {
  my $phenotypes = shift;
  my $variation_ids = shift;
  my $source_id = shift;
  my $db_adaptor = shift;
  
  # Prepared statements
  my $phen_check_stmt = qq{
    SELECT
      phenotype_id
    FROM
      phenotype
    WHERE
      description = ?
    LIMIT 1
  };
  my $phen_ins_stmt = qq{
    INSERT INTO phenotype ( name, description ) VALUES ( ?,? )
  };
  my $st_ins_stmt = qq{
    INSERT INTO study (
      source_id,
      external_reference,
      study_type,
      description
    )
    VALUES (
      $source_id,
      ?,
      ?,
      ?
    )
  };
  
  my $extra_cond = '';
  if ($source =~ m/uniprot/i) {
    $extra_cond = 'AND variation_names = ? ';
  }
  elsif ($source =~ m/omim/i) {
    $extra_cond = 'AND associated_variant_risk_allele = ? AND associated_gene = ?';
  }
  my $va_check_stmt = qq{
    SELECT
      variation_annotation_id
    FROM
      variation_annotation
    WHERE
      variation_id = ? AND
      phenotype_id = ? AND
      study_id = ? 
      $extra_cond
    LIMIT 1
  };
  
  my $va_ins_stmt = qq{
    INSERT INTO variation_annotation (
      variation_id,
      phenotype_id,
      study_id,
      associated_gene,
      associated_variant_risk_allele,
      variation_names,
      risk_allele_freq_in_controls,
      p_value
    )
    VALUES (
      ?,
      ?,
      ?,
      ?,
      ?,
      ?,
      ?,
      ?
    )
  };
  my $phen_check_sth = $db_adaptor->dbc->prepare($phen_check_stmt);
  my $phen_ins_sth   = $db_adaptor->dbc->prepare($phen_ins_stmt);
  my $st_ins_sth   = $db_adaptor->dbc->prepare($st_ins_stmt);
  my $va_check_sth   = $db_adaptor->dbc->prepare($va_check_stmt);
  my $va_ins_sth   = $db_adaptor->dbc->prepare($va_ins_stmt);

  # First, sort the array according to the phenotype description
  my @sorted = sort {$a->{"description"} cmp $b->{"description"}} @{$phenotypes};
  my $current = "";
  my $phenotype_id;
  my $study_count = 0;
  my $phenotype_count = 0;
  my $annotation_count = 0;
  
  while (my $phenotype = shift(@sorted)) {
  
    # If the rs could not be mapped to a variation id, skip it
    next if (!defined($variation_ids->{$phenotype->{"rsid"}}[0]));
    
    my $sql_study = '= ?';
    my $sql_type  = '= ?';
    
    # To avoid duplication of study entries
    if (!defined $phenotype->{"study"}) {$sql_study = 'IS NULL'; }
    if (!defined $phenotype->{"study_type"}) {$sql_type = 'IS NULL'; }  
    
    my $st_check_stmt = qq{
      SELECT study_id
      FROM study
      WHERE
      source_id = $source_id AND
      external_reference $sql_study AND
      study_type $sql_type
      LIMIT 1
    };
    my $st_check_sth = $db_adaptor->dbc->prepare($st_check_stmt);
    my $second_param_num = 1;
    
    my $study_id;
    if (defined $phenotype->{"study"}) {
      $st_check_sth->bind_param(1,$phenotype->{"study"},SQL_VARCHAR);
      $second_param_num = 2;
    }
    
    $st_check_sth->bind_param($second_param_num,$phenotype->{"study_type"},SQL_VARCHAR) if (defined $phenotype->{"study_type"});
    $st_check_sth->execute();
    $st_check_sth->bind_columns(\$study_id);
    $st_check_sth->fetch();
    
    if (!defined($study_id)) {
      $st_ins_sth->bind_param(1,$phenotype->{"study"},SQL_VARCHAR);
      $st_ins_sth->bind_param(2,$phenotype->{"study_type"},SQL_VARCHAR);
      $st_ins_sth->bind_param(3,$phenotype->{"study_description"},SQL_VARCHAR);
      $st_ins_sth->execute();
      $study_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
      $study_count++;
    }
    
    # If we have a new phenotype we need to see if it exists in the database and otherwise add it
    if (defined($phenotype->{"description"}) and $phenotype->{"description"} ne $current) {
      undef($phenotype_id);
      $phen_check_sth->bind_param(1,$phenotype->{"description"},SQL_VARCHAR);
      $phen_check_sth->execute();
      $phen_check_sth->bind_columns(\$phenotype_id);
      $phen_check_sth->fetch();
      
      # If no phenotype was found, we need to add it
      if (!defined($phenotype_id)) {
        $phen_ins_sth->bind_param(1,undef,SQL_VARCHAR);
        $phen_ins_sth->bind_param(2,$phenotype->{"description"},SQL_VARCHAR);
        $phen_ins_sth->execute();
        $phenotype_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
        $phenotype_count++;
      }
      $current = $phenotype->{"description"};
    }
    
    # Check if this phenotype already exists for this variation and source, in that case we probably want to skip it
    my $va_id;
    $va_check_sth->bind_param(1,$variation_ids->{$phenotype->{"rsid"}}[0],SQL_INTEGER);
    $va_check_sth->bind_param(2,$phenotype_id,SQL_INTEGER);
    $va_check_sth->bind_param(3,$study_id,SQL_INTEGER);
    # For uniprot data
    if ($source =~ m/uniprot/i) {
      $va_check_sth->bind_param(4,$phenotype->{"variation_names"},SQL_VARCHAR);
    }
    # For omim data
    if ($source =~ m/omim/i) {
      $va_check_sth->bind_param(4,$phenotype->{"associated_variant_risk_allele"},SQL_VARCHAR);
      $va_check_sth->bind_param(5,$phenotype->{"associated_gene"},SQL_VARCHAR);
    }
      
    $va_check_sth->execute();
    $va_check_sth->bind_columns(\$va_id);
    $va_check_sth->fetch();
    next if (defined($va_id));
    
    $phenotype->{"p_value"} = convert_p_value($phenotype->{"p_value"}) if (defined($phenotype->{"p_value"}));
    
    # Else, insert this variation annotation.
    $va_ins_sth->bind_param(1,$variation_ids->{$phenotype->{"rsid"}}[0],SQL_INTEGER);
    $va_ins_sth->bind_param(2,$phenotype_id,SQL_INTEGER);
    $va_ins_sth->bind_param(3,$study_id,SQL_INTEGER);
    $va_ins_sth->bind_param(4,$phenotype->{"associated_gene"},SQL_VARCHAR);
    $va_ins_sth->bind_param(5,$phenotype->{"associated_variant_risk_allele"},SQL_VARCHAR);
    $va_ins_sth->bind_param(6,$phenotype->{"variation_names"},SQL_VARCHAR);
    $va_ins_sth->bind_param(7,$phenotype->{"risk_allele_frequency_in_controls"},SQL_VARCHAR);
    $va_ins_sth->bind_param(8,$phenotype->{"p_value"},SQL_VARCHAR);
    $va_ins_sth->execute();
    $annotation_count++;
  }
  print STDOUT "$study_count new studies added\n" if ($verbose);
  print STDOUT "$phenotype_count new phenotypes added\n" if ($verbose);
  print STDOUT "$annotation_count variations were annotated with phenotypes\n" if ($verbose);
}

sub add_synonyms {
  my $synonyms = shift;
  my $variation_ids = shift;
  my $source_id = shift;
  my $db_adaptor = shift;
  
  # If we actually didn't get any synonyms, just return
  return if (!defined($synonyms) || !scalar(keys(%{$synonyms})));
  
  # Some prepeared statements needed for inserting the synonyms into database
  my $ins_stmt = qq{
    INSERT IGNORE INTO
      variation_synonym (
      variation_id,
      source_id,
      name
      )
    VALUES (
      ?,
      $source_id,
      ?
    )
  };
  my $ins_sth = $db_adaptor->dbc->prepare($ins_stmt);
  
  my $alt_count = 0;
  my $variation_count = 0;
  
  foreach my $rs_id (keys %{$variation_ids}) {
    
    my $var_id = $variation_ids->{$rs_id}[0];
    
    # If we have a variation id, we can proceed
    if (defined($var_id)) {
      
      $variation_count++;
      
      $ins_sth->bind_param(1,$var_id,SQL_INTEGER);
      
      # Handle all synonym ids for this rs_id
      while (my $alt_id = shift(@{$synonyms->{$rs_id}})) {
      
        # Add the id as synonym, if it is already present, it will just be ignored
        $ins_sth->bind_param(2,$alt_id,SQL_VARCHAR);
        $ins_sth->execute();
        $alt_count++;
      }
    }
  }
  
  print STDOUT "Added $alt_count synonyms for $variation_count rs-ids\n" if ($verbose);
}

sub add_set {
  my $set = shift;
  my $source_id = shift;
  my $db_adaptor = shift;
  
  return if (!defined($set));
  
  my $variation_set_id;
  
  # Get variation_set_id
  my $select_set_stmt = qq{
        SELECT v.variation_set_id
        FROM variation_set v, attrib a
        WHERE v.short_name_attrib_id=a.attrib_id 
          AND a.value = ?
  };
  my $sth1 = $db_adaptor->dbc->prepare($select_set_stmt);
  $sth1->bind_param(1,$set,SQL_VARCHAR);
  $sth1->execute();
  $sth1->bind_columns(\$variation_set_id);
  $sth1->fetch();
  return if (!defined($variation_set_id));
  
  # Insert into variation_set_variation
  my $insert_set_stmt = qq{ 
    INSERT IGNORE INTO variation_set_variation (variation_id,variation_set_id)
      SELECT distinct variation_id, ? 
      FROM variation_annotation WHERE study_id IN
        (SELECT study_id FROM study WHERE source_id=?)
  };
  my $sth2 = $db_adaptor->dbc->prepare($insert_set_stmt);
  $sth2->bind_param(1,$variation_set_id,SQL_INTEGER);
  $sth2->bind_param(2,$source_id,SQL_INTEGER);
  $sth2->execute();
}


sub convert_p_value {

  my $pval = shift;
  
  my $sci_pval = '';
  # If a scientific format is not found, then ...
  if ($pval !~ /^\d+.*e.+$/i) {  
    # If a range format is found (e.g. 10^-2 > p > 10^-3)
    if ($pval =~ /^\d+\^(-\d+)/) {
      if (length("$1")==1) { $1 = "0$1"; } 
      $sci_pval = "1.00e$1"; # e.g 10^-2 > p > 10^-3 => 1.00e-2
    }
    # If a decimal format is found (e.g. 0.0023)
    elsif ($pval =~ /^\d+/){
      $sci_pval = $pval;
    #$sci_pval = sprintf("%.2e",$pval); # e.g. 0.002 => 2,30e-3
    }
    elsif ($pval =~ /^\w+/) {
      $sci_pval = "NULL";
    }
    
  }
  else {
    $pval =~ tr/E/e/;
    if ($pval =~ /^(\d+)(e-?\d+)/) {
      $pval="$1.00$2";  
    }
    if ($pval =~ /^(\d+\.\d{1})(e-?\d+)/) {
      $pval="$1"."0$2";  
    }
    if ($pval =~ /^(\d+\.\d+e-?)(\d{1})$/) {
      $pval = "$1"."0$2";
    }
    $sci_pval = $pval;
  }
  return $sci_pval;
}


# Method to remove the carriage return character in each line of the input file
sub remove_carriage_return {
  my $infile = shift;
  my $replacement = shift;
  
  $replacement ||= "";
  
  my @ARGV_bak = @ARGV;
  @ARGV = ($infile);
  $^I = ".bak2";
  while (<>) {
    s/\r/$replacement/g;
    print;
  }
  # Restore the @ARGV variable
  @ARGV = @ARGV_bak;
}


sub usage {
  
  print qq{
  Usage: perl import_variation_annotation.pl [OPTION]
  
  Import variation annotation phenotype data into a Variation database
  
  Options:
    
      -verbose           Progress information is printed
      -help               Print this message
      
      -skip_phenotypes   Skip the study, variation_annotation and phenotype tables insertions.
      -skip_synonyms     Skip the variation_synonym table insertion.
      -skip_sets         Skip the variation_set_variation table insertion.
      
    Database credentials are specified on the command line
    
      -host      Variation database host name (Required)
      -dbname    Variation database name (Required)
      -user      Variation database user (Required)
      -pass      Variation database password (Required)
      -port      Variation database port (Default: 3306)
      
    An input file must be specified. This file contains the data that will be imported, typically tab-delimited
    and obtained from the UniProt or NHGRI GWAS catalog. If a new source is required, a method for parsing the
    source format into a "standard" data structure can be added.
    
      -infile          Typically a tab-delimited file (Required)
      
    The source of the data must be specified so that the correct parser and source table entry will be used. Currently
    supported sources are 'uniprot' and 'nhgri'.
    
      -source             String indicating the source of the data (Required)
      -version           Numerical version of the source (Required)
  } . "\n";
  exit(0);
}
