# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

# Script to generate an HTML page containing the variant sources of each species


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut


use Bio::EnsEMBL::Registry;
use DBI;
use strict;
use POSIX;
use Getopt::Long;
use XML::Simple;
use HTTP::Tiny;

###############
### Options ###
###############
my ($e_version,$user,$pass,$hlist,$help);

## EG options
my ($site, $etype);

usage() if (!scalar(@ARGV));
 
GetOptions(
     'v=s'        => \$e_version,
     'user|u=s'   => \$user,
     'pass|p=s'   => \$pass,
     'help!'      => \$help,
     'hlist=s'    => \$hlist
);

if ($help) {
  usage();
} elsif (!$e_version) {
  print STDERR "> Error! Please give an Ensembl version, using the option '-v' \n";
  usage();
} elsif (!$user) {
  print STDERR "> Error! Please give a database user name, using the option '-user' \n";
  usage();
} elsif (!$pass) {
  print STDERR "> Error! Please give a database password, using the option '-pass' \n";
  usage();
} elsif (!$hlist) {
  print STDERR "> Error! Please give the list of host names where the new databases are stored using the option '-hlist'\n";
  usage();
}



# Settings
my $database = "";
my @hostnames = split /,/, $hlist;

# Loop over hosts
foreach my $hostname (@hostnames) {

  my $sql = qq{SHOW DATABASES LIKE '%variation_$e_version%'};
  my $sth_h = get_connection_and_query($database, $hostname, $sql, 1);
  my $db_found = 0;
  
  # Loop over databases
  while (my ($dbname) = $sth_h->fetchrow_array) {
    next if ($dbname !~ /^[a-z]+_[a-z]+_variation_\d+_\d+$/i);
    next if ($dbname =~ /^master_schema/ || $dbname =~ /^homo_sapiens_variation_\d+_37$/ || $dbname =~ /private/);

    $dbname =~ /^(.+)_variation/;
    my $s_name = $1;
    print STDERR "# $s_name - [ $dbname ]\n";
    
    ## Source (data types)
    post_process_source($hostname,$dbname);
    
    ## Phenotype (related tables)
    post_process_phenotype($hostname,$dbname);

    ## Population (size)
    post_process_population($hostname,$dbname);
    
    ## Clinical significance (human only)
    if ($dbname =~ /^homo_sapiens_variation_$e_version/) {
      post_process_clin_sign($hostname,$dbname);
      post_process_publication($hostname,$dbname);
    }
  }
  $sth_h->finish;
}


# Source (data types)
sub post_process_source () {
  my $host = shift;
  my $db   = shift;

  print STDERR "\t Source:";

  my $table_tmp_1 = 'tmp_source_type';
  my $table_tmp_2 = 'tmp_source_type_grouped';

  my $sql_source_create_1 = "CREATE TABLE $table_tmp_1 ( source_id INT, source_type VARCHAR(255));";
  my $sql_source_create_2 = "CREATE TABLE $table_tmp_2 SELECT source_id, group_concat(source_type) AS data_types FROM $table_tmp_1 GROUP BY source_id";
  my @sql_source_ins = (
    "INSERT INTO $table_tmp_1 SELECT s.source_id, 'variation' FROM source s, variation v WHERE s.source_id = v.source_id GROUP BY s.source_id",
    "INSERT INTO $table_tmp_1 SELECT s.source_id, 'variation_synonym' FROM source s, variation_synonym v WHERE s.source_id = v.source_id GROUP BY s.source_id",
    "INSERT INTO $table_tmp_1 select s.source_id, 'structural_variation' FROM source s, structural_variation v WHERE s.source_id = v.source_id GROUP BY s.source_id",
    "INSERT INTO $table_tmp_1 SELECT s.source_id, 'phenotype_feature' FROM source s, phenotype_feature v WHERE s.source_id = v.source_id GROUP BY s.source_id",
    "INSERT INTO $table_tmp_1 SELECT s.source_id, 'study' FROM source s, study v WHERE s.source_id = v.source_id GROUP BY s.source_id"
  );
  my $sql_source_update_2 = "UPDATE source s, $table_tmp_2 t SET s.data_types = t.data_types WHERE s.source_id = t.source_id;";
  
  # Create tmp table #1
  get_connection_and_query($db, $host, $sql_source_create_1);

  # Update tmp table #1
  foreach my $sql_source_i (@sql_source_ins) {
    get_connection_and_query($db, $host, $sql_source_i);
  }
  # Create tmp table #2
  get_connection_and_query($db, $host, $sql_source_create_2);
  # Update tmp table #2
  get_connection_and_query($db, $host, $sql_source_update_2);

  # Delete tmp table #1
  get_connection_and_query($db, $host, "DROP TABLE $table_tmp_1");
  # Delete tmp table #2
  get_connection_and_query($db, $host, "DROP TABLE $table_tmp_2");

  print STDERR " update data types done";
  print STDERR "\n";
}


# Phenotype related tables
sub post_process_phenotype () {
  my $host = shift;
  my $db   = shift;

  print STDERR "\t Phenotype:";

  my $sql_phe_1 = "DELETE FROM phenotype WHERE phenotype_id NOT IN (SELECT phenotype_id FROM phenotype_feature)";
  my $sql_phe_2 = "DELETE FROM phenotype_ontology_accession WHERE phenotype_id NOT IN (SELECT phenotype_id FROM phenotype)";

  my @tables = ('phenotype','phenotype_feature','phenotype_ontology_accession','phenotype_feature_attrib');

  # Cleanup phenotype related tables
  get_connection_and_query($db, $host, $sql_phe_1);
  get_connection_and_query($db, $host, $sql_phe_2);
  print STDERR " Cleanup done";

  # Optimise phenotype related tables
  foreach my $table (@tables) {
    get_connection_and_query($db, $host, "OPTIMIZE TABLE ".$table);
  }
  print STDERR " | Optimizing: done";
  print STDERR "\n";
}


# Population table (size)
sub post_process_population () {
  my $host = shift;
  my $db   = shift;

  print STDERR "\t Population:";

  my $sql_pop_1 = "UPDATE population p SET p.size=(SELECT count(s.sample_id) from sample_population s WHERE s.population_id=p.population_id) WHERE p.size is null;";
  my $sql_pop_2 = "UPDATE population SET size=NULL WHERE size=0;";

  # Update population table
  get_connection_and_query($db, $host, $sql_pop_1);
  get_connection_and_query($db, $host, $sql_pop_2);
  print STDERR " size update done";
  print STDERR "\n";
}


# Clinical significance columns (only human)
sub post_process_clin_sign () {
  my $host = shift;
  my $db   = shift;

  print STDERR "\t Clinical significance:";

  my @tables = ('variation', 'variation_feature', 'structural_variation');

  # Update population table
  foreach my $table (@tables) {
    get_connection_and_query($db, $host, "UPDATE $table SET clinical_significance=NULL WHERE clinical_significance='';");
  }

  # Update population table
  print STDERR " clinical_significance columns update done";
  print STDERR "\n";
}

# Connects and execute a query
sub get_connection_and_query {
  my $dbname  = shift;
  my $hname   = shift;
  my $sql     = shift;
  my $rtn_sth = shift;
  
  my ($host, $port) = split /\:/, $hname;
  
  # DBI connection 
  my $dsn = "DBI:mysql:$dbname:$host:$port";
  my $dbh = DBI->connect($dsn, $user, $pass) or die "Connection failed";

  my $sth = $dbh->prepare($sql);
  $sth->execute;
  
  if ($rtn_sth) {
    return $sth;
  }
  else {
    $sth->finish;
  }
}

sub post_process_publication {
  my $host_port = shift;
  my $db        = shift;

  my ($host, $port) = split /\:/, $host_port;

  print STDERR "\t Publications from Phenotype feature:";

  my $reg = 'Bio::EnsEMBL::Registry';
  $reg->load_registry_from_db(
    -host => $host,
    -user => $user,
    -port => $port,
    -pass => $pass,
    # -db_version => $e_version # Testing
  );

  my $species = 'homo_sapiens';

  my $dba = $reg->get_DBAdaptor($species, 'variation') || die "Error getting db adaptor\n";

  ## extract all variants - cited variants failing QC are still displayed
  $dba->include_failed_variations(1);

  my $var_ad = $reg->get_adaptor($species, 'variation', 'variation');
  my $pub_ad = $reg->get_adaptor($species, 'variation', 'publication');
  my $source_ad = $reg->get_adaptor($species, 'variation', 'source');

  # Fetch all citation_source attribs
  my $citation_attribs = get_citation_attribs($dba);

  process_phenotype_feature($species, $dba, $var_ad, $pub_ad, $source_ad, $citation_attribs);
  process_phenotype_feature_attrib($species, $dba, $var_ad, $pub_ad, $source_ad, $citation_attribs);

}

sub process_phenotype_feature {
  my $species = shift;
  my $dba = shift;
  my $var_ad = shift;
  my $pub_ad = shift;
  my $source_ad = shift;
  my $citation_attribs = shift;

  ## Get studies from phenotype_feature
  my $attrib_ext_sth = $dba->dbc()->prepare(qq[ select s.study_id, s.source_id, s.external_reference, s.study_type
                                                from study s 
                                                inner join phenotype_feature p on s.study_id = p.study_id
                                                where p.type = 'variation' and p.study_id is not null and s.external_reference is not null
                                                group by s.study_id, s.source_id, s.external_reference, s.study_type ]);
  $attrib_ext_sth->execute()||die;
  my $data = $attrib_ext_sth->fetchall_arrayref();

  my $rsid_sth = $dba->dbc()->prepare(qq[ select object_id
                                          from phenotype_feature 
                                          where study_id = ? ]);

  foreach my $l (@{$data}){
    my $study_id = $l->[0];
    my $source_id = $l->[1];
    my $external_reference = $l->[2];
    $external_reference =~ s/PMID://;
    my $study_type = $l->[3];

    # Get publication with same PMID from study table
    my $publication = $pub_ad->fetch_by_pmid($external_reference);

    next unless !defined($publication);

    # Get attrib id for source - some are null 
    my $source_attrib_id;
    if(defined $study_type){
      $source_attrib_id = $citation_attribs->{$study_type};
      die "No attrib of type 'citation_source' was found for '$study_type'!\n" unless defined $source_attrib_id;
    }
    else{
      # Get source name from source table (dbGaP)
      my $source_obj = $source_ad->fetch_by_dbID($source_id);
      my $source_name = $source_obj->name();
      $source_attrib_id = $citation_attribs->{$source_name};
      die "No attrib of type 'citation_source' was found for '$source_name'!\n" unless defined $source_attrib_id;
    }

    # Get publication that is not in publication table 
    my $pub_data = get_epmc_data($external_reference);

    if(defined $pub_data->{resultList}->{result}->{title}){
      my $pub_title = $pub_data->{resultList}->{result}->{title};
      my $pub_pmcid = $pub_data->{resultList}->{result}->{pmcid};
      my $pub_author = $pub_data->{resultList}->{result}->{authorString};
      my $pub_year = $pub_data->{resultList}->{result}->{pubYear};
      my $pub_doi = $pub_data->{resultList}->{result}->{doi};
      
      my @variant_list;
      $rsid_sth->execute($study_id);
      my $rsids = $rsid_sth->fetchall_arrayref();
      
      foreach my $rsid (@{$rsids}){
        my $v = $var_ad->fetch_by_name($rsid->[0]);
        if (defined $v){
          push @variant_list, $v;
        }
      }

      ### create new object
      my $publication = Bio::EnsEMBL::Variation::Publication->new( 
                -title    => $pub_title,
                -authors  => $pub_author,
                -pmid     => $external_reference,
                -pmcid    => $pub_pmcid || undef,
                -year     => $pub_year,
                -doi      => $pub_doi,
                -ucsc_id  => undef,
                -variants => \@variant_list,
                -adaptor  => $pub_ad
                );

       $pub_ad->store($publication,$source_attrib_id);
    }
  }
}

sub process_phenotype_feature_attrib {
  my $species = shift;
  my $dba = shift;
  my $var_ad = shift;
  my $pub_ad = shift;
  my $source_ad = shift;
  my $citation_attribs = shift;

  my $attrib_type_sth = $dba->dbc()->prepare(qq[ select attrib_type_id
                                                 from attrib_type
                                                 where code = 'pubmed_id' ]);

  $attrib_type_sth->execute()||die;
  my $attrib_type_id = $attrib_type_sth->fetchall_arrayref();

  die "No attribute type 'pubmed_id' found\n" unless defined $attrib_type_id->[0]->[0];

  my $attrib_type_id_v = $attrib_type_id->[0]->[0];

  my $pheno_feature_sth = $dba->dbc()->prepare(qq[ select phenotype_feature_id, value 
                                                   from phenotype_feature_attrib 
                                                   where attrib_type_id = $attrib_type_id_v ]);

  my $pheno_sth = $dba->dbc()->prepare(qq[ select source_id, object_id 
                                           from phenotype_feature
                                           where phenotype_feature_id = ? ]);

  $pheno_feature_sth->execute()||die;
  my $pheno_feature_data = $pheno_feature_sth->fetchall_arrayref();

  # PMID -> list of phenotype_feature_id
  my %new_publications;

  # Create map to associate PMID with all phenotype features it's linked to
  foreach my $pheno_feat_data (@{$pheno_feature_data}){
    my $pheno_feat_id = $pheno_feat_data->[0];
    my $value_pubid = $pheno_feat_data->[1];
    
    my @value_pubid_splited = split /,/, $value_pubid;

    # Get publication with same PMID from study table
    foreach my $pmid (@value_pubid_splited){
      # PMID=25806920 is the same paper as PMID=25806919
      # PMID=25806919 is already in the publication table - faster to skip this PMID and avoid duplicated data
      next if $pmid eq '25806920';

      my $publication = $pub_ad->fetch_by_pmid($pmid);

      next unless !defined($publication);

      my @feature_id;
      if(defined $new_publications{$pmid}){
        push(@{$new_publications{$pmid}}, $pheno_feat_id);
      }
      else{
        push @feature_id, $pheno_feat_id;
        $new_publications{$pmid} = \@feature_id;
      }
    }
  }

  foreach my $pmid (keys(%new_publications)){
    my $feature_id_list = $new_publications{$pmid};

    my $pub_data = get_epmc_data($pmid);

    if(defined $pub_data->{resultList}->{result}->{title}){
      my $pub_title = $pub_data->{resultList}->{result}->{title};
      my $pub_pmcid = $pub_data->{resultList}->{result}->{pmcid};
      my $pub_author = $pub_data->{resultList}->{result}->{authorString};
      my $pub_year = $pub_data->{resultList}->{result}->{pubYear};
      my $pub_doi = $pub_data->{resultList}->{result}->{doi};

      my @variant_list;
      my $source_attrib_id;

      foreach my $feature_id (@{$feature_id_list}){
        # Get source and variant id
        $pheno_sth->execute($feature_id)||die;
        my $pheno_data = $pheno_sth->fetchall_arrayref();

        next unless defined $pheno_data->[0]->[0];

        my $source_id = $pheno_data->[0]->[0];
        my $var_name = $pheno_data->[0]->[1];

        my $v = $var_ad->fetch_by_name($var_name);
        # Add variant if it exists - avoid undef values
        if (defined $v){
          push @variant_list, $v;
        }

        # Get source name from source table (ClinVar)
        my $source_obj = $source_ad->fetch_by_dbID($source_id);
        my $source_name = $source_obj->name();
        $source_attrib_id = $citation_attribs->{$source_name};
        die "No attrib of type 'citation_source' was found for '$source_name'!\n" unless defined $source_attrib_id;
      }

      # Clean title before insertion
      if($pub_title =~ /^\[ ?[A-Za-z]{2}/){
          $pub_title =~ s/^\[//;
          $pub_title =~ s/\]//;
      }

      # Skip publication if title 'Not Available'
      next if ($pub_title =~ /Not Available/);

      ## create new object
      my $publication = Bio::EnsEMBL::Variation::Publication->new( 
                -title    => $pub_title,
                -authors  => $pub_author,
                -pmid     => $pmid,
                -pmcid    => $pub_pmcid || undef,
                -year     => $pub_year,
                -doi      => $pub_doi,
                -ucsc_id  => undef,
                -variants => \@variant_list,
                -adaptor  => $pub_ad
                );

      $pub_ad->store($publication,$source_attrib_id);
    }
  }
}

sub get_epmc_data{

  my $id = shift; ## specific part of URL including pmid or pmcid

  my $http = HTTP::Tiny->new();

  return undef unless defined $id && $id =~/\d+/;

  my $xs = XML::Simple->new();
  # Get from source medline
  my $request = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=ext_id:' . $id . '%20src:med';

  my $response = $http->get($request, {
      headers => { 'Content-type' => 'application/xml' }
                            });
  unless ($response->{success}){
      warn "Failed request: $request :$!\n" ;
      return;
  }
  return $xs->XMLin($response->{content} );
}

# Fetch all attribs that have type citation_source
sub get_citation_attribs {
  my $dba = shift;

  my %attribs;
  my $attrib_type_id;

  my $stm_1 = $dba->dbc()->prepare(qq[ SELECT attrib_type_id FROM attrib_type WHERE code = 'citation_source' ]);
  $stm_1->execute();
  my $attrib_type = $stm_1->fetchall_arrayref();

  die "No attribute of type 'citation_source' was found in attrib_type table!\n" unless defined $attrib_type->[0];

  $attrib_type_id = $attrib_type->[0]->[0];

  my $stm_2 = $dba->dbc()->prepare(qq[ SELECT attrib_id,value FROM attrib where attrib_type_id = $attrib_type_id ]);
  $stm_2->execute();
  my $citation_attribs = $stm_2->fetchall_arrayref();

  die "No attribute of type 'citation_source' was found in attrib table!\n" unless defined $citation_attribs->[0];

  foreach my $data (@{$citation_attribs}){
    my $id = $data->[0];
    my $name = $data->[1];
    $attribs{$name} = $id;
  }

  return \%attribs;
}

sub usage {
  
  print qq{
  Usage: perl post_process_databases_tables.pl [OPTION]
  
  Post processing of some of the Ensembl Variation tables.
  
  Options:

    -help       Print this message
      
    -v          Ensembl version, e.g. 90 (Required)
    -user|u     Database login user name (Required)
    -pass|p     Database user password name (Required)
    -hlist      The list of host names (with port) where the new databases are stored, separated by a coma,
                e.g. ensembldb.ensembl.org1:3334, ensembldb.ensembl.org2:1234 (Required)
  } . "\n";
  exit(0);
}
