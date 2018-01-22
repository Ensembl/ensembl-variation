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
use HTTP::Tiny;
use Getopt::Long;

our ($output_dir, $e_version, $species, $hlist, $user, $port, $help);

GetOptions('output_dir=s' => \$output_dir,
           'v=s'          => \$e_version,
           'species=s'    => \$species,
           'hlist=s'      => \$hlist,
           'user=s'       => \$user,
           'help!'        => \$help,
          );

usage() if ($help);

usage('-output_dir argument is required') if (!$output_dir);
die ("Output dir '$output_dir' doesn't exist") unless (-d $output_dir);

if (!$e_version) {
  print "> Error! Please give an Ensembl version, using the option '-v' \n";
  usage();
}
if (!$hlist) {
  print "> Error! Please give the list of host names where the new databases are stored using the option '-hlist'\n";
  usage();
}
if (!$user) {
  print "> Error! Please give user name using the option '-user'\n";
  usage();
}

my $ens_type    = 'Translation';
my $ens_linkage = 'IMP';
my $ens_prot_db = '%UniProtKB%';
my $pmid_col    = 'REFERENCE';
my $rest_url    = 'http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&goid=%s&evidence=%s&protein=%s';

my $database = "";
my $pswd = "";
my $db_core_type = 'core';
my $db_var_type  = 'variation';
my @hostnames = split /,/, $hlist;
my %species_list;
%species_list = map {$_ => 1} split /,/, $species if ($species);
my %species_found;


my $sql_var  = qq{SHOW DATABASES LIKE '%$db_var_type\_$e_version%'};
my $sql_core = qq{SHOW DATABASES LIKE '%$db_core_type\_$e_version%'};

my $sql1 = qq{
  SELECT t.translation_id, t.stable_id, t.transcript_id, x1.dbprimary_acc, x2.display_label, x2.description
  FROM translation t, ontology_xref o, object_xref obj, xref x1, xref x2, external_db e
  WHERE
    o.object_xref_id=obj.object_xref_id AND
    o.source_xref_id=x1.xref_id AND
    t.translation_id=obj.ensembl_id AND
    obj.xref_id=x2.xref_id AND
    x1.external_db_id=e.external_db_id AND
    e.db_display_name like "$ens_prot_db" AND
    o.linkage_type="$ens_linkage" AND
    obj.ensembl_object_type="$ens_type"         
};


foreach my $hostname (@hostnames) {
  
  my $db_type = (scalar(keys(%species_list))) ? $db_core_type : $db_var_type;

  my $sql = (scalar(keys(%species_list))) ? $sql_core : $sql_var;
  my $sth = get_connection_and_query($database, $hostname, $sql);
  
  # loop over databases
  while (my ($dbname) = $sth->fetchrow_array) {
    next if ($dbname =~ /^master_schema/);
    next if ($dbname =~ /sample$/);
    
    my %gene_list;
    my %go_uniprot;
    my %gene_tr;
    
    $dbname =~ /^(.+)_$db_type/;
    my $s_name = $1;
    
    print "# $s_name\n";
    print STDERR "# $s_name";

    if (%species_list) {
      if (!$species_list{$s_name}) {
        print STDERR " ... skipped\n";
        next;
      }
      $species_found{$s_name} = 1;
    }
    print STDERR "\n";

    my $output_file = $output_dir."/".$s_name."_imp.txt"; 
    
    my $core_dbname = $dbname;
       $core_dbname =~ s/$db_type/$db_core_type/i;

    my ($translation_id, $translation, $transcript_id, $uniprot_id, $go_term, $go_desc);
    my $sth1 = get_connection_and_query($core_dbname, $hostname, $sql1);
    $sth1->bind_columns(\$translation_id, \$translation, \$transcript_id, \$uniprot_id, \$go_term, \$go_desc);
    my $count_query_result = 0;
    
    while(my @res = $sth1->fetchrow_array()) {
  
      # Gene
      my $gene = ($gene_tr{$transcript_id}{'gene'}) ? $gene_tr{$transcript_id}{'gene'} : undef;
      my $transcript = ($gene_tr{$transcript_id}{'transcript'}) ? $gene_tr{$transcript_id}{'transcript'} : undef;
      if (!$gene || !$transcript) {
        my $sql2 = qq{
          SELECT g.stable_id, t.stable_id FROM gene g, transcript t WHERE t.gene_id=g.gene_id AND t.transcript_id=$transcript_id
        };      
        my $sth2 = get_connection_and_query($core_dbname, $hostname, $sql2);
        $sth2->bind_columns(\$gene, \$transcript);
        $sth2->fetch();
        $sth2->finish();
        $gene_tr{$transcript_id} = {'gene' => $gene, 'transcript' => $transcript};
      }
      
      next if (!$gene && !$transcript);

      $go_uniprot{$go_term}{'go_desc'} = $go_desc;
  
      # Ensembl transcript
      $gene_list{$gene}{$go_term}{'transcript_list'}{$transcript} = 1;
  
      # UniProt
      $gene_list{$gene}{$go_term}{'uniprot_list'}{$uniprot_id} = 1;
      $go_uniprot{$go_term}{'uniprot_list'}{$uniprot_id} = 1;

      $count_query_result++;
    }

    $sth1->finish();


    if (scalar(keys(%go_uniprot)) == 0 || !%go_uniprot) {
      print STDERR "No data\n";
      next;
    }
    
    print STDERR "Count query results = $count_query_result\n";
    
    ######
    # 1 - Run web services
    # 2 - Parse results (search results with PMIDs)
    # 3 - Export data in a tabulated datafile
    # 4 - Import the file using the script "import_phenotype_data.pl"
    ######
    my $http = HTTP::Tiny->new();

    foreach my $go (keys(%go_uniprot)) {
      foreach my $uniprot (keys(%{$go_uniprot{$go}{'uniprot_list'}})) {
        my $rest_url_with_params = sprintf($rest_url, $go, $ens_linkage, $uniprot);
        my $response = $http->get($rest_url_with_params);
       
        if ($response->{success}) {
    
          my $pmids = parse_result($response->{content});
    
          $go_uniprot{$go}{'uniprot_list'}{$uniprot} = join(',',@{$pmids});
        }
        else {
          delete $go_uniprot{$go}{'uniprot_list'}{$uniprot};
        }
      }
    }

    
    open OUT, "> $output_file" or die $!;
    print OUT "Gene_ID\tGO_term\tGO_description\tUniProt_IDs\tPMIDs\n";
    my $count_rows = 0;
    foreach my $gene (keys(%gene_list)) {
      foreach my $go_term (keys(%{$gene_list{$gene}})) {
        my %uniprot_ids;
        my %pmids;
        my $go_desc = $go_uniprot{$go_term}{'go_desc'};
        foreach my $uniprot (keys(%{$gene_list{$gene}{$go_term}{'uniprot_list'}})) {
      
          next if (!$go_uniprot{$go_term}{'uniprot_list'}{$uniprot});
          
          $uniprot_ids{$uniprot} = 1;
          
          if ($go_uniprot{$go_term}{'uniprot_list'}{$uniprot} != 1) {
            my @pubmeds_list = split(',',$go_uniprot{$go_term}{'uniprot_list'}{$uniprot});
            $pmids{$_} = 1 for @pubmeds_list;
          }
        }
        
        if (scalar keys(%uniprot_ids)) {
          print OUT "$gene\t$go_term\t$go_desc\t".join(',',keys(%uniprot_ids))."\t".join(',',keys(%pmids))."\n";
          $count_rows++;
          print STDERR "$count_rows entries done\n" if ($count_rows =~ /^\d+000$/);
        }
      }
    }
    close(OUT);

    print STDERR "Count rows = $count_rows\n";
  }
}

# List the species not found
if (%species_list) {
  if (scalar(keys(%species_list)) != scalar(keys(%species_found))) {
    foreach my $sp (keys(%species_list)) {
      print STDERR "Species '$sp' not found in $hlist\n" if (!$species_found{$sp});
    }
  }
}

sub parse_result {
  my $content = shift;
  
  my %headers;
  my %pmids;
  
  my @rows = split(/\n/, $content);
  
  
  my @header_row = split(/\t/,shift(@rows));
  $headers{uc($header_row[$_])} = $_ for 0..$#header_row;
  
  foreach my $row (@rows) {
    my @row_data = split(/\t/,$row);
    my %content;
    $content{$_} = $row_data[$headers{$_}] for keys %headers;
    
    next if ($content{$pmid_col} eq '-');
    
    $pmids{$content{$pmid_col}} = 1;
  }
  
  my @pmids_list = (scalar(keys(%pmids))) ? keys(%pmids) : ();
  
  return \@pmids_list;
}


# Connects and execute a query
sub get_connection_and_query {
  my $dbname = shift;
  my $hname  = shift;
  my $sql    = shift;
  
  my ($host, $port) = split /\:/, $hname;
  
  # DBI connection 
  my $dsn = "DBI:mysql:$dbname:$host:$port";
  my $dbh = DBI->connect($dsn, $user, $pswd) or die "Connection failed";

  my $sth = $dbh->prepare($sql);
  $sth->execute;
  return $sth;
}


sub usage {
  my $msg = shift;
     $msg ||= '';
  
  print qq{$msg
  
  Usage: perl get_all_goa-imp_data.pl [OPTION]
  
  Fetch and store the GO-IMP data into a tabulated text file. This file will be used to import the data in the Variation
  database, using the script "import_phenotype_data.pl" with the option "-source go-imp"
  
  Options:
    
      -help        Print this message

      -v           Ensembl version, e.g. 82 (Required)
      -output_dir  Path to the output files. (Required)
      -hlist       The list of host names (with port) where the new databases are stored, separated by a coma,
                    e.g. ensembldb.ensembl.org1:3334, ensembldb.ensembl.org2:1234 (Required)
      -user        MySQL user name (Required)
      -species     The list of species from which the GOA data will be retrieved, separated by a coma,
                   e.g. homo_sapiens,mus_musculus. By default the script retrieves GOA data from every species (optional)
  } . "\n";
  exit(0);
}

