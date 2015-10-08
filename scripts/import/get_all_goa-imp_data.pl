#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

our ($species, $output_file, $registry_file, $help);

GetOptions('species=s'     => \$species,
           'output_file=s' => \$output_file,
           'registry=s'    => \$registry_file,
           'help!'         => \$help,
          );
          
$registry_file ||= "./ensembl.registry";
$species ||= 'human';

usage() if ($help);

usage('-output_file argument is required') if (!$output_file);

my $ens_type    = 'Translation';
my $ens_linkage = 'IMP';
my $ens_prot_db = 'UniProtKB/Swiss-Prot';
my $pmid_col    = 'REFERENCE';
my $rest_url    = 'http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&goid=%s&evidence=%s&protein=%s';

# connect to databases
Bio::EnsEMBL::Registry->load_all( $registry_file );
my $db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $dbCore = $db_adaptor->dbc->db_handle;

my $sql1 = qq{
  SELECT t.translation_id, t.stable_id, t.transcript_id, x1.dbprimary_acc, x2.display_label, x2.description
  FROM translation t, ontology_xref o, object_xref obj, xref x1, xref x2
  WHERE
    o.object_xref_id=obj.object_xref_id AND
    o.source_xref_id=x1.xref_id AND
    t.translation_id=obj.ensembl_id AND
    obj.xref_id=x2.xref_id AND
    o.linkage_type="$ens_linkage" AND
    obj.ensembl_object_type="$ens_type"         
};

my $sql2 = qq{
  SELECT g.stable_id FROM gene g, transcript t WHERE t.gene_id=g.gene_id AND t.transcript_id=?
};


my %gene_list;
my %go_uniprot;

my $sth1 = $db_adaptor->dbc->prepare($sql1);
my $sth2 = $db_adaptor->dbc->prepare($sql2);

my ($translation_id, $translation, $transcript_id, $uniprot_id, $go_term, $go_desc);
$sth1->execute();
$sth1->bind_columns(\$translation_id, \$translation, \$transcript_id, \$uniprot_id, \$go_term, \$go_desc);
my $count_query_result = 0;
while(my @res = $sth1->fetchrow_array()) {
  
  # Gene
  my $gene;
  $sth2->execute($transcript_id);
  $sth2->bind_columns(\$gene);
  $sth2->fetch();

  $go_uniprot{$go_term}{'go_desc'} = $go_desc;
  
  # UniProt
  $gene_list{$gene}{$go_term}{'uniprot_list'}{$uniprot_id} = 1;
  $go_uniprot{$go_term}{'uniprot_list'}{$uniprot_id} = 1;

  $sth2->finish();
  $count_query_result++;
}
$sth1->finish();


######
# 1 - Run web services
# 2 - Parse results (search results with PMIDs)
# 3 - Export data in a tabulated datafile
# 4 - Import the file using the script "import_phenotype_data.pl"
######
my $http = HTTP::Tiny->new();
my $server = 'http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&goid=%s&evidence=%s&protein=%s';

foreach my $go (keys(%go_uniprot)) {
  foreach my $uniprot (keys(%{$go_uniprot{$go}{'uniprot_list'}})) {
    #my $rest_url_with_params = "goid=$go&evidence=$ens_linkage&protein=$uniprot";
    my $rest_url_with_params = sprintf($rest_url, $go, $ens_linkage, $uniprot);
    my $response = $http->get($rest_url_with_params);
       
    die "Failed for GO: $go and UniProt: $uniprot! ($server.$ext)\n" unless $response->{success};
    
    my $pmids = parse_result($response->{content});
    
    $go_uniprot{$go}{'uniprot_list'}{$uniprot} = join(',',@{$pmids});
  }
}

open OUT, "> $output_file" or die $!;
print OUT "Gene_ID\tGO_term\tGO_description\tUniProt_ID\tPMID(s)\n";
my $count_rows = 0;
foreach my $gene (keys(%gene_list)) {
  foreach my $go_term (keys(%{$gene_list{$gene}})) {
    foreach my $uniprot (keys(%{$gene_list{$gene}{$go_term}{'uniprot_list'}})) {
    
      next if (!$go_uniprot{$go_term}{'uniprot_list'}{$uniprot});
      
      my $go_desc = $go_uniprot{$go_term}{'go_desc'};
      my $pmids = ($go_uniprot{$go_term}{'uniprot_list'}{$uniprot} == 1) ? '' : $go_uniprot{$go_term}{'uniprot_list'}{$uniprot};
      
      print OUT "$gene\t$go_term\t$go_desc\t$uniprot\t$pmids\n";
      $count_rows++;
    }
  }
}
close(OUT);

print STDOUT "Count query results = $count_query_result\n";
print STDOUT "Count rows = $count_rows\n";


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

sub usage {
  my $msg = shift;
     $msg ||= '';
  
  print qq{$msg
  
  Usage: perl get_go-imp_data.pl [OPTION]
  
  Fetch and store the GO-IMP data into a tabulated text file. This file will be used to import the data in the Variation
  database, using the script "import_phenotype_data.pl" with the option "-source go-imp"
  
  Options:
    
      -help              Print this message
      
      -species       Species from which you want to get the data. (optional)
                     The default value is 'human'.
      -output_file   Path and name of the output text file. (Required)
      -registry      Path to the registry file. (optional)
                     The default value is './ensembl.registry'.
  } . "\n";
  exit(0);
}

