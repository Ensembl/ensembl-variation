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
use JSON;

# Mean to parse the schema version 1.2.3

our ($species, $input_file, $output_file, $registry_file, $help);

GetOptions('species=s'     => \$species,
           'input_file=s'  => \$input_file,
           'output_file=s' => \$output_file,
           'registry=s'    => \$registry_file,
           'help!'         => \$help,
          );
          
$registry_file ||= "./ensembl.registry";
$species ||= 'human';

usage() if ($help);

usage('-input_file argument is required') if (!$input_file);
usage('-output_file argument is required') if (!$output_file);

my %efos;
my %data;

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all( $registry_file );

my $ga = $registry->get_adaptor($species,'core','gene');
die("ERROR: Could not get gene adaptor") unless defined($ga);

my $ota = $registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );
die("ERROR: Could not get ontology term adaptor") unless defined($ota);

# Source: Cancer Gene Census (http://cancer.sanger.ac.uk/census/)

if($input_file =~ /gz$/) {
  open IN, "zcat $input_file |" or die ("Could not open $input_file for reading");
}
else {
  open(IN,'<',$input_file) or die ("Could not open $input_file for reading");
}

while(<IN>) {
  chomp $_;
  my $json_hash = decode_json($_);
  
  my $source = $json_hash->{'sourceID'};
  next unless ($source =~ /cancer_gene_census/);
  
  my $type          = $json_hash->{'type'};
  my $gene_symbol   = $json_hash->{'target'}{'gene_info'}{'symbol'};
  my $pmids         = parse_publications($json_hash->{'literature'}{'references'});
  my $phenotype_url = $json_hash->{'unique_association_fields'}{'disease_uri'};
  
  # Phenotype fetching and parsing
  $phenotype_url =~ /\/(\w+)$/;
  my $phenotype_id = $1;
  my $phenotype;
  if ($phenotype_id =~ /^EFO/) {
    $phenotype = ($efos{$phenotype_id}) ? $efos{$phenotype_id} : get_phenotype_desc($phenotype_id);
  }
  
  if (!$phenotype) {
    print STDERR "$gene_symbol: no phenotype desc found for $phenotype_id\n";
    $phenotype = 'ND';
  }
  else {
    $efos{$phenotype_id} = $phenotype;
  }
  
  my $genes = $ga->fetch_all_by_external_name($gene_symbol, 'HGNC');		
	# we don't want any LRG genes
  @$genes = grep {$_->stable_id !~ /^LRG_/} @$genes;

  my %stable_ids;
  if (scalar @$genes != 1) {
    my $gene_id = $json_hash->{'target'}{'id'};
    $gene_id =~ /(ENSG\d+)$/i;
    if ($1) {
      $stable_ids{$1} = 1;
    }
    print STDERR "WARNING: Found ".(scalar @$genes)." matching Ensembl genes for HGNC ID $gene_symbol\n";
  }

	next unless scalar @$genes;
			
  foreach my $gene(@$genes) {
    my $gene_id = $gene->stable_id;
    next if (%stable_ids && !$stable_ids{$gene_id});
    
    $data{$gene_symbol}{$gene_id}{$phenotype}{'type'} = $type;
    $data{$gene_symbol}{$gene_id}{$phenotype}{'source'} = $source;
    $data{$gene_symbol}{$gene_id}{$phenotype}{'phenotype_id'} = $phenotype_id;
    foreach my $pmid (@$pmids) {
      $data{$gene_symbol}{$gene_id}{$phenotype}{'pmids'}{$pmid} = 1; 
    }
  }
}
close(F);



open OUT, "> $output_file" || die $!;
foreach my $gene_symbol (sort(keys(%data))) {
  foreach my $gene_id (keys(%{$data{$gene_symbol}})) {
    foreach my $phenotype (keys(%{$data{$gene_symbol}{$gene_id}})) {
      my $type         = $data{$gene_symbol}{$gene_id}{$phenotype}{'type'};
      my $source       = $data{$gene_symbol}{$gene_id}{$phenotype}{'source'};
      my $phenotype_id = $data{$gene_symbol}{$gene_id}{$phenotype}{'phenotype_id'};
      my $pmids = join(',', keys(%{$data{$gene_symbol}{$gene_id}{$phenotype}{'pmids'}}));
      print OUT "$gene_symbol\t$gene_id\t$type\t$source\t$phenotype\t$phenotype_id\t$pmids\n";
    }
  }
}

close(OUT);


sub get_phenotype_desc {
  my $id = shift;
  
  $id =~ s/ //g;
  $id =~ s/_/:/g;
  
  my $phenotype;
  
  my $term = $ota->fetch_by_accession($id);

  $phenotype = $term->name if ($term);
  
  return $phenotype;
}

sub parse_publications {
  my $pubs = shift;
  
  my @pmids = ();
  foreach my $pub (@$pubs) {
    my $pub_id = $pub->{'lit_id'};
    $pub_id=~ /^http:\/\/europepmc.org\/abstract\/MED\/(\d+)$/i;
    if ($1) {
      my $pmid = "PMID:$1";
      push @pmids,$pmid;
    }
  }
  return \@pmids;
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
      -input_file    Path and name of the input text file. (Required)
      -output_file   Path and name of the output text file. (Required)
      -registry      Path to the registry file. (optional)
                     The default value is './ensembl.registry'.
  } . "\n";
  exit(0);
}

