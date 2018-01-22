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


## Import cross references between OMIM to phenotype and disease ontologies 
## provided by Orphanet and HPO
## Orphanet import requires OMIM data + studies to be available;


use strict;
use warnings;

use JSON;
use HTTP::Tiny;
use Getopt::Long;
use Data::Dumper;
use Bio::EnsEMBL::Registry;


my ( $registry_file, $data_file, $type);

GetOptions ("data_file=s"  => \$data_file,
            "registry=s"   => \$registry_file,
            "type=s"       => \$type,
            );

usage() unless defined $registry_file && ( $type eq 'Orphanet' ||( $type eq "HP" && defined $data_file)) ;

my $reg = 'Bio::EnsEMBL::Registry';
$reg->no_version_check(); 
$reg->load_all($registry_file);


my $pheno_adaptor  = $reg->get_adaptor('homo_sapiens', 'variation', 'phenotype');
my $phenof_adaptor = $reg->get_adaptor('homo_sapiens', 'variation', 'phenotypefeature');
my $study_adaptor  = $reg->get_adaptor('homo_sapiens', 'variation', 'study'); 

my $http = HTTP::Tiny->new();

my $phenos;

if( $type eq 'Orphanet'){
  $phenos = get_Orphanet();
}
elsif($type eq 'HP'){
  $phenos = get_HP($data_file);
  store_HP($phenos)

}
else{
  usage();
}


## take from OLS
sub get_Orphanet{
   get_ontol();
}

## Page through Orphanet data extracting OMIM xrefs
sub get_ontol{

  my $request = shift;

  $request = 'http://www.ebi.ac.uk/ols/api/ontologies/ordo/terms' unless defined $request;

  my @links;

  #warn "Looking for $request\n\n" ;
  my $response = $http->get($request, {
        headers => { 'Content-type' => 'application/xml' }
                              });
  unless ($response->{success}){
    warn "Failed request: $request :$!\n" ;
    die;
  }

   my $data =  JSON->new->decode($response->{content} );


  foreach my $term (@{$data->{_embedded}->{terms} }){
    foreach my $xref(@{ $term->{annotation}->{hasDbXref}}){
      push @links, [ $term->{obo_id}, $xref ] if $xref =~/OMIM/;;
    } 
  }

  store_links(\@links);

  my $next = $data->{_links}->{next}->{href} if defined $data->{_links}->{next};
  $data ='';

  ## move onto next page
  get_ontol($next) if defined $next;

}

## store Orphanet annotation
sub store_links{

  my $links = shift;

  my $pheno_adaptor  = $reg->get_adaptor('homo_sapiens', 'variation', 'phenotype');
  my $phenof_adaptor = $reg->get_adaptor('homo_sapiens', 'variation', 'phenotypefeature');
  my $study_adaptor  = $reg->get_adaptor('homo_sapiens', 'variation', 'study'); 

  foreach my $link (@{$links}){

    my ($orphanet, $omim) = @{$link};
    $omim =~ s/^O//;
warn "Storing link $orphanet, $omim\n";
    ## OMIM id held in the study record 
    my $study = $study_adaptor->fetch_all_by_external_reference($omim);
    next unless defined $study->[0];

    if (defined $study->[1]){
      print "Two studies found for $omim\n";
      next;
    }

    ## count the number of phenotypes attached to the study and don't proceed if there are more than one. 
    my $feats = $phenof_adaptor->fetch_all_by_Study($study->[0]);
    next unless $feats->[0];

    my %phen;
    foreach my $feat (@{$feats}){
      $phen{$feat->phenotype()->description()} =1;
    }
    if(scalar(keys %phen) >1){
      warn "Skipping - 2 pheno for $omim\n";
      print join(",", (keys %phen)). "\n";
      next;
    }

    ## if only one phenotype attach Orphanet accession
    my $omim_pheno = $feats->[0]->phenotype();
    $omim_pheno->add_ontology_accession( {accession      => $orphanet, 
                                          mapping_source => 'Orphanet', 
                                          mapping_type   => 'is'}
                                        );   
    $pheno_adaptor->store_ontology_accessions($omim_pheno);
  }
}


## take OMIM phenotype/disease descriptions and xrefs from HP annotation file
sub get_HP {

  my $infile = shift;

  my %phenos;

  open my $link, $infile || die "Failed to open annotations file $infile : $!\n";
  while(<$link>){
    next unless /HP:/;
    next if /MOVED TO/;

    chomp;
    my @a = split/\t/;

    ## remove OMIM id
    $a[2] =~ s/^\#\d+\s+|^\%\d+\s+//g;

    $a[2] =~ s/\'/\\'/g;
 
    ## divide synonyms
    my @d = split/\;/,$a[2];
    foreach my $d (@d){
      push @{$phenos{$d}},  $a[4];
    }
  }
  return \%phenos;
}



sub store_HP{

  my $phenos = shift;

  foreach my $desc (keys %$phenos){
    my $pheno;
    eval{
       $pheno = $pheno_adaptor->fetch_by_description($desc);
    };
    print "Error with $desc : $@\n" unless $@ eq '';

    next unless $pheno->[0];  ## we don't expect to hold them all
  
    print "$desc\t" . join(",",@{$phenos->{$desc}}) ."\t" . $pheno->[0]->description() . "\n";
  
    foreach my $acc(@{$phenos->{$desc}}){
      $pheno->[0]->add_ontology_accession( {accession      => $acc, 
                                            mapping_source => 'HPO',
                                            mapping_type   => 'involves'}
                                         );
    }
    $pheno_adaptor->store_ontology_accessions($pheno->[0]);
  }
}



sub usage{

  die "\n\n\timport_ontology_xref.pl -registry [registry file] -type [ Orphanet or HP]

                   -data_file [annotations file] required for HP\n\n";
}
