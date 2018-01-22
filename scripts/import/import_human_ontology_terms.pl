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


=head1

Seek ontology accessions for phenotype descriptions

  First pass:  Seek Zooma mappings curated by EVA for all descriptions
  Second pass: Seek OLS exact matches for all descriptions
  Third pass:  If no mapping, trim descriptions and seek OLS exact matches

  Supported ontologies: EFO, ORDO, HPO

=cut


use strict;
use warnings;
use HTTP::Tiny;
use JSON;
use Getopt::Long;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::PhenotypeAdaptor;

my ($registry_file);
GetOptions ( "registry=s"   => \$registry_file  );

die "Error registry file needed\n" unless defined $registry_file;

open my $out, ">import_phenotype_accessions.log" ||die "Failed to open log file :$!\n";

my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all($registry_file);
my $dba = $reg->get_DBAdaptor('homo_sapiens','variation');



## add exact matches where available for all phenotypes
my $all_phenos = get_all_phenos($dba->dbc->db_handle );

my $zooma_terms  = add_zooma_matches($all_phenos);
store_terms($reg, $zooma_terms);

my $ols_terms = add_ols_matches($all_phenos);
store_terms($reg, $ols_terms );



## seek matches for parent descriptions missing terms
## eg for 'Psoriasis 13' seek 'Psoriasis'
my $non_matched_phenos = get_termless_phenos($dba->dbc->db_handle);
my $ols_parent_terms   = add_ols_matches($non_matched_phenos, 'parent');
store_terms($reg, $ols_parent_terms );


=head2 add_OLS_matches

 find exact match terms from supported ontologies

=cut

sub add_ols_matches{

  my $phenos   = shift;
  my $truncate = shift;

  my %terms;

  foreach my $id (keys %{$phenos}){
    my $search_term = $phenos->{$id};

    ## modify term if no exact match available
    if (defined $truncate && $truncate eq 'parent'){

      $search_term = (split/\,/, $search_term)[0] if $search_term =~ /\,/;
      $search_term =~ s/\s*\((\w+\s*)+\)\s*//;         ## (one family) or (DDG2P abbreviation)
      $search_term =~ s/SUSCEPTIBILITY TO\s*|SUSCEPTIBILITY//i;
      $search_term =~ s/(\,+(\s*\w+\s*)+)+\s*$//;      ## remove ", with other  condition" ", one family" type qualifiers
      $search_term =~ s/\s+type\s+\d+\s*$//i;          ## remove " type 34"
      $search_term =~ s/\s+\d+\s*$//;                  ## remove just numeric subtype
      $search_term =~ s/\s+\d+\w\s*$//;                ## remove numeric+letter subtype 1c
      $search_term =~ s/primary\s*//i;
      $search_term =~ s/\s*\,\s*$//; 
      $search_term =~ s/\s+/ /g;
      $search_term =~ s/nonsyndromic //i;

      next if length($search_term) == length($phenos->{$id});
      print "Seeking $search_term from $phenos->{$id}\n";
    }

    my $ontol_data = get_ols_terms( $search_term );
    next unless defined $ontol_data;

    my @terms;
    foreach my $doc (@{$ontol_data->{response}->{docs}}){ 
      next unless $doc->{ontology_prefix} =~/EFO|Orphanet|ORDO|HP/;

      push @terms, iri2acc($doc->{iri});
    }
    $terms{$id}{terms} = \@terms;
    $terms{$id}{type} = (defined $truncate ? 'OLS partial' : 'OLS exact');
  }
  return \%terms;
}

=head2 add_zooma_matches

  use EVA curated mappings
=cut
sub add_zooma_matches{
  
  my $phenos = shift;

  my %terms;

  foreach my $id (keys %{$phenos}){
    my $desc = $phenos->{$id};
    $desc =~ s/\s+\(\w+\)$//; ## remove DDG2P abbreviation 

    my $eva_terms = get_eva_zooma_terms($desc);   
    $terms{$id}{terms} = $eva_terms if defined $eva_terms;
    $terms{$id}{type} = "Zooma exact";
  }
  return \%terms;
}


sub store_terms{

  my $reg   = shift;
  my $terms = shift;

  my $pheno_adaptor = $reg->get_adaptor('homo_sapiens','variation', 'Phenotype');

  foreach my $id (keys %{$terms}){

    my $pheno = $pheno_adaptor->fetch_by_dbID( $id );
    die "Not in db: $id\n" unless defined $pheno; ## this should not happen
    foreach my $accession (@{$terms->{$id}->{terms}}){
      next if $accession =~ /UBERON|NCBITaxon|NCIT|CHEBI|PR|MPATH|MA|PATO/; ## filter out some EFO term types 
      print $out "$id\t$accession\t$terms->{$id}->{type}\t" . $pheno->description() ."\n";

      $pheno->add_ontology_accession({ accession      => $accession, 
                                       mapping_source => $terms->{$id}->{type},
                                       mapping_type   => 'is'
                                      });
    }
    $pheno_adaptor->store_ontology_accessions($pheno);
  }
}



=head2 get_all_phenos

look up descriptions for all phenotypes to add exact matching terms
from the supported ontologies.

=cut
sub get_all_phenos{

  my $dbh = shift;

  my $desc_ext_sth = $dbh->prepare(qq [ select phenotype_id, description from phenotype ]);

  $desc_ext_sth->execute()||die "Problem extracting all phenotype descriptions\n";
  my $data = $desc_ext_sth->fetchall_arrayref();
  my %pheno;

  foreach my $l (@{$data}){
   $pheno{$l->[0]} = $l->[1];
  }

  return \%pheno;
}

=head2 get_termless_phenos

look up descriptions for phenotypes which don't have supplied or
exact match terms of 'is' type

=cut
sub get_termless_phenos{

  my $dbh = shift;

  my $desc_ext_sth =   $dbh->prepare(qq [ select phenotype_id, description
                                          from phenotype
                                          where phenotype.phenotype_id not in 
                                          (select phenotype_id from phenotype_ontology_accession where mapping_type ='is')  
                                        ]);

  $desc_ext_sth->execute()||die "Problem extracting termless phenotype descriptions\n";
  my $data = $desc_ext_sth->fetchall_arrayref();
  my %pheno;
  foreach my $l (@{$data}){
   $pheno{$l->[0]} = $l->[1];
  }

  return \%pheno;
}

=head2 get_high_quality_zooma_terms

return only high quality matches for a phenotype description

=cut

sub get_high_quality_zooma_terms{

  my $desc = shift;

  my $ontol_data = get_zooma_terms($desc);
  return undef unless defined $ontol_data;

  my @terms;

  foreach my $annot(@{$ontol_data}){
    next unless $annot->{confidence} eq 'HIGH';
    foreach my $term (@{$annot->{semanticTags}} ){
      next unless $term =~ /EFO|Orphanet|ORDO|HP/;       
      print $out "$term\t$annot->{confidence}\t$desc\n";
      push @terms, iri2acc($term) ;
    }
  }
  @terms ? return \@terms : return undef;
}


=head2 get_eva_zooma_terms

return only eva curated matches for a phenotype description

=cut

sub get_eva_zooma_terms{

  my $desc = shift;

  my $ontol_data = get_zooma_terms($desc);
  return undef unless defined $ontol_data;

  my @terms;

  foreach my $annot(@{$ontol_data}){
    next unless $annot->{derivedFrom}->{provenance}->{source}->{name} eq 'eva-clinvar'; 

    next unless $annot->{ontology_prefix} =~/EFO|Orphanet|ORDO|HP/;

    foreach my $term (@{$annot->{semanticTags}} ){
      print $out "$term\t$annot->{confidence}\t$desc\tEVA\n";
      push @terms, iri2acc($term) ;
    }
  }
  return undef unless defined $terms[0];
  return \@terms;
}


=head2 get_ols_terms

look up phenotype description in OLS seeking exact match

=cut

sub get_ols_terms{

  my $pheno = shift;

  $pheno =~ s/\s+/\%20/g;

  my $http = HTTP::Tiny->new();
  my $server = 'http://www.ebi.ac.uk/ols/api/search?q=';
  my $request  = $server . $pheno . "&queryFields=label,synonym&exact=1";


  my $response = $http->get($request, {
        headers => { 'Content-type' => 'application/xml' }
                              });
  unless ($response->{success}){
        warn "Failed request: $request :$!\n" ;
        return;
  }

  return JSON->new->decode($response->{content} );
}

=head2 get_zooma_terms

look up phenotype description in zooma

=cut
sub get_zooma_terms{

  my $pheno = shift;

  $pheno =~ s/\(\w+\)//; ## for DDG2P
  $pheno =~ s/\s+/\+/g;

  my $http = HTTP::Tiny->new();
  my $server = 'http://www.ebi.ac.uk/spot/zooma/v2/api/services/annotate?propertyValue=';
  my $request  = $server . '"'. $pheno .'"' ;

  #warn "Looking for $request\n\n" ;
  my $response = $http->get($request, {
      headers => { 'Content-type' => 'application/xml' }
                     });
  unless ($response->{success}){
    warn "Failed request: $request :$!\n" ;
    return;
  }

 return JSON->new->decode($response->{content} );

}

sub iri2acc{

  my $iri = shift;
  my @a = split/\//, $iri;
  my $acc = pop @a;
  $acc =~ s/\_/\:/;

  return $acc;
}
