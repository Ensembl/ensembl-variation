=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut


=head1 OntologyMapping

This module add ontology accessions for phenotype descriptions present in the database.

  First pass: Seek OLS exact matches for all descriptions
  Second pass: If no mapping, trim descriptions and seek OLS exact matches
  Third pass: For homo_sapiens only: If no mapping, seek Zooma mappings curated by EVA for all descriptions
  Supported ontologies: EFO, Orphanet, ORDO, HPO, VT

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::OntologyMapping;

use strict;
use warnings;

use File::Path qw(make_path);
use HTTP::Tiny;
use JSON;
use Data::Dumper;

use base qw(Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation);

my $DEFAULT_SPECIES = 'homo_sapiens';

sub fetch_input {
  my $self = shift;

  my $species = $self->required_param('species');
  my $workdir = $self->param('workdir');
  $workdir ||= $self->required_param('pipeline_dir')."/OntologyMap";
  unless (-d $workdir) {
    my $err;
    make_path($workdir, {error => \$err});
    die "make_path failed: ".Dumper($err) if $err && @$err;
  }

  $self->debug($self->param('debug_mode'));
  $self->workdir($workdir);

  open(my $logFH, ">", $workdir."/".'log_import_out_ontologyMapping_'.$species) || die ("Could not open file for writing: $!\n");
  open(my $errFH, ">", $workdir."/".'log_import_err_ontologyMapping_'.$species) || die ("Failed to open file: $!\n");
  open(my $pipelogFH, ">", $workdir."/".'log_import_debug_pipe_ontologyMapping_'.$species) || die ("Failed to open file: $!\n");
  $self->logFH($logFH);
  $self->errFH($errFH);
  $self->pipelogFH($pipelogFH);

}

sub run {
  my $self = shift;

  ## add exact matches where available for all phenotypes
  my $all_phenos = get_all_phenos($self->variation_db_adaptor->dbc->db_handle );

  my $ols_terms = $self->add_ols_matches($all_phenos);
  $self->store_terms( $ols_terms );

  ## seek matches for parent descriptions missing terms
  ## eg for 'Psoriasis 13' seek 'Psoriasis'
  my $non_matched_phenos = get_termless_phenos($self->variation_db_adaptor->dbc->db_handle);
  my $ols_parent_terms   = $self->add_ols_matches($non_matched_phenos, 'parent');
  $self->store_terms( $ols_parent_terms );

  if ($self->param('species') eq $DEFAULT_SPECIES) {
    my $non_matched_phenos = get_termless_phenos($self->variation_db_adaptor->dbc->db_handle);
    my $zooma_terms  = $self->add_zooma_matches($non_matched_phenos);
    $self->store_terms( $zooma_terms);
  }

  $self->param('output_ids', [{species => $self->param('species')}]);

}

sub write_output {
  my $self = shift;

  $self->print_pipelogFH("Passing ".$self->param('species').") for summary counts (finish_phenotype_annotation)\n") if ($self->debug);
  close($self->logFH) if defined $self->logFH ;
  close($self->errFH) if defined $self->errFH ;
  close($self->pipelogFH) if defined $self->pipelogFH ;

  $self->dataflow_output_id($self->param('output_ids'), 2);

}

sub store_terms{
  my ($self, $terms)   = @_;

  my $species = $self->param('species');
  my $phenotype_dba  = $self->variation_db_adaptor->get_PhenotypeAdaptor;

  foreach my $id (keys %{$terms}){

    my $pheno = $phenotype_dba->fetch_by_dbID( $id );
    die "Not in db: $id\n" unless defined $pheno; ## this should not happen
    foreach my $accession (@{$terms->{$id}->{terms}}){
      if ($species eq $DEFAULT_SPECIES) {
        next if $accession =~ /UBERON|NCBITaxon|NCIT|CHEBI|PR|MPATH|MA|PATO/; ## filter out some EFO term types
      }
      $self->print_logFH("$id\t$accession\t$terms->{$id}->{type}\t" . $pheno->description() ."\n");

      $pheno->add_ontology_accession({ accession      => $accession, 
                                       mapping_source => $terms->{$id}->{type},
                                       mapping_type   => 'is'
                                      });
    }
    $phenotype_dba->store_ontology_accessions($pheno);
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


=head2 get_eva_zooma_terms

return only eva curated matches for a phenotype description

=cut

sub get_eva_zooma_terms{
  my ($self, $desc, $confLevel) = @_;

  my $ontol_data = get_zooma_terms($desc);
  return undef unless defined $ontol_data;

  my @terms;

  foreach my $annot(@{$ontol_data}){
    next unless $annot->{derivedFrom}->{provenance}->{source}->{name} eq 'eva-clinvar';

    if (defined $confLevel) {
      next unless $annot->{confidence} eq $confLevel;
    }
    next unless grep(/EFO|Orphanet|ORDO|HP/, @{$annot->{semanticTags}});

    foreach my $term (@{$annot->{semanticTags}} ){
      $self->print_logFH("$term\t$annot->{confidence}\t$desc\tEVA\n");
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

  $pheno =~ s/\s+/+/g;

  my $http = HTTP::Tiny->new();
  my $server = 'https://www.ebi.ac.uk/ols/api/search?q=';
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
  my $server = 'https://www.ebi.ac.uk/spot/zooma/v2/api/services/annotate?propertyValue=';
  my $request  = $server . $pheno  ;

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


=head2 add_OLS_matches

 find exact match terms from supported ontologies

=cut

sub add_ols_matches{
  my ($self, $phenos, $truncate) = @_;

  my $species = $self->param('species');
  my %terms;

  foreach my $id (keys %{$phenos}){
    my $search_term = $phenos->{$id};

    ## modify term if no exact match available
    if (defined $truncate && $truncate eq 'parent'){

      $search_term = (split/\,/, $search_term)[0] if $search_term =~ /\,/;
      $search_term =~ s/\s*\((\w+\s*)+\)\s*//;         ## (one family) or (DDG2P abbreviation)
      $search_term =~ s/SUSCEPTIBILITY TO\s*|SUSCEPTIBILITY//i;
      $search_term =~ s/(\,+(\s*\w+\s*)+)+\s*$//;      ## remove ", with other  condition" ", one family" type qualifiers
      $search_term =~ s/\s+type\s*\w*\s*$//i;          ## remove " type 34"
      $search_term =~ s/\s+\d+\s*$//;                  ## remove just numeric subtype
      $search_term =~ s/\s+\d+\w\s*$//;                ## remove numeric+letter subtype 1c
      $search_term =~ s/primary\s*//i;
      $search_term =~ s/\s*\,\s*$//; 
      $search_term =~ s/\s+/ /g;
      $search_term =~ s/nonsyndromic //i;

      next if length($search_term) == length($phenos->{$id});
      $self->print_errFH("Seeking $search_term from $phenos->{$id}\n");
    }

    my $ontol_data = get_ols_terms( $search_term );
    next unless defined $ontol_data;

    my @terms;
    foreach my $doc (@{$ontol_data->{response}->{docs}}){ 
      if ($species eq $DEFAULT_SPECIES) {
        next unless $doc->{ontology_prefix} =~/EFO|Orphanet|ORDO|HP/;
      } elsif($species ne $DEFAULT_SPECIES) {
        next unless $doc->{ontology_prefix} =~/VT/;
      }

      push @terms, iri2acc($doc->{iri});
    }
    next unless scalar(@terms) > 0;
    $terms{$id}{terms} = \@terms;
    $terms{$id}{type} = (defined $truncate ? 'OLS partial' : 'OLS exact');
  }

  return \%terms;
}


=head2 add_zooma_matches

  use EVA curated mappings
=cut

sub add_zooma_matches{
  my ($self, $phenos) = @_;

  my %terms;

  foreach my $id (keys %{$phenos}){
    my $desc = $phenos->{$id};
    $desc =~ s/\s+\(\w+\)$//; ## remove DDG2P abbreviation 

    my $eva_terms = $self->get_eva_zooma_terms($desc, 'HIGH');
    $terms{$id}{terms} = $eva_terms if defined $eva_terms;
    $terms{$id}{type} = "Zooma exact" if defined $eva_terms;
  }

  return \%terms;
}

=head2 iri2acc

  extract accession number from iri identifiers (https://cloud.identifiers.org)

=cut

sub iri2acc{
  my $iri = shift;

  my @a = split/\//, $iri;
  my $acc = pop @a;
  $acc =~ s/\_/\:/;

  return $acc;
}

1;

