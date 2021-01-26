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


=head1 MouseBasePhenotypeAnnotation

General mouse phenotype import methods.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::MouseBasePhenotypeAnnotation;

use strict;
use warnings;
use HTTP::Tiny;
use JSON;

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

sub get_input_file {
  my ($self, $working_dir, $data_source, $url) = @_;

  my $phenotype_file = "$working_dir/$data_source\_phenotypes.txt";

  open(OUT, ">".$phenotype_file) || die ("Could not open $phenotype_file for writing\n");
  my $http = HTTP::Tiny->new();

  my $i      = 0;
  my $rows   = 0;
  my $server = 'https://www.ebi.ac.uk';
  my $ext     = "$url/select?q=*:*&rows=$i&wt=json";

  # get total number of rows
  while ($i == $rows) {
    $i += 1000;
    $ext = "$url/select?q=*:*&rows=$i&wt=json";
    my $response = $http->get($server.$ext, {
      headers => { 'Content-type' => 'application/json' }
    });
    my $hash = decode_json($response->{content});
    $rows = scalar @{$hash->{response}->{docs}};
  }
  my $start = 0;
  while ($start <= $rows) {
    $ext = "$url/select?q=*:*&rows=100&wt=json&start=$start";
    my $response = $http->get($server.$ext, {
      headers => { 'Content-type' => 'application/json' }
    });

    die "Failed!\n" unless $response->{success};
    if (length $response->{content}) {
      my $hash = decode_json($response->{content});
      foreach my $result (@{$hash->{response}->{docs}}) {
        my @pairs = ();
        foreach my $key (keys %$result) {
          my $value = $result->{$key};
          $value ||= '\N';
          push @pairs, $key . "=" . $value;
        }
        print OUT join("\t", @pairs), "\n";
      }
    }
    $start += 100;
  }
  close OUT;
  return $data_source."\_phenotypes.txt";
}

sub fetch_data_version {
  my $self = $_;

  my $http = HTTP::Tiny->new();

  my $url = 'https://www.mousephenotype.org/data/release.json';
  my $response = $http->get($url, {
    headers => { 'Content-type' => 'application/json' }
  });
  my $hash = decode_json($response->{content});
  my $data_release_date = $hash->{data_release_date};
  die "data_release_date not defined\n" unless $data_release_date;

 # e.g. 20 October 2015
  my ($date, $month, $year) = split(/[\s-]/, $hash->{data_release_date}) ;

  my $months = {
    'January' => '01',
    'February' => '02',
    'March' => '03',
    'April' => '04',
    'May' => '05',
    'June' => '06',
    'July' => '07',
    'August' => '08',
    'September' => '09',
    'October' => '10',
    'November' => '11',
    'December' => '12',
  };

  my $month_number = $months->{$month};
  die "month_number not defined\n" unless $month_number;

  my $version = sprintf '%d%02d%02d', $year, $month_number, $date;

  return $version;
}

sub get_marker_coords {
  my ($self, $phenotype_file, $coord_file) = @_;

  my $gene_adaptor = $self->core_db_adaptor->get_GeneAdaptor;

  my $markers;
  my $column_headers;
  # filter for marker accession ids (MGI or IMPC)
  open(IN, "<".$phenotype_file) || die ("Could not open $phenotype_file for reading\n");
  while (<IN>) {
    chomp;
    my @pairs = split("\t", $_);
    my $hash;
    foreach my $pair (@pairs) {
      my ($key, $value) = split("=", $pair);
      $hash->{$key} = $value;
      $column_headers->{$key} = 1;
    }
    $markers->{$hash->{marker_accession_id}} = 1;
  }
  close IN;

  open(IN, "<".$coord_file) || die ("Could not open $coord_file for reading\n");
  my $header;
  my $marker_coords = {};
  while (<IN>) {
    chomp;
    # test column names are still the same
    if (/^1\./) {
      my @values = split("\t", $_);
      foreach my $value (@values) {
        my ($column_number, $column_name) = split(/\s/, $value, 2);
        $header->{$column_number} = $column_name;
      }
      unless ($header->{'1.'} eq 'MGI Marker Accession ID') { die 'Header in MGI_MRK_Coord.rpt has changed.'};
      unless ($header->{'6.'} eq 'Chromosome') { die 'Header in MGI_MRK_Coord.rpt has changed.'};
      unless ($header->{'7.'} eq 'Start Coordinate') { die 'Header in MGI_MRK_Coord.rpt has changed.'};
      unless ($header->{'8.'} eq 'End Coordinate') { die 'Header in MGI_MRK_Coord.rpt has changed.'};
      unless ($header->{'9.'} eq 'Strand') { die 'Header in MGI_MRK_Coord.rpt has changed.'};
    }
   my @values = split("\t", $_);
    my $marker_acc = $values[0];
    if ($markers->{$marker_acc}) {
      unless ($marker_coords->{$marker_acc}) {
        my $genes = $gene_adaptor->fetch_all_by_external_name($marker_acc);
        if (scalar @$genes != 1) {
          my $number_of_genes = scalar @$genes;
          $self->print_errFH("WARNING: Found $number_of_genes matching Ensembl genes for gene $marker_acc\n");
          my $unique_seq_region_id = {};
          foreach my $gene (@$genes) {
            $unique_seq_region_id->{$gene->slice->get_seq_region_id} = 1;
          }
          $self->print_errFH("WARNING: genes $number_of_genes are on different slices\n") unless (scalar keys %$unique_seq_region_id == 1);
        }
        next unless scalar @$genes;

        my $marker_type = $values[1];
        $self->print_errFH("WARNING: type is not Gene (type = $marker_type)\n") unless ($marker_type eq 'Gene');

        my $chromosome = $values[5];
        my $start_coord = $values[6];
        my $end_coord = $values[7];
        my $strand = $values[8];

        $marker_coords->{$marker_acc}->{type} = $marker_type;
        $marker_coords->{$marker_acc}->{object_id} = join(";", map {$_->stable_id} @$genes);
        $marker_coords->{$marker_acc}->{chromosome} = $chromosome;
        $marker_coords->{$marker_acc}->{seq_region_id} = $genes->[0]->slice->get_seq_region_id;
        $marker_coords->{$marker_acc}->{seq_region_start} = $start_coord;
        $marker_coords->{$marker_acc}->{seq_region_end} = $end_coord;
        $marker_coords->{$marker_acc}->{seq_region_strand} = ($strand eq '+') ? 1 : -1;
      }
    }
  }

  close IN;
  return $marker_coords;
}

sub parse_input_file {
  my ($self, $infile, $marker_coords, $source_name, $source_id) = @_;

  my $individual_adaptor = $self->variation_db_adaptor->get_IndividualAdaptor;

  my @phenotypes = ();
  my $already_inserted = {};

  my @data_attrib = (
    'allele_symbol',
    'allele_accession_id',
    'marker_accession_id',);

  open(IN, "<".$infile) || die ("Could not open $infile for reading\n");
  while (<IN>) {
    chomp;
    my %data = ();
    my @pairs = split("\t", $_);
    my $hash;
    foreach my $pair (@pairs) {
      my ($key, $value) = split("=", $pair);
      if ($key eq 'allele_symbol') { # change e.g. Cdk5rap2<tm1a(EUCOMM)Wtsi> to Cdk5rap2_tm1a(EUCOMM)Wtsi
        $value =~ s/</_/g;
        $value =~ s/>//g;
      }
      $hash->{$key} = $value;
    }

    foreach my $value (@data_attrib) {
      $data{$value} = $hash->{$value};
    }
    my $marker_accession_id = $hash->{marker_accession_id};
    my $description = $hash->{mp_term_name};
    if (!$description) {
      $self->print_errFH("WARNING: No description/mp_term_name for $_\n");
      next;
    }
    $data{description} = $description;
    $data{accession} = $hash->{mp_term_id};

    my $source = '';
    if ($source_name eq 'IMPC') {
      $source = $hash->{resource_name};
    } else {
      $source = $hash->{project_name};
    }
    $data{source} = $source;
    $data{source_id} = $source_id;

    my $pf_data           = $marker_coords->{$marker_accession_id};
    my $type              = $pf_data->{type};
    my $seq_region_id     = $pf_data->{seq_region_id};
    my $seq_region_start  = $pf_data->{seq_region_start};
    my $seq_region_end    = $pf_data->{seq_region_end};
    my $seq_region_strand = $pf_data->{seq_region_strand};

    $data{seq_region_id} = $seq_region_id;
    $data{seq_region_start} = $seq_region_start;
    $data{seq_region_end} = $seq_region_end;
    $data{seq_region_strand} = $seq_region_strand; 

    my $strain_name       = $hash->{strain_name};
    my $gender            = $hash->{sex};
    if ($gender) {
      unless ($gender eq 'female' || $gender eq 'male') {
        $gender = 'unknown';
      }
    }

    my $strain_id = 0;
    my $create_new_individual = 0;
    if ($strain_name) {
      # attribs
      my $p_value;
      my $strains = $individual_adaptor->fetch_all_by_name($strain_name);
      if (scalar @$strains > 0) {
        my @unique = ();
        foreach my $strain (@$strains) {
          if (lc $strain->gender eq $gender) {
            push @unique, $strain;
          }
        }
        if (scalar @unique > 0) {
          $strain_id = $unique[0]->dbID();
        } else {
          $create_new_individual = 1;
        }
      }
      if ($create_new_individual) {
        my $individual = new Bio::EnsEMBL::Variation::Individual(
            -name => $strain_name,
            -gender => $gender,
            -type_individual => 'fully_inbred',
        );
        $individual_adaptor->store($individual);
        $strain_id = $individual_adaptor->last_insert_id();
      }
    }
    $data{associated_gene} = $hash->{marker_symbol};

    if ($hash->{p_value}) {
      $data{p_value} = convert_p_value($hash->{p_value});
    }
    $data{external_id} = $hash->{mp_term_id};

    if ($strain_id) {
      $data{strain_id} = $strain_id;
    }

    my @object_ids = ();
    if ($pf_data->{object_id}) {
      @object_ids = split(';', $pf_data->{object_id});
      foreach my $object_id (@object_ids) {
        my $accession = $data{accession};
        my $source = $data{source};
        my $strain_id = $data{strain_id};
        my $key = join('-', ($accession, $source, $object_id, $strain_id || ''));
        unless ($already_inserted->{$key}) {
          my %phenotype_data = ();
          foreach my $key (keys %data) {
            $phenotype_data{$key} = $data{$key};
          }
          $phenotype_data{id} = $object_id;   
          $phenotype_data{accessions} = [$accession];
          $phenotype_data{ontology_mapping_type} = 'is'; 
          push(@phenotypes, \%phenotype_data);
        }
        $already_inserted->{$key} = 1;
      }
    }
  }
  close IN;
  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
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
    elsif ($pval =~ /(\d+.*)/){
      $sci_pval = $1; #$sci_pval = sprintf("%.2e",$pval); # e.g. 0.002 => 2,30e-3
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
    $sci_pval = $pval;
  }
  return $sci_pval;
}

1;
