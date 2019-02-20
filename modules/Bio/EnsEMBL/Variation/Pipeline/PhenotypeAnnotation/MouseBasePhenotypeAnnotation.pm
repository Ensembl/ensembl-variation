=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::MouseBasePhenotypeAnnotation;

use strict;
use warnings;
use HTTP::Tiny;
use JSON;

#TODO: check that all libraries from import_script are imported. Any missing? Were they actually used?
use Data::Dumper;

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

sub get_mouse_phenotype_data {
  my $self = shift;
  my $working_dir = shift;
  my $data_source = shift;
  my $url = shift;
  my $phenotype_file = "$working_dir/$data_source\_phenotypes.txt";
  open(IN, ">".$phenotype_file) or die ("Could not open $phenotype_file for reading");
  my $http = HTTP::Tiny->new();

  my $i      = 0;
  my $rows   = 0;
  my $server = 'http://www.ebi.ac.uk';
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
  return $phenotype_file;
}

sub get_mouse_phenotype_data_source_ids {
  my $self = shift;
  my $phenotype_file = shift;
  my $source_info = shift;
  my $variation_dba = shift;

  my $data_source = $source_info->{source};
  my $version = $source_info->{source_version};
  
  my $dbh = $variation_dba->dbc->db_handle;
  open(IN, "<".$phenotype_file) or die ("Could not open $phenotype_file for reading");
  my $source_names = {};
  my $source_name2id = {};

  while (<IN>) {
    chomp;
    my @pairs = split("\t", $_);
    my $hash;
    foreach my $pair (@pairs) {
      my ($key, $value) = split('=', $pair);
      if ($key eq 'resource_name' && $data_source eq 'impc') {
        $source_names->{$value} = 1;
      }
      if ($key eq 'project_name' && $data_source eq 'mgi') {
        $source_names->{$value} = 1;
      }
    }
  }
  close IN;
  foreach my $source_name (keys %$source_names) {
    my $stmt = qq{ SELECT source_id FROM source WHERE name = '$source_name' LIMIT 1};
    my $sth = $dbh->prepare($stmt);
    $sth->execute();
    my $source_id;
    $sth->bind_columns(\$source_id);
    $sth->fetch();
    if (defined($source_id)) {
      $source_name2id->{$source_name} = $source_id;
    } else {
      warn "Could not fetch dbID for source name $source_name\n";
    }
    $sth->finish();
  }

  return $source_name2id;
}

sub update_mouse_phenotype_data_version {
  my $self = shift;
  my $source_name2ids = shift;
  my $variation_db = shift;

  my $http = HTTP::Tiny->new();

  my $url = 'http://www.mousephenotype.org/data/release.json';
  my $response = $http->get($url, {
    headers => { 'Content-type' => 'application/json' }
  });
  my $hash = decode_json($response->{content});
  my $data_release_date = $hash->{data_release_date};
  die "data_release_date not defined\n" unless $data_release_date;

 # e.g. 20 October 2015
  my ($date, $month, $year) = split(' ', $hash->{data_release_date}) ;

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

  my $dbh = $variation_db->dbc->db_handle;

  foreach my $source_name (keys %$source_name2ids) {
    $dbh->do(qq{UPDATE source SET version=$version WHERE name='$source_name';}) or die $dbh->errstr;
  }
  return $version;
}

sub clear_mouse_phenotype_data_from_last_release {
  my ($self,$source_name2id,$variation_db)  =  @_;

  my $dbh = $variation_db->dbc->db_handle;
  my $source_ids = join(',', values %$source_name2id);
  print $source_ids;
  $dbh->do(qq{ DELETE pfa FROM phenotype_feature_attrib pfa JOIN phenotype_feature pf ON pfa.phenotype_feature_id = pf.phenotype_feature_id AND pf.source_id IN ($source_ids);} );
  $dbh->do(qq{ DELETE FROM phenotype_feature WHERE source_id IN ($source_ids);} );
}

sub get_marker_coords {
  my ($self, $phenotype_file, $coord_file, $core_dba) = @_;

  my $gene_adaptor = $core_dba->get_GeneAdaptor;

  my $markers;
  my $column_headers;
  # filter for marker accession ids (MGI)
  open(IN, "<".$phenotype_file) or die ("Could not open $phenotype_file for reading");
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

  open(IN, "<".$coord_file) or die ("Could not open $coord_file for reading");
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
          print STDERR "WARNING: Found $number_of_genes matching Ensembl genes for gene $marker_acc\n";
          my $unique_seq_region_id = {};
          foreach my $gene (@$genes) {
            $unique_seq_region_id->{$gene->slice->get_seq_region_id} = 1;
          }
          print STDERR "WARNING: genes $number_of_genes are on different slices\n" unless (scalar keys %$unique_seq_region_id == 1);
        }
        next unless scalar @$genes;

        my $marker_type = $values[1];
        print STDERR "WARNING: type is not Gene (type = $marker_type)\n" unless ($marker_type eq 'Gene');

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

sub parse_mouse_phenotype_data {
  my $self = shift;
  my $infile = shift;
  my $marker_coords = shift;
  my $data_source = shift;
  my $source_name2ids = shift;
  my $variaiton_dba = shift;

  my $individual_adaptor = $variaiton_dba->get_IndividualAdaptor;

  my @phenotypes = ();
  my $already_inserted = {};

  my @data_attrib = (
    'allele_symbol',
    'allele_accession_id',
    'marker_accession_id',);

  open(IN, "<".$infile) or die ("Could not open $infile for reading");
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
      warn "No description/mp_term_name for $_\n";
      next;
    }
    $data{description} = $description;
    $data{accession} = $hash->{mp_term_id};

    my $source = '';
    if ($data_source eq 'impc') {
      $source = $hash->{resource_name};
    } else {
      $source = $hash->{project_name};
    }
    $data{source} = $source;
    $data{source_id} = $source_name2ids->{$source};

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
            -individual_type_id => 1,
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
        my $key = join('-', ($accession, $source, $object_id, $strain_id));
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
