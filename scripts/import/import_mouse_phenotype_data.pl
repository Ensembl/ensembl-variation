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

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use DBI qw(:sql_types);
use FileHandle;
use Getopt::Long;
use HTTP::Tiny;
use JSON;

usage() if (!scalar(@ARGV));

my $config = {};

GetOptions(
  $config,
  'working_dir=s',
  'coord_file=s',
  'phenotype_file=s',
  'registry=s',
  'host=s',
  'dbname=s',
  'user=s',
  'pass=s',
  'port=i',
  'chost=s',
  'cdbname=s',
  'cuser=s',
  'cpass=s',
  'cport=i',
  'help!',
) or die "Error: Failed to parse command line arguments\n";

usage() if ($config->{help});

die ("MGI_MRK_Coord.rpt file (--coord_file) is required") unless (defined($config->{coord_file}));
die ("A working directory (--working_dir) is required") unless (defined($config->{working_dir}));
my $variation_credentials = defined($config->{host}) && defined($config->{port}) && defined($config->{dbname}) && defined($config->{user}) && defined($config->{pass});
my $core_credentials      = defined($config->{chost}) && defined($config->{port}) && defined($config->{cdbname});
my $registry = defined($config->{registry});
die ("Database credentials or a registry file are required (try --help)") unless (($variation_credentials && $core_credentials) || $registry);
die ("Variation database credentials (--host, --dbname, --user, --pass, --port) are required") if (!$variation_credentials && !$registry);
die ("Core database credentials (--chost, --cdbname, --cport) are required") if (!$core_credentials && !$registry);

#http://www.informatics.jax.org/marker/MGI:101877
$config->{urls} = {
  impc => '/mi/impc/solr/genotype-phenotype',
  mgi => '/mi/impc/solr/mgi-phenotype',
};

set_up_db_connections($config);
pre_cleanup($config);
get_version($config);
foreach my $data_source (qw/mgi impc/) {
  get_data($config, $data_source);
  get_marker_coords($config, $data_source);
  clear_data_from_last_release($config, $data_source); # update source version
  populate_phenotype_table($config, $data_source);
  import_phenotype_features($config, $data_source);
}

sub set_up_db_connections {
  my $config = shift;
  my $registry = 'Bio::EnsEMBL::Registry';
  my $species = 'mouse';

  my ($gene_adaptor, $individual_adaptor, $phenotype_adaptor, $dbh);

  if ($config->{registry}) {
    $registry->load_all($config->{registry});
    $gene_adaptor       = $registry->get_adaptor($species, 'core', 'gene');
    $phenotype_adaptor  = $registry->get_adaptor($species, 'variation', 'phenotype');
    $individual_adaptor = $registry->get_adaptor($species, 'variation', 'individual');
    my $vdba            = $registry->get_DBAdaptor($species, 'variation');
    $dbh                = $vdba->dbc->db_handle;
  } else {
    # Connect to the variation database
    my $db_adaptor = new Bio::EnsEMBL::Variation::DBSQL::DBAdaptor(
        -host => $config->{host},
        -user => $config->{user},
        -pass => $config->{pass},
        -port => $config->{port},
        -dbname => $config->{dbname},
        ) or die("Could not get a database adaptor for $config->{dbname} on $config->{host}:$config->{port}");

    # connect to core DB for genes
    my $core_db_adaptor;

    $core_db_adaptor = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host => $config->{chost},
        -user => 'ensro',
        -pass => '',
        -port => $config->{cport},
        -dbname => $config->{cdbname},
        ) or die("Could not get a database adaptor for $config->{cdbname} on $config->{chost}:$config->{cport}");
    $gene_adaptor       = $core_db_adaptor->get_GeneAdaptor;
    $individual_adaptor = $db_adaptor->get_IndividualAdaptor;
    $phenotype_adaptor  = $db_adaptor->get_PhenotypeAdaptor;

    my $dbname = $config->{dbname};
    my $host   = $config->{host};
    my $port   = $config->{port};
    my $user   = $config->{user};
    my $pass   = $config->{pass};
    $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;port=$port;user=$user;password=$pass", {RaiseError => 1});
  }
  $config->{gene_adaptor}       = $gene_adaptor;
  $config->{phenotype_adaptor}  = $phenotype_adaptor;
  $config->{individual_adaptor} = $individual_adaptor;
  $config->{dbh} = $dbh;
}

sub pre_cleanup {
  my $config = shift;
  my $dbh = $config->{dbh};
  $dbh->do(qq{CREATE TABLE IF NOT EXISTS TMP_phenotype_feature LIKE phenotype_feature;}) or die $dbh->errstr;
  $dbh->do(qq{TRUNCATE TABLE TMP_phenotype_feature;}) or die $dbh->errstr;
  $dbh->do(qq{INSERT INTO TMP_phenotype_feature SELECT * FROM phenotype_feature;}) or die $dbh->errstr;

  $dbh->do(qq{CREATE TABLE IF NOT EXISTS TMP_phenotype_feature_attrib LIKE phenotype_feature_attrib;}) or die $dbh->errstr;
  $dbh->do(qq{TRUNCATE TABLE TMP_phenotype_feature_attrib;}) or die $dbh->errstr;
  $dbh->do(qq{INSERT INTO TMP_phenotype_feature_attrib SELECT * FROM phenotype_feature_attrib;}) or die $dbh->errstr;
}

sub clear_data_from_last_release {
  my $config = shift;
  my $data_source = shift;
  my $dbh = $config->{dbh};
  # create backup tables for: individual, phenotype_feature_attrib, phenotype, phenotype_feature --> not needed,
  # if something goes wrong, tables are available on ens-livemirror

  my $phenotype_file = $config->{"$data_source\_phenotype_file"};
  my $fh = FileHandle->new($phenotype_file, 'r');
  my $source_names = {};
  my $source_name2id = {};

  while (<$fh>) {
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
  $fh->close();
  foreach my $source_name (keys %$source_names) {  
    my $stmt = qq{ SELECT source_id FROM source WHERE name = '$source_name' LIMIT 1};
    my $sth = $dbh->prepare($stmt);
    $sth->execute();
    my $source_id;
    $sth->bind_columns(\$source_id);
    $sth->fetch();
    if (defined($source_id)) {
      $source_name2id->{$source_name} = $source_id;  
    }
    $sth->finish();
  }

  # update version
  my $version = $config->{version};
  foreach my $source_name (keys %$source_names) {
    $dbh->do(qq{UPDATE source SET version=$version WHERE name='$source_name';}) or die $dbh->errstr;
  }

  my $source_ids = join(',', values %$source_name2id);

  $dbh->do(qq{ DELETE pfa FROM phenotype_feature_attrib pfa JOIN phenotype_feature pf ON pfa.phenotype_feature_id = pf.phenotype_feature_id AND pf.source_id IN ($source_ids);} );
  $dbh->do(qq{ DELETE FROM phenotype_feature WHERE source_id IN ($source_ids);} );

}

sub get_version {
  my $config = shift;
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

  $config->{version} = $version;

}

sub get_data {
  my $config = shift;
  my $data_source = shift;

  my $working_dir = $config->{working_dir};

  my $url = $config->{urls}->{$data_source};
  my $phenotype_file = "$working_dir/$data_source\_phenotypes.txt"; 	
  $config->{"$data_source\_phenotype_file"} = $phenotype_file;
  my $fh = FileHandle->new($phenotype_file, 'w');
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
        print $fh join("\t", @pairs), "\n";
      }
    }
    $start += 100;
  }
  $fh->close();
}

sub populate_phenotype_table {
  my $config = shift;	
  my $data_source = shift;

  # filter mp_term_ids from phenotype_file
  my $phenotype_file = $config->{"$data_source\_phenotype_file"};

  my $fh = FileHandle->new($phenotype_file, 'r');
  my $mp_term_names = {};

  while (<$fh>) {
    chomp;
    my @pairs = split("\t", $_);
    my $hash;
    foreach my $pair (@pairs) {
      my ($key, $value) = split('=', $pair);
      if ($key eq 'mp_term_name') {
        $mp_term_names->{$value} = 1;
      }
    }
  }
  $fh->close();

  # populate phenotype table
  my $dbh = $config->{dbh};
  
  my $phenotype_descriptions = {};
  my $stmt = qq{ SELECT description FROM phenotype;};
  my $sth = $dbh->prepare($stmt);
  $sth->execute();
  my $description;
  $sth->bind_columns(\$description);
  while ($sth->fetch()) {
    $phenotype_descriptions->{$description} = 1;
  }
  $sth->finish();

  my $sql = "INSERT INTO phenotype (description) values (?)";
  $sth = $dbh->prepare($sql);
  foreach my $description (sort keys %$mp_term_names) {
    if (!$phenotype_descriptions->{$description}) {
      $sth->bind_param(1, $description, SQL_VARCHAR);
      $sth->execute();
    }
  }
  $sth->finish();
}

sub import_phenotype_features {
  my $config = shift;
  my $data_source = shift;

  my $dbh                = $config->{dbh};
  my $individual_adaptor = $config->{individual_adaptor};
  my $phenotype_adaptor  = $config->{phenotype_adaptor};
  my $marker_coords      = $config->{marker_coords};

  my @attrib_types = @{get_attrib_types($config)};

  my @data_attrib = (
      'allele_symbol', 
      'allele_accession_id', 
      'marker_accession_id',);

  my $pf_ins_sth = qq{
    INSERT INTO phenotype_feature (
        phenotype_id, source_id, type, object_id, is_significant, seq_region_id, seq_region_start, seq_region_end, seq_region_strand
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
  };

  my $attrib_ins_sth = qq{
    INSERT INTO phenotype_feature_attrib (
        phenotype_feature_id, attrib_type_id, value
        ) SELECT ?, at.attrib_type_id, ?
      FROM attrib_type at
      WHERE at.code = ?
  };

  my $attrib_ins_cast_sth = qq{
    INSERT INTO phenotype_feature_attrib (
        phenotype_feature_id, attrib_type_id, value
        ) SELECT ?, at.attrib_type_id, CAST(? AS CHAR)
      FROM attrib_type at
      WHERE at.code = ?
  };

  $pf_ins_sth          = $dbh->prepare($pf_ins_sth);
  $attrib_ins_sth      = $dbh->prepare($attrib_ins_sth);
  $attrib_ins_cast_sth = $dbh->prepare($attrib_ins_cast_sth);

  my $attribs;
  my $already_inserted;

  my $phenotype_file = $config->{"$data_source\_phenotype_file"};
  my $fh = FileHandle->new($phenotype_file, 'r');
  while (<$fh>) {
    chomp; 
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
    my $attribs = {};

    foreach my $value (@data_attrib) {
      $attribs->{$value} = $hash->{$value};
    }
    my $marker_accession_id = $hash->{marker_accession_id};
    my $phenotypes = $phenotype_adaptor->fetch_by_description($hash->{mp_term_name});
    my $phenotype = @{$phenotypes}[0];
    if (! $hash->{mp_term_name}) {
      print STDERR "No mp_term_name:\n$_\n";
      next;  
    } 
    unless ($phenotype) {
      die("No phenotype description for term: ", $hash->{mp_term_name});
    }
    my $phenotype_id      = $phenotype->dbID();

    my $source_id;
    if ($data_source eq 'impc') {
      $source_id = source($config, $hash->{resource_name}, $hash->{resource_fullname});
    } else {
      $source_id = source($config, $hash->{project_name}, $hash->{project_fullname});
    }

    my $pf_data           = $marker_coords->{$marker_accession_id};
    my $type              = $pf_data->{type};
    my $seq_region_id     = $pf_data->{seq_region_id};
    my $seq_region_start  = $pf_data->{seq_region_start};
    my $seq_region_end    = $pf_data->{seq_region_end};
    my $seq_region_strand = $pf_data->{seq_region_strand};
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
      if (scalar @$strains == 1) {
        my $strain = pop @$strains;
        if (lc $strain->gender ne $gender) {
          $create_new_individual = 1;
        } else {
          $strain_id = $strain->dbID();
        }
      } elsif (scalar @$strains > 1) {
        my @unique = ();
        foreach my $strain (@$strains) {
          if (lc $strain->gender eq $gender) {
            push @unique, $strain;
          }
        }
        if (scalar @unique != 1) {
          $create_new_individual = 1; 
        } else {
          $strain_id = $unique[0]->dbID();
        }
      } else {
        $create_new_individual = 1; 
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
    $attribs->{associated_gene} = $hash->{marker_symbol};
    if ($strain_id) {
      $attribs->{strain_id} = $strain_id;
    }
    if ($hash->{p_value}) {
      $attribs->{p_value} = convert_p_value($hash->{p_value});
    }
    $attribs->{external_id} = $hash->{mp_term_id};

    my @object_ids = ();
    if ($pf_data->{object_id}) {
      @object_ids = split(';', $pf_data->{object_id});
      my $is_significant = 0;
      if ($hash->{resource_name} && $hash->{resource_name} eq 'MGP') {
        $is_significant = 1;
      } else {
        if ($attribs->{p_value} && ($attribs->{p_value} <= 0.0001)) {
          $is_significant = 1;
        }
      } 
      foreach my $object_id (@object_ids) {
        my $key = join('-', ($phenotype_id, $source_id, $object_id, $strain_id));
        unless ($already_inserted->{$key}) {			
          $pf_ins_sth->bind_param(1,$phenotype_id,SQL_INTEGER);
          $pf_ins_sth->bind_param(2,$source_id,SQL_INTEGER);
          $pf_ins_sth->bind_param(3,$type,SQL_VARCHAR);
          $pf_ins_sth->bind_param(4,$object_id,SQL_VARCHAR);
          $pf_ins_sth->bind_param(5,$is_significant,SQL_INTEGER);
          $pf_ins_sth->bind_param(6,$seq_region_id,SQL_INTEGER);
          $pf_ins_sth->bind_param(7,$seq_region_start,SQL_INTEGER);
          $pf_ins_sth->bind_param(8,$seq_region_end,SQL_INTEGER);
          $pf_ins_sth->bind_param(9,$seq_region_strand,SQL_INTEGER);
          $pf_ins_sth->execute();
          $already_inserted->{$key} = 1;
          # get inserted ID
          my $pf_id = $dbh->{'mysql_insertid'};

          # add attribs
          foreach my $attrib_type(grep {defined($attribs->{$_}) && $attribs->{$_} ne ''} @attrib_types) {
            my $value = $attribs->{$attrib_type};
            my $sth = $value =~ m/^\d+(\.\d+)?$/ ? $attrib_ins_cast_sth : $attrib_ins_sth;
            $sth->bind_param(1,$pf_id,SQL_INTEGER);
            $sth->bind_param(2,$value,SQL_VARCHAR);
            $sth->bind_param(3,$attrib_type,SQL_VARCHAR);
            $sth->execute();
          }
        }
      } # end foreach object_id
    } else {
      print STDERR "No object_id for $marker_accession_id\n";
    }
  } # end read phenotype feature file
  $fh->close();
}

sub get_marker_coords {
  my $config = shift;
  my $data_source = shift;
  my $gene_adaptor = $config->{gene_adaptor};
  my $phenotype_file = $config->{"$data_source\_phenotype_file"};
  my $markers;
  my $column_headers;
  # filter for marker accession ids (MGI)
  my $fh = FileHandle->new($phenotype_file, 'r');
  while (<$fh>) {
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
  $fh->close();

  my $coord_file = $config->{coord_file};	
  $fh = FileHandle->new($coord_file, 'r');
  my $header;
  my $marker_coords = $config->{marker_coords} || {};
  while (<$fh>) {
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

  $fh->close();
  $config->{marker_coords} = $marker_coords;
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
  } else {
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

sub get_attrib_types {
  my $config = shift;
  my $dbh = $config->{dbh};
  my $sth = $dbh->prepare(qq{ SELECT code FROM attrib_type });
  $sth->execute();
  my ($attrib_type, @tmp_types);
  $sth->bind_columns(\$attrib_type);
  push @tmp_types, $attrib_type while $sth->fetch();
  $sth->finish;
  return \@tmp_types;
}

sub source {
  my $config = shift;
  my $source_name = shift;
  my $source_description = shift;

  my $dbh = $config->{dbh};

  my $stmt = qq{ SELECT source_id FROM source WHERE name = '$source_name' LIMIT 1};
  my $sth = $dbh->prepare($stmt);
  $sth->execute();
  my $source_id;
  $sth->bind_columns(\$source_id);
  $sth->fetch();

  if (!defined($source_id)) {
    $stmt = qq{
      INSERT INTO source (name, description)
        VALUES ('$source_name', '$source_description')
    };
    $stmt = $dbh->prepare($stmt);
    $stmt->execute();
    $source_id = $dbh->{'mysql_insertid'};
    $stmt->finish();

    my $version = $config->{version};
    if ($version) {
      $dbh->do(qq{UPDATE source SET version=$version WHERE name='$source_name';}) or die $dbh->errstr;
    }
    print STDERR "Added source for $source_name (source_id = $source_id)\n";
  }
  return $source_id;
}

sub usage {

  print qq{
  Usage: perl import_mouse_phenotype_data.pl [OPTION]

  Import mouse phenotype data into a Variation database

  Options:

  -help               Print this message

  Either provide all the necessary database credentials or a registry file:	
  Database credentials are specified on the command line

  -host      Variation database host name (Required)
  -dbname    Variation database name (Required)
  -user      Variation database user (Required)
  -pass      Variation database password (Required)
  -port      Variation database port (Default: 3306)

  -chost     Core database host name (Required for gene-type sources)
  -cdbname   Core database name (Required for gene-type sources)
  -cuser     Core database user (Required for gene-type sources)
  -cpass     Core database password (Default: )
  -cport     Core database port (Default: 3306)

  -registry  Location of registry file with database connections to mouse variation and core databases

  Location of MGI_MRK_Coord.rpt file. The file contains coordinates for markers/genes. The file can be
  downloaded from ftp://ftp.informatics.jax.org/pub/reports/index.html.

  -coord_file      Location of MGI_MRK_Coord.rpt file (Reguired)

  Working directory. The directory is needed for temporary files created during the import. After running the
  script the directory will still contain the file impc_phenotypes.txt which contains all the phenotypes downloaded
  from the IMPC webserver. The file could be relevant if there have been problems during the import.

  -working_dir     Location of a directory to save temporary results

  } . "\n";
  exit(0);
}
