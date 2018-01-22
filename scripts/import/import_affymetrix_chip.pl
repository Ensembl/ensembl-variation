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
use warnings;
use File::Basename qw(fileparse);
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

my ( $infile, $registry_file, $species, $chip_name, $chip_desc, $chip_url, $short_name, $version, $debug, $help );

my $usage = "Usage:\n\timport_illumina_chip.pl -input [Affymetrix annotation file] -chip_name [Affymetrix chip name] 
-description [Affymetrix chip desc] -version [chip version] -chip_url [chip url (optional)]  -short_name [Variation set short name] 
-registry [Path to the registry file] -species [Species name] -debug [debug mode]\n";

GetOptions(
  "input|i=s"     => \$infile,
  "registry|r=s"  => \$registry_file,
  "species=s"     => \$species,
  "chip_name=s"   => \$chip_name,
  "description=s" => \$chip_desc,
  "chip_url=s"    => \$chip_url,
  "short_name=s"  => \$short_name,
  "version=s"     => \$version,
  "debug!"        => \$debug,
  "help|h"        => \$help,
);

die $usage if ($help);

die "\nError - no registry \n\n $usage"    unless defined $registry_file;
die "\nError - no species \n\n $usage"     unless defined $species;
die "\nError - no chip_name \n\n $usage"   unless defined $chip_name;
die "\nError - no description \n\n $usage" unless defined $chip_desc;
die "\nError - no input file \n\n $usage"  unless defined $infile;
die "\nError - no short name \n\n $usage"  unless defined $short_name;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_all($registry_file);

my $slice_adaptor = $registry->get_adaptor($species, 'core', 'slice');
my $var_adaptor   = $registry->get_adaptor($species, 'variation', 'variation');
my $dbh = $var_adaptor->dbc->db_handle;

my $prefix = $chip_name;
   $prefix =~ s/ /_/g;

# Source
my $source_id = source();

# Class attrib ID
my $class_attrib_id = get_snp_class_attrib_id();

## get allele codes to convert bases for new schema
my $al_codes = get_allele_code();

## get sequence ids 
my $seq_region_ids = get_seq_region_ids();


#####################
# Parsing data file #
#####################
my %var_id_list;
my %headers;

debug("\n>>>> Parsing input file <<<<") if ($debug);

open INPUT, "< $infile "||die "Failed to open $infile:$!\n";
while(<INPUT>) {
  chomp $_;
  
  next if ($_ =~ /^#/);
  
  $_ =~ s/"//g;
  
  my @row_data = split(',',$_);
    
  # header
  if(/^Probe\sSet\sID/) { # Headers
    @row_data = map { s/ /_/g; $_} @row_data;
    $headers{uc($row_data[$_])} = $_ for 0..$#row_data;
  }
  elsif (!%headers){
    next;
  }
  else {
    die "ERROR: Couldn't find header data\n" unless %headers;
    
    my %data;
    $data{$_} = $row_data[$headers{$_}] for keys %headers;
    
    my $custom_name = $data{'CUST_SNPID'};
    
    debug("#### Variant $custom_name ####") if ($debug);
    
    my ($var_id, $var_name) = @{get_variant_id(\%data)};
 
    # Check by location
    if (!$var_id) {
      ($var_id, $var_name) = @{check_by_location(\%data, $seq_region_ids)};
    }

    debug("Existing variant $var_name") if ($debug && $var_name);
    
    # Add synonym
    if ($var_id) {
      debug("Variant synonym $custom_name") if ($debug && $var_name ne $custom_name);
      add_variation_synonym($var_id, $custom_name) if ($var_name ne $custom_name);
    }
    # Add new entry
    else {
      # 1 - Add a variation entry
      $var_name = $custom_name;
      $var_id   = add_variation($var_name);
      
      # 2 - Check alleles (flip ?) => check reference using location
      my $seq_region_id = $seq_region_ids->{$data{'CHROMOSOME'}};
      my ($alleles, $match) = @{check_alleles(\%data, $seq_region_id)};

      # 3 - Add alleles
      add_alleles($var_id, $alleles, \%data);
      
      if ($seq_region_id) {
        # 4 - Add VF if alleles and location OK (QC ?) 
        add_variation_feature($var_id, $var_name, $seq_region_id, $data{'PHYSICAL_POSITION'}, $alleles)
      }
      
      if (!$seq_region_id || $match == 0) {
        # 5 - Failed variants
        my $failed_type = (!$seq_region_id) ? 'seq_region' : 'allele';
        add_failed_variation($var_id,$failed_type);
      }
    }
    
    $var_id_list{$var_id} = {'var' => $var_name, 'syn' => $custom_name} if ($var_id);
  }
}

## Variation set + output file ##

debug("\n>>>> Variation set + output file <<<<") if ($debug);

# Get the variation set ID
my $var_set_id = get_variation_set_id($short_name);

# Print variants in a file
my ($filename,$dir) = fileparse($infile);
open OUT, "> $dir/$prefix.tsv" or die $!;

foreach my $var_id (sort { $var_id_list{$a}{'var'} cmp $var_id_list{$b}{'var'} } keys(%var_id_list)) {
  my $name = $var_id_list{$var_id}{'var'};
  my $syn  = $var_id_list{$var_id}{'syn'};
  print OUT "$var_id\t$name\t$syn\n";
  # Add entries in variation_set
  if ($var_set_id) {
    add_variation_set_variation($var_id, $var_set_id);
  }
}
close(OUT);



###########
# Methods #
###########
sub source {

	$chip_url ||= 'http://www.affymetrix.com/';
  # Check if the source already exists, else it create the entry
  if ($dbh->selectrow_arrayref(qq{SELECT source_id FROM source WHERE name='$chip_name';})) {
    $dbh->do(qq{UPDATE IGNORE source SET description='$chip_desc',url='$chip_url',version=$version where name='$chip_name';});
  }
  else {
    $dbh->do(qq{INSERT INTO source (name,description,url,version,type) VALUES ('$chip_name','$chip_desc','$chip_url',$version,'chip');});
  }
  my @source_id = @{$dbh->selectrow_arrayref(qq{SELECT source_id FROM source WHERE name='$chip_name';})};
  return $source_id[0];
}

sub get_variant_id {
  my $data = shift;
  
  my $var_name = ($data->{'DBSNP_RS_ID'} =~ /^rs\d+$/) ? $data->{'DBSNP_RS_ID'} : $data->{'CUST_SNPID'};
  
  my $var = get_variant_by_name($var_name);
  
  $var = get_variant_by_name($data->{'CUST_SNPID'}) if (!$var && $var_name ne $data->{'CUST_SNPID'});
  
  return ($var) ? [$var->dbID,$var->name] : [undef,undef];
}

sub get_variant_by_name {
  my $name = shift;
  
  my $var = $var_adaptor->fetch_by_name($name);
  
  return $var;
}


sub get_snp_class_attrib_id {

  my $sth = $dbh->prepare(qq{SELECT attrib_id FROM attrib WHERE value="SNV"});
  $sth->execute()|| die "Problem getting class attrib ID\n";
  my $res = $sth->fetchall_arrayref();

  return (defined($res)) ? $res->[0]->[0] : 2;
}

=head2  get_allele_code
  Look up allele codes to use to enter alleles if using new variation schema 
  A/T/C/G/- expected for chips
=cut
sub get_allele_code{

  my %code;

  my $allele_code_ext_sth = $dbh->prepare(qq{SELECT allele_code_id, allele FROM allele_code});
  $allele_code_ext_sth->execute()|| die "Problem getting codes\n";
  my $al_data = $allele_code_ext_sth->fetchall_arrayref();

  foreach my $l(@{$al_data}){
    $code{$l->[1]} = $l->[0];
  }
  return \%code;
}

=head2  get_seq_region_ids

  Look up hash of seq_region_ids to use to enter VARIATION_FEATUREs

=cut
sub get_seq_region_ids{

  my %seq_region_ids;

  my $seq_region_ids_ext_sth = $dbh->prepare(qq{SELECT seq_region_id, name FROM seq_region});

  $seq_region_ids_ext_sth->execute() || die "Error extracting seq_region_ids\n";
  my $data = $seq_region_ids_ext_sth->fetchall_arrayref();

  unless(defined $data->[0]->[0]){die "No seq region ids available\n";} 
    
  foreach my $l(@{$data}){
    $seq_region_ids{$l->[1]} = $l->[0];
  }
  return (\%seq_region_ids);

}

=head2 check_by_location
 Look for variant at the same location and with same alleles as chip variant. Stick everything on forward strand to run check
=cut
sub check_by_location {

    my ($data_row, $seq_region_ids )  = @_;
 
    my $chr     = $data_row->{'CHROMOSOME'};
    my $start   = $data_row->{'PHYSICAL_POSITION'};
    my @alleles = ($data_row->{'ALLELE_A'}, $data_row->{'ALLELE_B'});
    
    my $allele_string = join('/',@alleles);
 
    my $seq_region_id = $seq_region_ids->{$chr};
    unless (defined $seq_region_id) {
      debug("check_by_location: Can't find a seq region ID for the chromosome $chr") if ($debug);
      return [undef,undef];
    }
    
    debug("check_by_location: Seeking ($chr, $start, $allele_string)") if ($debug);

    my $var_id_chrom_pos_ext_sth = $dbh->prepare(qq{SELECT 
                                                      vf.variation_id,
                                                      vf.variation_name,
                                                      vf.allele_string,
                                                      vf.seq_region_strand
                                                    FROM variation_feature vf 
                                                    WHERE vf.seq_region_start =?
                                                      AND vf.seq_region_id = ?                                       
                                                   });

    $var_id_chrom_pos_ext_sth->execute($start,$seq_region_id) || die "Error extracting var by location for $start,$chr,$allele_string\n";
    my $data = $var_id_chrom_pos_ext_sth->fetchall_arrayref();
      
    return [undef,undef] unless defined $data->[0]->[0];
    
    debug("check_by_location: Name: $data->[0]->[1] | Alleles: $data->[0]->[2] | Strand: $data->[0]->[3]") if ($debug);

    my $check_allele = $data->[0]->[2];
    my $same = 1;
    if ($check_allele ne $allele_string) {
      debug("check_by_location: Found existing variant at this location, with different alleles");
      
      my %exp;
      my @exp_al = split/\//,$check_allele;
      foreach my $exp_al(@exp_al){
        $exp{$exp_al} = 1;
      }

      foreach my $al(@alleles){
	      unless (defined $exp{$al}){$same = 0;}
      }
      
      if( $same == 0 && $data->[0]->[3] eq "-1"){
        $same = 1; # Reintialise the "$same" variable
        
        $check_allele = get_reverse_comp($check_allele);
        debug("check_by_location: Try to reverse complement the existing variant alleles on the FORWARD strand") if ($debug);
        
        %exp = ();
        @exp_al = split/\//,$check_allele;
        foreach my $exp_al(@exp_al){
          $exp{$exp_al} = 1;
        }

        foreach my $al(@alleles){
	        unless (defined $exp{$al}){$same = 0;}
        }
      }
    }
    
    if( $same == 1 ){        
      ## all alleles found
      debug("check_by_location: Existing variant found: $data->[0]->[2]") if ($debug);
      return [$data->[0]->[0],$data->[0]->[1]];
    }
    else {
      debug("check_by_location: Strand & allele string incompatibile! See: $allele_string (data) vs $check_allele (ref)  | $data->[0]->[1] & $data->[0]->[2] (Strand $data->[0]->[3])") if ($debug);
      return [undef,undef];
    }
}


sub check_alleles {
  my ($data, $seq_region_id)  = @_;
 
  my $chr      = $data->{'CHROMOSOME'};
  my $position = $data->{'PHYSICAL_POSITION'};
  my @alleles  = ($data->{'ALLELE_A'},$data->{'ALLELE_B'});

  my $allele_string = join('/',@alleles);

  my $match = 0;
  
  return [$allele_string, $match] unless $seq_region_id;

  my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $position, $position, 1);
  return [$allele_string, $match] unless $slice;
  
  my $ref_allele = $slice->seq;
  
  # Sort by putting the reference allele first, if possible
  @alleles = sort { ($a !~ /$ref_allele/ cmp $b !~ /$ref_allele/) || $a cmp $b } @alleles;
  $allele_string = join('/',@alleles);
  
  debug("check_alleles: Seeking ($chr, $position, $allele_string)") if ($debug);
  
  foreach my $al (@alleles){
	  $match = 1 if ($ref_allele eq $al);
  }
    
  ## all alleles found
  return [$allele_string, $match];
}


sub add_variation {
  my $name = shift;

  debug("Add variation $name") if ($debug);
  $dbh->do(qq{INSERT INTO variation (name,source_id,class_attrib_id) VALUES ('$name', $source_id, $class_attrib_id);});
  
  my $var_id = $dbh->{'mysql_insertid'};
  
  return $var_id;
}

sub add_variation_synonym {
  my $var_id   = shift;
  my $syn_name = shift;
  
  debug("Add variation synonym $syn_name") if ($debug);
  $dbh->do(qq{INSERT INTO variation_synonym (variation_id,name,source_id) VALUES ($var_id, '$syn_name', $source_id);});
}

sub add_alleles {
  my $var_id  = shift;
  my $alleles = shift;
  my $data    = shift;
  
  my @alleles_list;
  
  if ($alleles && $alleles ne '') {
    @alleles_list = split('/',$alleles);
  }
  else {
    @alleles_list = ($data->{'ALLELE_A'},$data->{'ALLELE_B'});
  }

  @alleles_list = map{ $al_codes->{$_} } @alleles_list;
  
  foreach my $al_id (@alleles_list) {
    next if (!$al_id || $al_id eq '' || $al_id !~ /^\d+$/);
    $dbh->do(qq{INSERT INTO allele (variation_id,allele_code_id) VALUES ($var_id, $al_id);});
    debug("Add allele $al_id") if ($debug);
  }
}

sub add_variation_feature {
  my $var_id        = shift;
  my $var_name      = shift;
  my $seq_region_id = shift;
  my $position      = shift;
  my $alleles       = shift;
  
  debug("Add variation_feature") if ($debug);
  $dbh->do(qq{ INSERT INTO variation_feature (
                 variation_id,
                 variation_name,
                 seq_region_id,
                 seq_region_start,
                 seq_region_end,
                 seq_region_strand,
                 allele_string,
                 map_weight,
                 class_attrib_id,
                 source_id
               ) VALUES (
                 $var_id, '$var_name',
                 $seq_region_id, $position, $position, 1,
                 '$alleles', 1, $class_attrib_id, $source_id)
             });
}

sub add_failed_variation {
  my $var_id = shift;
  my $type   = shift;
  
  my $failed_id = ($type eq 'seq_region') ? 5 : 2;
  
  debug("Add failed_variation") if ($debug);
  $dbh->do(qq{ INSERT INTO failed_variation (variation_id,failed_description_id) VALUES ($var_id,$failed_id)});
  $dbh->do(qq{ UPDATE variation SET display=0 WHERE variation_id=$var_id});
  $dbh->do(qq{ UPDATE variation_feature SET display=0 WHERE variation_id=$var_id});
}

sub get_variation_set_id {
  my $short_name = shift;
  
  my $sth1 = $dbh->prepare(qq{ SELECT vs.variation_set_id
                               FROM variation_set vs, attrib a
                               WHERE vs.short_name_attrib_id=a.attrib_id
                                 AND a.value=?
                             });

  $sth1->execute($short_name) || die "Error extracting variation set ID from the short name $short_name\n";
  my $data = $sth1->fetchall_arrayref();
  
  if (defined $data->[0]->[0]) {
    return $data->[0]->[0];
  }
  else {
    print STDERR "ERROR: Variation set with the short name '$short_name' not found!\n";
    return;
  }
}

sub add_variation_set_variation {
  my $var_id     = shift;
  my $var_set_id = shift;
  
  $dbh->do(qq{ INSERT INTO variation_set_variation (variation_id,variation_set_id) VALUES ($var_id,$var_set_id) }); 
}

sub get_reverse_comp {
  my $al = shift;
  
  my @alleles = split('/', $al);

  foreach my $allele (@alleles) {
    reverse_comp(\$allele);
  }
  return join('/',@alleles);
}

sub debug {
  my $msg = shift;
  print STDERR "$msg\n";
}

