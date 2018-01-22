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

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::VariationFeature;

use Getopt::Long;
use FileHandle;

my ($host, $port, $user, $pass, $chosen_species, $version, $dir, $formats);

GetOptions(
  'host=s'   => \$host,
  'user=s'   => \$user,
  'pass=s'   => \$pass,
  'port=i'   => \$port,
  'version=i' => \$version,
  'species=s' => \$chosen_species,
  'dir=s'      => \$dir,
  'formats=s'  => \$formats,
);

if(defined($host) && $host =~ /staging|variation|livemirror/) {
  $port ||= 3306;
  $user ||= 'ensro';
  $pass ||= '';
}

else {
  $host ||= 'ensembldb.ensembl.org';
  $port ||= 5306;
  $user ||= 'anonymous';
  $pass ||= '';
}

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
  -host       => $host,
  -user       => $user,
  -pass       => $pass,
  -port       => $port,
  -db_version => $version,
);

$dir ||= '.';
$formats ||= 'ensembl,vcf,id,hgvs,pileup';

my @formats = split(',', $formats);

# special case ID
my $doing_id = 0;
if(grep {$_ eq 'id'} @formats) {
  $doing_id = 1;
  @formats = grep {$_ ne 'id'} @formats;
}

my @all_species = sort @{$reg->get_all_species()};

@all_species = grep {$_ eq $chosen_species} @all_species if defined($chosen_species);

my @alts = qw(A C G T);
my $do = 0;

SPECIES: foreach my $species(@all_species) {

  my %web_data;
  
  my $sa = $reg->get_adaptor($species, 'core', 'slice');
  my $ta = $reg->get_adaptor($species, 'core', 'transcript');
  my ($highest_cs) = @{$reg->get_adaptor($species, 'core', 'coordsystem')->fetch_all()};
  my $assembly = $highest_cs->version();
  my $vfa = Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor->new_fake($species);
  my $real_vfa = $reg->get_adaptor($species, 'variation', 'variationfeature');
  
  print STDERR "Doing $species, assembly $assembly\n";
  
  if($real_vfa && $doing_id) {
    my $name;
    my $sth = $real_vfa->db->dbc->prepare(qq{
      SELECT variation_name
      FROM variation_feature
      WHERE consequence_types LIKE ?
      LIMIT 1
    });

    my $fn = $dir.'/'.$species.'_'.$assembly.'.id';
    open OUT, ">$fn" or die("ERROR: Could not open file $fn\n");
    
    foreach my $type(qw(missense intron frameshift)) {
      $sth->execute('%'.$type.'%');
      $sth->bind_columns(\$name);
      $sth->fetch();
      print OUT "$name\n" if $name;
      push @{$web_data{"VEP_ID"}}, $name;
    }
    
    close OUT;
    
    $sth->finish(); 
  }
  
  my %files;
  foreach my $format(@formats) {
    $files{$format} = FileHandle->new;
    my $fn = $dir.'/'.$species.'_'.$assembly.'.'.$format;
    $files{$format}->open('>'.$fn) or die("ERROR: Could not open file $fn\n");
  }
  
  my @slices = sort {
    $b->seq_region_name =~ /^[\d+]$/ <=> $a->seq_region_name =~ /^[\d+]$/ ||
    $b->seq_region_name =~ /^[MXY+]$/ <=> $a->seq_region_name =~ /^[MXY+]$/ ||
    $b->length <=> $a->length} @{$sa->fetch_all('toplevel')
  };
  
  SLICE:
  my ($slice, $trs);
  
  # create a missense
  my $tr = undef;
  while(!defined($tr)) {
    next SPECIES unless scalar @slices;
    ($slice, $trs) = @{select_slice(\@slices)};
    $tr = select_transcript($trs);
  }
  
  my $pos = $tr->coding_region_start + 3;
  my $sub_slice = $slice->sub_Slice($pos, $pos);
  my $ref_seq = $sub_slice->seq;
  
  my $con = '';
  my $tmp_vf;
  my @tmp_alts = sort {rand() <=> rand()} grep {$_ ne $ref_seq} @alts;
  my $alt;
  
  while($con !~ /missense/ && scalar @tmp_alts) {
    $alt = shift @tmp_alts;
    
    $tmp_vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
      start          => $pos,
      end            => $pos,
      allele_string  => $ref_seq.'/'.$alt,
      strand         => 1,
      map_weight     => 1,
      adaptor        => $vfa,
      slice          => $slice
    });
    
    $con = join ",", @{$tmp_vf->consequence_type};
  }
  
  if($con =~ /missense/) {
    dump_vf($tmp_vf, \%files, \%web_data);
    # printf("%s %s %i %i %s\/%s 1\n", $species, $slice->seq_region_name, $pos, $pos, $ref_seq, $alt);
  }
  else {
    goto SLICE;
  }
  
  # create an intron_variant
  $tr = undef;
  while(!defined($tr)) {
    next SPECIES unless scalar @slices;
    ($slice, $trs) = @{select_slice(\@slices)};
    $tr = select_transcript($trs);
  }
  
  my @introns = @{$tr->get_all_Introns};
  $con = '';
  
  while($con !~ /intron/ && scalar @introns) {
    my $intron = shift @introns;
    $pos = $intron->start + 15;
    $sub_slice = $slice->sub_Slice($pos, $pos);
    $ref_seq = $sub_slice->seq;
    
    ($alt) = sort {rand() <=> rand()} grep {$_ ne $ref_seq} @alts;
    
    $tmp_vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
      start          => $pos,
      end            => $pos,
      allele_string  => $ref_seq.'/'.$alt,
      strand         => 1,
      map_weight     => 1,
      adaptor        => $vfa,
      slice          => $slice
    });
    
    $con = join ",", @{$tmp_vf->consequence_type};
  }
  
  if($con =~ /intron/) {
    dump_vf($tmp_vf, \%files, \%web_data);
    # printf("%s %s %i %i %s\/%s 1\n", $species, $slice->seq_region_name, $pos, $pos, $ref_seq, $alt);
  }
  
  # create a frameshift
  $tr = undef;
  while(!defined($tr)) {
    next SPECIES unless scalar @slices;
    ($slice, $trs) = @{select_slice(\@slices)};
    $tr = select_transcript($trs);
  }
  
  $pos = $tr->coding_region_start + 3;
  $con = '';
  
  while($con !~ /frameshift/ && scalar @tmp_alts) {
    $pos++;
    $sub_slice = $slice->sub_Slice($pos, $pos);
    $ref_seq = $sub_slice->seq;
    
    $tmp_vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
      start          => $pos,
      end            => $pos,
      allele_string  => $ref_seq.'/-',
      strand         => 1,
      map_weight     => 1,
      adaptor        => $vfa,
      slice          => $slice
    });
    
    $con = join ",", @{$tmp_vf->consequence_type};
  }
  
  if($con =~ /frameshift/) {
    dump_vf($tmp_vf, \%files, \%web_data);
    # printf("%s %s %i %i %s\/%s 1\n", $species, $slice->seq_region_name, $pos, $pos, $ref_seq, '-');
  }

  my $fn = $dir.'/'.ucfirst($species).'.txt';
  open OUT, ">$fn";
  foreach my $key(keys %web_data) {
    print OUT sprintf("%-17s = %s\n", $key, join('\n', @{$web_data{$key}}));
  }
  close OUT;
  # exit 0;
}

sub dump_vf {
  my $vf = shift;
  my $files = shift;
  my $web_data = shift;
  
  foreach my $format(keys %$files) {
    my $method = 'convert_to_'.lc($format);
    my $method_ref = \&$method;
    my $file = $files->{$format};
    my @ret = grep {defined($_)} @{&$method_ref({}, $vf)};
    print $file join("\t", @ret)."\n";

    push @{$web_data->{"VEP_".uc($format)}}, join(" ", @ret);
  }
  
  return;
}

sub select_transcript {
  my $trs = shift;
  
  # we want a transcript on the fwd strand with an intron that's protein coding
  my $biotype = '';
  my $intron_count = 0;
  my $strand = 0;
  my $crs = undef;
  my $tr;
  
  while($biotype ne 'protein_coding' || $intron_count == 0 || $strand ne 1 || !defined($crs)) {
    return undef unless scalar @$trs;
    
    $tr = shift @$trs;
    $biotype = $tr->biotype;
    $strand = $tr->seq_region_strand;
    $crs = $tr->coding_region_start;
    $intron_count = scalar @{$tr->get_all_Introns};
  }
  
  print STDERR "Chose transcript ".$tr->stable_id."\n";
  
  return $tr;
}

sub select_slice {
  my $slices = shift;
  my $trs = [];
  my $slice;
  
  my $ta = $slices->[0]->adaptor->db->get_TranscriptAdaptor;
  
  # we want a slice with transcripts on
  while(scalar @$trs == 0) {
    $slice = shift @$slices;
    $trs = $ta->fetch_all_by_Slice($slice);
  }
  
  print STDERR "Chose slice ".$slice->seq_region_name."\n";
  
  return [$slice, $trs];
}

# converts to Ensembl format
sub convert_to_ensembl {
  my $config = shift;
  my $vf = shift;
    
  return [
    $vf->{chr} || $vf->seq_region_name,
    $vf->start,
    $vf->end,
    $vf->allele_string,
    $vf->strand,
    $vf->variation_name
  ];
}


# converts to pileup format
sub convert_to_pileup {
  my $config = shift;
  my $vf = shift;
    
  # look for imbalance in the allele string
  my %allele_lengths;
  my @alleles = split /\//, $vf->allele_string;
    
  foreach my $allele(@alleles) {
    $allele =~ s/\-//g;
    $allele_lengths{length($allele)} = 1;
  }
    
  # in/del
  if(scalar keys %allele_lengths > 1) {
        
    if($vf->allele_string =~ /\-/) {
            
      # insertion?
      if($alleles[0] eq '-') {
        shift @alleles;
            
        for my $i(0..$#alleles) {
          $alleles[$i] =~ s/\-//g;
          $alleles[$i] = '+'.$alleles[$i];
        }
      }
            
      else {
        @alleles = grep {$_ ne '-'} @alleles;
                
        for my $i(0..$#alleles) {
          $alleles[$i] =~ s/\-//g;
          $alleles[$i] = '-'.$alleles[$i];
        }
      }
            
      @alleles = grep {$_ ne '-' && $_ ne '+'} @alleles;
            
      return [
        $vf->{chr} || $vf->seq_region_name,
        $vf->start - 1,
        '*',
        (join "/", @alleles),
      ];
    }
        
    else {
      warn "WARNING: Unable to convert variant to pileup format on line number ", $config->{line_number} unless defined($config->{quiet});
      return [];
    }
        
  }
    
  # balanced sub
  else {
    return [
      $vf->{chr} || $vf->seq_region_name,
      $vf->start,
      shift @alleles,
      (join "/", @alleles),
    ];
  }
}

# converts to HGVS (hackily returns many lines)
sub convert_to_hgvs {
  my $config = shift;
  my $vf = shift;
    
  # ensure we have a slice
  # $vf->{slice} ||= get_slice($config, $vf->{chr}, undef, 1);
    
  my $tvs = $vf->get_all_TranscriptVariations;
    
  my @return;# = values %{$vf->get_all_hgvs_notations()};
    
  if(defined($tvs)) {
    push @return, map {values %{$vf->get_all_hgvs_notations($_->transcript, 'c')}} @$tvs;
    # push @return, map {values %{$vf->get_all_hgvs_notations($_->transcript, 'p')}} @$tvs;
  }
  
  @return = grep {defined($_)} @return;
    
  return [$return[0] || undef];
}

sub convert_to_vcf {
  my $config = shift;
  my $vf = shift;
    
  # look for imbalance in the allele string
  if($vf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
    my %allele_lengths;
    my @alleles = split /\//, $vf->allele_string;
        
    map {reverse_comp(\$_)} @alleles if $vf->strand < 0;
        
    foreach my $allele(@alleles) {
      $allele =~ s/\-//g;
      $allele_lengths{length($allele)} = 1;
    }
        
    # in/del/unbalanced
    if(scalar keys %allele_lengths > 1) {
            
      # we need the ref base before the variation
      # default to N in case we can't get it
      my $prev_base = 'N';
            
      if(defined($vf->slice)) {
        my $slice = $vf->slice->sub_Slice($vf->start - 1, $vf->start - 1);
        $prev_base = $slice->seq if defined($slice);
      }
            
      for my $i(0..$#alleles) {
        $alleles[$i] =~ s/\-//g;
        $alleles[$i] = $prev_base.$alleles[$i];
      }
            
      return [
        $vf->{chr} || $vf->seq_region_name,
        $vf->start - 1,
        $vf->variation_name || '.',
        shift @alleles,
        (join ",", @alleles),
        '.', '.', '.'
      ];
            
    }
        
    # balanced sub
    else {
      return [
        $vf->{chr} || $vf->seq_region_name,
        $vf->start,
        $vf->variation_name || '.',
        shift @alleles,
        (join ",", @alleles),
        '.', '.', '.'
      ];
    }
  }
    
  # SV
  else {
        
    # convert to SO term
    my %terms = (
      'insertion' => 'INS',
      'deletion' => 'DEL',
      'tandem_duplication' => 'TDUP',
      'duplication' => 'DUP'
    );
        
    my $alt = '<'.($terms{$vf->class_SO_term} || $vf->class_SO_term).'>';
        
    return [
      $vf->{chr} || $vf->seq_region_name,
      $vf->start,
      $vf->variation_name || '.',
      '.',
      $alt,
      '.', '.', '.'
    ];
  }
}