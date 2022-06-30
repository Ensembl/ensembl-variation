#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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
  <http://www.ensembl.org/Help/Contact>.

=cut

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::VariationFeature;

use Getopt::Long;
use FileHandle;

# Print instructions if run without parameters
usage() unless (scalar(@ARGV));

my ($host, $port, $user, $pass, $chosen_species, $version, $dir, $formats, $write_to_db, $help);

GetOptions(
  'host=s'   => \$host,
  'user=s'   => \$user,
  'pass=s'   => \$pass,
  'port=i'   => \$port,
  'version=i' => \$version,
  'species=s' => \$chosen_species,
  'dir=s'      => \$dir,
  'formats=s'  => \$formats,
  'write_to_db' => \$write_to_db,
  'help!' => \$help,
);
usage() if $help;
die "Error: provide the Ensembl version with option '-v'\n" if !$version;

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
$formats ||= 'ensembl,vcf,id,hgvs,spdi';

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
  
  my $div_bacteria = $sa->dbc->dbname() =~ /^bacteria.*/ ? 1 : 0;

  print STDERR "Generating examples for $species, assembly $assembly...\n";
  
  if($real_vfa && $doing_id) {
    my $name;
    my $sth = $real_vfa->db->dbc->prepare(qq{
      SELECT variation_name
      FROM variation_feature vf
        LEFT JOIN failed_variation fv
        ON vf.variation_id = fv.variation_id
      WHERE consequence_types LIKE ?
      AND fv.variation_id IS NULL
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
  
  # special address of bacteria species due to:
  # a) undef seq_region_start/end issue unless explicitly defined in variationFeature and
  # b) VariationFeature with seq_region_start/end defined for non-bacteria species returning different results than with them left out
  # c) transcriptVariation consequence apprearing intergenic when missense expected possibly due to real/fake adaptors mix
  if ($div_bacteria) {
    my ($slice, $trs);

    # create a missense
    my $tr = undef;
    while(!defined($tr)) {
      next SPECIES unless scalar @slices;
      ($slice, $trs) = @{select_slice(\@slices)};
      $tr = select_transcript($trs, $div_bacteria);
    }
    my $pos_snp = $tr->coding_region_start + 3;
    my $sub_slice_snp = $slice->sub_Slice($pos_snp, $pos_snp);
    my $ref_seq_snp = $sub_slice_snp->seq;

    my @tmp_alts = sort {rand() <=> rand()} grep {$_ ne $ref_seq_snp} @alts;
    my $alt = shift @tmp_alts;

    #VariationFeature objects for bacteria have to be created explicitly specifying seq_region_start,seq_region_end, this will prevent these two being undef in subsequent uses, this is known in eg bacteria
    my $tmp_vf_snp = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
      start          => $pos_snp,
      end            => $pos_snp,
      allele_string  => $ref_seq_snp.'/'.$alt,
      strand         => 1,
      map_weight     => 1,
      adaptor        => $vfa,
      slice          => $slice,
      seq_region_start => $slice->seq_region_start,
      seq_region_end   => $slice->seq_region_end
    });
    dump_vf($tmp_vf_snp, \%files, \%web_data);

    # create a frameshift in a different transcript on the same slice (for bacteria most of the time there is one single top level slice)
    $tr = undef;
    while(!defined($tr)) {
      $tr = select_transcript($trs,$div_bacteria);
    }
    my $pos_fs = $tr->coding_region_start + 3;
    my $sub_slice_fs = $slice->sub_Slice($pos_fs, $pos_fs);
    my $ref_seq_fs = $sub_slice_fs->seq;
    my $tmp_vf_fs = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
      start          => $pos_fs,
      end            => $pos_fs,
      allele_string  => $ref_seq_fs.'/-',
      strand         => 1,
      map_weight     => 1,
      adaptor        => $vfa,
      slice          => $slice,
      seq_region_start => $slice->seq_region_start,
      seq_region_end   => $slice->seq_region_end
    });

    dump_vf($tmp_vf_fs, \%files, \%web_data);
  } else {
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

  if ($write_to_db) {
    write_to_db($species, \%web_data);
  }

  # exit 0;
  }
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
  my $div_bacteria = shift;

  $div_bacteria ||= 0;

  # we want a transcript on the fwd strand with an intron that's protein coding
  my $biotype = '';
  my $intron_count = ($div_bacteria == 1 ? 1 :0 );
  my $strand = 0;
  my $crs = undef;
  my $tr;
  
  while($biotype ne 'protein_coding' || $intron_count == 0 || $strand ne 1 || !defined($crs)) {
    return undef unless scalar @$trs;
    
    $tr = shift @$trs;
    $biotype = $tr->biotype;
    $strand = $tr->seq_region_strand;
    $crs = $tr->coding_region_start;
    $intron_count = ($div_bacteria == 1 ? 1 : scalar @{$tr->get_all_Introns});
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
    return [$slice, $trs] if (! $slice);
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

sub convert_to_spdi {
  my $config = shift;
  my $vf = shift;

  my $tvs = $vf->get_all_TranscriptVariations;
  my @return;
  if(defined($tvs)) {
    push @return, map {values %{$vf->spdi_genomic()}} @$tvs;
  }

  @return = grep {defined($_)} @return;

  return [$return[0] || undef];
}

sub write_to_db {
  my $species = shift;
  my $web_data = shift;

  my $mca = $reg->get_adaptor($species, 'core', 'MetaContainer');

  foreach my $key(keys %$web_data) {
    my $meta_key = 'sample.'.lc($key);
    # Note that the newline is intentionally not interpreted;
    # we want the literal '\n' embedded in the database value,
    # the webcode will subsequently treat that as a proper newline.
    my $meta_value = join('\n', @{$$web_data{$key}});

    $mca->delete_key($meta_key);
    $mca->store_key_value($meta_key, $meta_value);
  }
}

sub usage {
  print qq{
  Usage: perl generate_vep_examples.pl -v [VERSION] [OPTIONS]
  
  Generate VEP examples for one or more species. If no databases are set up, it
  automatically connects to the public database (ensembldb.ensembl.org).
  
  Mandatory arguments:
  
    -v            Ensembl version, e.g. 108 (Required)
  
  Optional arguments:
  
    -dir          Directory where to write output (default: current directory)
    -write_to_db  Write generated VEP examples to 'meta' table from the core
                  database (turned off by default)
    -formats      Comma-separated list of formats to output; prints all formats
                  by default: 'ensembl,vcf,id,hgvs,spdi'
    -species      Filter by species (deafult: uses all species in database)
    -help         Print this message
    
  Load database from parameters:
      
    -host         Host
    -port         Port number
    -user         MySQL user name
    -pass         MySQL password
  } . "\n";
  exit(0);
}
