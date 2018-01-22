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
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Mapper::RangeRegistry;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(ambiguity_code);

my $species;
my $dump_file;
my $region;
my $new_db_version;

GetOptions('dump_file=s' => \$dump_file,
	   'species=s'   => \$species,
           'new_db_version=i' => \$new_db_version,
	   'region=i'    => \$region
	   );

$region = $ENV{LSB_JOBINDEX} if (!defined $region); #if there is no region as an argument, try to get it from LSF

$species ||= 'mouse'; #by default, dump mouse data
usage('You need to enter the file name where you want to dump the data') if (!defined $dump_file); 
usage('You need to give the region you want to dump data') if(!defined $region);

#we need to transform the region name to chromosome names if taken from LSF: X, Y or MT
if ($species eq 'mouse' && defined $region){
    if($region == 20){
	$region = 'X';
    }
    elsif ($region == 21){
	$region = 'Y';
    }
    elsif ($region == 22){
	$region = 'MT';
    }
}
elsif ($species eq 'rat' && defined $region){
    if ($region == 21){
	$region = 'X';
    }
    elsif ($region == 22){
	$region = 'MT';
    }
}
elsif ($species eq 'human' && defined $region){
    if($region == 23){
	$region = 'X';
    }
    elsif ($region == 24){
	$region = 'Y';
    }
    elsif ($region == 25){
	$region = 'MT';
    }
}
else{
    die "Species $species not supported to dump data\n\n";
}

Bio::EnsEMBL::Registry->load_registry_from_db( -host => 'ens-staging',
                                               -db_version => $new_db_version,
                                               -user => 'ensro',
					      );

my $dbVar = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');

print "dbCore is ",ref($dbCore),"\n";
my $slice_adaptor = $dbCore->get_SliceAdaptor();
my $rc_adaptor = $dbVar->get_ReadCoverageAdaptor();
my $ind_adaptor = $dbVar->get_IndividualAdaptor();
#the hashs contains all base conversion to avoid comparisons when printing the sequences
my %read_code = qw(M 2 V 2 N 2 H 2 R 2 D 2 W 2 S 2 B 2 Y 2 K 2 C 2 A 2 T 2 G 2 * 2
~ 0 
m 1 v 1 n 1 h 1 r 1 d 1 w 1 s 1 b 1 y 1 k 1 c 1 a 1 t 1 g 1 - 1);

#find out MAX_LEVEL in this specie
my $levels = $rc_adaptor->get_coverage_levels();
die("The dumping process only works with 2 different levels: 1 and a second one") if (@{$levels} != 2);
my $MAX_LEVEL;
$MAX_LEVEL =$levels->[0] if ($levels->[0] > $levels->[1]);
$MAX_LEVEL =$levels->[1] if ($levels->[0] < $levels->[1]);

open DUMP, ">$dump_file" . $region . ".emf"|| die "Could not open file to dump data: $!\n";

my $strains = $ind_adaptor->fetch_all_strains_with_coverage();   #get strains with coverage information, all the columns in the file
if ($species eq 'rat'){
#let's start with the exceptions... we only dump CELERA strain in rat
    @{$strains} = grep {$_->name eq 'SD'} @{$strains};
}
#we only want to dump Celera data
if ($species eq 'human'){
    @{$strains} = grep {$_->name =~ /Hu\w\w|Venter|Watson/i} @{$strains};
    #print $$strains[0]->name,"\n";
}
#my $slice = $slice_adaptor->fetch_by_region('chromosome','8',109_700_000,110_000_000); #dump this region to find problem 108213779-108237682
#my $slice = $slice_adaptor->fetch_by_region('chromosome','1',1_900_900,2_300_300);
my $slice = $slice_adaptor->fetch_by_region('chromosome',$region);
my $subSlice;
#for each chromosome, get the union of all coverages for all strains
my $range_registry = &get_chromosome_coverage($rc_adaptor,$slice,$strains);
#print "size regions_covered: ",total_size($regions_covered)," and length array ", scalar(@{$regions_covered}),"\n";
&create_file_header() if (defined $range_registry->get_ranges(1)); #create the file header, some regions might not have coverage at all
foreach my $region (@{$range_registry->get_ranges(1)}){
    $region->[0] = 1 if ($region->[0] < 0); #if the region lies outside the boundaries of the slice
    $region->[1] = ($slice->end - $slice->start + 1) if ($region->[1] + $slice->start > $slice->end); 
    $subSlice = $slice->sub_Slice($region->[0],$region->[1],1);
    #foreach of the subSlices with coverage for one strain, print the header, and the base information
    &print_seq_header($subSlice,$strains); #method to print the SEQ and SCORE block
    #and print the sequence information, one base per row
    &print_base_info($rc_adaptor,$subSlice,$strains);
    print DUMP "//\n"; #end of block
}

close DUMP || die "Could not close file to dump data: $!\n";

#method to print all the bases in the region, one per line
sub print_base_info{
    my $rc_adaptor = shift;
    my $slice = shift;
    my $strains = shift;
    my %strain_seq;
    #my $rcs;

    print DUMP "DATA\n";
    foreach my $strain (@{$strains}){
	my $seq='';
	foreach my $level (1,$MAX_LEVEL){
	    my $rcs = $rc_adaptor->fetch_all_by_Slice_Sample_depth($slice,$strain,$level);	    
	    &get_strain_seq($slice,$rcs,\$seq); #with the seq and the coverage info, make the seq
	    #@{$rcs} = ();
	}	
	#we need to apply AlleleFeature to the sequence, adding ambiguity codes and indels
	&apply_AF_to_seq($slice,$strain->name,\$seq);
	$strain_seq{$strain->name} = $seq;
    }
    &print_sequences($slice,$strains,\%strain_seq); #method that prints the actual sequence
    
}

#gets the strain and the AF compared to the Slice, and applies them to the sequence
sub apply_AF_to_seq{
    my $slice = shift;
    my $strain = shift;
    my $ref_seq = shift;

    my $strainSlice = $slice->get_by_strain($strain);
    my $allele;
    my $afs = $strainSlice->get_all_AlleleFeatures_Slice(1);
    foreach my $af (@{$afs}){
	my $format = "@" . ($af->start-1) . "A" . ($af->end - $af->start +1);
	my $base = unpack($format,$$ref_seq);
	$allele = ambiguity_code($af->allele_string);
	#if it is a deletion, use the * to mean it is covered by $max_level reads
	if ($base =~ /[A-Z]/){
	    #we have to consider deletions of more than 1 base in the reference !!!
	    substr($$ref_seq,$af->start-1,$af->end - $af->start + 1,uc($allele) . '*' x ($af->end - $af->start)) if ($allele ne '-');
	    substr($$ref_seq,$af->start-1,$af->end - $af->start + 1,'*' x ($af->end - $af->start +1 )) if ($allele eq '-');
	}
	if ($base =~ /[a-z-]/){
	    #if it is a deletion, use the '-' to mean it is covered by 1 read
	    substr($$ref_seq,$af->start-1,$af->end - $af->start + 1,lc($allele) . '-' x ($af->end - $af->start)) if ($allele ne '-');
	    substr($$ref_seq,$af->start-1,$af->end - $af->start + 1,'-' x ($af->end - $af->start + 1)) if ($allele eq '-');
	}
    }
}

#mehtod to print each line at a time
sub print_sequences{
    my $slice = shift;
    my $strains = shift;
    my $strain_seq = shift;

    my $strain_reads;
    my @strain_array;
    my @ref_seq = split//,$slice->seq;
    my $index_strain = 0;
    foreach my $strain (@{$strains}){
	push @{$strain_array[$index_strain]},split //,$strain_seq->{$strain->name};
	$index_strain++;
    }
    for (my $i=0;$i<@ref_seq;$i++){
	print DUMP join(" ",$ref_seq[$i],map {$strain_reads .= $read_code{$_->[$i]}. ' ';($_->[$i] eq '*') ? '-' : uc($_->[$i]) } @strain_array), " ";
	chomp $strain_reads;
	print DUMP $strain_reads,"\n";
	$strain_reads = '';
   }
}
#with a slice and the regions covered for a particular strain, make the strain_seq
sub get_strain_seq{
    my $slice = shift;
    my $rcs = shift;
    my $ref_seq = shift; #ref to sequence

    my $end = 0;
    my $end_level1 = 0;
    my $format;
    foreach my $rc (@{$rcs}){
	$rc->start(1) if ($rc->start <= 0); #if the region lies outside the boundaries of the slice
	$rc->end($slice->end - $slice->start + 1) if ($rc->end + $slice->start > $slice->end); 
	$$ref_seq .= '~' x ($rc->start - 1 - $end_level1) if ($rc->level == 1);
	$format = '@' . ($rc->start - 1) . 'A' . ($rc->end - $rc->start + 1);
	$$ref_seq .= lc(unpack($format,$slice->seq)) if ($rc->level == 1);
	substr($$ref_seq,$rc->start-1,$rc->end-$rc->start+1,uc(unpack($format,$$ref_seq))) if ($rc->level == $MAX_LEVEL);
		
	$end = $rc->end;
	$end_level1 = $rc->end if ($rc->level == 1);	
    }
    $$ref_seq .= '~' x ($slice->length - $end);
  
}

sub print_seq_header{
    my $slice = shift;
    my $strains = shift;

    print DUMP "SEQ $species reference ", $slice->seq_region_name," ",$slice->start," ",$slice->end," ",$slice->strand,"\n";
    #print the SEQ
    foreach my $strain (@{$strains}){
	print DUMP "SEQ $species ",$strain->name, " WGS\n";
    }
    #and print the SCORE
    foreach my $strain (@{$strains}){
	print DUMP "SCORE aligned ",$strain->name, " reads\n";
    }

}

sub get_chromosome_coverage{
    my $rc_adaptor = shift;
    my $slice = shift;
    my $strains = shift;

    #need to overlap the regions using the RangeRegistry module
    my $range_registry = Bio::EnsEMBL::Mapper::RangeRegistry->new();
    
    #we have to do it for strain, since we might not want to dump all strains....
    foreach my $strain (@{$strains}){
      my $rcs = $rc_adaptor->fetch_all_by_Slice_Sample_depth($slice,$strain,1);

      #get all regions covered in the chromsome for all strains
      foreach my $rc (@{$rcs}){	
	#insert a new region, without bothering about the strain
	$range_registry->check_and_register(1,$rc->start,$rc->end);
      }
    }
    return $range_registry;
}





sub create_file_header{
    
   
    print DUMP "##FORMAT (resequencing)\n";
    print DUMP "##DATE ",scalar(localtime),"\n";
    print DUMP "##RELEASE $new_db_version\n\n";
    #print DUMP "##RELEASE ",Bio::EnsEMBL::Registry->software_version(),"\n\n";

}

sub usage{
    my $msg = shift;

    print STDERR <<EOF;

usage: perl dump_strain_seq.pl <options>

options:
    -dump_file <filename>    file where you want to dump the resequencing data
    -species   <species>     species you want to dump data (default = mouse)
    -region    <region_name> region to dump the data

EOF
   
die ("\n$msg\n\n");
}
