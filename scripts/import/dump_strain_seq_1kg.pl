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
use FindBin qw( $Bin );
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Mapper::RangeRegistry;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(ambiguity_code);

my $species;
my $dump_file;
my $region;
my $new_db_version;
my $chunk_size;

GetOptions('dump_file=s' => \$dump_file,
	   'species=s'   => \$species,
           'new_db_version=i' => \$new_db_version,
	   'region=i'    => \$region,
	   'chunk_size=i' => \$chunk_size
	   );

$region = $ENV{LSB_JOBINDEX} if (!defined $region); #if there is no region as an argument, try to get it from LSF
$chunk_size = 100000 if (!defined $chunk_size); #set a default chunk size of 100kb;

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

# IMPORTANT - MAKE SURE YOUR REGISTRY IS POINTING TO THE CORRECT DATABASES!!!

#Bio::EnsEMBL::Registry->load_registry_from_db( -host => 'ens-staging',
#                                               -db_version => $new_db_version,
#                                               -user => 'ensro',
#					      );

my $registry_file ||= $Bin . "/ensembl.registry.forEMFscript";

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $dbVar = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');


my $slice_adaptor = $dbCore->get_SliceAdaptor();
my $rc_adaptor = $dbVar->get_ReadCoverageAdaptor();
my $ind_adaptor = $dbVar->get_IndividualAdaptor();

#this hash contains all base conversion to avoid comparisons when printing the sequences
my %read_code = qw(M 2 V 2 N 2 H 2 R 2 D 2 W 2 S 2 B 2 Y 2 K 2 C 2 A 2 T 2 G 2 * 2
~ 0 
m 1 v 1 n 1 h 1 r 1 d 1 w 1 s 1 b 1 y 1 k 1 c 1 a 1 t 1 g 1 - 1);

#find out MAX_LEVEL in this species
my $levels = $rc_adaptor->get_coverage_levels();

open DUMP, ">$dump_file" . $region . ".emf"|| die "Could not open file to dump data: $!\n";

#get strains with coverage information, all the columns in the file
my $strains = $ind_adaptor->fetch_all_strains_with_coverage();

if ($species eq 'rat'){
    #let's start with the exceptions... we only dump CELERA strain in rat
    @{$strains} = grep {$_->name eq 'SD'} @{$strains};
}

#we only want to dump Celera data
if ($species eq 'human'){
    @{$strains} = grep {$_->name =~ /1KG\_NA\d+|Venter|Watson/i} @{$strains};
}

# get a slice for the specified region
my $slice = $slice_adaptor->fetch_by_region('chromosome',$region);

&create_file_header();# create the file header

# define initial start and end coordinates
my $start = 1;
my $end   = $start + $chunk_size -1;

# fix the end if it is beyond the end of our slice
$end = $slice->length if($end > $slice->length);

# iterate through the slice, taking a subslice of size $chunk_size each time
while($start <= $slice->length) {
    my $subSlice = $slice->sub_Slice($start, $end);

    #foreach of the subSlices with coverage for one strain, print the header, and the base information
    &print_seq_header($subSlice,$strains); #method to print the SEQ and SCORE block

    #and print the sequence information, one base per row
    &print_base_info($rc_adaptor,$subSlice,$strains);

    print DUMP "//\n"; #end of block
    
    # adjust start and end for next chunk
    $start = $end + 1;
    $end = $start + $chunk_size -1;
    $end = $slice->length if($end > $slice->length);
}

close DUMP || die "Could not close file to dump data: $!\n";

#method to print all the bases in the region, one per line
sub print_base_info{
    my $rc_adaptor = shift;
    my $slice = shift;
    my $strains = shift;
    my %strain_seq;
    
    # start the data chunk in the dump file
    print DUMP "DATA\n";

    foreach my $strain (@{$strains}){
	
	# initialise a blank string to represent this strain's sequence
	my $seq='';
	
	# iterate through each coverage level (NB may include some that don't apply to this strain)
	foreach my $level (sort {$a<=>$b} @$levels){

	    my $rcs;

	    # allow "normal" i.e. non-1KG individual
	    if ($strain->name !~ /^1KG/ and $level <3) {
		$rcs = $rc_adaptor->fetch_all_by_Slice_Sample_depth($slice,$strain,$level);
		get_strain_seq($slice,$rcs,\$seq); #with the seq and the coverage info, make the seq
	    }
	    
	    # for 1KG individual
	    elsif ($strain->name =~ /^1KG/ and $level >2) {
		$rcs = $rc_adaptor->fetch_all_by_Slice_Sample_depth($slice,$strain,$level);
	  
		my $pos_level = get_strain_seq($slice,$rcs,\$seq); #with the seq and the coverage info, make the seq
		
		my %pos_level = %$pos_level;
		
		# 1KG individuals return the coverage level as a hash of level for each pos in this slice
		# we need to copy this to the %strain_seq hash so it can be used when writing out the data
		foreach my $pos (keys %{$pos_level}) {
		    $strain_seq{$strain->name}{$pos} = $pos_level{$pos};
		}
	    }
	}
	
	#we need to apply AlleleFeature to the sequence, adding ambiguity codes and indels
	&apply_AF_to_seq($slice,$strain->name,\$seq);
	
	# put the finished sequence in the %strain_seq hash
	$strain_seq{$strain->name}{'seq'} = $seq;
    }
    
    # method that prints the actual sequence to the dump file
    &print_sequences($slice,$strains,\%strain_seq);

}

#gets the strain and the AF compared to the Slice, and applies them to the sequence
sub apply_AF_to_seq{
    my $slice = shift;
    my $strain = shift;
    my $ref_seq = shift;

    my $strainSlice = $slice->get_by_strain($strain);

    my $allele;
    my $afs = $strainSlice->get_all_AlleleFeatures_Slice();#get rid of 1 as level

    foreach my $af (@{$afs}){
	my $format = "@" . ($af->start-1) . "A" . ($af->end - $af->start +1);
	
	my $base = unpack($format,$$ref_seq);
	
	$allele = ambiguity_code($af->allele_string);

	#if it is a deletion, use the * to mean it is covered by $max_level reads
	if ($base =~ /[A-Z]/) {
	    #we have to consider deletions of more than 1 base in the reference !!!
	    substr($$ref_seq,$af->start-1,$af->end - $af->start + 1,uc($allele) . '*' x ($af->end - $af->start)) if ($allele ne '-');
	    substr($$ref_seq,$af->start-1,$af->end - $af->start + 1,'*' x ($af->end - $af->start +1 )) if ($allele eq '-');
	}

	if ($base =~ /[a-z]/) {
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
    my (@strain_array,%rec);
    my @ref_seq = split//,$slice->seq;
    my $index_strain = 0;

    foreach my $strain (@{$strains}){
      push @{$strain_array[$index_strain]},split //,$strain_seq->{$strain->name}{'seq'};
      
      $rec{$strain_array[$index_strain]} = $strain->name;#keep same name-seq order by use array rather than hash
      $index_strain++;
    }

    for (my $i=0;$i<@ref_seq;$i++){
	#for 1KG individual,use level in hash strain_seq,for normal individual use read_code to get level
	#NB we have to adjust $i to $i+1 for the 1KG individuals as the RC-level data positions are recorded 1-indexed
	#instead of 0-indexed as the sequence etc is
	print DUMP join(
			" ",
			$ref_seq[$i],
			map {
			    $strain_reads .= ($_->[$i] eq '~') ? '0'.' ' : 
				($strain_seq->{$rec{$_}}{$i+1}) ? $strain_seq->{$rec{$_}}{$i+1}.' ' : $read_code{$_->[$i]}.' ';
			    ($_->[$i] eq '*') ? '-' : uc($_->[$i])
			    } @strain_array
			), " ";
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
    my %pos_level;

    my $format;
    
    # if ref_seq is uninitialized, make it a string of ~'s the length of the slice
    $$ref_seq .= '~' x ($slice->end - $slice->start + 1) if length($$ref_seq) == 0;

    while(my $rc = shift @{$rcs}){
      	# fix the coordinates if they are outside the slice boundaries
	$rc->start(1) if ($rc->start <= 0);
	$rc->end($slice->end - $slice->start + 1) if ($rc->end + $slice->start > $slice->end);
	
	# make the unpack format based on the coordinates of the rc relative to the slice
	$format = '@' . ($rc->start - 1) . 'A' . ($rc->end - $rc->start + 1);
	
	#padding level 1 with small letters
	if($rc->level == 1 or $rc->level >2) {
	    # this replaces the appropriate chunk in $$ref_seq with a 
	    # lower-case formatted chunk from the slice sequence
	    substr($$ref_seq, $rc->start - 1, $rc->end - $rc->start + 1) = lc(unpack($format,$slice->seq));
	}

	#padding max_level with big letters
	if($rc->level == 2 or $rc->level > 3) {
	    # as above, but uc() makes upper case for other read levels
	    substr($$ref_seq, $rc->start - 1, $rc->end - $rc->start + 1) = uc(unpack($format,$slice->seq));
	}
	
	# record each position covered at this level - NB these are 1-indexed, not 0-indexed
	# hence the adjustment used in the print_sequences subroutine
	for my $i($rc->start..$rc->end) {
	    $pos_level{$i} = $rc->level;
	}
    }

    return (\%pos_level);
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
	my $rcs = $rc_adaptor->fetch_all_by_Slice_Sample_depth($slice,$strain);#delete 1 here to get data for all levels

	#get all regions covered in the chromsome for all strains
	while(my $rc = shift @{$rcs}){	
	    #insert a new region, without bothering about the strain
	    $range_registry->check_and_register(1,$rc->start,$rc->end);
	}
    }
    
    return $range_registry;
}





sub create_file_header{
    
   
    print DUMP "##FORMAT (resequencing)\n";
    print DUMP "##DATE ",scalar(localtime),"\n";
    print DUMP "##RELEASE ",Bio::EnsEMBL::Registry->software_version(),"\n\n";

    print DUMP "##READ COVERAGE LEVELS\n";
    print DUMP "# 0: No coverage        (Watson, Venter)\n";
    print DUMP "# 1: 1x                 (Watson, Venter)\n";
    print DUMP "# 2: 2x and above       (Watson, Venter)\n";
    print DUMP "# 3: 0-5x               (1000 Genomes)\n";
    print DUMP "# 4: 6-25x              (1000 Genomes)\n";
    print DUMP "# 5: 26x and above      (1000 Genomes)\n\n";

}

sub usage{
    my $msg = shift;

    print STDERR <<EOF;

usage: perl dump_strain_seq.pl <options>

options:
    -dump_file <filename>    file where you want to dump the resequencing data
    -species   <species>     species you want to dump data (default = mouse)
    -region    <region_name> region to dump the data
    -chunk     <chunk_size>  size of chunks (default = 100kb)

EOF
   
die ("\n$msg\n\n");
}
