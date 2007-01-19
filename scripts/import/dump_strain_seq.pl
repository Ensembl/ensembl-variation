#!/usr/local/ensembl/bin/perl -w

use strict;
use warnings;
use FindBin qw( $Bin );
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Mapper::RangeRegistry;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(ambiguity_code);

use Data::Dumper;

my $species;
my $dump_file;

GetOptions('dump_file=s' => \$dump_file,
	   'species=s'   => \$species,
	   );

warn("Make sure you have a updated ensembl.registry file!\n");

$species ||= 'mouse'; #by default, dump mouse data
usage('You need to enter the file name where you want to dump the data') if (!defined $dump_file); 
my $registry_file ||= $Bin . "/ensembl.registry";

Bio::EnsEMBL::Registry->load_all( $registry_file );


my $dbVar = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');

my $slice_adaptor = $dbCore->get_SliceAdaptor();
my $rc_adaptor = $dbVar->get_ReadCoverageAdaptor();
my $ind_adaptor = $dbVar->get_IndividualAdaptor();

#find out MAX_LEVEL in this specie
my $levels = $rc_adaptor->get_coverage_levels();
die("The dumping process only works with 2 different levels: 1 and a second one") if (@{$levels} != 2);
my $MAX_LEVEL;
$MAX_LEVEL =$levels->[0] if ($levels->[0] > $levels->[1]);
$MAX_LEVEL =$levels->[1] if ($levels->[0] < $levels->[1]);

open DUMP, ">$dump_file" || die "Could not open file to dump data: $!\n";

&create_file_header(); #create the file header
my $strains = $ind_adaptor->fetch_all_strains_with_coverage();#get strains with coverage information, all the columns in the file

#my $slices = $slice_adaptor->fetch_all('chromosome');
#foreach my $slice (@{$slices}){
#my $slice = $slice_adaptor->fetch_by_region('chromosome','12',36_085_000,36_085_300);
my $slice = $slice_adaptor->fetch_by_region('chromosome','19');
print "Processing chromosome ", $slice->seq_region_name,"\n";
#for each chromosome, get the union of all coverages for all strains
my $regions_covered = &get_chromosome_coverage($rc_adaptor,$slice);
foreach my $region (@{$regions_covered}){
#foreach of the subSlices with coverage for one strain, print the header, and the base information
    &print_seq_header($region,$strains); #method to print the SEQ and SCORE block
    #and print the sequence information, one base per row
    &print_base_info($rc_adaptor,$region,$strains);
    print DUMP "//\n"; #end of block
}
close DUMP || die "Could not close file to dump data: $!\n";
#}

#method to print all the bases in the region, one per line
sub print_base_info{
    my $rc_adaptor = shift;
    my $slice = shift;
    my $strains = shift;
    my %strain_seq;
    my $rcs;
    my $seq;

    print DUMP "DATA\n";
    foreach my $strain (@{$strains}){
	$rcs = $rc_adaptor->fetch_all_by_Slice_Sample_depth($slice,$strain);
	$seq =  &get_strain_seq($slice,$rcs); #with the seq and the coverage info, make the seq
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

    my $afs = $strainSlice->get_all_AlleleFeatures_Slice();
    foreach my $af (@{$afs}){
	my $base = substr($$ref_seq,$af->start-1,$af->end - $af->start + 1);
	substr($$ref_seq,$af->start-1,$af->end - $af->start + 1,uc(ambiguity_code($af->allele_string))) if ($base =~ /[A-Z]/);
	substr($$ref_seq,$af->start-1,$af->end - $af->start + 1,lc(ambiguity_code($af->allele_string))) if ($base =~ /[a-z]/);
    }
}

#mehtod to print each line at a time
sub print_sequences{
    my $slice = shift;
    my $strains = shift;
    my $strain_seq = shift;

    my %strain_reads;
    my $base;
    for (my $i = 0;$i<$slice->length - 1;$i++){
	#print the reference base
	print DUMP substr($slice->seq,$i,1), " ";
	foreach my $strain (@{$strains}){
	    #get the base and the number read information
	    $base = substr($strain_seq->{$strain->name},$i,1);
	    $strain_reads{$strain->name} = 0 if ($base eq '.');
	    $strain_reads{$strain->name} = 1 if ($base =~ /[a-z]/);
	    $strain_reads{$strain->name} = 2 if ($base =~ /[A-Z]/);
	    #and print the base 
	    print DUMP uc($base)," ";
	}
	#and print the score information
	foreach my $strain (@{$strains}){
	    print DUMP $strain_reads{$strain->name}, " ";
	}
	print DUMP "\n";
    }
}
#with a slice and the regions covered for a particular strain, make the strain_seq
sub get_strain_seq{
    my $slice = shift;
    my $rcs = shift;
    
    my $seq;
    my $end = 0;
    my $end_level1 = 0;
    foreach my $rc (@{$rcs}){
	$rc->start(1) if ($rc->start < 0); #if the region lies outside the boundaries of the slice
	$rc->end($slice->end - $slice->start + 1) if ($rc->end + $slice->start > $slice->end); 
	$seq .= '.' x ($rc->start - 1 - $end_level1) if ($rc->level == 1);
	$seq .= lc(substr($slice->seq,$rc->start-1,$rc->end - $rc->start +1)) if ($rc->level == 1);
	substr($seq,$rc->start-1,$rc->end-$rc->start+1,uc(substr($seq,$rc->start-1,$rc->end - $rc->start +1))) if ($rc->level == $MAX_LEVEL);
	$end = $rc->end;
	$end_level1 = $rc->end if ($rc->level == 1);	
    }
    $seq .= '.' x ($slice->length - $end);
    return $seq;
}

sub print_seq_header{
    my $slice = shift;
    my $strains = shift;

    print DUMP "SEQ mouse reference ", $slice->seq_region_name," ",$slice->start," ",$slice->end," ",$slice->strand,"\n";
    #print the SEQ
    foreach my $strain (@{$strains}){
	print DUMP "SEQ mouse ",$strain->name, " WGS\n";
    }
    #and print the SCORE
    foreach my $strain (@{$strains}){
	print DUMP "SCORE aligned ",$strain->name, " reads\n";
    }

}

sub get_chromosome_coverage{
    my $rc_adaptor = shift;
    my $slice = shift;

    my $rcs = $rc_adaptor->fetch_all_by_Slice_Sample_depth($slice);
    #need to overlap the regions using the RangeRegistry module
    my $range_registry = Bio::EnsEMBL::Mapper::RangeRegistry->new();

    #get all regions covered in the chromsome for all strains
    foreach my $rc (@{$rcs}){	
	#insert a new region, without bothering about the strain
	$range_registry->check_and_register(1,$rc->start,$rc->end);
    }

    #and return slices for all the regions covered
    my @sub_Slices;
    foreach my $region (@{$range_registry->get_ranges(1)}){
	$region->[0] = 1 if ($region->[0] < 0); #if the region lies outside the boundaries of the slice
	$region->[1] = ($slice->end - $slice->start + 1) if ($region->[1] + $slice->start > $slice->end); 
	push @sub_Slices, $slice->sub_Slice($region->[0],$region->[1],1);#create the subSlice
    }
    return \@sub_Slices;   
}





sub create_file_header{
    
   
    print DUMP "##FORMAT (resequencing)\n";
    print DUMP "##DATE ",scalar(localtime),"\n";
    print DUMP "##RELEASE ",Bio::EnsEMBL::Registry->software_version(),"\n\n";

}

sub usage{
    my $msg = shift;

    print STDERR <<EOF;

usage: perl dump_strain_seq.pl <options>

options:
    -dump_file <filename>    file where you want to dump the resequencing data
    -species   <species>     species you want to dump data (default = mouse)

EOF
   
die ("\n$msg\n\n");
}
