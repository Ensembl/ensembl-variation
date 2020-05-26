#!/usr/bin/env perl

## extract genotypes relative to the reference
## from http://data.nextstrain.org/ncov.json
use strict;
use warnings;
use JSON;

die "\nJSON input needed\n\n" unless scalar(@ARGV) ==1;

local $/ = undef;

open my $in, $ARGV[0] ||die "Failed to open file to read $!\n";
open my $out, ">", "$ARGV[0]\_variants.txt" ||die "Failed to open file to write $!\n";

while(<$in>) {
  chomp $_;
  my $h = decode_json($_);

  my $tr = $h->{tree}->{children};

  foreach my $ch (@{$tr}){

    ## reset
    my %var;
    ## seed with a reference match to save ref-matching samples
    $var{149} = "nuc\t149\tG149G\tG\tG\tORF14\tT39I";

    parse_children($ch, "top", \%var );
  }
}

# print each variant to a line with the sample name
# check for other children & resubmit to check for more children
sub parse_children{

  my $ch = shift;
  my $parent = shift;      ## keeping for checking purposes
  my $inherit_var = shift; ## may be undef - hash of variant locations from parent

  #for readability
  my $clade = $ch->{node_attrs}->{clade_membership}->{value};
  my $muts = $ch->{branch_attrs}->{mutations};

  # get the gene name 
  # used as a key for protein level changes 
  my $gene = ".";
  foreach my $loc (keys %{$muts}){   ##nuc & ORF 
    $gene = $loc unless $loc =~ /nuc/;
  }
  
  # build list of variants to pass on
  my %changed_pos;

  if($muts->{nuc}){
    ## loop through mutations linking nuc & gene level names
    for( my $n =0; $n < scalar(@{$muts->{nuc}}) ; $n++ ){

      my $var  = $muts->{nuc}->[$n];
      my $prot = $muts->{$gene}->[$n] || "."; 
      my $pos = $var;
      $pos =~ s/\D+//g;
      my ($ref, $alt) = split/\d+/,$var;

      ## save for children
      ## by position as reference may change (is relative to parent)
      $changed_pos{$pos} = "nuc\t$pos\t$var\t$ref\t$alt\t$gene\t$prot";
 
      ##write out for this sample
      unless ($ch->{name} =~/travel_history|NODE/){
        print $out "$changed_pos{$pos}\t$ch->{name}\t$clade\t$parent\n";
      }
    }
  }

  ## write out parental genotypes unless site already seen
  foreach my $pos(keys %{$inherit_var}){
    ## don't write or pass on parent genotype if changed in child
    next if $changed_pos{$pos};

    ## save to pass to child
    $changed_pos{$pos} = $inherit_var->{$pos};

    ## write out, excluding predicted nodes
    unless( $ch->{name} =~/travel_history|NODE/){
      print $out $inherit_var->{$pos} . "\t$ch->{name}\t$clade\t$parent\n";
    }
  }


  ## there may be children   
  foreach my $child (@{$ch->{children}}){
    parse_children($child, $ch->{name}, \%changed_pos);
  }
}
