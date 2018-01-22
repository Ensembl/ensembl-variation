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


=head1 Docs

From Lon Cardon

The EM algorithm is used to estimate the recombination fraction (parameter theta)
between two genetic markers under the assumption of Hardy-Weinberg equilibrium.

For SNPs, consider the 3 x 3 table of marker 1 with alleles 'A' and 'a' 
and marker 2 with alleles 'B' and 'b'.  Number the cells as follows:

      AA  |  Aa  |  aa
    -------------------
 BB |  0  |   1  |  2
 Bb |  3  |   4  |  5
 bb |  6  |   7  |  8
  
Note that all cells arise from unique haplotypes (e.g., cell 0 can only comprise
two 'AB' haplotypes; cell 1 has one 'AB' and one 'aB', etc), _except_ the double
heterozygote cell 4, which can have either AB/ab or Ab/aB.  We need to estimate the
probability of these two events.  We assume that pair AB/ab arises when no
recombination occurs (1 - theta) whereas pair Ab/aB arises in the presence of recombination
(probability theta).  Thus, for the four haplotype possibilities at two markers 
(AB,ab,Ab,aB), we have:

 prAB=prAb=praB=prab=.25;  
  
  nAB=(float)(2*cells[0]+cells[1]+cells[3]); 
  nab=(float)(2*cells[8]+cells[7]+cells[5]);
  nAb=(float)(2*cells[6]+cells[7]+cells[3]);
  naB=(float)(2*cells[2]+cells[1]+cells[5]);

  N = nAB + nab + nAb + naB + 2*cell4;

  while(theta-thetaprev > CONVERGENCE_CRITERIA) 
   {
    thetaprev=theta;
    prAB=(nAB + (1-theta)*cells[4])/N;
    prab=(nab + (1-theta)*cells[4])/N;
    prAb=(nAb + theta*cells[4])/N;
    praB=(naB + theta*cells[4])/N;
    theta=(prAb*praB)/(prAB*prab + prAb*praB); 
  }


i.e., D = prAB - frq(A)*frq(B), r^2 = D^2/(frq(A)*frq(B)*frq(a)*frq(b)),
     Dprime = D/Dmax, etc.

=cut

# expect marker format to be
# locus chr position person-id genotype variation_feature_id population_id
# Genotype must be AA, Aa, aa

use strict;
use warnings;

use constant WINDOW_SIZE => 100_000;
use Bio::EnsEMBL::Utils::Cache;

my %ld_cache;
tie(%ld_cache, 'Bio::EnsEMBL::Utils::Cache', 1_000_000 );

my %genotypes;
my %snps_ordered; #hash $snps_ordered{$chr}{$position} = $snp
my %additional_info; #to store additional information to have in the database: variation_feature_id
my @positions_ordered; #to keep the order of the markers
#assuming the script is called calc_genotypes.pl genotypes.txt genotypes_output.txt


 open FH, ">>$ARGV[1]"
    or die("Could not open tmp file: $ARGV[1]\n");
open IN, "$ARGV[0]"
 or die("Could not open input file: $ARGV[0]\n");

my $seq_region_id_previous = -1;
my $last_position = -1;

my $position_removed; #position where it has been deleted the marker. Remove it from the rest of the hashes
my %iterations_stats;
my $iterations = 0;

while( <IN> ) {
  chomp;
  /^\#/ && next;
  my($locus,$position,$personid,$genotype,$variation_feature_id,
     $seq_region_id,$seq_region_end,$population_id) = split;
  
  if( $position !~ /\d+/ ) {
    die "bad position $position";
  }

  if( $genotype eq 'aA' ) {
    # unphased
    $genotype = 'Aa';
  }

  if( $genotype !~ /AA|Aa|aa/ ) {
    die "Genotype must be AA, Aa or aa, not $genotype";
  }
    

  #new region, flush all structures and calculate the ld information for the remaining markers
  if ($seq_region_id != $seq_region_id_previous ){
    #calculate the LD for the remaining markers in the region
    while (@positions_ordered){
      &calculate_ld(\%snps_ordered,\%genotypes,
		    \%additional_info,\@positions_ordered,
		    $seq_region_id_previous);  
    }
    $seq_region_id_previous = $seq_region_id;
    $last_position = -1;
    %snps_ordered = ();
    %genotypes = ();
    %additional_info = ();
    @positions_ordered = ();
  }

  my $h = {'locus' => $locus,
	   'personid' => $personid,
	   'genotype' => $genotype };
    
  $genotypes{$locus} ||= [];
  push(@{$genotypes{$locus}},$h);    

  if( $position == $last_position ) {
    next;
  }

  $last_position = $position;

  # check if the position is farther than the limit.
  # if so, calculate the ld information for the values in the array

  while ( @positions_ordered && (abs($positions_ordered[0]-$position) > WINDOW_SIZE())){
    $position_removed = &calculate_ld(\%snps_ordered,\%genotypes,
				      \%additional_info,\@positions_ordered,
				      $seq_region_id_previous);
    #necessary to remove the information from the different hashes already calculated
    delete $genotypes{$snps_ordered{$position_removed}};
    delete $additional_info{$snps_ordered{$position_removed}};
    delete $snps_ordered{$position_removed};
  }

  #add the position of the marker to the array if a new position
  push( @positions_ordered, $position );

  #store the information necessary to retrieve all data from the marker
  $snps_ordered{$position} = $locus;

  #this is to store additional informaiton
  $additional_info{$locus}->{position} = $position;
  $additional_info{$locus}->{variation_feature_id} = $variation_feature_id;
  $additional_info{$locus}->{seq_region_end} = $seq_region_end;
  $additional_info{$locus}->{population_id} = $population_id;
}


#calculate LD for the last region
while (@positions_ordered){
   &calculate_ld(\%snps_ordered,\%genotypes,
		 \%additional_info,\@positions_ordered,
		 $seq_region_id_previous);  
}

#print STDERR "Iterations: $iterations\n";

close IN;
close FH;


sub calculate_ld{

  my $snps_ordered = shift;
  my $genotypes = shift;
  my $additional_info = shift;
  my $positions_ordered = shift;
  my $seq_region_id = shift;
  
#  my $snp_count = 0; #to count the number of snps between the 2  in the region

  my $position = shift @{$positions_ordered}; #remove first element from the array, it has already been compared
  #calculate LD against all SNPs in the window
  my $locus1 = $snps_ordered->{$position};
  my $seq_region_start = $additional_info->{$locus1}->{position};
  
  my $locus2;
  my $seq_region_end;

  # format of the output file
  #  variation_feature_id_1 variation_feature_id_2 population_id seq_region_id seq_region_start 
  #   seq_region_end snp_distance_count r2 d_prime samplesize

  foreach my $position2 (@{$positions_ordered}){
    
    $locus2 = $snps_ordered->{$position2};
    $seq_region_end = $additional_info->{$locus2}->{seq_region_end};
    
    my $stats_hash = &calculate_pairwise_stats( $genotypes->{$snps_ordered->{$position}},
						$genotypes->{$snps_ordered->{$position2}} );
    

    next unless ($stats_hash->{'r2'} >= 0.05 && $stats_hash->{'N'} >= 40); #cut-off in 5% and at least 20 individuals genotyped
#    $snp_count++;

    print FH join("\t",
		  $additional_info->{$locus1}->{variation_feature_id},
		  $additional_info->{$locus2}->{variation_feature_id},
		  $additional_info->{$locus1}->{population_id},
		  $seq_region_id, $seq_region_start, $seq_region_end,
#		  $snp_count,
		  $stats_hash->{'r2'}, 
		  abs($stats_hash->{'d_prime'}),
		  $stats_hash->{'N'}),"\n";	
  }

  return $position;
}


sub calculate_pairwise_stats {
    my $first  = shift;
    my $second = shift;


    my %genotypes = ( 'AABB' => 0,
		      'AaBB' => 0,
		      'aaBB' => 0,

		      'AABb' => 0,
		      'AaBb' => 0,
		      'aaBb' => 0,

		      'AAbb' => 0,
		      'Aabb' => 0,
		      'aabb' => 0 );


    my %people;

    foreach my $g ( @{$first} ) {
      $people{$g->{'personid'}} = $g->{'genotype'};
    }

    foreach my $h ( @{$second} ) {
      if( exists $people{$h->{'personid'}} ) {
	my $hg = $h->{'genotype'};
	$hg =~ tr/Aa/Bb/;
	$people{$h->{'personid'}} .= $hg;
      }
    }

    for my $pers_id ( keys %people ) {
      my $genotype = $people{ $pers_id };
      if( length( $genotype ) < 4 ) {
	delete $people{$pers_id};
      } else {
	$genotypes{ $genotype }++;
      }
    }


    # my ($pAB,$pAb,$paB,$pab);
    my ($nAB,$nAb,$naB,$nab);

    #    print  "Table: $AABB  $AaBB $aaBB\n";
#    print  "Table: $AABb  $AaBb $aaBb\n";
#    print  "Table: $AAbb  $Aabb $aabb\n";


    $nAB = 2*$genotypes{'AABB'}+$genotypes{'AaBB'}+$genotypes{'AABb'};
    $nab = 2*$genotypes{'aabb'}+$genotypes{'Aabb'}+$genotypes{'aaBb'};
    $nAb = 2*$genotypes{'AAbb'}+$genotypes{'Aabb'}+$genotypes{'AABb'};
    $naB = 2*$genotypes{'aaBB'}+$genotypes{'AaBB'}+$genotypes{'aaBb'};
    my $AaBb = $genotypes{'AaBb'};

    my $ld_cache_key = join( "-", $nab,$naB, $nAb, $nAB, $AaBb );
    if( exists $ld_cache{ $ld_cache_key } ) {
      return $ld_cache{$ld_cache_key};
    }

    my $theta = 0.5;
    my $thetaprev = 2;

#    print  "Table: $AABB  $AaBB $aaBB\n";
#    print  "Table: $AABb  $AaBb $aaBb\n";
#    print  "Table: $AAbb  $Aabb $aabb\n";

#    $pAB = $pAb = $paB = $pab = 0.25;

    my $N = $nAB + $nab + $nAb + $naB + 2*$AaBb;
    
   
    while( abs($theta-$thetaprev) > 0.0001 ) {	
	$thetaprev = $theta;
#	if ($N != 0){
	  # restimate probabilities of haplotypes using
	  # theta to split out the double het
	  # $pAB = ($nAB + (1-$theta)*$AaBb)/$N;
	  # $pab = ($nab + (1-$theta)*$AaBb)/$N;
	  # $pAb = ($nAb + $theta*$AaBb)/$N;
	  # $paB = ($naB + $theta*$AaBb)/$N;
#	}
 #       print "...estimate is $pAB $pab $pAb $paB\n";
	eval{
	  #  $theta = (($pAb*$paB)/($pAB*$pab + $pAb*$paB));
	  $theta = (($nAb + $theta*$AaBb)*($naB + $theta*$AaBb))/
		    (($nAB + (1-$theta)*$AaBb)*($nab + (1-$theta)*$AaBb) + 
		     ($nAb + $theta*$AaBb)*($naB + $theta*$AaBb));
	};
	if ($@){$theta = 0.5;} #included to avoid division by 0
#	print "New theta $theta\n";
	$iterations++;
    }

    # now calculate stats
    my ( $f_A, $f_B ) = major_freqs( \%people );
    
    # my $D  = $pAB - ($f_A*$f_B);
    my $D;
    my $r2;
    eval{
      $D  = ($nAB+(1-$theta)*$AaBb)/$N - ($f_A*$f_B);
      $r2 = $D*$D/($f_A*$f_B*(1-$f_A)*(1-$f_B)); 
    };

    if ($@){
      $D = 0;
      $r2 = 0; #for some cases is not possible to calculate the r2 due to a 0 in the divisor
    }
    
    my $Dmax = 0;
    my $d_prime;
    
    if ($D < 0){
	$Dmax = $f_A*$f_B if ($f_A*$f_B < ((1-$f_A)*(1-$f_B)));
	$Dmax = (1-$f_A)*(1-$f_B) if ($f_A*$f_B >= ((1-$f_A)*(1-$f_B)));	
    }
    if ($D > 0){
	$Dmax = $f_A*(1-$f_B) if ($f_A*(1-$f_B) < (1-$f_A)*$f_B);
	$Dmax = (1-$f_A)*$f_B if ($f_A*(1-$f_B) >= (1-$f_A)*$f_B);
    }
    eval{
	$d_prime = $D/$Dmax;
    };
    if ($@){
	$d_prime = 0;
    }

    my $o = { 'D'=> $D,
	      'r2' => $r2,
	      'theta' => $theta,
	      'N' => $N,
	      'd_prime' => $d_prime,
	      'people' => scalar(keys %people)};

    $ld_cache{ $ld_cache_key } = $o;

    return $o;

}


sub major_freqs {
  my $people = shift;
  my ( $f_a, $f_A, $f_b, $f_B, $total );

  for my $genotype ( values %$people ) {
    $f_a += ( $genotype =~ tr/a/a/ );
    $f_A += ( $genotype =~ tr/A/A/ );
    $f_b += ( $genotype =~ tr/b/b/ );
    $f_B += ( $genotype =~ tr/B/B/ );
    $total += 2;
  }

  if( ! $total ) {
    return ( 0,0 );
  }

  #just use $f_A and $f_B, no need to find bigger value
  #if( $f_a > $f_A ) {
  #  ( $f_a, $f_A ) = ( $f_A, $f_a );
  #}

  #if( $f_b > $f_B ) {
  #  ( $f_b, $f_B ) = ( $f_B, $f_b );
  #}

  return ( $f_A/$total, $f_B/$total );
}
