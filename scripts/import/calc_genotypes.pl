#!/usr/local/bin/perl

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
two 'AB' haplotypes; cell 1 has one 'AB' and one 'aB, etc), _except_ the double
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
     D' = D/Dmax, etc.

=cut

# expect marker format to be
# locus chr position person-id genotype variation_feature_id population_id
# Genotype must be AA, Aa, aa

use strict;
use Data::Dumper;
use Benchmark;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);


my %genotypes;
my %snps_ordered; #hash $snps_ordered{$chr}{$position} = $snp
my %additional_info; #to store additional information to have in the database: variation_feature_id

#assuming the script is called calc_genotypes.pl genotypes.txt genotypes_output.txt
 open FH, ">>$ARGV[1]"
    or throw("Could not open tmp file: $ARGV[1]\n");
open IN, "$ARGV[0]"
 or throw("Could not open input file: $ARGV[0]\n");

my $seq_region_id_previous = 0;

while( <IN> ) {
    chomp;
    /^#/ && next;
    my($locus,$position,$personid,$genotype,$variation_feature_id,$seq_region_id,$seq_region_end,$population_id) = split;
    if( !defined $genotype) {
      warn "bad format for locus, $_, skipping";
    }

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
    if ($seq_region_id != $seq_region_id_previous && $seq_region_id_previous != 0){
       &calculate_ld(\%snps_ordered,\%genotypes,\%additional_info,$seq_region_id_previous);  
       $seq_region_id_previous = $seq_region_id;
       %snps_ordered = ();
       %genotypes = ();
       %additional_info = ();
    }
    if ($seq_region_id_previous == 0){$seq_region_id_previous = $seq_region_id};
    #improves version: substituting 2 hashes, one for snps in chromosome and the other for snps in position for 1 hash: $snps_ordered{chr}{position} = sn
     if (defined $snps_ordered{$seq_region_id}{$position} && $snps_ordered{$seq_region_id}{$position} != $locus){
           warn "position $seq_region_id,$position for $locus disagrees with earlier position";
     }
     else{
         $snps_ordered{$seq_region_id}{$position} = $locus;
         #this is to store additional informaiton
         $additional_info{$locus}->{position} = $position;
         $additional_info{$locus}->{variation_feature_id} = $variation_feature_id;
         $additional_info{$locus}->{seq_region_end} = $seq_region_end;
         $additional_info{$locus}->{population_id} = $population_id;


         my $h = {};
         $h->{'locus'} = $locus;
         $h->{'personid'} = $personid;
         $h->{'genotype'} = $genotype;
    
        if( !defined $genotypes{$locus} ) {
          $genotypes{$locus} = [];
         }

        push(@{$genotypes{$locus}},$h);    
    }
}

close IN;
close FH;

sub calculate_ld{
    my $snps_ordered = shift;
    my $genotypes = shift;
    my $additional_info = shift;
    my $seq_region_id = shift;

    my %output;
    my $snp_count = 0; #to count the number of snps between the 2  in the region
    #my approach to the problem: have a hash by chromosome and position and compare things    
    #then, sort all the positions
    my @positions_ordered = (sort {$a <=> $b} keys %{$snps_ordered{$seq_region_id}});
#    while ( @positions_ordered ) { 
#     my $position = shift( @positions_ordered );
#     foreach my $position2 ( @positions_ordered ) {
#       compare position with position2

    foreach my $position (sort {$a <=> $b} keys %{$snps_ordered->{$seq_region_id}}){
	#delete the own snp to not copy it again
	shift @positions_ordered;#remove first element from the array, it has already been compared
	    $snp_count = 0;
	#compare against all the SNPs, until we find one outside the window
	foreach my $position2 (@positions_ordered){
	    last if (abs($position - $position2) > 100_000);
	    my $key = "$snps_ordered->{$seq_region_id}{$position}:$snps_ordered->{$seq_region_id}{$position2}";
	    $snp_count++;
	    if( &shared_people($genotypes->{$snps_ordered->{$seq_region_id}{$position}},$genotypes->{$snps_ordered->{$seq_region_id}{$position2}}) < 20 ) {
#               print "skipping $key as not enough shared people\n";
#		next;
	    }
	    
	    my $stats_hash = &calculate_pairwise_stats($genotypes->{$snps_ordered->{$seq_region_id}{$position}},$genotypes->{$snps_ordered->{$seq_region_id}{$position2}},$snp_count);
	    
	    $output{$key} = $stats_hash;

	}
    }
    
    #and finally, print the ld calculation for the chromosome into a file
    my $locus1;
    my $locus2;
    my $seq_region_start;
    my $seq_region_end;
    foreach my $key ( keys %output ) {
	($locus1,$locus2) = split /:/,$key;
#format of the output file
#  variation_feature_id_1 variation_feature_id_2 population_id seq_region_id seq_region_start seq_region_end snp_distance_count r2 Dprime samplesize
#it will be necessary to find out the start and end of the region, that will correspond to the start of the region for the variation in the lowest position and the end of the region for the variation in the higher position
	if ($additional_info->{$locus1}->{position} < $additional_info->{$locus2}->{position}){
	    $seq_region_start = $additional_info->{$locus1}->{position};
	}
	else{
	    $seq_region_start = $additional_info->{$locus2}->{position};
	}
	if ($additional_info->{$locus1}->{seq_region_end} < $additional_info->{$locus2}->{seq_region_end}){
	    $seq_region_end = $additional_info->{$locus2}->{seq_region_end};
	}
	else{
	    $seq_region_end = $additional_info->{$locus1}->{seq_region_end};
	}
	print FH join("\t",$additional_info->{$locus1}->{variation_feature_id},$additional_info->{$locus2}->{variation_feature_id},$additional_info->{$locus1}->{population_id},$seq_region_id,$seq_region_start,$seq_region_end,$output{$key}->{'snp_count'},$output{$key}->{'r2'}, abs($output{$key}->{'Dprime'}),$output{$key}->{'N'}),"\n";	
    }   
}


sub calculate_pairwise_stats {
    my $first  = shift;
    my $second = shift;
    my $snp_count = shift;

    my ($AABB,$AABb,$AAbb,$AaBB,$AaBb,$Aabb,$aaBB,$aaBb,$aabb);
    $AABB = $AABb = $AAbb = $AaBB = $AaBb = $Aabb = $aaBB = $aaBb = $aabb = 0;

    my %people;

    # double loop potentially expensive. Could set up reverse lookup
    # hashes reading in, but would probably cause memory loops.
    my $count = 0;
    foreach my $g ( @{$first} ) {
	foreach my $h ( @{$second} ) {
	    if( $g->{'personid'} ne $h->{'personid'} ) {
		next;
	    }
#	    print "Doing a new person ",$g->{'personid'},"$count\n";
	    $count++;
	    $people{$g->{'personid'}}  = 1;
	    # now assign to double genotype.
	    if( $g->{'genotype'} eq 'AA' ) {
		if( $h->{'genotype'} eq 'AA' ) {
		    $AABB++;
		} elsif ( $h->{'genotype'} eq 'Aa' ) {
		    $AABb++;
		} else {
		    $AAbb++;
		}
	    } elsif ( $g->{'genotype'} eq 'Aa' ) {
		if( $h->{'genotype'} eq 'AA' ) {
		    $AaBB++;
		} elsif ( $h->{'genotype'} eq 'Aa' ) {
		    $AaBb++;
		} else {
		    $Aabb++;
		}		
	    } else { 
		# is aa
		if( $h->{'genotype'} eq 'AA' ) {
		    $aaBB++;
		} elsif ( $h->{'genotype'} eq 'Aa' ) {
		    $aaBb++;
		} else {
		    $aabb++;
		}
	    }

	}
    }

    my $theta = 0.5;
    my $thetaprev = 2;
    my ($pAB,$pAb,$paB,$pab);
    my ($nAB,$nAb,$naB,$nab);

#    print  "Table: $AABB  $AaBB $aaBB\n";
#    print  "Table: $AABb  $AaBb $aaBb\n";
#    print  "Table: $AAbb  $Aabb $aabb\n";

    $pAB = $pAb = $paB = $pab = 0.25;

    $nAB = 2*$AABB+$AaBB+$AABb;
    $nab = 2*$aabb+$Aabb+$aaBb;
    $nAb = 2*$AAbb+$Aabb+$AABb;
    $naB = 2*$aaBB+$AaBB+$aaBb;

#    my $N = $nAB + $nab + $nAB + $naB + 2*$AaBb;
    my $N = $nAB + $nab + $nAb + $naB + 2*$AaBb;

    while( abs($theta-$thetaprev) > 0.0001 ) {
	$thetaprev = $theta;
	if ($N != 0){
	    # restimate probabilities of haplotypes using
	    # theta to split out the double het
	    $pAB = ($nAB + (1-$theta)*$AaBb)/$N;
	    $pab = ($nab + (1-$theta)*$AaBb)/$N;
	    $pAb = ($nAb + $theta*$AaBb)/$N;
	    $paB = ($naB + $theta*$AaBb)/$N;
	}
 #       print "...estimate is $pAB $pab $pAb $paB\n";
	eval{
	    $theta = (($pAb*$paB)/($pAB*$pab + $pAb*$paB));
	};
	if ($@){$theta = 0;} #included to avoid division by 0
#	print "New theta $theta\n";
    }
    
    

    # now calculate stats

    my $f_A = &major_frequency($first,\%people);
    my $f_B = &major_frequency($second,\%people);
    my $D  = $pAB - ($f_A*$f_B);
    my $r2;
    eval{
	$r2 = $D*$D/($f_A*$f_B*(1-$f_A)*(1-$f_B)); 
    };
    if ($@){
	$r2 = 0; #for some cases is not possible to calculate the r2 due to a 0 in the divisor
    }
    
    my $Dmax = 0;
    my $Dprima;
    
    if ($D < 0){
	$Dmax = $f_A*$f_B if ($f_A*$f_B < ((1-$f_A)*(1-$f_B)));
	$Dmax = (1-$f_A)*(1-$f_B) if ($f_A*$f_B >= ((1-$f_A)*(1-$f_B)));	
    }
    if ($D > 0){
	$Dmax = $f_A*(1-$f_B) if ($f_A*(1-$f_B) < (1-$f_A)*$f_B);
	$Dmax = (1-$f_A)*$f_B if ($f_A*(1-$f_B) >= (1-$f_A)*$f_B);
    }
    eval{
	$Dprima = $D/$Dmax;
    };
    if ($@){
	$Dprima = 0;
    }

    my $o = {};

    $o->{'D'}= $D;
    $o->{'r2'} = $r2;
    $o->{'theta'} = $theta;
    $o->{'N'} = $N;
    $o->{'Dprime'} = $Dprima;
    $o->{'snp_count'} = $snp_count;
    $o->{'people'} = scalar(keys %people);

    return $o;
}

sub major_frequency {
    my ($genotypes,$people) = @_;
    
    my $major = 0;
    my $total = 0;
    foreach my $g ( @{$genotypes} ) {
	if( exists $people->{$g->{'personid'}} && $people->{$g->{'personid'}} == 1 ) {
	    if( $g->{'genotype'} eq 'AA' ) {
		$major += 2;
	    } elsif ( $g->{'genotype'} eq 'Aa' ) {
		$major++;
	    }

	    $total += 2;
	}
    }

    return $major/$total if ($total !=0);
    return 0 if ($total == 0);

}


sub shared_people {
    my $first  = shift;
    my $second = shift;
    my %h;
    my $c = 0;
    #new way to find number repeated elements in 2 arrays
    foreach my $elem (@{$first},@{$second}){
	$h{$elem->{'personid'}}++;
    }
    foreach my $elem (keys %h){
	if ($h{$elem} > 1){
	    $c++;
	}
    }
    return $c;
}
