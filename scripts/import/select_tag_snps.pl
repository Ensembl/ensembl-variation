#!/usr/local/ensembl/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;

my ($TMP_DIR);

my ($vhost, $vport, $vdbname, $vuser, $vpass,
    $population_id);

GetOptions('vhost=s'   => \$vhost,
	   'vuser=s'   => \$vuser,
	   'vpass=s'   => \$vpass,
	   'vport=i'   => \$vport,
	   'vdbname=s' => \$vdbname,
	   'tmpdir=s'  => \$TMP_DIR,
	   'population_id=i' => \$population_id);


my $dbVariation = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new
    (-host => $vhost,
     -user => $vuser,
     -pass => $vpass,
     -port => $vport,
     -dbname => $vdbname
     );

#first, calculate the MAF for the SNPs in the chromosome
my $file = $ARGV[0] if (defined @ARGV);
die "Not possible to calculate SNP tagging without file with SNPs" if (!defined @ARGV);

open IN, "$file" or die("Could not open input file: $ARGV[0]\n");

my $r2 = 0.99; #default value for r2

my $previous_variation_feature = 0;
my $previous_seq_region_start;
my ($variation_feature_id, $allele_1, $allele_2, $seq_region_start);
my $genotypes_snp = {};
my $last_snp = 0; #to indicate if there is a last snp that needs the MAF calculation
my $MAF_snps = {};
while (<IN>){
    chomp;
    /^\#/ && next;
    ($variation_feature_id, $allele_1, $allele_2, $seq_region_start) = split;
    #when we change variation_feature, calculate the MAF
    if (($previous_variation_feature != $variation_feature_id) && $previous_variation_feature != 0){
	$MAF_snps->{$previous_variation_feature} = &calculate_MAF($genotypes_snp) . '-' . $previous_seq_region_start;
	$genotypes_snp = {};
	$previous_variation_feature = $variation_feature_id;
	$previous_seq_region_start = $seq_region_start;
	$last_snp = 1;

    }
    if ($previous_variation_feature == 0){ #initialize the variables
	$previous_variation_feature = $variation_feature_id;
	$previous_seq_region_start = $seq_region_start;
    }
    $genotypes_snp->{$allele_1}++ if ($allele_1 ne 'N');
    $genotypes_snp->{$allele_2}++ if ($allele_2 ne 'N');    
    $last_snp = 0;
}
close IN;
#calculate the MAF for the last SNP, if necessary
$MAF_snps->{$previous_variation_feature} = &calculate_MAF($genotypes_snp) . '-' . $previous_seq_region_start if ($last_snp == 0);



#get LD values for the chromosome
$file =~/.*:(\d+)\.txt/; #extract the seq_region_id from the name of the file
my $LD_values = &get_LD_chromosome($dbVariation,$1,$r2,$population_id);
#do the algorithm
my $remove_snps = {}; #hash containing the snps that must be removed from the entry, they have been ruled out
foreach $variation_feature_id (sort {$MAF_snps->{$a} cmp $MAF_snps->{$b}} keys %{$MAF_snps}){
    if (!defined $remove_snps->{$variation_feature_id}){
	#add the SNPs that should be removed in future iterations
	#and delete from the hash with the MAF_snps, the ones that have a r2 greater than r2 with $variation_feature_id
	map {$remove_snps->{$_}++;delete $MAF_snps->{$_}} @{$LD_values->{$variation_feature_id}};
    }
}
my $host = `hostname`;
chop $host;
open OUT, ">$TMP_DIR/snps_tagged_$population_id\_$host\:$$\.txt" or die ("Could not open output file");
map {print OUT join("\t",$_,$population_id),"\n"} keys %{$MAF_snps};
close OUT or die ("Could not close output file with tagged SNPs");


#for a list of genotypes, get the MAF
sub calculate_MAF{
    my $genotypes_snp = shift;
    my $MAF;
    my $total_genotypes = 0;
    my $allele_freq;

    if (keys %{$genotypes_snp} == 2){
	foreach (values %{$genotypes_snp}){
	    $total_genotypes += $_;
	    $allele_freq = ($_ / $total_genotypes) if ($total_genotypes != $_);
	}
	return $allele_freq if ($allele_freq <= (1 - $allele_freq));
	return 1 - $allele_freq if ($allele_freq > (1 - $allele_freq));
    }    
    return 0;
#    die "genotype with more than 2 alleles!!";
}

#creates a hash with all the variation_features in the chromosome with a r2 greater than r2
sub get_LD_chromosome{
    my $dbVariation = shift;
    my $seq_region_id = shift;
    my $r2 = shift;
    my $population_id = shift;
    
    my $variation_features = {};
    my $sth = $dbVariation->dbc->prepare(qq{SELECT variation_feature_id_1, variation_feature_id_2
						FROM pairwise_ld
						WHERE seq_region_id = ?
						AND r2 > ?
						AND population_id = ?
					    },{mysql_use_result =>1});
    $sth->execute($seq_region_id,$r2,$population_id);
    my ($variation_feature_id_1, $variation_feature_id_2);
    $sth->bind_columns(\$variation_feature_id_1, \$variation_feature_id_2);
    while ($sth->fetch()){
	push @{$variation_features->{$variation_feature_id_1}}, $variation_feature_id_2;
	push @{$variation_features->{$variation_feature_id_2}}, $variation_feature_id_1;
    }
    $sth->finish();
    return $variation_features;
}
