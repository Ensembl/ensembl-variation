#!/usr/local/ensembl/bin/perl

use strict;
use warnings;

use Getopt::Long;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use ImportUtils qw(load);
use FindBin qw( $Bin );

use constant MAX_SIZE => 30_000_000;

my ($TMP_DIR, $TMP_FILE, $species);


GetOptions('species=s' => \$species,
	   'tmpdir=s'  => \$ImportUtils::TMP_DIR,
	   'tmpfile=s' => \$ImportUtils::TMP_FILE);

warn("Make sure you have a updated ensembl.registry file!\n");

my $registry_file ||= $Bin . "/ensembl.registry";

##
## bsub -q long -W12:00 -o [tmpdir]/output_tag.txt perl tag_snps.pl -tmpdir [tmpdir]/tag_snps -tmpfile tag_snps.txt
##
Bio::EnsEMBL::Registry->load_all( $registry_file );

#added default options
$species ||= 'human';

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

my $dbVariation = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');

my $population_id;
my $population_name;
#get all populations to be tagged (HapMap and PerlEgen)
my $sth = $dbVariation->dbc->prepare(qq{SELECT s.sample_id, s.name
					    FROM population p, sample s
					    WHERE (s.name like 'PERLEGEN:AFD%'
					    OR s.name like 'CSHL-HAPMAP%')
					    AND s.sample_id = p.sample_id
					});

$sth->execute();
$sth->bind_columns(\$population_id,\$population_name);
#get all the children that we do not want in the genotypes
my $siblings = {}; # hash {$individual_id} ,where the individual is sibling of another one
my @pops;
while($sth->fetch){
    if($population_name =~ /CEU|YRI|MEX/){
	&get_siblings($dbVariation,$population_id,$siblings);
    }
    push @pops, $population_id;
}

my $in_str = " IN (" . join(',', @pops). ")";

#select all SNPs in the population given with genotypes, and order them by variation
$sth = $dbVariation->dbc->prepare(qq{SELECT STRAIGHT_JOIN c.genotypes, c.seq_region_id, 
				            c.seq_region_start, ip.individual_sample_id, ip.population_sample_id
						FROM compressed_genotype_single_bp c FORCE INDEX(pos_idx), individual_population ip
						WHERE c.sample_id = ip.individual_sample_id
						AND   ip.population_sample_id $in_str
						ORDER BY c.seq_region_id,c.seq_region_start
#						LIMIT 10000
					    },{mysql_use_result=>1});

my ($genotypes, $seq_region_id, $seq_region_start,$individual_id);
my $buffer = {}; #buffer with data to write to file
my $blob; # contains the genotype information compressed
my @genotypes = ();
my $seq_region = {};
$sth->execute();
$sth->bind_columns(\$genotypes, \$seq_region_id, \$seq_region_start,\$individual_id,\$population_id);
while ($sth->fetch()){
    $seq_region->{$seq_region_id}++; #store the seq_region_id
    #only print genotypes without parents genotyped
    if (!exists $siblings->{$population_id.'-'.$individual_id}){
	#uncompress the genotypes
	$blob = substr($genotypes,2);
	@genotypes = unpack("naa" x (length($blob)/4),$blob);

	unshift @genotypes, substr($genotypes,1,1); #add the second allele of the first genotype
	unshift @genotypes, substr($genotypes,0,1); #add the first allele of the first genotype
	unshift @genotypes, 0; #the first SNP is in the position indicated by the seq_region1
	my ($snp_start, $allele_1, $allele_2);
	#print the information to the file
	for (my $i=0; $i < @genotypes -1;$i+=3){
	    if ($genotypes[$i] == 0){
		$snp_start = $seq_region_start; #first SNP is in the beginning of the region
	    }
	    else{
		$snp_start += $genotypes[$i] +1;
	    }
	    #genotype
	    $allele_1 = $genotypes[$i+1];
	    $allele_2 = $genotypes[$i+2];
	    
	    #print the data to the file with the information in the chromosome
	    print_buffered($buffer,"$TMP_DIR/tag_snps_$population_id\-$seq_region_id\.txt",
			   join("\t", $allele_1, $allele_2, $snp_start) . "\n");
	}
    }
}
$sth->finish();
print_buffered($buffer);
my $call;
my @stats;
#when the data is selected in the files (one per chromosome), send the jobs to the farm
my $files_size = {};
foreach my $file (keys %{$buffer}){
    @stats = stat($file);
    $files_size->{$file} = $stats[7];
}
foreach my $file (sort {$files_size->{$b} <=> $files_size->{$a}} keys %{$files_size}){
    $file =~ /tag_snps_(\d+)\-.*/; #extract the population_id from the name of the file
    if ($files_size->{$file} > MAX_SIZE){
	$call = "bsub -W2:00 -J tag_snps_$1 -R'select[myens_genomics1 < 200] rusage[myens_genomics1=10]' ";
    }
    $call .= "/usr/local/ensembl/bin/perl select_tag_snps.pl $file ";
    $call .= " -tmpdir $TMP_DIR -population_id $1 -species $species";
    
#   print $call,"\n";
    system($call); #send the job to one queue
    $call = '';
    
}

 $call = "bsub -K -w 'done(tag_snps*)' -J waiting_process sleep 1"; #waits until all snp_tagging have finished to continue
 system($call);

&import_tagg_snps($dbVariation);    

#for a given population, gets all individuals that are children (have father or mother)
sub get_siblings{
    my $dbVariation = shift;
    my $population_id = shift;
    my $siblings = shift;

    my $sth_individual = $dbVariation->dbc()->prepare(qq{SELECT i.sample_id
							     FROM individual i, individual_population ip
							     WHERE ip.individual_sample_id = i.sample_id
							     AND ip.population_sample_id = ? 
							     AND i.father_individual_sample_id IS NOT NULL
							     AND i.mother_individual_sample_id IS NOT NULL
							 });
    my ($individual_id);
    $sth_individual->execute($population_id);
    $sth_individual->bind_columns(\$individual_id);
    while ($sth_individual->fetch){
	$siblings->{$population_id.'-'.$individual_id}++; #necessary to have in the key the population, since some individuals are shared between
	                                                   #populations
    }
    return $siblings;
}

sub print_buffered {
    my $buffer = shift;
    my $filename = shift;
    my $text = shift;

    local *FH;

    if( ! $filename ) {
	# flush the buffer
	foreach my $file (keys %{$buffer}){
	    open( FH, ">>$file" ) or die;
	    print FH $buffer->{ $file };
	    close FH;
	}

    } else {
	$buffer->{ $filename } .= $text;
	if( length( $buffer->{ $filename } ) > 10_000 ) {
	    open( FH, ">>$filename" ) or die;
	    print FH $buffer->{ $filename };
	    close FH;
	    $buffer->{ $filename } = '';
	}
    }
}

sub import_tagg_snps{
    my $dbVariation = shift;
    
    #group all the fragments in 1 file
    my $call = "cat $TMP_DIR/snps_tagged_*.txt > $TMP_DIR/$TMP_FILE";
    system($call);

    unlink(<$TMP_DIR/snps_tagged_*>);
    #and finally, load the information
    load($dbVariation->dbc(), qw(tagged_variation_feature variation_feature_id sample_id));
#    unlink(<$TMP_DIR/tag_snps_*>);
   
}
sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: perl tag_snps.pl <options>

options:
    -species <name>      species name (default:human)
    -tmpdir <dir>        temp directory to use (with lots of space!)
    -tmpfile <filename>  name of temp file to use

EOF

  die("\n$msg\n\n");
}
