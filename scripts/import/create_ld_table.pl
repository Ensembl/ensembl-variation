#!/usr/local/ensembl/bin/perl
use strict;
use warnings;

use Getopt::Long;
use ImportUtils qw(load);
use Bio::EnsEMBL::Registry;
use FindBin qw( $Bin );


##
##  run with bsub -q long -W14:00 -o [tmpdir]/output_tag.txt perl create_ld_table.pl -tmpdir [tmpdir]
##  -tmpfile tag_snps.txt
use constant MAX_SIZE => 500_000_000;
my ($TMP_DIR, $TMP_FILE, $species);


GetOptions('species=s' => \$species,
	   'tmpdir=s'  => \$ImportUtils::TMP_DIR,
	   'tmpfile=s' => \$ImportUtils::TMP_FILE);

warn("Make sure you have a updated ensembl.registry file!\n");

my $registry_file ||= $Bin . "/ensembl.registry";

$TMP_DIR = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

#added default options
$species ||= 'human';
Bio::EnsEMBL::Registry->load_all( $registry_file );
my $dbVariation = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');

my $siblings;
create_LD_table($dbVariation); #create the LD table if it is not there yet
my $in_str = get_LD_populations($dbVariation,$siblings);

my $sth = $dbVariation->dbc->prepare(qq{SELECT c.sample_id,c.seq_region_id,c.seq_region_start,c.seq_region_end,c.genotypes,ip.population_sample_id
				 FROM compressed_genotype_single_bp c FORCE INDEX(pos_idx), individual_population ip
				 WHERE  ip.individual_sample_id = c.sample_id
				 AND   ip.population_sample_id $in_str
#				 LIMIT 100000
			     },{mysql_use_result => 1});

$sth->execute();
create_files($sth,$siblings);
$sth->finish();
my @files = glob("$TMP_DIR/*");
my $call = '';
my @stats;
foreach my $file (@files){
    $file =~ /.*_(\d+)/;
    @stats = stat($file);
    if ($stats[7] > MAX_SIZE){
#once the files are created, we have to calculate the ld
	$call = "bsub -W2:00 -M3500000 -R'select[mem>3500] rusage[mem=3500]' -q normal -J ld_calculation_$1 ";
    }
    $call .= "perl calculate_ld_table.pl -tmpdir $TMP_DIR -tmpfile $TMP_FILE -ldfile $file ";
  #  print $call,"\n";
    system($call);
    $call = '';
}
sleep(60);
$call = "bsub -K -W2:00 -w 'done(ld_calculation*)' -J waiting_ld date";
system($call);
#finally, import the data in the LD table
$call = "cat $TMP_DIR/tag_snps*.out > $TMP_DIR/$TMP_FILE";
system($call);
load($dbVariation->dbc,qw(pairwise_ld sample_id seq_region_id seq_region_start seq_region_end r2));
unlink(<$TMP_DIR/tag_snps*.out>);

sub store_file{
    my $buffer = shift;
    my $individual_id = shift;
    my $seq_region_start = shift;
    my $genotype = shift;
    my $population_id = shift;
    my $seq_region_id = shift;

    #get the first byte of the string, and unpack it (the genotype, without the gaps)
    my $blob = substr($genotype,2);
    #the array contains the uncompressed value of genotype, always in the format number_gaps . genotype		  
    my @genotypes = unpack("naa" x (length($blob)/4),$blob);
    unshift @genotypes, substr($genotype,1,1); #add the second allele of the first genotype
    unshift @genotypes, substr($genotype,0,1); #add the first allele of the first genotype
    unshift @genotypes, 0; #the first SNP is in the position indicated by the seq_region1
    my $snp_start;
    my $allele_1;
    my $allele_2;
    for (my $i=0; $i < @genotypes -1;$i+=3){
	#number of gaps
	if ($i == 0){
	    $snp_start = $seq_region_start; #first SNP is in the beginning of the region
	}
	else{
            #ignore when there is more than 1 genotype in the same position
	    if ($genotypes[$i] == 0){
		$snp_start += 1;
		next; 
	    }
	    $snp_start += $genotypes[$i] +1;
	}
	#genotype
	$allele_1 = $genotypes[$i+1];
	$allele_2 = $genotypes[$i+2];
	print_buffered($buffer,"$TMP_DIR/dump_data_$seq_region_id.txt",join("\t",$snp_start,$individual_id,$population_id,$allele_1,$allele_2)."\n");
    }
}


#creates all files, one per chromosome, to run the LD calculation
sub create_files{
    my $sth = shift;
    my $siblings = shift;
    
    my $buffer = {};
   
    my ($individual_id, $seq_region_id, $seq_region_start,$seq_region_end,$genotypes, $population_id);
    $sth->bind_columns(\$individual_id, \$seq_region_id, \$seq_region_start, \$seq_region_end, \$genotypes, \$population_id);
    while($sth->fetch()) {
	if (!exists $siblings->{$population_id . '-' . $individual_id}){ #necessary to use the population_id
	    store_file($buffer,$individual_id,$seq_region_start,$genotypes,$population_id,$seq_region_id);
	}
    }
    print_buffered($buffer);
    
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


sub get_LD_populations{
    my $dbVariation = shift;
    my $siblings = shift;
    my ($pop_id,$population_name);
    my $sth = $dbVariation->dbc->prepare(qq{SELECT s.sample_id, s.name
				     FROM population p, sample s
				     WHERE (s.name like 'PERLEGEN:AFD%'
				     OR s.name like 'CSHL-HAPMAP%')
				     AND s.sample_id = p.sample_id});

    $sth->execute();
    $sth->bind_columns(\$pop_id,\$population_name);
    #get all the children that we do not want in the genotypes
    my @pops;
    while($sth->fetch){
	if($population_name =~ /CEU|YRI|MEX/){
	    get_siblings($dbVariation,$pop_id,$siblings);
	}
	push @pops, $pop_id;
    }
    
    my $in_str = " IN (" . join(',', @pops). ")";
    
    return $in_str if (defined $pops[0]);
    return '' if (!defined $pops[0]);

}

#for a given population, gets all individuals that are children (have father or mother)
sub get_siblings{
    my $dbVariation = shift;
    my $population_id = shift;
    my $siblings = shift;

    my $sth_individual = $dbVariation->dbc->prepare(qq{SELECT i.sample_id
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
}


sub create_LD_table{
    my $dbVariation = shift;
    
    #need to create ld table if it has not been create yet
    $dbVariation->dbc->do(qq{
	CREATE TABLE IF NOT EXISTS pairwise_ld(
					       sample_id int not null,
					       seq_region_id int not null,
					       seq_region_start int not null,
					       seq_region_end int not null,
					       r2 float not null,
					       
					       key seq_region_idx(seq_region_id,seq_region_start)
					       )
	});

}
