#!/usr/local/ensembl/bin/perl

use strict;
use warnings;

use Getopt::Long;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use ImportUtils qw(load);

my ($TMP_DIR, $TMP_FILE);

my ($vhost, $vport, $vdbname, $vuser, $vpass);



GetOptions('vhost=s'   => \$vhost,
	   'vuser=s'   => \$vuser,
	   'vpass=s'   => \$vpass,
	   'vport=i'   => \$vport,
	   'vdbname=s' => \$vdbname,
	   'tmpdir=s'  => \$ImportUtils::TMP_DIR,
	   'tmpfile=s' => \$ImportUtils::TMP_FILE);

$vhost    ||='ia64g';
$vport    ||= 3306;
$vuser    ||= 'ensadmin';

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

usage('-vdbname argument is required') if(!$vdbname);

my $dbVariation = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new
    (-host => $vhost,
     -user => $vuser,
     -pass => $vpass,
     -port => $vport,
     -dbname => $vdbname
     );


my $population_id;
my $population_name;
#get all populations to be tagged (HapMap and PerlEgen)
my $sth = $dbVariation->dbc->prepare(qq{SELECT population_id, name
					    FROM population
					    WHERE name like 'perlegen:afd%'
					    OR name like 'cshl-hapmap%'
					});

$sth->execute();
$sth->bind_columns(\$population_id,\$population_name);
#get all the children that we do not want in the genotypes
my $siblings = {}; # hash {$individual_id} ,where the individual is sibling of another one
my @pops;
while($sth->fetch){
    if($population_name =~ /CEU|YRI/){
	&get_siblings($dbVariation,$population_id,$siblings);
    }
    push @pops, $population_id;
}

my $in_str = " IN (" . join(',', @pops). ")";

#select all SNPs in the population given with genotypes, and order them by variation
$sth = $dbVariation->dbc->prepare(qq{SELECT vf.variation_feature_id, ig.allele_1, ig.allele_2, vf.seq_region_id, 
				            vf.seq_region_start, ip.individual_id, ip.population_id
				     FROM individual_genotype_single_bp ig, variation_feature vf, individual_population ip
				     WHERE ig.variation_id = vf.variation_id
				     AND ig.individual_id = ip.individual_id
				     AND vf.map_weight = 1
				     AND ip.population_id $in_str
				 },{mysql_use_result =>1});
my ($variation_feature_id, $allele_1, $allele_2, $seq_region_id, $seq_region_start,$individual_id);
my $buffer = {}; #buffer with data to write to file
$sth->execute();
$sth->bind_columns(\$variation_feature_id, \$allele_1, \$allele_2, \$seq_region_id, \$seq_region_start,\$individual_id,\$population_id);
while ($sth->fetch()){
    #only print genotypes without parents genotyped
    if (!exists $siblings->{$population_id.'-'.$individual_id}){
	#print the data to the file with the information in the chromosome
	print_buffered($buffer,"$TMP_DIR/tag_snps_$population_id\:$seq_region_id\.txt",
		       join("\t", $variation_feature_id, $allele_1, $allele_2, $seq_region_start) . "\n");
    }
}
$sth->finish();
print_buffered($buffer);
my $call;
#when the data is selected in the files (one per chromosome), send the jobs to the farm
foreach my $file (keys %{$buffer}){
    $file =~ /tag_snps_(\d+)\:.*/; #extract the population_id from the name of the file
    $call = "bsub -q normal -m 'bc_hosts ecs4_hosts' -J tag_snps_$1 ";
    $call .= "-o $TMP_DIR/output_tag.txt /usr/local/ensembl/bin/perl select_tag_snps.pl $file ";
    $call .= " -vhost $vhost -vuser $vuser -vport $vport -vdbname $vdbname -tmpdir $TMP_DIR ";
    $call .= " -population_id $1 ";
    
    $call .= "-vpass $vpass " if ($vpass);
    system($call); #send the job to one queue
    
}
$call = "bsub -K -w 'done(tag_snps*)' -J waiting_process sleep 1"; #waits until all snp_tagging have finished to continue
system($call);

&import_tagg_snps($dbVariation);    



#for a given population, gets all individuals that are children (have father or mother)
sub get_siblings{
    my $dbVariation = shift;
    my $population_id = shift;
    my $siblings = shift;

    my $sth_individual = $dbVariation->dbc()->prepare(qq{SELECT i.individual_id
							     FROM individual i, individual_population ip
							     WHERE ip.individual_id = i.individual_id
							     AND ip.population_id = ? 
							     AND i.father_individual_id IS NOT NULL
							     AND i.mother_individual_id IS NOT NULL
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
    load($dbVariation->dbc(), qw(tagged_variation_feature variation_feature_id population_id));
    unlink(<$TMP_DIR/tag_snps_*>);
   
}
sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: perl tag_snps.pl <options>

options:
    -vhost <hostname>    hostname of variation MySQL database to write to
    -vuser <user>        username of variation MySQL database to write to (default = ensadmin)
    -vpass <pass>        password of variation MySQL database to write to
    -vport <port>        TCP port of variation MySQL database to write to (default = 3306)
    -vdbname <dbname>    dbname of variation MySQL database to write to
    -tmpdir <dir>        temp directory to use (with lots of space!)
    -tmpfile <filename>  name of temp file to use

EOF

  die("\n$msg\n\n");
}
