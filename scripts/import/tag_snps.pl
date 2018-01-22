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
use ImportUtils qw(load);
use Bio::EnsEMBL::Registry;
use FindBin qw( $Bin );


##
##  run with bsub -q long -o [tmpdir]/output_tag.txt perl tag_snps.pl -tmpdir [tmpdir] -tmpfile [tmpfile]
use constant MAX_SIZE => 500_000_000;
my ($TMP_DIR, $TMP_FILE, $species, $registry_file, $selected_seq_region);


GetOptions('species=s' => \$species,
	   'tmpdir=s'  => \$ImportUtils::TMP_DIR,
	   'tmpfile=s' => \$ImportUtils::TMP_FILE,
	   'registry_file=s' => \$registry_file,
	   'seq_region=i' => \$selected_seq_region);

warn("Make sure you have a updated ensembl.registry file!\n");

$selected_seq_region ||= $ENV{LSB_JOBINDEX} if defined($ENV{LSB_JOBINDEX});

$registry_file ||= $Bin . "/ensembl.registry";

$TMP_DIR = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

#added default options
$species ||= 'human';
my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all( $registry_file );
my $dbVariation = $reg->get_DBAdaptor($species,'variation');

my $siblings;
create_LD_table($dbVariation); #create the LD table if it is not there yet
my $in_str = get_LD_populations($dbVariation,$siblings);

$in_str .= ' AND c.seq_region_id = '.$selected_seq_region if defined($selected_seq_region);

# genotype codes
my $gca = $reg->get_adaptor($species, 'variation', 'genotypecode');
my %codes = map {$_->dbID => (join "|", @{$_->genotype})} @{$gca->fetch_all()};

my $sth = $dbVariation->dbc->prepare(qq{SELECT c.sample_id,c.seq_region_id,c.seq_region_start,c.seq_region_end,c.genotypes,ip.population_sample_id
				 FROM compressed_genotype_region c FORCE INDEX(pos_idx), individual_population ip
				 WHERE  ip.individual_sample_id = c.sample_id
				 AND   ip.population_sample_id $in_str
			     },{mysql_use_result => 1});

$sth->execute();
debug("Dumping genotype data");
create_files($sth,$siblings);
$sth->finish();

my @files = glob("$TMP_DIR/dump*");
my $call = '';
my @stats;
debug("Submitting calculate LD jobs");


my $files_size = {};
foreach my $file (@files){
    @stats = stat($file);
    $files_size->{$file} = $stats[7];
}
foreach my $file (sort {$files_size->{$b} <=> $files_size->{$a}} keys %{$files_size}){
	
    $file =~ /.*_(\d+)_(\d+)/;
	$call = "bsub -J ld_calculation_$1\_$2 -o $TMP_DIR\/$1\_$2\_farm.out ";
	
    @stats = stat($file);
    if ($stats[7] > MAX_SIZE){
		$call .= " -q hugemem  -R'select[mem>28000] rusage[mem=28000]' -M28000000 ";
    }
	else {
		$call .= " -q basement  -R'select[mem>8000] rusage[mem=8000]' -M8000000 ";
	}
    $call .= "perl calculate_ld_table.pl -tmpdir $TMP_DIR -tmpfile $TMP_FILE -ldfile $file ";
    print $call,"\n";
    system($call);
    $call = '';
}

sleep(60);
$call = "bsub -K -W2:00 -w -q basement 'done(ld_calculation*)' -J waiting_ld date";
system($call);

# import the data in the LD table
debug("Importing LD values to DB");


my @ld_files = glob("$TMP_DIR/$TMP_FILE*.out");
foreach my $file(@ld_files) {
	system("mv $file $TMP_DIR/$TMP_FILE");
	load($dbVariation->dbc,qw(pairwise_ld sample_id seq_region_id seq_region_start seq_region_end r2));
}

debug("Submitting find tag SNPs jobs");
$call = '';

foreach my $file (sort {$files_size->{$b} <=> $files_size->{$a}} keys %{$files_size}){
    $file =~ /.*_(\d+)_(\d+)/;
	
	$call = "bsub -J tag_snps_$1\_$2 -o $TMP_DIR\/$1\_$2\_farm.out ";
	
    if ($files_size->{$file} > MAX_SIZE){
		$call .= " -q hugemem  -R'select[mem>28000] rusage[mem=28000]' -M28000000 ";
    }
	else {
		$call .= " -q basement  -R'select[mem>8000] rusage[mem=8000]' -M8000000 ";
	}
    $call .= "/software/bin/env perl select_tag_snps.pl $file ";
    $call .= " -tmpdir $TMP_DIR -species $species";
	$call .= " -registry_file $registry_file" if defined($registry_file);
    
    print $call,"\n";
    system($call); #send the job to one queue
    $call = '';
    
}

$call = "bsub -K -w 'done(tag_snps*)' -q basement -J waiting_process sleep 1"; #waits until all snp_tagging have finished to continue
system($call);

debug("Importing tag SNPs");
&import_tagg_snps($dbVariation);


unlink(<$TMP_DIR/tag_snps*.out>);

sub store_file{
    my $buffer = shift;
    my $individual_id = shift;
    my $seq_region_start = shift;
    my $genotype = shift;
    my $population_id = shift;
    my $seq_region_id = shift;
	
	my @genotypes = unpack '(www)*', $genotype;
	my $snp_start = $seq_region_start;
	
	while( my( $variation_id, $gt_code, $gap ) = splice @genotypes, 0, 3 ) {
		
		my ($allele_1, $allele_2) = split /\|/, $codes{$gt_code};
		
		$allele_1 ||= '-';
		$allele_2 ||= '-';
		
		print_buffered($buffer,"$TMP_DIR/dump_data_$population_id\_$seq_region_id.txt",join("\t",$snp_start,$individual_id,$allele_1,$allele_2)."\n");
		
		$snp_start += $gap + 1 if defined $gap;
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
	    open( FH, ">>$file" ) or die "Could not write to file $file";
	    print FH $buffer->{ $file };
	    close FH;
	}

    } else {
	$buffer->{ $filename } .= $text;
	if( length( $buffer->{ $filename } ) > 10_000 ) {
	    open( FH, ">>$filename" ) or die "Could not write to file $filename";
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
    my $sth = $dbVariation->dbc->prepare(qq{
		SELECT s.sample_id, s.name
		FROM population p, sample s
		WHERE (s.name like 'PERLEGEN:AFD%'
		OR s.name like 'CSHL-HAPMAP%'
		OR s.name like '1000GENOMES%'
		OR s.display = 'LD')
		AND s.sample_id = p.sample_id
	});

    $sth->execute();
    $sth->bind_columns(\$pop_id,\$population_name);
    #get all the children that we do not want in the genotypes
    my @pops;
    while($sth->fetch){
	#if($population_name =~ /CEU|YRI|MEX/){
	    get_siblings($dbVariation,$pop_id,$siblings);
	#}
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

    my $sth_individual = $dbVariation->dbc->prepare(qq{
		SELECT i.sample_id
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
		# necessary to have in the key the population,
		# since some individuals are shared between populations
		$siblings->{$population_id.'-'.$individual_id}++;
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
		);
	});

    $dbVariation->dbc->do(qq{TRUNCATE pairwise_ld});

}


sub import_tagg_snps{
    my $dbVariation = shift;
    
	my @ld_files = glob("$TMP_DIR/snps_tagged_*.txt");
	foreach my $file(@ld_files) {
		system("mv $file $TMP_DIR/$TMP_FILE");
		load($dbVariation->dbc(), qw(tagged_variation_feature variation_feature_id tagged_variation_feature_id sample_id));
	}
}

# gets time
sub get_time() {
	my @time = localtime(time());

	# increment the month (Jan = 0)
	$time[4]++;

	# add leading zeroes as required
	for my $i(0..4) {
		$time[$i] = "0".$time[$i] if $time[$i] < 10;
	}

	# put the components together in a string
	my $time =
 		($time[5] + 1900)."-".
 		$time[4]."-".
 		$time[3]." ".
		$time[2].":".
		$time[1].":".
		$time[0];

	return $time;
}



# prints debug output with time
sub debug {
	my $text = (@_ ? (join "", @_) : "No message");
	my $time = get_time;
	
	print $time." - ".$text.($text =~ /\n$/ ? "" : "\n");
}
