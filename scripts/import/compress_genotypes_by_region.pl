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
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use DBH;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Registry;
use FindBin qw( $Bin );
use POSIX;
use Data::Dumper;

use constant DISTANCE => 100_000;
use constant MAX_SHORT => 2**16 -1;

$| = 1;

## Backslash Escaping


my %Printable = ( "\\"=>'\\', "\r"=>'r', "\n"=>'n', "\t"=>'t', "\""=>'"' );

my ($TMP_DIR, $TMP_FILE, $species, $selected_seq_region, $registry_file, $genotype_table, $monoploid, $no_flip, $jump, $has_proxy, $restart);

GetOptions(
	'tmpdir=s'        => \$TMP_DIR,
	'tmpfile=s'       => \$TMP_FILE,
	'species=s'       => \$species,
	'seq_region=i'    => \$selected_seq_region,
	'registry_file=s' => \$registry_file,
	'table=s'         => \$genotype_table,
	'monoploid'       => \$monoploid,
	'no_flip'         => \$no_flip,
	'jump=i'          => \$jump,
	'proxy=i'         => \$has_proxy,
	'start=i'         => \$restart,
);


$genotype_table ||= 'tmp_sample_genotype_single_bp';
$monoploid = 0 unless defined($monoploid);
$jump ||= 10000000;

#$selected_seq_region ||= $ENV{LSB_JOBINDEX} if defined($ENV{LSB_JOBINDEX});
undef $selected_seq_region unless $selected_seq_region;

$TMP_FILE .= "_".$selected_seq_region if defined($selected_seq_region);

my $dump_file = 'compressed_genotype_region.txt';
$dump_file .= "_".$selected_seq_region if defined($selected_seq_region);

warn("Make sure you have an updated ensembl.registry file!\n");

usage('-TMP_DIR argument is required') if(!$TMP_DIR);
usage('-TMP_FILE argument is required') if(!$TMP_FILE);
usage('-species argument is required') if(!$species);

$registry_file ||= $Bin . "/ensembl.registry";

my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all( $registry_file );
$reg->set_reconnect_when_lost();

my $vdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');

my $dbVar = $vdba->dbc->db_handle;

$ImportUtils::TMP_DIR = $TMP_DIR;
$ImportUtils::TMP_FILE = $TMP_FILE;

# genotype codes
my $gca = $reg->get_adaptor($species, 'variation', 'genotypecode');
my %codes = map {(join "|", @{$_->genotype}) => $_->dbID} @{$gca->fetch_all()};
my $max_code = (sort {$a <=> $b} values %codes)[-1] || 0;

my $flipped_status = defined($no_flip) ? {} : flipped_status($dbVar, $genotype_table);

$has_proxy = 0;

if(!defined($has_proxy)) {
	# find out if this table has subsnp_proxy_id or variation_id+subsnp_id
	my $sth = $dbVar->prepare(qq{DESCRIBE $genotype_table});
	$sth->execute();
	
	while(my $ref = $sth->fetchrow_arrayref) {
		$has_proxy = 1 if $ref->[0] eq 'subsnp_proxy_id';
	}
	$sth->finish();
}

compress_genotypes($dbCore,$dbVar);
update_meta_coord($dbCore,$dbVar,"compressed_genotype_region");

#reads the genotype and variation_feature data and compress it in one table with blob field
sub compress_genotypes{
    my $dbCore = shift;
    my $dbVar = shift;
    my $buffer;
    my $blob = '';
    my $count = 0;
	
	# get sample_ids
	my $sth = $dbVar->prepare(qq{
		SELECT distinct(s.sample_id)
		FROM sample s, population p, sample_population sp
		WHERE s.sample_id = sp.sample_id
		AND sp.population_id = p.population_id
		AND (
		   s.display IN ('DEFAULT','DISPLAYABLE','REFERENCE')
		   OR p.display = 'LD'
		) 
	});
	$sth->execute;
	my $sample_id;
	$sth->bind_columns(\$sample_id);
	
	my @sample_ids;
	push @sample_ids, $sample_id while $sth->fetch();
	$sth->finish;

	print "Will compress genotypes from ", scalar @sample_ids, " sample_ids\n";
	die "No samples found to compress" unless scalar @sample_ids;
	
	my $sample_string = (scalar @sample_ids == 1 ? ' = '.$sample_ids[0] : ' IN ('.(join ",", @sample_ids).') ');
	
	
	# get seq_regions
	my @seq_regions;
	
	if(defined($selected_seq_region)) {
		foreach my $val(split /\,/, $selected_seq_region) {
            my @nnn = split /\-/, $val;
            
            foreach my $sr($nnn[0]..$nnn[-1]) {
                push @seq_regions, $sr;
            }
        }
	}
	else {
		$sth = $dbVar->prepare(qq{SELECT distinct(seq_region_id) FROM variation_feature});
		$sth->execute;
		
		my $sr;
		$sth->bind_columns(\$sr);
		push @seq_regions, $sr while $sth->fetch;
		$sth->finish();
	}
	
	
	# create a temporary table to hold VF entries
	my $tmp_table = 'variation_feature_'.$$;
	$dbVar->do(qq{CREATE TABLE IF NOT EXISTS $tmp_table LIKE variation_feature});
	$dbVar->do(qq{TRUNCATE $tmp_table});
	
	my $sr_counter = 0;
	
	foreach my $seq_region(@seq_regions) {
		
		$sr_counter++;
		
		debug("Dumping genotypes from seq_region $seq_region ($sr_counter\/".(scalar @seq_regions).")");
		
		$sth = $dbVar->prepare(qq{
			SELECT min(seq_region_start), max(seq_region_start)
			FROM variation_feature
			WHERE seq_region_id = ?
		});
		$sth->execute($seq_region);
		
		my ($min_pos, $max_pos);
		$sth->bind_columns(\$min_pos, \$max_pos);
		$sth->fetch;
		$sth->finish;
		
		my $origin = $min_pos;
		
		# restart
		if(defined($restart)) {
			debug("Resuming from position $restart");
			$min_pos = $restart;
			undef $restart;
		}
		
		while($min_pos <= $max_pos) {
			
			&progress($min_pos - $origin, $max_pos - $origin);
			
			# populate tmp vf table
			$dbVar->do(qq{TRUNCATE $tmp_table});
			$sth = $dbVar->prepare(qq{
				INSERT INTO $tmp_table
				SELECT * FROM variation_feature
				WHERE seq_region_id = ?
				AND seq_region_start >= ?
				AND seq_region_start < ?
			});
			$sth->execute($seq_region, $min_pos, $min_pos + $jump);
			$sth->finish();
			
			$dbVar->do(qq{ALTER TABLE $tmp_table ORDER BY seq_region_start ASC});
			
			my %new_gt_codes;
			
			if($has_proxy) {
				$sth = $dbVar->prepare(qq{
					SELECT vf.seq_region_id, vf.seq_region_start, vf.seq_region_end, vf.seq_region_strand, gt.allele_1, gt.allele_2, gt.sample_id, vf.variation_id
					FROM $tmp_table vf, $genotype_table as gt, subsnp_proxy_sp
					WHERE sp.variation_id = vf.variation_id
					AND sp.subsnp_proxy_id = gt.subsnp_proxy_id
					AND gt.sample_id $sample_string
					ORDER BY vf.seq_region_start
				}, {mysql_use_result => 1});
			}
			else {
				$sth = $dbVar->prepare(qq{
					SELECT vf.seq_region_id, vf.seq_region_start, vf.seq_region_end, vf.seq_region_strand, gt.allele_1, gt.allele_2, gt.sample_id, vf.variation_id
					FROM $tmp_table vf, $genotype_table as gt
					WHERE gt.variation_id = vf.variation_id
					AND gt.sample_id $sample_string
					ORDER BY vf.seq_region_start
				}, {mysql_use_result => 1});
			}
			
			#print "Time starting to dump data from seq_region_id $seq_region: ",scalar(localtime(time)),"\n";
			$sth->execute();
			
			#information dumping from database
			my ($seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $allele_1, $allele_2, $variation_id);
			my $previous_seq_region_id = 0;
		
			$sth->bind_columns(\$seq_region_id, \$seq_region_start, \$seq_region_end, \$seq_region_strand, \$allele_1, \$allele_2, \$sample_id, \$variation_id);
			
			my $genotypes = {}; #hash containing all the genotypes
			
			while ($sth->fetch) {
				next unless defined($sample_id);
				
				#unless(defined($no_flip)) {
				my $flipped = (defined($flipped_status->{$variation_id}) ? $flipped_status->{$variation_id} : undef);
				reverse_alleles($seq_region_strand, $flipped, \$allele_1, \$allele_2);
				#}
				
				$allele_1 ||= '';
				$allele_2 ||= '';
				
				my $genotype = ($monoploid ? $allele_1 : $allele_1.'|'.$allele_2);
				my $genotype_code = $codes{$genotype};
				
				if(!defined($genotype_code)) {
					$genotype_code = ++$max_code;
					$new_gt_codes{$genotype_code} = $genotype;
					$codes{$genotype} = $genotype_code;
				}
				
				#new chromosome, print all remaining genotypes and upload the file
				if ($previous_seq_region_id != $seq_region_id && $previous_seq_region_id != 0){
					print_file("$TMP_DIR/$dump_file",$genotypes, $previous_seq_region_id);
					$genotypes = {}; #and flush the hash
					#need to fork for upload the data
					my $pid = fork;
					if (! defined $pid){
					throw("Not possible to fork: $!\n");
					}
					elsif ($pid == 0){
					#you are the child.....
					my $dbVar_write = $vdba->dbc->db_handle;
					&import_genotypes($dbVar_write);
					  POSIX:_exit(0);
					}
					else{
					#the parent waits for the child to finish;
					waitpid($pid,0);
					}
				}
				next if (defined $genotypes->{$sample_id}->{region_end} && $seq_region_start == $genotypes->{$sample_id}->{region_end}); #same variation and sample but different genotype !!!
				#new region for the sample, update coordinates
				if (!defined $genotypes->{$sample_id}->{region_start}){
					$genotypes->{$sample_id}->{region_start} = $seq_region_start;
					$genotypes->{$sample_id}->{region_end} = $seq_region_end;
				}
				#compare with the beginning of the region if it is within the DISTANCE of compression
				if ((abs($genotypes->{$sample_id}->{region_start} - $seq_region_start) > DISTANCE()) || (abs($seq_region_start - $genotypes->{$sample_id}->{region_end}) > MAX_SHORT)){
					#snp outside the region, print the region for the sample we have already visited and start a new one
					print_file("$TMP_DIR/$dump_file",$genotypes, $seq_region_id, $sample_id);
					delete $genotypes->{$sample_id}; #and remove the printed entry
					$genotypes->{$sample_id}->{region_start} = $seq_region_start;
				}
				#escape characters (tab, new line)
				#and write it to the buffer
				if ($seq_region_start != $genotypes->{$sample_id}->{region_start}){					
					#compress information
					$blob = pack ("w",$seq_region_start - $genotypes->{$sample_id}->{region_end} - 1);
					$genotypes->{$sample_id}->{genotypes} .= escape($blob) .escape(pack("w", $variation_id)). escape(pack("w", $genotype_code));
				}
				else{
					#first genotype starts in the region_start, not necessary the number
					$genotypes->{$sample_id}->{genotypes} = escape(pack("w", $variation_id)).escape(pack("w", $genotype_code));
				}
				#and update the new region_end
				$genotypes->{$sample_id}->{region_end} = $seq_region_start;    #to avoid nasty effects of indels coordinates
				$previous_seq_region_id = $seq_region_id;
			
			}
			$sth->finish();
			
			create_genotype_code($dbVar, $new_gt_codes{$_}, $_) for sort {$a <=> $b} keys %new_gt_codes;
			
			#print "Time finishing dumping data: ",scalar(localtime(time)),"\n";
			#print last region
			print_file("$TMP_DIR/$dump_file",$genotypes, $previous_seq_region_id);
			#and import remainig genotypes
			&import_genotypes($dbVar);
			
			$min_pos += $jump;
		}
		
		&end_progress();
	}
	
	# clean up
	$dbVar->do(qq{DROP TABLE $tmp_table});
	$dbVar->do(qq{DELETE FROM compressed_genotype_region WHERE seq_region_id = 0;});
}


sub print_file{
    my $file = shift;
    my $genotypes = shift;
    my $seq_region_id = shift;
    my $sample_id = shift;

    open( FH, ">>$file") or die "Could not add compressed information: $!\n";
    if (!defined $sample_id){
	#new chromosome, print all the genotypes and flush the hash
	foreach my $sample_id (keys %{$genotypes}){
	    print FH join("\t",$sample_id,$seq_region_id, $genotypes->{$sample_id}->{region_start}, $genotypes->{$sample_id}->{region_end}, 1, $genotypes->{$sample_id}->{genotypes}) . "\n";
	}
    }
    else{
	#only print the region corresponding to sample_id
	print FH join("\t",$sample_id,$seq_region_id, $genotypes->{$sample_id}->{region_start}, $genotypes->{$sample_id}->{region_end}, 1, $genotypes->{$sample_id}->{genotypes}) . "\n";
    }
    close FH;
}

sub import_genotypes{
    my $dbVar = shift;

    #debug("Importing compressed genotype data");
	#debug("mv $TMP_DIR/$dump_file $TMP_DIR/$TMP_FILE");
	
	$dbVar->do(qq{LOAD DATA LOCAL INFILE "$TMP_DIR/$dump_file" INTO TABLE compressed_genotype_region});
	unlink("$TMP_DIR/$dump_file");
	
    #my $call = "mv $TMP_DIR/$dump_file $TMP_DIR/$TMP_FILE";
    #system($call);
    #load($dbVar,qw(compressed_genotype_region sample_id seq_region_id seq_region_start seq_region_end seq_region_strand genotypes));
}

#
# updates the meta coord table
#
sub update_meta_coord {
    my $dbCore = shift;
    my $dbVar  = shift;
    my $table_name = shift;
    my $csname = shift || 'chromosome';
    
    my $csa = $dbCore->get_CoordSystemAdaptor();

    my $cs = $csa->fetch_by_name($csname);

    my $sth = $dbVar->prepare
	('INSERT IGNORE INTO meta_coord (table_name,coord_system_id,max_length) VALUES (?,?,?)');

    $sth->execute($table_name, $cs->dbID(),DISTANCE+1);
    
    $sth->finish();
    
  return;
}


#checks if the genotype is in the opposite strand. If so, reverse the alleles
sub reverse_alleles{
    my $seq_region_strand = shift;
    my $flipped = shift;
    my $ref_allele_1 = shift;
    my $ref_allele_2 = shift;
    
	# work out whether we're flipping or not
	my $flip = 0;
	
	# flipped column has been set
	if(defined($flipped)) {
		if($flipped == 1 && $seq_region_strand == 1) {
			$flip = 1;
		}
		elsif($flipped == 0 && $seq_region_strand == -1) {
			$flip = 1;
		}
	}
	
	# flipped column has not been set - use old behaviour
	else {
		if($seq_region_strand == -1) {
			$flip = 1;
		}
	}

    if ($flip){
		reverse_comp($ref_allele_1);
		reverse_comp($ref_allele_2);
    }
}

# $special_characters_escaped = printable( $source_string );
sub escape ($) {
  local $_ = ( defined $_[0] ? $_[0] : '' );
  s/([\r\n\t\\\"])/\\$Printable{$1}/sg;
  return $_;
}

sub flipped_status {
	my $dbVar = shift;
	my $table = shift;
	
	my %flipped_status;
	
	my ($var_id, $flipped);
	
	my $sth = $dbVar->prepare("SELECT variation.variation_id, variation.flipped FROM variation, $table 
                                  WHERE variation.flipped =1 and variation.variation_id = $table.variation_id",
                                  {mysql_use_result => 1});
	$sth->execute;
	$sth->bind_columns(\$var_id, \$flipped);
		
	while($sth->fetch()) {
		$flipped_status{$var_id} = $flipped;
	}
	
	$sth->finish();
	
	return \%flipped_status;
}

# creates a new genotype_code entry
sub create_genotype_code {
	my $dbVar = shift;
	my $genotype = shift;
	my $gt_code = shift;
	my $hap_id = 1;
	
	my $sth = $dbVar->prepare("INSERT INTO genotype_code (genotype_code_id, allele_code_id, haplotype_id) VALUES(?, ?, ?)");
	
	foreach my $allele(split /\|/, $genotype) {
		my $allele_code = allele_code($dbVar, $allele);
		$sth->execute($gt_code, $allele_code, $hap_id++);
	}
	
	$sth->finish;
	return $gt_code;
}

# getter/setter for allele code
sub allele_code {
	my $dbVar = shift;
	my $allele = shift;
	
	my $sth = $dbVar->prepare("SELECT allele_code_id FROM allele_code WHERE allele = ?");
	$sth->execute($allele);
	
	my $allele_code;
	$sth->bind_columns(\$allele_code);
	$sth->fetch();
	$sth->finish();
	
	# create if doesn't exist
	if(!defined($allele_code)) {
		$sth = $dbVar->prepare("INSERT INTO allele_code(allele) VALUES(?)");
		$sth->execute($allele);
		$sth->finish();
		
		$allele_code = $dbVar->last_insert_id(undef, undef, qw(allele_code allele_code_id));
	}
	
	return $allele_code;
}


sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: perl compress_genotypes_by_region.pl <options>

options:
  -tmpdir <dir>            temp directory to use (with lots of space!)
  -tmpfile <filename>      name of temp file to use
  -registry <file>         registry file
  -species <species_name>  name of the species
  -table <table>           table name to read genotypes from
  -no_flip                 don't flip genotypes
  -monoploid               indicate this species is monoploid
EOF

  die("\n$msg\n\n");
}

# update or initiate progress bar
sub progress {
    my ($i, $total) = @_;
    
    my $width = 60;
    my $percent = $total > 0 ? int(($i/$total) * 100) : 0;
    my $numblobs = $total > 0 ? (($i/$total) * $width) - 2 : 0;
	
    printf("\r% -${width}s% 1s% 10s", '['.('=' x $numblobs).($numblobs == $width - 2 ? '=' : '>'), ']', "[ " . $percent . "% ]");
}

# end progress bar
sub end_progress {
    progress(1,1);
    print "\n";
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
