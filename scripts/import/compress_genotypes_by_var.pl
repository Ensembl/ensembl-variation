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
use ImportUtils qw(debug load);
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use DBH;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use FindBin qw( $Bin );
use POSIX;
use Data::Dumper;

$| = 1;

my %Printable = ( "\\"=>'\\', "\r"=>'r', "\n"=>'n', "\t"=>'t', "\""=>'"' );


my ($TMP_DIR, $TMP_FILE, $species, $registry_file, $genotype_table, $start_sp_id, $monoploid, $jump, $flip, $straight, $no_sort);

GetOptions(
	'tmpdir=s'        => \$TMP_DIR,
	'tmpfile=s'       => \$TMP_FILE,
	'table=s'         => \$genotype_table,
	'species=s'       => \$species,
	'registry_file=s' => \$registry_file,
	'variation_id=s'  => \$start_sp_id,
	'monoploid'       => \$monoploid,
	'jump=s'          => \$jump,
	'flip'            => \$flip,
	'straight'        => \$straight,
	'no_sort'         => \$no_sort,
);


my $dump_file = 'compressed_genotype_var.txt';
$genotype_table ||= 'tmp_sample_genotype_single_bp';
$start_sp_id ||= 1;
$jump ||= 100000;
$monoploid = 0 unless defined($monoploid);
$no_sort ||= 0;

warn("Make sure you have an updated ensembl.registry file!\n");

usage('-TMP_DIR argument is required') if(!$TMP_DIR);
usage('-TMP_FILE argument is required') if(!$TMP_FILE);
usage('-species argument is required') if(!$species);

$registry_file ||= $Bin . "/ensembl.registry";

my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all( $registry_file );
#$reg->set_reconnect_when_lost();

my $vdba = $reg->get_DBAdaptor($species,'variation');
my $dbCore = $reg->get_DBAdaptor($species,'core');

my $dbVar = $vdba->dbc->db_handle;

$ImportUtils::TMP_DIR = $TMP_DIR;
$ImportUtils::TMP_FILE = $TMP_FILE;

debug("Examining table types");

# find out if this is a MRG_MyISAM table
my $desc_sth = $dbVar->prepare("SHOW CREATE TABLE $genotype_table");
$desc_sth->execute();

my ($tmp_table, $create_sql);
$desc_sth->bind_columns(\$tmp_table, \$create_sql);
$desc_sth->fetch();

if($create_sql =~ /mrg_myisam/i) {
	debug("Specified genotype table looks like a MRG_MyISAM table - forcing --straight");
	$straight = 1;
	
	unless($no_sort) {
		debug("Sorting tables - this may take a while!");
		
		if($create_sql =~ m/UNION\=\(([\w\`\,]+?)\)/i) {
			my @sub_tables = split /\,/, $1;
			$_ =~ s/\`//g for @sub_tables;
			
			foreach my $sub_table(@sub_tables) {
				debug("Sorting $sub_table");
				$dbVar->do("ALTER TABLE $sub_table ORDER BY variation_id, subsnp_id");
			}
		}
	}
}

#exit(0);

debug("Loading genotype codes");

# genotype codes
my $gca = $reg->get_adaptor($species, 'variation', 'genotypecode');
my %codes = map {(join "|", @{$_->genotype}) => $_->dbID} @{$gca->fetch_all()};
my $max_code = (sort {$a <=> $b} values %codes)[-1] || 0;

debug("Getting flipped status");

# flipped status
my $flipped_status = defined($flip) ? flipped_status($dbVar) : undef;

# set value for progress
my $prev_percent;

compress_genotypes($dbCore,$dbVar);

#reads the genotype and variation_feature data and compress it in one table with blob field
sub compress_genotypes{
	my $dbCore = shift;
	my $dbVar = shift;
	my %genotypes = (); #hash containing genotypes
	
	# find out the maximum variation_id in the genotype table
	my ($max_var_id, $row_count, $sth);
	my $sth1 = $dbVar->prepare(qq{SELECT max(variation_id) from $genotype_table});
	$sth1->execute;
	$sth1->bind_columns(\$max_var_id);
	$sth1->fetch;
	$sth1->finish;
	
	if(!defined($max_var_id)) {
		die("ERROR: No genotypes found in table $genotype_table\n");
	}
	
	my $low = $start_sp_id;
	my $tmp_var_table = 'tmp_var_'.$$;
	
	# straight, no ordering (for human)
	if(defined($straight)) {
		$sth = $dbVar->prepare(qq{
			SELECT gt.variation_id, gt.subsnp_id, gt.allele_1, gt.allele_2, gt.sample_id FROM
			$genotype_table gt
		}, {mysql_use_result => 1});
		
		# count rows
		my $sth3 = $dbVar->prepare(qq{SELECT count(*) from $genotype_table});
		$sth3->execute;
		$sth3->bind_columns(\$row_count);
		$sth3->fetch;
		$sth3->finish;
		
		open OUT, ">$TMP_DIR/$dump_file";
	}
	
	else {
		# create a tmp_var table
		$dbVar->do(qq{
			CREATE TABLE $tmp_var_table (
				variation_id int(11) unsigned NOT NULL,
				PRIMARY KEY (variation_id)
			);
		});
		
		$sth = $dbVar->prepare(qq{
			SELECT STRAIGHT_JOIN gt.variation_id, gt.subsnp_id, gt.allele_1, gt.allele_2, gt.sample_id FROM
			$tmp_var_table tv, $genotype_table gt 
			WHERE gt.variation_id = tv.variation_id
		}, {mysql_use_result => 1});
	}
	
	my $got_gts = 0;
	
	while($low <= $max_var_id) {
		
		unless(defined($straight)) {
			# we do a bit of trickery here to help slow old MySQL
			# first write out a file containing a range of var IDs
			$dbVar->do(qq{TRUNCATE $tmp_var_table});
			open TMP_VAR, ">$TMP_DIR/$tmp_var_table\.txt" or die "Could not write to temporary file $TMP_DIR/$tmp_var_table\.txt\n";
			print TMP_VAR "$_\n" for ($low..($low+$jump)-1);
			close TMP_VAR;
			
			# load that into our tmp_var table
			my $call = "mv $TMP_DIR/$tmp_var_table\.txt $TMP_DIR/$TMP_FILE";
			system($call);
			load($dbVar, ($tmp_var_table, "variation_id"));
			
			&progress($low, $max_var_id);
			open OUT, ">$TMP_DIR/$dump_file";
		}
		
		# now execute statement that joins the genotype table to the tmp_var table
		$sth->execute();#$low, $low+$jump);
		$low += defined($straight) ? $max_var_id : $jump;
		
		my ($var_id, $ss_id, $allele_1, $allele_2, $sample_id);
		my $previous_var_id = 0;
		
		$sth->bind_columns(\$var_id, \$ss_id, \$allele_1, \$allele_2, \$sample_id);
		
		my %new_gt_codes;
		%genotypes = ();
		
		my $i = 0;
		
		while ($sth->fetch){
			
			$i++;
			
			&progress($i, $row_count) if defined($straight);
			
			$got_gts = 1;
			
			# fill blanks
			$ss_id    ||= '\N';
			$allele_1 ||= "";
			$allele_2 ||= "";
			
			next unless defined($sample_id);
			
			my $combo_id = $var_id."\t".$ss_id;
			my $use_id = defined($straight) ? $combo_id : $var_id;
			
			# new var, print genotypes
			if(
			   $previous_var_id ne 0 && $previous_var_id ne $use_id && scalar keys %genotypes
			){
				print OUT "$_\t".$genotypes{$_}."\n" for keys %genotypes;
				%genotypes = ();
			}
			
			if(defined($flip)) {
				my $flipped = $flipped_status->{$var_id};
				reverse_alleles(\$allele_1, \$allele_2) if defined($flipped) && $flipped;
			}
			
			my $genotype = ($monoploid ? $allele_1 : $allele_1.'|'.$allele_2);
			my $genotype_code = $codes{$genotype};
			
			if(!defined($genotype_code)) {
				$genotype_code = ++$max_code;
				$new_gt_codes{$genotype_code} = $genotype;
				$codes{$genotype} = $genotype_code;
			}
			
			$genotypes{$combo_id} .= escape(pack("ww", $sample_id, $genotype_code));
			
			#print "Encoding $genotype as $sample_id $genotype_code\n";
			#exit(0);
			
			$previous_var_id = $use_id;
		}
		
		$sth->finish();
		
		if(scalar keys %genotypes && $previous_var_id ne 0) {
			print OUT "$_\t".$genotypes{$_}."\n" for keys %genotypes;
			%genotypes = ();
		}
		close OUT;
	
		#and import remaining genotypes
		&import_genotypes($dbVar) if $got_gts;
		
		create_genotype_code($dbVar, $new_gt_codes{$_}, $_) for sort {$a <=> $b} keys %new_gt_codes;
	}
	
	$dbVar->do(qq{DROP TABLE $tmp_var_table}) unless defined($straight);
	
	&end_progress();
}

sub import_genotypes{
    my $dbVar = shift;
    my $call = "mv $TMP_DIR/$dump_file $TMP_DIR/$TMP_FILE";
    system($call);
	
	load($dbVar,qw(compressed_genotype_var variation_id subsnp_id genotypes));
}

# creates a new genotype_code entry
sub create_genotype_code {
	my $dbVar = shift;
	my $genotype = shift;
	my $gt_code = shift;
	my $hap_id = 1;
	
	my $sth = $dbVar->prepare("INSERT INTO genotype_code (genotype_code_id, allele_code_id,  haplotype_id) VALUES(?, ?, ?)");
	
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




#checks if the genotype is in the opposite strand. If so, reverse the alleles
sub reverse_alleles{
    my $ref_allele_1 = shift;
    my $ref_allele_2 = shift;

	reverse_comp($ref_allele_1);
	reverse_comp($ref_allele_2);
}

# $special_characters_escaped = printable( $source_string );
sub escape ($) {
  local $_ = ( defined $_[0] ? $_[0] : '' );
  s/([\r\n\t\\\"])/\\$Printable{$1}/sg;
  return $_;
}

sub flipped_status {
	my $dbVar = shift;
	
	my %flipped_status;
	
	my ($var_id, $flipped);
	
	my $sth = $dbVar->prepare(qq{
		SELECT variation_id, flipped
		FROM variation v
		WHERE flipped = 1
	}, {mysql_use_result => 1});
	$sth->execute;
	$sth->bind_columns(\$var_id, \$flipped);
		
	while($sth->fetch()) {
		$flipped_status{$var_id} = $flipped;
	}
	
	$sth->finish();
	
	return \%flipped_status;
}

sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: perl compress_genotypes_by_var.pl <options>

options:
	-tmpdir <dir>            temp directory to use (with lots of space!)
	-tmpfile <filename>      name of temp file to use
	-registry <file>         registry file
	-species <species_name>  name of species
	-table                   name of table to get genotypes from
	-variation_id            variation_id to restart broken process from
	-flip                    flip genotypes according to flipped column in variation table
	-monoploid               indicate this species is monoploid
EOF

  die("\n$msg\n\n");
}

# update or initiate progress bar
sub progress {
    my ($i, $total) = @_;
    
    my $width = 100;
    my $percent = sprintf("%.1f", ($i/$total) * 100);
    my $numblobs = (($i/$total) * $width) - 2;
	
	return if defined($prev_percent) && $percent == $prev_percent;
	$prev_percent = $percent;
	
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
