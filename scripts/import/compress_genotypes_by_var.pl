#!/usr/local/ensembl/bin/perl
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


my ($TMP_DIR, $TMP_FILE, $species, $registry_file, $genotype_table, $start_sp_id, $monoploid, $jump, $no_flip);

GetOptions(
	'tmpdir=s'        => \$TMP_DIR,
	'tmpfile=s'       => \$TMP_FILE,
	'table=s'         => \$genotype_table,
	'species=s'       => \$species,
	'registry_file=s' => \$registry_file,
	'subsnp_proxy=s'  => \$start_sp_id,
	'monoploid'       => \$monoploid,
	'jump=s'          => \$jump,
	'no_flip'         => \$no_flip,
);


my $dump_file = 'compressed_genotype_var.txt';
$genotype_table ||= 'tmp_individual_genotype_single_bp';
$start_sp_id ||= 1;
$jump ||= 1000000;
$monoploid = 0 unless defined($monoploid);

warn("Make sure you have an updated ensembl.registry file!\n");

usage('-TMP_DIR argument is required') if(!$TMP_DIR);
usage('-TMP_FILE argument is required') if(!$TMP_FILE);
usage('-species argument is required') if(!$species);

$registry_file ||= $Bin . "/ensembl.registry";

my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all( $registry_file );
$reg->set_reconnect_when_lost();

my $vdba = $reg->get_DBAdaptor($species,'variation');
my $dbCore = $reg->get_DBAdaptor($species,'core');

my $dbVar = $vdba->dbc->db_handle;

$ImportUtils::TMP_DIR = $TMP_DIR;
$ImportUtils::TMP_FILE = $TMP_FILE;

print "Loading genotype codes\n";

# genotype codes
my $gca = $reg->get_adaptor($species, 'variation', 'genotypecode');
my %codes = map {(join "|", @{$_->genotype}) => $_->dbID} @{$gca->fetch_all()};
my $max_code = (sort {$a <=> $b} values %codes)[-1] || 0;

print "Getting flipped status\n";

# flipped status
my $flipped_status = defined($no_flip) ? undef : flipped_status($dbVar);


	
# find out if the genotype table has subsnp_proxy_id or variation_id+subsnp_id
my $sth1 = $dbVar->prepare(qq{DESCRIBE $genotype_table});
$sth1->execute();
my $has_proxy = 0;

while(my $ref = $sth1->fetchrow_arrayref) {
	$has_proxy = 1 if $ref->[0] eq 'subsnp_proxy_id';
}
$sth1->finish();

unless($has_proxy) {
	$dbVar->do(qq{
		CREATE TABLE IF NOT EXISTS tmp_compressed_genotype_var (
			variation_id int(11) unsigned NOT NULL,
			subsnp_id int(11) unsigned NOT NULL,
			genotypes blob,
			PRIMARY KEY (variation_id, subsnp_id)
		) ENGINE=MyISAM DEFAULT CHARSET=latin1;
	});
	$dbVar->do(qq{TRUNCATE tmp_compressed_genotype_var;});
}


compress_genotypes($dbCore,$dbVar);

#reads the genotype and variation_feature data and compress it in one table with blob field
sub compress_genotypes{
    my $dbCore = shift;
    my $dbVar = shift;
    my $buffer;
    my $genotypes = ''; #string containing genotypes
    my $blob = '';
    my $count = 0;
	
	my $max_var_id;
	my $tmp_type = $has_proxy ? 'subsnp_proxy' : 'variation';
	my $sth1 = $dbVar->prepare(qq{select max($tmp_type\_id) FROM $tmp_type});
	$sth1->execute;
	$sth1->bind_columns(\$max_var_id);
	$sth1->fetch;
	$sth1->finish;
	
	print "$genotype_table ".($has_proxy ? "has" : "does not have")." subsnp_proxy_id\n";
	
	my $low = $start_sp_id;
	my $high = $jump;
	
	while($low <= $max_var_id) {
		&progress($low, $max_var_id);
		
		open OUT, ">$TMP_DIR/$dump_file";
		
		my $sth;
		
		if($has_proxy) {
			$sth = $dbVar->prepare(qq{
				SELECT subsnp_proxy_id, allele_1, allele_2, sample_id FROM
				$genotype_table
				WHERE subsnp_proxy_id >= ?
				AND subsnp_proxy_id < ?
				ORDER BY subsnp_proxy_id, sample_id
			}, {mysql_use_result => 1});
		}
		
		else {
			$sth = $dbVar->prepare(qq{
				SELECT variation_id, subsnp_id, allele_1, allele_2, sample_id FROM
				$genotype_table gt
				WHERE variation_id >= ?
				AND variation_id < ?
				ORDER BY variation_id, subsnp_id, sample_id
			}, {mysql_use_result => 1});
		}
		
		$sth->execute($low, $high);
		$low += $jump;
		$high += $jump;
		
		#information dumping from database
		my ($var_id, $ss_id, $allele_1, $allele_2, $sample_id);
		my $previous_var_id = 0;
		
		if($has_proxy) {
			$sth->bind_columns(\$var_id, \$allele_1, \$allele_2, \$sample_id);
		}
		else {
			$sth->bind_columns(\$var_id, \$ss_id, \$allele_1, \$allele_2, \$sample_id);
		}
		
		my %new_gt_codes;
		
		while ($sth->fetch){
			
			# fill blanks
			$ss_id    ||= 0;
			$allele_1 ||= "";
			$allele_2 ||= "";
			
			# new var, print genotypes
			if(
			   $previous_var_id ne 0 &&
				(($has_proxy && $previous_var_id != $var_id) || $previous_var_id ne $var_id."\t".$ss_id)
			){ 
				print OUT "$previous_var_id\t$genotypes\n";
				$genotypes = '';
			}
			
			unless(defined($no_flip)) {
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
			
			$genotypes .= escape(pack("ww", $sample_id, $genotype_code));
			
			$previous_var_id = $has_proxy ? $var_id : $var_id."\t".$ss_id;
		}
		$sth->finish();
		print OUT "$previous_var_id\t$genotypes\n";
		
		create_genotype_code($dbVar, $new_gt_codes{$_}, $_) for sort {$a <=> $b} keys %new_gt_codes;
		
		close OUT;
		
		#and import remaining genotypes
		&import_genotypes($dbVar);
	}
	
	&end_progress();
	
	unless($has_proxy) {
		print "Converting table\n";
		
		$dbVar->do(qq{
			INSERT IGNORE INTO compressed_genotype_var
			SELECT sp.subsnp_proxy_id, gt.genotypes
			FROM subsnp_proxy sp FORCE INDEX (variation_idx), tmp_compressed_genotype_var gt
			WHERE sp.variation_id = gt.variation_id
			AND sp.subsnp_id = gt.subsnp_id;
		});
	}
}

sub import_genotypes{
    my $dbVar = shift;
    my $call = "mv $TMP_DIR/$dump_file $TMP_DIR/$TMP_FILE";
    system($call);
	
	my $table = $has_proxy ? load($dbVar,qw(compressed_genotype_var subsnp_proxy_id genotypes)) : load($dbVar,qw(tmp_compressed_genotype_var variation_id subsnp_id genotypes));
}

# creates a new genotype_code entry
sub create_genotype_code {
	my $dbVar = shift;
	my $genotype = shift;
	my $gt_code = shift;
	my $hap_id = 1;
	
	my $sth = $dbVar->prepare("INSERT INTO genotype_code VALUES(?, ?, ?)");
	
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
		SELECT sp.subsnp_proxy_id, v.flipped
		FROM variation v, subsnp_proxy sp
		WHERE v.variation_id = sp.variation_id
		AND v.flipped = 1
	});
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
	-subsnp_proxy            subsnp_proxy_id to restart broken process from
	-no_flip                 don't flip genotypes
	-monoploid               indicate this species is monoploid
EOF

  die("\n$msg\n\n");
}

# update or initiate progress bar
sub progress {
    my ($i, $total) = @_;
    
    my $width = 60;
    my $percent = int(($i/$total) * 100);
    my $numblobs = (($i/$total) * $width) - 2;
	
    printf("\r% -${width}s% 1s% 10s", '['.('=' x $numblobs).($numblobs == $width - 2 ? '=' : '>'), ']', "[ " . $percent . "% ]");
}

# end progress bar
sub end_progress {
    progress(1,1);
    print "\n";
}
