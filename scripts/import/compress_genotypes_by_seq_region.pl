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
use ImportUtils qw(dumpSQL debug load);
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use DBH;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use FindBin qw( $Bin );
use POSIX;
use Data::Dumper;

use constant DISTANCE => 100_000;
use constant MAX_SHORT => 2**16 -1;

## Backslash Escaping


my %Printable = ( "\\"=>'\\', "\r"=>'r', "\n"=>'n', "\t"=>'t', "\""=>'"' );


my ($TMP_DIR, $TMP_FILE, $species, $allow_nulls);


GetOptions(   'tmpdir=s'  => \$TMP_DIR,
	      'tmpfile=s' => \$TMP_FILE,
	      'species=s' => \$species,
	      'allownulls' => \$allow_nulls,
		   );

warn("Make sure you have an updated ensembl.registry file!\n");

usage('-TMP_DIR argument is required') if(!$TMP_DIR);
usage('-TMP_FILE argument is required') if(!$TMP_FILE);
usage('-species argument is required') if(!$species);

my $registry_file ||= $Bin . "/ensembl.registry";

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $vdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');

my $dbVar = $vdba->dbc->db_handle;

$ImportUtils::TMP_DIR = $TMP_DIR;
$ImportUtils::TMP_FILE = $TMP_FILE;

compress_genotypes($dbCore,$dbVar, $allow_nulls);
update_meta_coord($dbCore,$dbVar,"compressed_genotype_single_bp");

#reads the genotype and variation_feature data and compress it in one table with blob field
sub compress_genotypes{
  my $dbCore = shift;
  my $dbVar = shift;
  my $allow_nulls = shift;

  my $buffer;
  my $genotypes = {};		#hash containing all the genotypes
  my $blob = '';
  my $count = 0;
  my %rec_seq_region;
    
  my $sth = $dbCore->dbc->prepare(qq{SELECT seq_region_id 
                                  FROM seq_region_attrib sa, attrib_type at
                                  WHERE sa.attrib_type_id=at.attrib_type_id
                                  AND at.code = "toplevel"});
  $sth->execute();

  while (my ($seq_region_id) = $sth->fetchrow_array()) {
    $rec_seq_region{$seq_region_id} = 1;
  }

  debug("Checking whether compressed_genotype table exist or not");
  my $table = "compressed_genotype_single_bp";
  my $res_ref = $dbVar->selectall_arrayref(qq{SELECT * FROM $table limit 10});
  my $table_rec = $res_ref->[0][0] ;
  if ($table_rec) {
    debug("TABLE $table have records exists, truncate the $table?");
    $dbVar->do(qq{TRUNCATE TABLE $table});
  }
  print "Time starting to dump data from database: ",scalar(localtime(time)),"\n";

  my $sql;

  foreach my $seq_region_id (keys %rec_seq_region) {
    if (! -e "$TMP_DIR/genotype_dump\_$seq_region_id") {
      ##for small amount of SNPs and a lot of seq_region_ids, such as zfish, it's slower to use this method
      debug("Dumping seq_region_id $seq_region_id...");

      $sql = qq{SELECT vf.seq_region_id, vf.seq_region_start, vf.seq_region_end, vf.seq_region_strand, ig.allele_1, ig.allele_2, ig.sample_id, vf.allele_string
	     FROM variation_feature vf, tmp_individual_genotype_single_bp ig
	     WHERE ig.variation_id = vf.variation_id
	     AND vf.map_weight = 1
             };

      if (!$allow_nulls)
      {
        $sql .= qq{AND ig.allele_1 <> 'N'
	         #AND ig.allele_2 <> 'N'
                 };
      }

      $sql .= qq{AND vf.seq_region_id = $seq_region_id
              #AND ig.sample_id=11346
              ORDER BY vf.seq_region_start
	      #ORDER BY vf.seq_region_id, vf.seq_region_start
              };

      dumpSQL($dbVar, $sql);

      
      #information dumping from database
      #system("sort -k 2 -g -o $TMP_DIR/genotype_dump\_$seq_region_id $TMP_DIR/$TMP_FILE");
      system("mv $TMP_DIR/$TMP_FILE $TMP_DIR/genotype_dump\_$seq_region_id");
    }
    open IN, "$TMP_DIR/genotype_dump\_$seq_region_id" or die "Can't open tmp_file : $!";

    while (<IN>) {
      my ($seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $allele_1, $allele_2, $sample_id, $allele_string) = split;
      #same variation and individual but different genotype !!!
      next if (defined $genotypes->{$sample_id}->{region_end} && $seq_region_start == $genotypes->{$sample_id}->{region_end}); 
      #new region for the individual, update coordinates
      if (!defined $genotypes->{$sample_id}->{region_start}) {
	$genotypes->{$sample_id}->{region_start} = $seq_region_start;
	$genotypes->{$sample_id}->{region_end} = $seq_region_end;
      }
      #compare with the beginning of the region if it is within the DISTANCE of compression
      if ((abs($genotypes->{$sample_id}->{region_start} - $seq_region_start) > DISTANCE()) || (abs($seq_region_start - $genotypes->{$sample_id}->{region_end}) > MAX_SHORT)) {
	#snp outside the region, print the region for the sample we have already visited and start a new one
	print_file("$TMP_DIR/compressed_genotype\_$seq_region_id",$genotypes, $seq_region_id, $sample_id);
	delete $genotypes->{$sample_id}; #and remove the printed entry
	$genotypes->{$sample_id}->{region_start} = $seq_region_start;
      }
      #escape characters (tab, new line)
      reverse_alleles($seq_region_strand, $allele_string, \$allele_1, \$allele_2);
      #and write it to the buffer
      if ($seq_region_start != $genotypes->{$sample_id}->{region_start}) {
	#compress information
	$blob = pack ("n",$seq_region_start - $genotypes->{$sample_id}->{region_end} - 1);
	$genotypes->{$sample_id}->{genotypes} .= escape($blob) . $allele_1 . $allele_2;
      } else {
	#first genotype starts in the region_start, not necessary the number
	$genotypes->{$sample_id}->{genotypes} = $allele_1 . $allele_2;
      }
      #and update the new region_end
      $genotypes->{$sample_id}->{region_end} = $seq_region_start; #to avoid nasty effects of indels coordinates
    }
    #print all remaining genotypes and upload the file
    print_file("$TMP_DIR/compressed_genotype\_$seq_region_id",$genotypes, $seq_region_id);
    $genotypes = {};		#and flush the hash
    #need to fork for upload the data
    my $pid = fork;
    if (! defined $pid) {
      throw("Not possible to fork: $!\n");
    } elsif ($pid == 0) {
      #you are the child.....
      my $dbVar_write = $vdba->dbc;
      debug("Import for seq_region_id $seq_region_id");
      &import_genotypes($dbVar_write,"$TMP_DIR/compressed_genotype\_$seq_region_id");
      #unlink "$TMP_DIR/genotype_dump\_$seq_region_id";
    POSIX:_exit(0);
    } else {
      #the parent waits for the child to finish;
      waitpid($pid,0);
    }    
  }

  print "Time finishing dumping data: ",scalar(localtime(time)),"\n";
  #unlink "$TMP_DIR/genotype_dump_*";
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
    my $file = shift;
    
    debug("Importing compressed genotype data");
    my $call = "mv $file $TMP_DIR/$TMP_FILE";
    system($call);
    load($dbVar,qw(compressed_genotype_single_bp sample_id seq_region_id seq_region_start seq_region_end seq_region_strand genotypes));
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
	('INSERT INTO meta_coord (table_name,coord_system_id,max_length) VALUES (?,?,?)');

    $sth->execute($table_name, $cs->dbID(),DISTANCE+1);
    
    $sth->finish();
    
  return;
}


#checks if the genotype is in the opposite strand. If so, reverse the alleles
sub reverse_alleles{
    my $seq_region_strand = shift;
    my $allele_string = shift;
    my $ref_allele_1 = shift;
    my $ref_allele_2 = shift;
    
    my @alleles = split("/",$allele_string);
    #the variation_feature is in the opposite strand, reverse the genotypes

    if ($seq_region_strand == -1){
	reverse_comp($ref_allele_1);
	reverse_comp($ref_allele_2);
#	reverse_comp(\$alleles[0]);
#	reverse_comp(\$alleles[1]);
    }
    #Only reverse genotypes with valid base
#    if ($$ref_allele_1 ne 'N'){
#	#compare the genotype with the forward alleles, if they are different, reverse them
#	if ((@alleles == 2) && (($$ref_allele_1 ne $alleles[0]) && ($$ref_allele_1 ne $alleles[1]))){
#	    reverse_comp($ref_allele_1);
#	    reverse_comp($ref_allele_2);
#	}	
#    }
}

# $special_characters_escaped = printable( $source_string );
sub escape ($) {
  local $_ = ( defined $_[0] ? $_[0] : '' );
  s/([\r\n\t\\\"])/\\$Printable{$1}/sg;
  return $_;
}

sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: perl compress_genotypes.pl <options>

options:
    -tmpdir <dir>        temp directory to use (with lots of space!)
    -tmpfile <filename>   name of temp file to use
    -species <species_name> name of the specie you want to compress the genotypes
EOF

  die("\n$msg\n\n");
}
