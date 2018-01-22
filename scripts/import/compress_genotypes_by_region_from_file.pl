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

  Use this script to re-import files dumped by
  ensembl-variation/scripts/export/dump_compressed_genotypes.pl

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

use strict;
use warnings;

use Getopt::Long;
use ImportUtils qw(loadfile);
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use DBH;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use FindBin qw( $Bin );
use POSIX;
use Data::Dumper;

use constant DISTANCE => 100_000;
use constant MAX_SHORT => 2**16 -1;

$| = 1;

## Backslash Escaping
my %Printable = ( "\\"=>'\\', "\r"=>'r', "\n"=>'n', "\t"=>'t', "\""=>'"' );

my ($TMP_DIR, $TMP_FILE, $species, $registry_file);

GetOptions(
  'tmpdir=s'        => \$TMP_DIR,
  'tmpfile=s'       => \$TMP_FILE,
  'species=s'       => \$species,
  'registry_file=s' => \$registry_file,
);


$TMP_FILE .= "_".$ENV{HOST}.'_'.$$;

my $dump_file = 'compressed_genotype_region.txt_'.$ENV{HOST}.'_'.$$;

warn("Make sure you have an updated ensembl.registry file!\n");

usage('-TMP_DIR argument is required') if(!$TMP_DIR);
usage('-TMP_FILE argument is required') if(!$TMP_FILE);
usage('-species argument is required') if(!$species);

$registry_file ||= $Bin . "/ensembl.registry";

my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all( $registry_file );
$reg->set_reconnect_when_lost();

my $vdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');

my $dbVar = $vdba->dbc;

$ImportUtils::TMP_DIR = $TMP_DIR;
$ImportUtils::TMP_FILE = $TMP_FILE;


compress_genotypes($dbVar);

#reads the genotype and variation_feature data and compress it in one table with blob field
sub compress_genotypes{
  my $dbVar = shift;
  my $buffer;
  my $blob = '';
  my $count = 0;
  
  #information dumping from database
  my $previous_seq_region_id = 0;
  
  my $genotypes = {}; #hash containing all the genotypes
  
  while (<>) {
    my ($variation_name, $seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $variation_id, $subsnp_id, $sample_id, $genotype_code) = split("\t", $_);
    
    next unless defined($sample_id) && defined($genotype_code) && defined($variation_id) && $variation_id > 0 && $genotype_code > 0;
    
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
    next if (defined $genotypes->{$sample_id}->{region_end} && $seq_region_start == $genotypes->{$sample_id}->{region_end}); #same variation and individual but different genotype !!!
    #new region for the individual, update coordinates
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
  
  #print "Time finishing dumping data: ",scalar(localtime(time)),"\n";
  #print last region
  print_file("$TMP_DIR/$dump_file",$genotypes, $previous_seq_region_id);
  #and import remainig genotypes
  &import_genotypes($dbVar);
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
  
  # check for lock file
	my $sleeps = 0;
	
	my $lock_file = $ENV{HOME}.'/.compress_genotypes_lockfile';
	
	while(-e $lock_file) {
		$sleeps++;
		sleep(1);
		
		if($sleeps % 60 == 0) {
			print STDERR
				"WARNING: Process ".
				"has been waiting to insert for ".
				($sleeps / 60)." min".($sleeps > 60 ? 's' : '').
				", you may need to delete the lock file ".$lock_file.
				" if you are sure no LOAD DATA MySQL process is still running\n";
		}
	}
	
	# create lock file
	open TMP, '> '.$lock_file or die "Could not write to lock file\n";
	print TMP '1';
	close TMP;

  #debug("Importing compressed genotype data");
  #debug("mv $TMP_DIR/$dump_file $TMP_DIR/$TMP_FILE");
  
  #$dbVar->do(qq{LOAD DATA INFILE "$TMP_DIR/$dump_file" INTO TABLE compressed_genotype_region});
  
  #my $call = "mv $TMP_DIR/$dump_file $TMP_DIR/$TMP_FILE";
  #system($call);
  loadfile("$TMP_DIR/$dump_file", $dbVar,qw(compressed_genotype_region individual_id seq_region_id seq_region_start seq_region_end seq_region_strand genotypes));
  
  #unlink("$TMP_DIR/$dump_file");
  
  # remove lock file
  unlink($lock_file);
}



# $special_characters_escaped = printable( $source_string );
sub escape ($) {
  local $_ = ( defined $_[0] ? $_[0] : '' );
  s/([\r\n\t\\\"])/\\$Printable{$1}/sg;
  return $_;
}
