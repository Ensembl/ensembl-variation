#!/usr/local/ensembl/bin/perl -w

use strict;
use warnings;
use FindBin qw( $Bin );
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);

use Bio::EnsEMBL::Registry;

#use Data::Dumper;

my $species;
my $dump_file;

GetOptions('dump_file=s' => \$dump_file,
	   'species=s'   => \$species,
	   );

warn("Make sure you have a updated ensembl.registry file!\n");

$species ||= 'mouse'; #by default, dump mouse data
usage('You need to enter the file name where you want to dump the data') if (!defined $dump_file); 
my $registry_file ||= $Bin . "/ensembl.registry";

Bio::EnsEMBL::Registry->load_all( $registry_file );


my $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');

my $slice_adaptor = $dbCore->get_SliceAdaptor();
my $slices = $slice_adaptor->fetch_all('chromosome');
#find out all possible chromosomes we want to dump and create the different job arrays
print "Time starting to dump data: ", scalar(localtime),"\n";
my $call = "bsub  -J dump_strain'[1-" . @{$slices} . "]' -q long /usr/local/bin/perl dump_strain_seq.pl -dump_file $dump_file -species $species";
system($call);
#print $call,"\n";    

#wait for all the process to go to LSF
sleep(120);
$call = "bsub  -o finish_dump.txt -w 'done(dump_strain)' -J waiting_process date";
#print $call,"\n";
system($call);

sub usage{
    my $msg = shift;

    print STDERR <<EOF;

usage: perl dump_strain_seq.pl <options>

options:
    -dump_file <filename>    file where you want to dump the resequencing data (including path)
    -species   <species>     species you want to dump data (default = mouse)

EOF
   
die ("\n$msg\n\n");
}
