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
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);

use Bio::EnsEMBL::Registry;

#use Data::Dumper;

my $species;
my $dump_file;
my $dump_dir;
my $new_db_version;

GetOptions('dump_file=s' => \$dump_file,
	   'dump_dir=s' => \$dump_dir,
	   'species=s'   => \$species,
           'new_db_version=n' => \$new_db_version,
	   );

print "species is $species and new_db_version is $new_db_version and dump_file is $dump_file\n";

usage('You need to enter species,new_db_version as well as the file name where you want to dump the data') if (!defined $dump_file or !$species or !$new_db_version); 


Bio::EnsEMBL::Registry->load_registry_from_db( -host => 'ens-staging',
                                               -db_version => $new_db_version,
                                               -user => 'ensro',
					      );
my $queue = 'normal';
my $memory = "'select[mem>4000] rusage[mem=4000]' -M4000000";

$queue = 'long' if ($species eq 'human');
$memory = "'select[mem>5000] rusage[mem=5000]' -M5000000" if ($species eq 'human');

my $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');

print "dbCore is ",ref($dbCore),"\n";
my $slice_adaptor = $dbCore->get_SliceAdaptor();
my $slices = $slice_adaptor->fetch_all('chromosome');
#find out all possible chromosomes we want to dump and create the different job arrays
print "Time starting to dump data: ", scalar(localtime),"\n";
my $call = "bsub -q $queue -R$memory  -o $dump_dir/out_dump_strain_$species\_$new_db_version -J dump_strain_$species'[1-" . @{$slices} . "]' ./dump_strain_seq.pl -dump_file \"$dump_dir/$dump_file\" -species $species -new_db_version $new_db_version";
system($call);
#print $call,"\n";    

#wait for all the process to go to LSF
sleep(30);
$call = "bsub  -o finish_dump.txt -w 'done(dump_strain_" . $species . ")' -J waiting_process gzip $dump_dir/$dump_file*";
#print $call,"\n";
system($call);

sub usage{
    my $msg = shift;

    print STDERR <<EOF;

usage: perl dump_strain_seq.pl <options>

options:
    -dump_file <filename>    file where you want to dump the resequencing data
    -dump_dir <path>         path of the dump_file
    -species   <species>     species you want to dump data (default = mouse)
    -new_db_version <version number>  release version number, such as 56

EOF
   
die ("\n$msg\n\n");
}
