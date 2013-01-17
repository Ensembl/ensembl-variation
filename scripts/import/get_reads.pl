use strict;

=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

use warnings;


use Bio::Index::Fastq;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Cache;

my ($individual_name,$index_name,$database_name,$encode_region, $reads_file, $reads_names);

{
    GetOptions('-index_name=s' => \$index_name,
#	       '-individual_name=s' => \$individual_name,
#	       '-encode_region=s' => \$encode_region,
	       '-reads_names=s' => \$reads_names,
	       '-reads_file=s' => \$reads_file,
	       '-database_name=s' => \$database_name);

}

#we need to use cache to speed up fetch read sequence and quality data
my $MAX_READS_CACHE_SIZE = 500;
my %CACHE;
tie(%CACHE, 'Bio::EnsEMBL::Utils::Cache', $MAX_READS_CACHE_SIZE);
my $seq;
my $new_index_name = "/tmp/index_name";
# #get read name from database
# my ($name,$start,$end) = split /\:/,$encode_region;
my $sql = qq{select s.query_name from ssahaSNP_feature s where substring_index(substring_index(target_name,'-',2),'-',1) = 'X' limit 10};
## my $sql = qq{select s.query_name from ssahaSNP_feature s, query_match_length_strain q where target_start > $start and substring_index(substring_index(target_name,'-',2),'-',1) = '$name' and target_end < $end and target_start < $end and s.query_name = q.query_name and q.strain_name = 'HuBB'};
 #my $call = qq{mysql -h ens-genomics1 -u ensro $database_name -e "$sql" > [reads_dir]/$encode_region\.$individual_name};
 my $call = qq{mysql -h ens-genomics1 -u ensro $database_name -e "$sql" > $reads_names};
print $call,"\n";
system($call);
system("cp $index_name $new_index_name");
# #open IN, "<[reads_dir]/$encode_region\.$individual_name" || die "Could not open file with reads:$!\n";
# open IN, "<[reads_dir]/$encode_region\.$individual_name\_yuan" || die "Could not open file with reads:$!\n";
# #open OUT , ">[reads_dir]/$encode_region\.$individual_name.out" || die "Could not open output file $!\n";
# open OUT , ">/[reads_dir]/$encode_region\.$individual_name\_yuan.out" || die "Could not open output file $!\n";
open IN, "<$reads_names" || die "Could not open file with reads:$!\n";
open OUT, ">$reads_file"|| die "Could not open output file:$!\n";
<IN>; #get rid of first line, the column name
while (<IN>){
    chomp;
#    next if ($_ eq 'gnl|ti|1823655463'); #this read is not present in the files: removed manually from fasta file
    $seq = get_pfetch_sequence_and_quality($index_name,$_);
    print OUT ">$_\n";
    print OUT $seq->seq,"\n";
}
close OUT;
close IN;
#unlink "$new_index_name";



sub get_pfetch_sequence_and_quality {
    my $index_filename = shift;
    my $query_name = shift;
    my $fastq_index = Bio::Index::Fastq->new(-filename => $index_filename);

    if (! defined($fastq_index)) {
	$fastq_index = Bio::Index::Fastq->new(-filename => $index_filename);
    }

    if (exists($CACHE{$query_name})) {
        return $CACHE{$query_name};
    }
    
    my $seq = $fastq_index->fetch($query_name);
    die "could not fetch $query_name" unless ($seq);
    
    $CACHE{$query_name} = $seq;

    return $seq;
  }
