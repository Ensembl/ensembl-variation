
use strict;
use warnings;

use lib '/nfs/acari/dr2/projects/src/branch-normalized-alleles/ensembl-variation/modules';
use lib '/nfs/acari/dr2/projects/src/branch-normalized-alleles/ensembl/modules';


use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Data::Dumper;

my $dbVariation = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new
    (-host => 'ecs2',
     -user => 'ensro',
     -pass => '',
     -port => 3365,
     -dbname => 'mus_musculus_variation_30_33f'
     );

my $dbCore = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (-host => 'ecs2',
     -user => 'ensro',
     -pass => '',
     -port => 3365,
     -dbname => 'mus_musculus_core_30_33f'
     );
	
my $slice_adaptor = $dbCore->get_SliceAdaptor();
#my $slice = $slice_adaptor->fetch_by_region("chromosome", '14',113798260,113798272); #deletion
my $slice = $slice_adaptor->fetch_by_region("chromosome", '2',71417810,71417820); #insertion
#my $slice = $slice_adaptor->fetch_by_region("chromosome", '11',53512843,53512863); #a couple of SNPs
print "Sequence for the slice: ",$slice->seq,"\n";
#my $sliceStrain = $slice->get_by_strain('CAST/Ei'); #deletion and the couple of SNPs
my $sliceStrain = $slice->get_by_strain('FVB'); #insertion
#my $sliceStrain2 = $slice->get_by_strain('SPRET/Ei'); #deletion
my $sliceStrain2 = $slice->get_by_strain('BALB/cJ'); #insertion and the snps
print "strain 1: ",$sliceStrain->seq,"\n";
print "strain 2: ",$sliceStrain2->seq,"\n";
my $differences = $sliceStrain->get_all_differences_StrainSlice($sliceStrain2);
foreach my $difference (@{$differences}){
    print "Difference in: ",$difference->start,"-",$difference->end," with allele ",$difference->allele_string,"\n";
}
my $subSlice_strain = $sliceStrain->sub_Slice(2,9,1);
my $subSlice = $slice->sub_Slice(2,9,1);
print "SubSequence for the slice: ",$subSlice->seq,"\n";
print "SuSequence for the strain: ",$subSlice_strain->seq,"\n";
print "subseq:                    ",$sliceStrain->subseq(2,9,1),"\n";




