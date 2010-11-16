package Bio::EnsEMBL::Variation::Pipeline::InitJobs;

use strict;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::EnsEMBL::Hive::Process');

my $DEBUG = 0;

sub fetch_input {
   
    my $self = shift;

    my $pph_dir = $self->param('pph_dir')
        or die "pph_dir is a required parameter";

    my $reg_file = $self->param('ensembl_registry')
        or die "ensembl_registry is a required parameter";

    my $species = $self->param('species')
        or die "species is a required parameter";
    
    my $reg = 'Bio::EnsEMBL::Registry';
    $reg->load_all($reg_file, 0, 1);
  
    my $core_dba = $reg->get_DBAdaptor($species, 'core')
        or die "failed to get core DBA for $species";
   
    my $var_dba = $reg->get_DBAdaptor($species, 'variation')
        or die "failed to get variation DBA for $species";

    # drop the indexes on tables we're going to insert into as
    # this significantly speeds up the TranscriptEffect process
    
    my $dbh = $var_dba->dbc->db_handle;

    $dbh->do(q{
        DROP INDEX feature_type_idx         ON vf_overlap;
        DROP INDEX variation_feature_idx    ON vf_overlap;
        DROP INDEX feature_stable_idx       ON vf_overlap;
        DROP INDEX vf_overlap_idx           ON vf_overlap_allele;
        DROP INDEX type_val_idx             ON vf_overlap_attrib;
        DROP INDEX val_only_idx             ON vf_overlap_attrib;
        DROP INDEX overlap_idx              ON vf_overlap_attrib;
        DROP INDEX type_val_idx             ON vf_overlap_allele_attrib;
        DROP INDEX val_only_idx             ON vf_overlap_allele_attrib;
        DROP INDEX allele_idx               ON vf_overlap_allele_attrib;
    }) or die "Failed to drop indexes: ".$dbh->errstr;

    my $sa = $core_dba->get_SliceAdaptor
        or die "Failed to get slice adaptor";

    my @gene_stable_ids;
    
    my $gene_count = 0;
    
    CHR : for my $chr (@{ $sa->fetch_all('chromosome') }) {
        for my $gene (@{ $chr->get_all_Genes }) {
            $gene_count++;
            push @gene_stable_ids, $gene->stable_id;
            
            if ($DEBUG) {
                last CHR if $gene_count >= 100;
            }
        }
    }
    
    my @output_ids = map { { 
            gene_stable_id => $_, 
            pph_dir => $pph_dir,
            ensembl_registry => $reg_file,
        } 
    } @gene_stable_ids;

    #@output_ids = ({ gene_stable_id => 'ENSG00000139618', pph_dir => $pph_dir, ensembl_registry => $reg_file });
    #push @output_ids, {gene_stable_id => 'ENSG00000171368', pph_dir => $pph_dir};

    $self->param('output_ids', \@output_ids);
}

sub write_output {
    my $self = shift;
    
    my $output_ids = $self->param('output_ids');

    $self->dataflow_output_id($output_ids, 2);
}

1;
