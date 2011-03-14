package Bio::EnsEMBL::Variation::Pipeline::InitVariationClass;

use strict;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Hive::Process;

use POSIX qw(ceil);

use base ('Bio::EnsEMBL::Hive::Process');

my $DEBUG = 0;

sub fetch_input {
    
    my $self = shift;

    my $reg_file = $self->param('ensembl_registry')
        or die "ensembl_registry is a required parameter";

    my $species = $self->param('species')
        or die "species is a required parameter";
    
    my $num_chunks = $self->param('num_chunks')
        or die "num_chunks is a required parameter";
    
    my $reg = 'Bio::EnsEMBL::Registry';
    
    $reg->load_all($reg_file, 0, 1);

    my $var_dba = $reg->get_DBAdaptor($species, 'variation')
        or die "failed to get variation DBA for $species";
    
    my $aa = $var_dba->get_AttributeAdaptor;
    
    my $dbh = $var_dba->dbc->db_handle;
    
    # first set everything in variation to 'sequence_alteration' by default
    # because sometimes we miss them because there is no variation_feature
    # or any alleles (though this should become unnecessary as we move to the
    # new approach to failing for all species)

    my $default_attrib_id = $aa->attrib_id_for_type_value('SO_term', 'sequence_alteration');

    die "No attrib_id for 'sequence_alteration'" unless defined $default_attrib_id;

    $dbh->do(qq{UPDATE variation SET class_attrib_id = $default_attrib_id});
    
    # now create some temp tables to store the class attribs

    my $temp_var_table = 'temp_variation_class';
    my $temp_var_feat_table = 'temp_variation_feature_class';
    
    $dbh->do(qq{DROP TABLE IF EXISTS $temp_var_table});
    $dbh->do(qq{DROP TABLE IF EXISTS $temp_var_feat_table});
    
    $dbh->do(qq{CREATE TABLE $temp_var_table LIKE variation});
    $dbh->do(qq{CREATE TABLE $temp_var_feat_table LIKE variation_feature});

    $dbh->do(qq{ALTER TABLE $temp_var_table DISABLE KEYS});
    $dbh->do(qq{ALTER TABLE $temp_var_feat_table DISABLE KEYS});

    # now get an ordered list of all the variation_ids
        
    my $get_var_ids_sth = $dbh->prepare(qq{
        SELECT variation_id FROM variation ORDER BY variation_id
    });

    $get_var_ids_sth->execute;

    my @var_ids;

    while (my ($var_id) = $get_var_ids_sth->fetchrow_array) {
        push @var_ids, $var_id;
    }

    # and split them up into as many chunks as requested

    my $num_vars = scalar @var_ids;

    my $chunk_size = ceil($num_vars / $num_chunks);

    my @output_ids;

    while (@var_ids) {

        my $start = $var_ids[0];
        my $stop  = $chunk_size <= $#var_ids ? $var_ids[$chunk_size - 1] : $var_ids[$#var_ids];

        splice(@var_ids, 0, $chunk_size);

        push @output_ids, {
            ensembl_registry    => $reg_file,
            species             => $species,
            variation_id_start  => $start,
            variation_id_stop   => $stop,
            temp_var_table      => $temp_var_table,
            temp_var_feat_table => $temp_var_feat_table,
        };
    }

    $self->param('chunk_output_ids', \@output_ids);

    $self->param(
        'finish_var_class', [{
            ensembl_registry    => $reg_file,
            species             => $species,
            temp_var_table      => $temp_var_table,
            temp_var_feat_table => $temp_var_feat_table,
        }]
    );
}

sub write_output {
    my $self = shift;

    $self->dataflow_output_id($self->param('finish_var_class'), 1);
    $self->dataflow_output_id($self->param('chunk_output_ids'), 2);
}

1;

