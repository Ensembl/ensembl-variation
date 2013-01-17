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
 <helpdesk@ensembl.org>.

=cut

package Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Pipeline::TranscriptFileAdaptor;

use Bio::EnsEMBL::Hive::AnalysisJob;

use base qw(Bio::EnsEMBL::Hive::Process);

sub param {
    my $self = shift;
    
    unless ($self->input_job) {
        # if we don't have an input job, add a dummy one (used when we're not 
        # running as part of a pipeline proper)
        $self->input_job(Bio::EnsEMBL::Hive::AnalysisJob->new);
    }

    return $self->SUPER::param(@_);
}

sub required_param {
    my $self        = shift;
    my $param_name  = shift;
   
    my $param_value = $self->param($param_name, @_);

    die "$param_name is a required parameter" unless defined $param_value;
    
    return $param_value;
}

sub get_transcript_file_adaptor {
    my $self = shift;
    my $transcripts = shift;
    
    unless ($self->{tfa}) {
        $self->{tfa} = Bio::EnsEMBL::Variation::Pipeline::TranscriptFileAdaptor->new(
            fasta_file  => $self->param('fasta_file'),
            transcripts => $transcripts,
        );
    }

    return $self->{tfa};
}

sub get_species_adaptor {
    my ($self, $group) = @_;

    my $species = $self->required_param('species');
    
    return $self->get_adaptor($species, $group);
}

sub get_adaptor {
    my ($self, $species, $group) = @_;
   
    my $dba;

    eval {
        $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $group);
    };

    unless (defined $dba) {
        $self->_load_registry();
        $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $group);
    }
    
    unless (defined $dba) {
        die "Failed to a get DBA for $species and group $group";
    }
    
    return $dba;
}

sub _load_registry {
    my ($self) = @_;
    
    my $reg_file = $self->required_param('ensembl_registry');
    
    Bio::EnsEMBL::Registry->load_all($reg_file, 0, 1);
    
    return;
}

1;
