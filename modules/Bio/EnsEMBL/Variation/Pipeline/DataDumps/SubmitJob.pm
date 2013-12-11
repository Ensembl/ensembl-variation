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
package Bio::EnsEMBL::Variation::Pipeline::DataDumps::SubmitJob;

use strict;
use FileHandle;
use base ('Bio::EnsEMBL::Hive::Process');

=begin
Call script for dumping with all necessaary parameters using bsub.
	- Check return code
	- grep output and error files for result messages: 'Successfully', Exit with out of memory?...
	
=end
=cut

sub fetch_input {}

sub run {
    my $self = shift;
    my $script          = $self->param('script'); 
    my $connection_args = $self->param('connection_args');
    my $species         = '--species ' . $self->param('species');
    my $seq_regions     = $self->param('seq_region_range');
    my $script_args     = $self->param('script_args');
    my $output_file     = $self->param('output_file');

	my $gvf_file        = $self->param('gvf_file');	
	my $vcf_file        = $self->param('vcf_file');

	my $err             = $self->param('err');
	my $out             = $self->param('out');
    
    my $debug = $self->param('debug') ? '--debug' : '';

	my $cmd = 'perl ' . join(' ', 
		$script, 
		$connection_args, 
		$species, 
		$seq_regions, 
		$output_file, 
		$script_args, 
		$gvf_file,
		$vcf_file,			
		$debug);
	if ($self->param('debug')) {
		my $fh = FileHandle->new($out, 'w');	
		print $fh $cmd, "\n";
		$fh->close();			
	} else {	
		$self->run_cmd("$cmd 1>$out 2>$err");					
	}
    return 1;
}

sub write_output {
    my $self = shift;
}

sub run_cmd {
	my $self = shift;
	my $cmd = shift;
	if (my $return_value = system($cmd)) {
		$return_value >>= 8;
		die "system($cmd) failed: $return_value";
	}
}


1;
