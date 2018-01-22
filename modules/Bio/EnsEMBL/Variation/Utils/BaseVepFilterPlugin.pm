=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Variation::Utils::BaseVepFilterPlugin

=head1 SYNOPSIS

    # a simple example filter plugin that excludes all
    # lines that are not non-synonymous changes (defined
    # as those which have alternate peptides in their
    # pep_allele_string)

    package NonSynonymousFilter;

    use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepFilterPlugin);

    sub feature_types {
        return ['Transcript'];
    }

    sub include_line {
        my ($self, $tva) = @_;
        
        if (my $pep_alleles = $tva->pep_allele_string) {
            return $pep_alleles =~ /\//;
        }

        return 0;
    }

    1;

=head1 DESCRIPTION

This is a subclass of BaseVepPlugin aimed to make plugins that act like
filters very straightforward to write. Users should subclass this module
and then need only override the include_line method which should return
a true value if the current line should be included, or a false value
to filter the line out. A filter can then be used with the VEP just as
any other plugin by using the --plugin <module name> command line option.

=cut

package Bio::EnsEMBL::Variation::Utils::BaseVepFilterPlugin;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

=head2 include_line

  Arg[1]     : An instance of a subclass of Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele
  Arg[2]     : A hashref containing all the data that will be printed on this line, keyed by column name
  Description: This method should return true if the plugin wants to filter out this line and
               false otherwise.
  Returntype : boolean
  Status     : Experimental

=cut

sub include_line {
    my ($self, $vfoa, $line_hash) = @_;
    warn "VEP filter plugins should implement the 'include_line' method\n";
    return 1;
}

# the following methods just override the BaseVepPlugin methods
# providing default configuration for filter plugins

sub run {
    my ($self, $vfoa, $line_hash) = @_;

    # all run does for filters is check if the plugin wants to 
    # include this line and then returns the appropriate type: 
    # an empty hashref to include the line, or undef to filter
    # out the line
    
    return $self->include_line($vfoa, $line_hash) ? {} : undef;
}

sub get_header_info {
    # we don't want to include anything in the header
    return {};
}

sub feature_types {
    # default to running for all feature types
    return ['Transcript', 'RegulatoryFeature', 'MotifFeature'];
}

1;

