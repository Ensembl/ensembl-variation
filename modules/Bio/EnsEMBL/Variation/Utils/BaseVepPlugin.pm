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

Bio::EnsEMBL::Variation::Utils::BaseVepPlugin

=head1 SYNOPSIS

    package FunkyPlugin;

    use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
    
    sub feature_types {
        return ['Transcript'];
    }

    sub get_header_info {
        return {
            FUNKY_PLUGIN => "Description of funky plugin"
        };
    }

    sub run {
        my ($self, $transcript_variation_allele) = @_;

        my $results = ... # do analysis
        
        return {
            FUNKY_PLUGIN => $results
        };
    }

    1;

=head1 DESCRIPTION

To make writing plugin modules for the VEP easier, get 
your plugin to inherit from this class, override (at least)
the feature_types, get_header_info and run methods to behave
according to the documentation below, and then run the VEP
with your plugin using the --plugin <module name> command
line option.

=cut

package Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use strict;
use warnings;

=head2 new

  Arg [1]    : a VEP configuration hashref
  Arg [>1]   : any parameters passed on the VEP command line, will be stored as a listref in $self->{params}
  Description: Creates and returns a new instance of this plugin
  Returntype : Bio::EnsEMBL::Variation::Utils::BaseVepPlugin instance (most likely a user-defined subclass)
  Status     : Experimental

=cut

sub new {
    my ($class, $config, @params) = @_;

    # default to the current VEP version, and analysing VariationFeatures and
    # Transcripts (which we expect to be the most common usage, this means that
    # the run method will be called with a TranscriptVariationAllele as the 
    # first argument)

    return bless {
        version                 => '2.3',
        feature_types           => ['Transcript'],
        variant_feature_types   => ['VariationFeature'],
        config                  => $config,
        params                  => \@params,
    }, $class;
}

=head2 params_to_hash

  Description: Returns a hashref of parameters as specified on the command line
               in the style "--plugin MyPlugin,param1=value1,param2=value2"
  Returntype : hashref
  Status     : Experimental

=cut

sub params_to_hash {
  my $self = shift;

  if(!exists($self->{_params_hash})) {

    my $params = $self->params;

    my %hash = ();
    return \%hash unless grep {/\=/} @$params;

    foreach my $param(@$params) {
      my ($key, $val) = split('=', $param);
      die("ERROR: Failed to parse parameter $param\n") unless defined($key) && defined($val);
      $hash{$key} = $val;
    }

    $self->{_params_hash} = \%hash;
  }

  return $self->{_params_hash};
}

=head2 version

  Arg [1]    : (optional) a version number string in the form N.N.N
  Description: Get/set the version of this plugin. The version should
               match the version of the VEP that this plugin works with
               (at least in the major version number). This is used to
               detect compatibility between the VEP and plugins.
  Returntype : string
  Status     : Experimental

=cut

sub version {
    my ($self, $version) = @_;
    $self->{version} = $version if $version;
    return $self->{version};
}

=head2 config

  Arg [1]    : a VEP configuration hashref
  Description: Get/set the VEP configuration hashref
  Returntype : hashref 
  Status     : Experimental

=cut

sub config {
    my ($self, $config) = @_;
    $self->{config} = $config if $config;
    return $self->{config};
}

=head2 params

  Arg [1]    : (optional) a listref of plugin parameters
  Description: Get/set the parameters of this plugin, typically as passed on the VEP command line.
  Returntype : listref
  Status     : Experimental

=cut

sub params {
    my ($self, $params) = @_;
    $self->{params} = $params if $params;
    return $self->{params} || [];
}

=head2 get_header_info

  Description: Return a hashref with any Extra column keys as keys and a description
               of the data as a value, this will be included in the VEP output file 
               header to help explain the results of this plugin. Plugins that do
               not want to include anything in the header should return undef.
  Returntype : hashref or undef
  Status     : Experimental

=cut

sub get_header_info {
    my ($self) = @_;
    return undef;
}

=head2 variant_feature_types

  Description: To indicate which types of variation features a plugin is interested
               in, plugins should return a listref of the types of variation feature 
               they can deal with. Currently this list should include one or more of: 
               'VariationFeature' or 'StructuralVariationFeature' 
  Returntype : listref
  Status     : Experimental

=cut

sub variant_feature_types {
    my ($self, $types) = @_;
    $self->{variant_feature_types} = $types if $types;
    return $self->{variant_feature_types};
}

=head2 feature_types

  Description: To indicate which types of genomic features a plugin is interested
               in, plugins should return a listref of the types of feature they can deal 
               with. Currently this list should include one or more of: 'Transcript', 
               'RegulatoryFeature' and 'MotifFeature'
  Returntype : listref
  Status     : Experimental

=cut

sub feature_types {
    my ($self, $types) = @_;
    $self->{feature_types} = $types if $types;
    return $self->{feature_types};
}

sub _check_types {
    my ($self, $type_type, $type) = @_;

    # if we're passed an object instead of a type string
    # get the type of reference and use that
    if (ref $type) {
        $type = ref $type;
    }

    # $type_type will either be 'variant_feature' or 'feature'
    # so construct the method to call and the hash key to
    # store the cached results under
    
    my $method   = $type_type.'_types';
    my $hash_key = $method.'_wanted';

    unless (defined $self->{$hash_key}->{$type}) {

        # cache the result so we don't have to loop each time

        $self->{$hash_key}->{$type} = 0;

        for my $wanted (@{ $self->$method }) {
           
            # special case the intergenic class
            
            if ($wanted eq 'Intergenic') {
                if ($wanted eq $type) {        
                    $self->{$hash_key}->{$type} = 1;
                    last;
                }
            }
            else {

                # we are fairly relaxed about how the user describes features, it can
                # be the fully qualified class name, or just the specific module name, 
                # (i.e. the text after the last '::') in which case we automatically 
                # fully qualify it

                if ($wanted !~ /::/) {

                    if ($type_type eq 'feature') {

                        if ($wanted eq 'RegulatoryFeature' || $wanted eq 'MotifFeature') {
                            $wanted = "Bio::EnsEMBL::Funcgen::$wanted";
                        }
                        else {
                            $wanted = "Bio::EnsEMBL::$wanted";
                        }
                    }
                    elsif ($type_type eq 'variant_feature') {
                        $wanted = "Bio::EnsEMBL::Variation::$wanted";
                    }
                }

                if ($type->isa($wanted)) {
                    $self->{$hash_key}->{$type} = 1;
                    last;
                }
            }
        }
    }

    return $self->{$hash_key}->{$type};
}

=head2 check_feature_type
  
  Arg[1]     : the feature type as a string (or a reference to an object of the desired type)
  Description: Check if this plugin is interested in a particular feature type 
  Returntype : boolean
  Status     : Experimental

=cut

sub check_feature_type {
    my ($self, $type) = @_;
    return $self->_check_types('feature', $type);
}

=head2 check_variant_feature_type
  
  Arg[1]     : the variant feature type as a string (or a reference to an object of the desired type)
  Description: Check if this plugin is interested in a particular variant feature type 
  Returntype : boolean
  Status     : Experimental

=cut

sub check_variant_feature_type {
    my ($self, $type) = @_;
    return $self->_check_types('variant_feature', $type);
}

=head2 run

  Arg[1]     : An instance of a subclass of Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele
  Arg[2]     : A hashref containing all the data that will be printed on this line, keyed by column name
  Description: Run this plugin, this is where most of the plugin logic should reside. 
               When the VEP is about to finish one line of output (for a given variation-allele-feature 
               combination) it will call this method, passing it a relevant subclass of a
               Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele object according to 
               feature types it is interested in, as specified by the plugin's feature_types method:
               
                feature type         argument type

                'Transcript'         Bio::EnsEMBL::Variation::TranscriptVariationAllele
                'RegulatoryFeature'  Bio::EnsEMBL::Variation::RegulatoryFeatureVariationAllele
                'MotifFeature'       Bio::EnsEMBL::Variation::MotifFeatureVariationAllele

               Once the plugin has done its analysis it should return the results as a hashref
               with a key for each type of data (which should match the keys described in 
               get_header_info) and a corresponding value for this particular object, or an empty 
               hash (*not* undef) if this plugin does not produce any annotation for this object. 
               Any edata will then be included in the Extra column in the VEP output file. 
               Please refer to the variation API documentation to see what methods are available 
               on each of the possible classes, bearing in mind that common functionality can be 
               found in the BaseVariationFeatureOverlapAllele superclass. 
               
               If the plugin wants to filter this line out of the VEP output it can indicate
               this by returning undef rather than a hashref. Using this mechanism a plugin
               can act as a filter on the VEP output to limit lines to particular events
               of interest. If you are writing a plugin to act as filter, consider subclassing
               Bio::EnsEMBL::Variation::Utils::BaseVepFilterPlugin, a subclass of this class
               which provides some convenient functionality for filter plugins.

  Returntype : hashref or undef if the plugin wants to filter out this line
  Status     : Experimental

=cut

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    warn "VEP plugins should implement the 'run' method\n";
    return {};
}

1;

