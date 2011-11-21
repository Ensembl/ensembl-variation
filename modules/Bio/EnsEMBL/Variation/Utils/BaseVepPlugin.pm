=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Variation::Utils::BaseVepPlugin

=head1 SYNOPSIS

    package FunkyPlugin;

    use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

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
  Description: Creates and returns a new instance of this plugin
  Returntype : Bio::EnsEMBL::Variation::Utils::BaseVepPlugin instance (most likely a user-defined subclass)
  Status     : Experimental

=cut

sub new {
    my ($class, $config, @params) = @_;

    # default to analysing Transcripts

    return bless {
        version         => '2.3',
        config          => $config,
        feature_types   => ['Transcript'],
        params          => \@params,
    }, $class;
}

=head2 version

  Arg [1]    : (optional) a version number in the form N.N.N
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
               header to help explain the results of this plugin.
  Returntype : hashref 
  Status     : Experimental

=cut

sub get_header_info {
    my ($self) = @_;
    return undef;
}

sub prefetch {
    my ($self) = @_;
    return undef;
}

=head2 feature_types

  Description: To indicate which types of genomic features a plugin is interested
               in, plugins should return a listref of the types of feature they can deal 
               with. Currently this list can only include 'Transcript', 'RegulatoryFeature'
               and 'MotifFeature'
  Returntype : listref
  Status     : Experimental

=cut

sub feature_types {
    my ($self, $types) = @_;
    $self->{feature_types} = $types if $types;
    return $self->{feature_types};
}

=head2 check_feature_type
  
  Arg[1]     : the feature type as a string (or a reference to an object of the desired type)
  Description: Check if this plugin is interested in a particular feature type 
  Returntype : boolean
  Status     : Experimental

=cut

sub check_feature_type {
    my ($self, $type) = @_;

    # if we're passed an object instead of a type string
    # get the type of reference and use that
    if (ref $type) {
        $type = ref $type;
    }

    unless (defined $self->{feature_types_wanted}->{$type}) {

        # cache the result so we don't have to loop each time

        $self->{feature_types_wanted}->{$type} = 0;

        for my $wanted (@{ $self->feature_types }) {
            if ($type =~ /$wanted/i) {
                $self->{feature_types_wanted}->{$type} = 1;
                last;
            }
        }
    }

    return $self->{feature_types_wanted}->{$type};
}

=head2 run

  Arg[1]     : An instance of a subclass of Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele
  Arg[2]     : A hashref containing all the data that will be printed on this line, keyed by column name
  Description: Run this plugin, this is where most of the plugin logic should reside. 
               When the VEP is about to finish one line of output (for a given variant-feature-allele 
               combination) it will call this method, passing it a relevant subclass of a
               Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele object according to 
               feature types it is interested in, as returned by the feature_types method:
               
                feature type         argument type

                'Transcript'         Bio::EnsEMBL::Variation::TranscriptVariationAllele
                'RegulatoryFeature'  Bio::EnsEMBL::Variation::RegulatoryFeatureVariationAllele
                'MotifFeature'       Bio::EnsEMBL::Variation::MotifFeatureVariationAllele

               Once the plugin has done its analysis it should return the results as a hashref
               with a key for each type of data (which should match the keys described in 
               get_header_info) and a corresponding value for this particular object, or an empty hash if 
               this plugin does not produce any annotation for this object. Any extra data will 
               then be included in the Extra column in the VEP output file. Please refer to the variation 
               API documentation to see what methods are available on each of the possible classes, 
               bearing in mind that common functionality can be found in the 
               VariationFeatureOverlapAllele superclass. 
               
               If the plugin wants to filter this line out of the VEP output it can indicate
               this by returning undef rather than a hashref. Using this mechanism a plugin
               can act as a filter on the VEP output to limit lines to particular events
               of interest.

  Returntype : hashref or undef if the plugin wants to filter out this line
  Status     : Experimental

=cut

sub run {
    my ($self, $vfoa, $line_hash) = @_;
    warn "VEP plugins should implement the 'run' method\n";
    return {};
}

1;

