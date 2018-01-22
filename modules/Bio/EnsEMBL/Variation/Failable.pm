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

# Ensembl module for Bio::EnsEMBL::Variation::BaseStructuralVariation
#
#

=head1 NAME

Bio::EnsEMBL::Variation::Failable

=head1 DESCRIPTION

This object returns failed information.

=cut


use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

package Bio::EnsEMBL::Variation::Failable;


=head2 failed_description

  Arg [1]    : $failed_description (optional)
               The new value to set the failed_description attribute to. Should 
               be a reference to a list of strings, alternatively a string can
               be passed. If multiple failed descriptions are specified, they should
               be separated with semi-colons.  
  Example    : $failed_str = $sv->failed_description();
  Description: Get/Sets the failed description for this variant. The failed
               descriptions are lazy-loaded from the database.
  Returntype : String (semi-colon separated)
  Exceptions : Thrown on illegal argument.
  Caller     : general
  Status     : Stable

=cut

sub failed_description {
    my $self = shift;
    my $description = shift;
  
    # Update the description if necessary
    if (defined($description)) {
        
        # If the description is a string, split it by semi-colon and take the reference
        if (check_ref($description,'STRING')) {
            my @pcs = split(/;/,$description);
            $description = \@pcs;
        }
        # Throw an error if the description is not an arrayref
        assert_ref($description.'ARRAY');
        
        # Update the cached failed_description
        $self->{'failed_description'} = $description;
    }
    # Else, fetch it from the db if it's not cached
    elsif (!defined($self->{'failed_description'})) {
        $self->{'failed_description'} = $self->get_all_failed_descriptions();
    }
    
    # Return a semi-colon separated string of descriptions
    return join(";",@{$self->{'failed_description'}});
}


=head2 get_all_failed_descriptions

  Example    :  
                if ($sv->is_failed()) {
                    my $descriptions = $sv->get_all_failed_descriptions();
                    print "Structural variant " . $sv->variation_name() . " has been flagged as failed because '";
                    print join("' and '",@{$descriptions}) . "'\n";
                }
                
  Description: Gets all failed descriptions associated with this Structural variant.
  Returntype : Listref of strings
  Exceptions : Thrown if an adaptor is not attached to this object.
  Caller     : general
  Status     : At risk

=cut

sub get_all_failed_descriptions {
  my $self = shift;
  
    # If the failed descriptions haven't been cached yet, load them from db
    unless (defined($self->{'failed_description'})) {
        
        # Check that this allele has an adaptor attached
        # unless () {
        #     throw('An adaptor must be attached to the ' . ref($self)  . ' object');
        # }
    
        $self->{'failed_description'} = defined($self->adaptor()) ? $self->adaptor->get_all_failed_descriptions($self) : [];
    }
    
    return $self->{'failed_description'};
}


=head2 is_failed

  Example    : print "Structural variant '" . $sv->variation_name() . "' has " . ($sv->is_failed() ? "" : "not ") . "been flagged as failed\n";
  Description: Gets the failed attribute for this structural variant. The failed attribute
               is lazy-loaded from the database.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub is_failed {
  my $self = shift;
  
  return (length($self->failed_description()) > 0);
}
1;
