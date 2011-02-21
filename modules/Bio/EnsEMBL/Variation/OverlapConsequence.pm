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

package Bio::EnsEMBL::Variation::OverlapConsequence;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::VariationEffect;

sub new {
    my ($class) = @_;
    return bless {}, $class;
}

sub new_fast {
    my ($class, $hashref) = @_;
    return bless $hashref, $class;
}

sub dbID {
    my ($self, $dbID) = @_;
    $self->{dbID} = $dbID if $dbID;
    return $self->{dbID};
}

sub SO_term {
    my ($self, $SO_term) = @_;
    $self->{SO_term} = $SO_term if $SO_term;
    return $self->{SO_term};
}

sub feature_SO_term {
    my ($self, $feature_SO_term) = @_;
    $self->{feature_SO_term} = $feature_SO_term if $feature_SO_term;
    return $self->{feature_SO_term};
}

sub predicate {
    my ($self, $predicate) = @_;
    
    $self->{predicate} = $predicate if $predicate;
    
    if ($self->{predicate} && ref $self->{predicate} ne 'CODE') {
        my $name = $self->{predicate};

        if (defined &$name && $name =~ /^Bio::EnsEMBL::Variation::Utils::VariationEffect/) {
            $self->{predicate} = \&$name;
        }
        else {
            die "Can't find a subroutine called $name in the VariationEffect module?";
        }
    }
    
    return $self->{predicate};
}

sub rank {
    my ($self, $rank) = @_;
    $self->{rank} = $rank if $rank;
    return $self->{rank};
}

sub ensembl_term {
    my ($self, $ensembl_term) = @_;
    $self->{ensembl_term} = $ensembl_term if $ensembl_term;
    return $self->{ensembl_term} || $self->SO_term;
}

sub SO_accession {
    my ($self, $SO_accession) = @_;
    $self->{SO_accession} = $SO_accession if $SO_accession;
    return $self->{SO_accession};
}

sub NCBI_term {
    my ($self, $NCBI_term) = @_;
    $self->{NCBI_term} = $NCBI_term if $NCBI_term;
    return $self->{NCBI_term};
}

sub is_definitive {
    my ($self, $is_definitive) = @_;
    $self->{is_definitive} = $is_definitive if defined $is_definitive;
    return $self->{is_definitive};
}

sub get_all_parent_SO_terms {
    my ($self) = @_;
    
    if (my $adap = $self->{adaptor}) {
        if (my $goa = $adap->db->get_SOTermAdaptor) {
            
        }
    }
}

1;
