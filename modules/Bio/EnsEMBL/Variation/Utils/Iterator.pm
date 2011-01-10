=head1 NAME

  Bio::EnsEMBL::Variation::Utils::Iterator

=head1 SYNOPSIS

  my $variation_iterator = $variation_adaptor->fetch_iterator_by_VariationSet($1kg_set);

  while (my $var = $variation_iterator->next) {
    # operate on variation object
  }


=head1 DESCRIPTION

  Some adaptor methods may return more objects than can fit in memory at once, in these cases 
  you can fetch an iterator object instead of the usual list reference. The iterator object 
  allows you to iterate over the set of objects (using the next() method) without loading the
  entire set into memory at once. The number of objects that are fetched and stored in memory
  can be controlled with the cache_size parameter to the the constructor. 

=head1 LICENSE

 Copyright (c) 1999-2010 The European Bioinformatics Institute and
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

package Bio::EnsEMBL::Variation::Utils::Iterator;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);

my $DEFAULT_CACHE_SIZE = 1000;

=head2 new

  Arg [dbIDs] :
    a listref of dbIDs of the objects you want this iterator to fetch

  Arg [adaptor] :
    Bio::EnsEMBL::Variation::DBSQL::*Adaptor

  Arg [cache_size] :
     The number of objects to fetch and store in memory at once

  Example    :

    my $iterator = Bio::EnsEMBL::Variation::Utils::Iterator->new(
        -dbIDS      => [1234, 2345, 3456],
        -adaptor    => $variation_adaptor,
        -cache_size => 1000,
    );

  Description: Constructor, creates a new iterator object
  Returntype : Bio::EnsEMBL::Variation::Utils::Iterator
  Exceptions : dies if the supplied adaptor does not have a fetch_all_by_dbID_list method
  Caller     : general
  Status     : Experimental

=cut

sub new {
    my $class = shift;

    my %args = @_;

    my ($dbids, $adaptor, $cache_size) = rearrange([qw(dbids adaptor cache_size)], @_);

    unless ($adaptor->can('fetch_all_by_dbID_list')) {
        die "The supplied adaptor does not implement the required fetch_all_by_dbID_list method";
    }

    my $self = {
        dbids       => $dbids,
        adaptor     => $adaptor,
        cache_size  => $cache_size || $DEFAULT_CACHE_SIZE,
        objs        => []
    };

    return bless $self, $class;
}

=head2 next

  Example    : $obj = $iterator->next()
  Description: returns the next object from this iterator, or undef if the iterator is exhausted
  Returntype : object reference (the type will depend on what this iterator is iterating over)
  Exceptions : none
  Caller     : general
  Status     : Experimental

=cut


sub next {
    my $self = shift;

    my @objs  = @{ $self->{objs} };
    my @dbids = @{ $self->{dbids} };

    if (@{ $self->{objs} } == 0 && @{ $self->{dbids} } > 0 ) {
        my @dbids = splice @{ $self->{dbids} }, 0, $self->{cache_size};
        $self->{objs} = $self->{adaptor}->fetch_all_by_dbID_list(\@dbids);
    }

    return shift @{ $self->{objs} };
}

=head2 has_next

  Example    : if ($iterator->has_next) { my $obj = $iterator->next }
  Description: returns true if this iterator has more objects to fetch, false when it is exhausted
  Returntype : boolean
  Exceptions : none
  Caller     : general
  Status     : Experimental

=cut

sub has_next {
    my $self = shift;
    return $self->num_remaining > 0;
}

1;

