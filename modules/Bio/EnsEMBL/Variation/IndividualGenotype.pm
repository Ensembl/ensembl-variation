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

# Ensembl module for Bio::EnsEMBL::Variation::IndividualGenotype
#
# Copyright (c) 2004 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::IndividualGenotype- Module representing the genotype
of a single individual at a single position

=head1 SYNOPSIS

    print $genotype->variation()->name(), "\n";
    print $genotype->allele1(), '/', $genotype->allele2(), "\n";
    print $genotype->individual()->name(), "\n";

=head1 DESCRIPTION

This is a class representing the genotype of a single diploid individual at
a specific position

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::IndividualGenotype;

use Bio::EnsEMBL::Variation::Genotype;
use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Variation::Genotype Bio::EnsEMBL::Feature);



=head2 new

  Arg [-adaptor] :
    Bio::EnsEMBL::Variation::DBSQL::IndividualAdaptor
  Arg [-START] :
    see superclass constructor
  Arg [-END] :
    see superclass constructor
  Arg [-STRAND] :
    see superclass constructor
  Arg [-SLICE] :
    see superclass constructor
  Arg [-allele1] :
    string - One of the two alleles defining this genotype
  Arg [-allele2] :
    string - One of the two alleles defining this genotype
  Arg [-variation] :
    Bio::EnsEMBL::Variation::Variation - The variation associated with this
    genotype
  Arg [-individual] :
    Bio::EnsEMBL::Individual - The individual this genotype is for.
  Example    : $ind_genotype = Bio:EnsEMBL::Variation::IndividualGenotype->new
                   (-start   => 100,
		    -end     => 100,
		    -strand  => 1,
		    -slice   => $slice,
		    -allele1 => 'A',
                    -allele2 => 'T',
                    -variation => $variation,
                    -individual => $ind);
  Description: Constructor.  Instantiates an IndividualGenotype object.
  Returntype : Bio::EnsEMBL::Variation::IndividualGenotype
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

=cut

sub new {
    my $caller = shift;
    my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my ($adaptor, $allele1, $allele2, $var, $ind) =
    rearrange([qw(adaptor allele1 allele2 variation individual)],@_);

  if(defined($var) &&
     (!ref($var) || !$var->isa('Bio::EnsEMBL::Variation::Variation'))) {
    throw("Bio::EnsEMBL::Variation::Variation argument expected");
  }

  if(defined($ind) &&
     (!ref($ind) || !$ind->isa('Bio::EnsEMBL::Variation::Individual'))) {
    throw("Bio::EnsEMBL::Variation::Individual argument expected");
  }

    $self->{'adaptor'} = $adaptor;
    $self->{'allele1'} = $allele1;
    $self->{'allele2'} = $allele2;
    $self->{'individual'} = $ind;
    $self->{'variation'} = $var;
    
    return $self;

}



sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}


=head2 individual

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Individual $ind
  Example    : $ind = $ind_genotype->individual();
  Description: Getter/Setter for the individual associated with this genotype
  Returntype : Bio::EnsEMBL::Variation::Individual
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

=cut


sub individual {
  my $self = shift;
  if(@_) {
    my $ind = shift;
    if(defined($ind) &&
       (!ref($ind) || !$ind->isa('Bio::EnsEMBL::Variation::Individual'))) {
      throw('Bio::EnsEMBL::Variation::Individual argument expected');
    }
    return $self->{'individual'} = $ind;
  }
  return $self->{'individual'};
}

=head2 variation

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Variation $var
  Example    : $var = $genotype->variation();
  Description: Getter/Setter for the Variation as
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

=cut

sub variation {
  my $self = shift;

  if(@_) {
      #Setter: check wether it is a variation and return it
      my $v = shift;
      if(defined($v) &&
	 (!ref($v) || !$v->isa('Bio::EnsEMBL::Variation::Variation'))) {
	  throw('Bio::EnsEMBL::Variation::Variation argument expected.');
      }
      $self->{'variation'} = $v;
  }
  else{
      if(!defined($self->{'variation'}) && $self->{'adaptor'})    {
		
		#lazy-load from database on demand
		my $vfa = $self->{'adaptor'}->db()->get_VariationFeatureAdaptor();
		%{$vfa->{'_slice_feature_cache'}} = (); #this is ugly, but I have no clue why the cache is not being properly stored...
		
		#print "FS: ", $self->feature_Slice->start, "-", $self->feature_Slice->end, "\n";
		
		# get all VFs on the feature slice
		my $vfs = $vfa->fetch_all_by_Slice($self->feature_Slice());
		
		# if there's only one that must be it
		if(scalar @$vfs == 1) {
		  $self->{'variation'} = $vfs->[0]->variation;
		}
		
		# otherwise we need to check co-ordinates match
		else {
		  foreach my $vf(@$vfs) {
			  #print "VF: ", $vf->variation_name, " ", $vf->seq_region_start, "-", $vf->seq_region_end, "\n";
			  
			  if(defined($self->{_table} && ($self->{_table} eq 'compressed'))) {
				next unless $vf->var_class =~ /mixed|sn[p|v]/;
			  }
			  
			  # only attach if the seq_region_start/end match the feature slice's
			  if($vf->seq_region_start == $self->feature_Slice->start and $vf->seq_region_end == $self->feature_Slice->end) {
				$self->{'variation'} = $vf->variation;
				last;
			  }
		  }
		}
		
		# if we still haven't got one
		if(!defined $self->{'variation'}) {
		  # try getting a bigger slice to find in-dels
		  my $new_slice = $self->feature_Slice->expand(1,1);
		  
		  #print "expanded FS: ", $new_slice->start, "-", $new_slice->end, "\n";
		  
		  # get VFs on the expanded slice
		  my $new_vfs = $vfa->fetch_all_by_Slice($new_slice);
		  
		  # if there's only one that must be it
		  if(scalar @$new_vfs == 1) {
			$self->{'variation'} = $new_vfs->[0]->variation;
		  }
		  
		  # otherwise we need to check start coord matches start coord of original feature slice
		  else {
			foreach my $vf(@$new_vfs) {
				if($vf->seq_region_start == $self->feature_Slice->start) {
					$self->{'variation'} = $vf->variation;
					last;
				}
			}
		  }
		}
		
		# old code just shifts off first variation found!!!
		#my $vf = shift @{$vfa->fetch_all_by_Slice($self->feature_Slice())};
		#$self->{'variation'} = $vf->variation;
      }
  }
  return $self->{'variation'};
}

1;
