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

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::CompressedGenotypeAdaptor
#
# Copyright (c) 2005 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::CompressedGenotypeAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $iga = $reg->get_adaptor("human","variation","individualgenotype");

  #returns all genotypes in a certain Slice

  $genotypes = $iga->fetch_by_Slice($slice);



=head1 DESCRIPTION

This adaptor provides database connectivity for IndividualGenotype objects.
IndividualGenotypes may be retrieved from the Ensembl variation database by
several means using this module.

=head1 METHODS

=cut

package Bio::EnsEMBL::Variation::DBSQL::CompressedGenotypeAdaptor;

use strict;
use warnings;

use vars qw(@ISA);

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Variation::IndividualGenotype;
use Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor;

@ISA = qw(Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor);

use Data::Dumper;

=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL::Variation $variation
  Example    : my $var = $variation_adaptor->fetch_by_name( "rs1121" )
               $igtypes = $igtype_adaptor->fetch_all_by_Variation( $var )
  Description: Retrieves a list of individual genotypes for the given Variation.
               If none are available an empty listref is returned.
  Returntype : listref Bio::EnsEMBL::Variation::IndividualGenotype 
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut


sub fetch_all_by_Variation {
    my $self = shift;
    my $variation = shift;
	my $individual = shift;

    my $res;
    if(!ref($variation) || !$variation->isa('Bio::EnsEMBL::Variation::Variation')) {
	throw('Bio::EnsEMBL::Variation::Variation argument expected');
    }

    if(!defined($variation->dbID())) {
	warning("Cannot retrieve genotypes for variation without set dbID");
	return [];
    }	
    my $vfa = $self->db->get_VariationFeatureAdaptor();
    if (!$vfa){
	throw("Cannot retrieve genotypes for variation without adaptor set");
	return [];
    }
	
	# only get compressed genotypes if this is a SNP
	if($variation->var_class =~ /snp/i) {
	
		my $variation_features = $vfa->fetch_all_by_Variation($variation);
		#foreach of the hitting variation Features, get the Genotype information
		foreach my $vf (@{$variation_features}){
			
			# skip this VF if the start and end are >1 apart
			# they should not be in the compressed table
			next if abs($vf->end - $vf->start) > 1;
			
			# get the feature slice for this VF
			my $fs = $vf->feature_Slice();
			
			# if the feature slice is start > end
			if($fs->start > $fs->end) {
				
				# get a new slice with the start and end the right way round
				# otherwise the call won't pick any variations up
				my $new_fs = $fs->{'adaptor'}->fetch_by_region($fs->coord_system->name,$fs->seq_region_name,$fs->end,$fs->start);
				$fs = $new_fs;
			}
			
			my $include_multi = 1;
			$include_multi = 0 if $variation->var_class =~ /snp/i;
			
			# get the IGs
			#my @igs = @{$self->fetch_all_by_Slice($fs, $individual, $include_multi)};
			
			#print "fS: ", $fs->start, " ", $fs->end, "\n";
			
			# iterate through to check
			#foreach my $ig(@igs) {
			#	#print $ig->variation->dbID, " ", $variation->dbID, "\n";
			#	
			#	# skip this if the variation attached to the IG does not match the query
			#	#next unless $ig->variation->dbID == $variation->dbID;
			#	
			#	# get the alleles
			#	my ($a1, $a2) = ($ig->allele1, $ig->allele2);
			#	
			#	# skip if the returned alleles are not in the allele_string for the VF
			#	#next unless $vf->allele_string =~ /^$a1\/|\/$a1\/|\/$a1$|^$a1$/ and $vf->allele_string =~ /^$a2\/|\/$a2\/|\/$a2$|^$a2$/;
			#	
			#	#$ig->variation($variation);
			#	push @{$res}, $ig;
			#}
			
			# old code without checks
			map {$_->variation($variation); push @{$res}, $_} @{$self->fetch_all_by_Slice($fs, $individual, $include_multi)};
		}
	}
	
	#print "Got ", (defined $res ? scalar @{$res} : 0), " from single bp\n";
	
    #and include the genotypes from the multiple genotype table
    $self->_multiple(1);
    push @{$res}, @{$self->SUPER::fetch_all_by_Variation($variation, $individual)} unless $variation->var_class =~ /snp/i;
    $self->_multiple(0);
	
	#print "Now have ", scalar @{$res}, " including multiple bp\n";
	
    return $res;
}

=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL:Slice $slice
  Arg [2]    : (optional) Bio::EnsEMBL::Variation::Individual $individual
  Example    : my @IndividualGenotypesFeatures = @{$ca->fetch_all_by_Slice($slice)};
  Description: Retrieves all IndividualGenotypeFeature features for a given slice for
               a certain individual (if provided). 
  Returntype : reference to list Bio::EnsEMBL::Variation::IndividualGenotype
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Slice{
    my $self = shift;
    my $slice = shift;
    my $individual = shift;
	my $include_multi = shift;
    my @results;
    my $features;
    my $constraint;
    if (!$self->_multiple){
	#if passed inividual or population, add constraint
	if (defined $individual && defined $individual->dbID){
	  my $instr;
	  
	  if($individual->isa("Bio::EnsEMBL::Variation::Population")) {
		my $inds = $individual->get_all_Individuals;
		my @list;
		push @list, $_->dbID foreach @$inds;
		$instr = (@list > 1)  ? " IN (".join(',',@list).")"   :   ' = \''.$list[0].'\'';
		$constraint = " c.sample_id $instr";
	  }
	  else {
		$constraint = ' c.sample_id = ' . $individual->dbID;
	  }
	  
	 # $constraint = ' c.sample_id = ?';
	 # $self->bind_param_generic_fetch($individual->dbID,SQL_INTEGER);
	  $features = $self->SUPER::fetch_all_by_Slice_constraint($slice,$constraint);
	}
	else{
	    $features = $self->SUPER::fetch_all_by_Slice($slice);
	}
	#need to check the feature is within the Slice
	
	my $seq_region_slice = $slice->seq_region_Slice;

	foreach my $indFeature (@{$features}){
	  #print "feature_start ",$indFeature->start," slice_end ",$slice->end," slice_start ",$slice->start," feature_end ",$indFeature->end, " a: ", $indFeature->allele1, "|", $indFeature->allele2, " in ", $indFeature->individual->name, "\n" if ($indFeature->end==1);
	    if ($indFeature->start > 0 && ($slice->end-$slice->start +1) >= $indFeature->end){
		
		# not sure we need this check now???
		#next unless defined $indFeature->variation;
			
		if ($indFeature->slice->strand == -1){ #ignore the different strand transformation

		  # Position will change if the strand is negative so change the strand to 1 temporarily
		    $indFeature->slice->{'strand'} = 1;
		    my $newFeature = $indFeature->transfer($seq_region_slice); 
		    $indFeature->slice->{'strand'} = -1;
		    $newFeature->slice->{'strand'} = -1;
                    $newFeature->variation($indFeature->variation);
		    push @results, $newFeature;
		}
		else{
		    push @results,$indFeature->transfer($seq_region_slice);
		}
	    }
		#else {
		#	print "ignored\n";
		#}
	}
	
	$self->_multiple(1);
	push @results, @{$self->fetch_all_by_Slice($slice,$individual)} if $include_multi;
	$self->_multiple(0);
	
    }
    else{
	#if passed inividual, add constraint
	if (defined $individual && defined $individual->dbID){
	  $constraint = ' ig.sample_id = ' . $individual->dbID;
	 # $constraint = ' c.sample_id = ?';
	 # $self->bind_param_generic_fetch($individual->dbID,SQL_INTEGER);
	  $features = $self->SUPER::fetch_all_by_Slice_constraint($slice,$constraint);
	}
	else{
	    $features = $self->SUPER::fetch_all_by_Slice($slice);
	}
	#and include the genotypes from the multiple genotype table
	push @results, @$features;
    }
	
    return \@results;
    
}

sub _tables{
    my $self = shift;
    return (['compressed_genotype_single_bp','c']) if (!$self->_multiple);
    return $self->SUPER::_tables if ($self->_multiple);
}

sub _columns{
    my $self = shift;

    return qw(sample_id seq_region_id seq_region_start seq_region_end seq_region_strand genotypes) if (!$self->_multiple);
    return $self->SUPER::_columns if ($self->_multiple);
}

sub _objs_from_sth{
    my ($self, $sth, $mapper, $dest_slice) = @_;
    
    return $self->SUPER::_objs_from_sth($sth,$mapper,$dest_slice) if ($self->_multiple);
    #
    # This code is ugly because an attempt has been made to remove as many
    # function calls as possible for speed purposes.  Thus many caches and
    # a fair bit of gymnastics is used.
    #
    
    my $sa = $self->db()->dnadb()->get_SliceAdaptor();
    
    my @results;
    my %slice_hash;
    my %sr_name_hash;
    my %sr_cs_hash;
    my %individual_hash;

  my ($sample_id, $seq_region_id, $seq_region_start,
      $seq_region_end, $seq_region_strand, $genotypes);

  $sth->bind_columns(\$sample_id, \$seq_region_id,
                     \$seq_region_start, \$seq_region_end, \$seq_region_strand,
                     \$genotypes);

  my $asm_cs;
  my $cmp_cs;
  my $asm_cs_vers;
  my $asm_cs_name;
  my $cmp_cs_vers;
  my $cmp_cs_name;
  if($mapper) {
    $asm_cs = $mapper->assembled_CoordSystem();
    $cmp_cs = $mapper->component_CoordSystem();
    $asm_cs_name = $asm_cs->name();
    $asm_cs_vers = $asm_cs->version();
    $cmp_cs_name = $cmp_cs->name();
    $cmp_cs_vers = $cmp_cs->version();
  }

  my $dest_slice_start;
  my $dest_slice_end;
  my $dest_slice_strand;
  my $dest_slice_length;
  if($dest_slice) {
    $dest_slice_start  = $dest_slice->start();
    $dest_slice_end    = $dest_slice->end();
    $dest_slice_strand = $dest_slice->strand();
    $dest_slice_length = $dest_slice->length();
  }

  FEATURE: while($sth->fetch()) {
    #get the slice object
    my $slice = $slice_hash{"ID:".$seq_region_id};
    if(!$slice) {
      $slice = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{"ID:".$seq_region_id} = $slice;
      $sr_name_hash{$seq_region_id} = $slice->seq_region_name();
      $sr_cs_hash{$seq_region_id} = $slice->coord_system();
    }
	
    #
    # remap the feature coordinates to another coord system
    # if a mapper was provided
    #
    if($mapper) {
      my $sr_name = $sr_name_hash{$seq_region_id};
      my $sr_cs   = $sr_cs_hash{$seq_region_id};

      ($sr_name,$seq_region_start,$seq_region_end,$seq_region_strand) =
        $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end,
                          $seq_region_strand, $sr_cs);

      #skip features that map to gaps or coord system boundaries
      next FEATURE if(!defined($sr_name));

      #get a slice in the coord system we just mapped to
      if($asm_cs == $sr_cs || ($cmp_cs != $sr_cs && $asm_cs->equals($sr_cs))) {
        $slice = $slice_hash{"NAME:$sr_name:$cmp_cs_name:$cmp_cs_vers"} ||=
          $sa->fetch_by_region($cmp_cs_name, $sr_name,undef, undef, undef,
                               $cmp_cs_vers);
      } else {
        $slice = $slice_hash{"NAME:$sr_name:$asm_cs_name:$asm_cs_vers"} ||=
          $sa->fetch_by_region($asm_cs_name, $sr_name, undef, undef, undef,
                               $asm_cs_vers);
      }
  }

    #
    # If a destination slice was provided convert the coords
    # If the dest_slice starts at 1 and is foward strand, nothing needs doing
    #
    if($dest_slice) {
	if($dest_slice_start != 1 || $dest_slice_strand != 1) {
	    if($dest_slice_strand == 1) {
		$seq_region_start = $seq_region_start - $dest_slice_start + 1;
		$seq_region_end   = $seq_region_end   - $dest_slice_start + 1;
	    } else {
		my $tmp_seq_region_start = $seq_region_start;
#		$seq_region_start = $dest_slice_end - $seq_region_end + 1;
		$seq_region_start = $seq_region_start - $dest_slice_start + 1;
		$seq_region_end   = $dest_slice_end - $tmp_seq_region_start + 1;
#		$seq_region_strand *= -1;
	    }
	    
	    #throw away features off the end of the requested slice
	    if($seq_region_end < 1 || $seq_region_start > $dest_slice_length) {
		next FEATURE;
	    }
	}
	$slice = $dest_slice;
    }	
	
	my @genotypes = unpack '(aan)*', $genotypes;
	my $snp_start = $seq_region_start;
  
	while( my( $allele_1, $allele_2, $gap ) = splice @genotypes, 0, 3 ) {
		my $igtype  = Bio::EnsEMBL::Variation::IndividualGenotype->new_fast({
			'start'   => $snp_start,
			'end'     => $snp_start,
			'strand'  => $seq_region_strand,
			'slice'   => $slice,
			'allele1' => $allele_1,
			'allele2' => $allele_2,
			'adaptor' => $self
		});
		
		$igtype->{_table} = 'compressed';
		
		$individual_hash{$sample_id} ||= [];
		push @{$individual_hash{$sample_id}}, $igtype;
		push @results, $igtype;
		$snp_start += $gap + 1 if defined $gap;

    }

	}
     # get all individual in one query (faster)
     # and add to already created genotypes
     my @ind_ids = keys %individual_hash;
     
     my $ia = $self->db()->get_IndividualAdaptor();
     my $inds = $ia->fetch_all_by_dbID_list(\@ind_ids);
     
     foreach my $i (@$inds) {
	 foreach my $igty (@{$individual_hash{$i->dbID()}}) {
	     $igty->individual($i);
	 }
     }

     return \@results;
}

sub _default_where_clause  {
    my $self = shift;
    return '' if (!$self->_multiple);
    return $self->SUPER::_default_where_clause if ($self->_multiple);
}

1;
