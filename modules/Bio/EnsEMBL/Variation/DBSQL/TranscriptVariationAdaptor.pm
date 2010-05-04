#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor

=head1 SYNOPSIS

  # connect to variation database
  $vdb = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(...);

  # connect to core database
  $db  = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);


  # tell the variation database where to obtain core database data
  $vdb->dnadb($db);

  $trva = $vdb->get_TranscriptVariationAdaptor();
  $vfa  = $vdb->get_VariationFeatureAdaptor();
  $tra  = $db->get_TranscriptAdaptor();


  # retrieve a TranscriptVariation by its internal identifier

  $trv = $tra->fetch_by_dbID(552);
  print $trv->transcript()->stable_id(), ' ',
        $trv->variation_feature()->variation_name(), ' ',
        $trv->consequence_type(), "\n";

  # retrieve all TranscriptVariations associated with a Transcript

  $tr = $tra->fetch_by_stable_id('ENST00000278995');
  foreach $trv (@{$trva->fetch_all_by_Transcripts([$tr])}) {
    print $trv->variation_feature->variation_name(), ' ', $trv->consequence_type(), "\n";
  }

  # retrieve all TranscriptVariations associated with a VariationFeature

  $vf = $vfa->fetch_by_dbID(99123);
  foreach $trv (@{$trva->fetch_all_by_VariationFeatures($vf)}) {
    print $trva->transcript->stable_id(), ' ', $trva->consequence_type(), "\n";
  }


=head1 DESCRIPTION

This adaptor provides database connectivity for TranscriptVariation objects.
TranscriptVariations which represent an association between a variation and
a Transcript may be retrieved from the Ensembl variation database via several
means using this module.

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);
use Bio::EnsEMBL::Utils::TranscriptAlleles qw(type_variation);

use Bio::EnsEMBL::Variation::TranscriptVariation;

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');

# size that we look either side of a variation for transcripts
our $UP_DOWN_SIZE = 5000;

=head2 fetch_all_by_Transcripts

    Arg[1]      : listref of Bio::EnsEMBL::Transcript
    Example     : $tr = $ta->fetch_by_stable_id('ENST00000278995');
                  @tr_vars = @{$tr_var_adaptor->fetch_all_by_Transcripts([$tr])};
    Description : Retrieves all TranscriptVariation objects associated with
                  provided Ensembl Transcript. Attaches them to the TranscriptVariation
    ReturnType  : ref to list of Bio::EnsEMBL::Variation::TranscriptVariations
    Exceptions  : throw on bad argument
    Caller      : general
    Status      : Stable

=cut

sub fetch_all_by_Transcripts{
    my $self = shift;
    my $transcript_ref = shift;

    if (ref($transcript_ref) ne 'ARRAY' or ! $transcript_ref->[0]->isa('Bio::EnsEMBL::Transcript')){
      throw('Array Bio::EnsEMBL::Transcript expected');
    }
    
    my %tr_by_id;

    foreach my $tr (@{$transcript_ref}) {
      if (!$tr->isa('Bio::EnsEMBL::Transcript')){
	throw('Bio::EnsEMBL::Transcript is expected');
      }
      $tr_by_id{$tr->dbID()} = $tr;
    }
    my $instr = join (",", keys( %tr_by_id));
    my $transcript_variations = $self->generic_fetch( "tv.transcript_id in ( $instr )");
    for my $tv (@{$transcript_variations}){
	#add to the TranscriptVariation object all the Transcripts
	$tv->{'transcript'} = $tr_by_id{$tv->{'_transcript_id'}};
	delete $tv->{'_transcript_id'}; #remove the transcript_id from the transcript_variation object
    }
    return $transcript_variations;
}

=head2 fetch_all_by_VariationFeatures

  Arg [1]    : Bio::EnsEMBL::Variation::VariationFeature $vf
  Example    :
     $vf = $vf_adaptor->fetch_by_dbID(1234);
     @tr_vars = @{$tr_var_adaptor->fetch_all_by_VariationFeatures([$vf])});
  Description: Retrieves all TranscriptVariation objects associated with
               provided Ensembl variation features. Attaches them to the given variation
               features.
  Returntype : ref to list of Bio::EnsEMBL::Variation::TranscriptVariations
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_VariationFeatures {
  my $self = shift;
  my $vf_ref = shift;

  if(ref($vf_ref) ne 'ARRAY') {
    throw('Array Bio::EnsEMBL::Variation::VariationFeature expected');
  }

  my %vf_by_id;
  my $tvs;
  
  my @have_dbID;
  my @no_dbID;
  
  foreach my $vf(@$vf_ref) {
	if(defined $vf->dbID) {
		push @have_dbID, $vf;
	}
	
	else {
		push @no_dbID, $vf;
	}
  }
  
  # already existing VFs - get TVs from the database
  if(scalar @have_dbID) {

	%vf_by_id = map {$_->dbID(), $_ } @have_dbID;
	my $instr = join (",",keys( %vf_by_id));
	$tvs = $self->generic_fetch( "tv.variation_feature_id in ( $instr )" );
	for my $tv ( @$tvs ) {
	  #add to the variation feature object all the transcript variations
	  $vf_by_id{ $tv->{'_vf_id'} }->add_TranscriptVariation( $tv );
	  #delete $tv->{'_vf_id'}; #remove the variation_feature_id from the transcript_variation object    
	}
  }
  
  # if we have some de novo VFs, we need to calculate
  if(scalar @no_dbID) {
	
	my $species;
	
	if(defined $self->db) {
		$species = $self->db()->species;
	}
	else {
		$species = $self->{'species'};
	}
	
	unless(defined $species) {
		warn("No species defined in adaptor");
		return [];
	}
	
	# get functional genomics adaptors
	my $rf_adaptor = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(-species => $species, -type => "RegulatoryFeature");
	my $ef_adaptor = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(-species => $species, -type => "ExternalFeature");
	
	unless(defined $rf_adaptor && defined $ef_adaptor) {
		warn("Must have functional genomics database attached to consider regulatory features");
	}
	
	# get a gene and transcript adaptor
	my $gene_adaptor = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(-species => $species, -type => "Gene");
	my $transcript_adaptor = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(-species => $species, -type => "Transcript");
	
	unless(defined $gene_adaptor && defined $transcript_adaptor) {
		warn("Must have core database attached to consider regulatory features");
	}
	
	# iterate through the VFs
	foreach my $vf(@no_dbID) {
		
		# sanity check the alleles
		if($vf->allele_string =~ /INS|DEL|LARGE|CNV|PROBE|\(\)/) {
			warn("Can't calculate consequences for alleles ".($vf->allele_string));
			next;
		}
		
		my @this_vf_tvs;
			
		## REGULATORY FEATURES
		######################
		
		if(defined $rf_adaptor && defined $ef_adaptor && defined $gene_adaptor && defined $transcript_adaptor) {
			
			# get the feature slice
			my $slice = $vf->feature_Slice;
			
			# invert it if needed
			if($slice->strand < 0) {
				$slice = $slice->invert();
			}
			
			# hash for storing IDs so we don't create the same TV twice
			my %done;
			
			# array for storing RFs
			my @rf;
			
			# first get external features
			foreach my $f  (@{$ef_adaptor->fetch_all_by_Slice($slice)}) {
				if ($f->feature_set->name =~ /miRanda/ or $f->feature_set->name =~ /VISTA\s+enhancer\s+set/i or $f->feature_set->name =~ /cisRED\s+motifs/i) {
					push @rf, $f;
				}
			}
			
			# get all the regulatory features aswell into the same array
			push @rf, @{$rf_adaptor->fetch_all_by_Slice($slice)};
			
			# now iterate through them all
			foreach my $rf(@rf) {
				
				next unless defined $rf;
				
				# go via gene first
				foreach my $dbEntry (@{$rf->get_all_DBEntries("$species\_core_Gene")}) {
					my $gene;
					
					# get the gene for the stable_id
					foreach my $g(@{$gene_adaptor->fetch_by_stable_id($dbEntry->primary_id)}) {
						if(defined $g && $g->stable_id eq $dbEntry->primary_id) {
							$gene = $g;
							last;
						}
					}
					
					# skip it if no gene found
					next unless defined $gene;
					
					# now get all of this gene's transcripts
					foreach my $tr (@{$gene->get_all_Transcripts()}) {
						
						# skip if we've already seen this transcript
						next if $done{$tr->dbID};
						
						# create a new TV
						my $trv = Bio::EnsEMBL::Variation::TranscriptVariation->new_fast( {
							'adaptor' 			=> $self,
							'consequence_type'	=> ['REGULATORY_REGION'],
							'transcript' 		=> $tr,
							'variation_feature' => $vf,
						} );
						
						# add it to the list we're returning
						push @this_vf_tvs, $trv;
						
						# add it to the VF object
						$vf->add_TranscriptVariation($trv);
						
						# record this transcript as done so we don't do it twice
						$done{$tr->dbID} = 1;
					}
				}
				
				
				# now go via transcript
				foreach my $dbEntry (@{$rf->get_all_DBEntries("$species\_core_Transcript")}) {
					my $tr;
					
					# get the gene for the stable_id
					foreach my $t(@{$transcript_adaptor->fetch_by_stable_id($dbEntry->primary_id)}) {
						if(defined $t && $t->stable_id eq $dbEntry->primary_id) {
							$tr = $t;
							last;
						}
					}
					
					next unless defined $tr;
					
					next if $done{$tr->dbID};
						
					# create a new TV
					my $trv = Bio::EnsEMBL::Variation::TranscriptVariation->new_fast( {
						'adaptor' 			=> $self,
						'consequence_type'	=> ['REGULATORY_REGION'],
						'transcript' 		=> $tr,
						'variation_feature' => $vf,
					} );
					
					# add it to the list we're returning
					push @this_vf_tvs, $trv;
					
					# add it to the VF object
					$vf->add_TranscriptVariation($trv);
						
					# record this in done so we don't do it twice
					$done{$tr->dbID} = 1;
				}
			}
		}
		
		
		
		## TRANSCRIPTS
		##############
		
		# get another slice, expanded to include up/down-stream regions
		my $expanded_slice = $vf->feature_Slice->expand($UP_DOWN_SIZE,$UP_DOWN_SIZE);
		
		# invert it if needed
		if($expanded_slice->strand < 0) {
			$expanded_slice = $expanded_slice->invert();
		}
		
		# get all the transcripts
		my @transcripts = @{$expanded_slice->get_all_Transcripts()};
		
		while(my $transcript = shift @transcripts) {
			
			# expand the allele string
			my $allele_string = $vf->allele_string;
			expand(\$allele_string);
			
			my @alleles = split /\//, $allele_string;
			my $strand = $vf->strand;
			
			# if we need to flip strand to match the transcript
			if ($strand != $transcript->strand()) {
				
				# flip feature onto same strand as transcript
				for (my $i = 0; $i < @alleles; $i++) {
				  reverse_comp(\$alleles[$i]);
				}
				
				$strand = $transcript->strand();
			}
			
			# shift off the reference allele
			shift @alleles;
			
			# convert the start and end to slice coords
			my $start = $vf->start - $expanded_slice->start + 1;
			my $end = $vf->end - $expanded_slice->start + 1;
			
			# create a consequence type object
			my $consequence_type = Bio::EnsEMBL::Variation::ConsequenceType->new(
				$transcript->dbID,
				undef,				# this is the VF dbID, we don't have one
				$start,
				$end,
				$strand,
				\@alleles
			);
			
			my $consequences;
			
			if($allele_string !~ /\+/) {
				$consequences = type_variation($transcript, "", $consequence_type);
			}
			
			foreach my $con(@$consequences) {
				my $trv = Bio::EnsEMBL::Variation::TranscriptVariation->new_fast( {
					'dbID' 				=> undef,
					'adaptor' 			=> $self,
					'cdna_start'		=> $con->cdna_start,
					'cdna_end'			=> $con->cdna_end,
					'translation_start' => $con->aa_start,
					'translation_end'	=> $con->aa_end,
					'pep_allele_string' => join("/", @{$con->aa_alleles || []}),
					'consequence_type'	=> $con->type,
					'transcript' 		=> $transcript,
					'variation_feature' => $vf,
				} );
				
				push @this_vf_tvs, $trv;
				$vf->add_TranscriptVariation($trv);
			}
		}
		
		# if we didn't get any consequences, make an INTERGENIC
		if(!(scalar @this_vf_tvs)) {
			my $trv = Bio::EnsEMBL::Variation::TranscriptVariation->new_fast( {
				'adaptor'			=> $self,
				'consequence_type'	=> ['INTERGENIC'],
			} );
			
			push @this_vf_tvs, $trv;
			$vf->add_TranscriptVariation($trv);
		}
		
		# add all TVs to the total list
		push @{$tvs}, @this_vf_tvs;
	}
  }
  
  return $tvs;
}


=head2 new_fake

  Arg [1]    : string $species
  Example    :
	$tva = Bio::EnsEMBL::Variation::TranscriptVariationAdaptor->new_fake('human');
  Description: Creates a TranscriptVariationAdaptor with no underlying database
			   attached. Should be used only when getting consequence types for
			   species with no variation database available.
  Returntype : Bio::EnsEMBL::Variation::TranscriptVariationAdaptor
  Exceptions : throw if no species given
  Caller     : called from Bio::EnsEMBL::VariationFeatureAdaptor
  Status     : Stable

=cut

sub new_fake {
  my $class = shift;
  my $species = shift;
  
  throw("No species defined") unless defined $species;
  
  my $self = bless {}, $class;
  
  $self->{'species'} = $species;
  
  return $self;
}


#
# internal method responsible for constructing transcript variation objects
# from an executed statement handle.  Ordering of columns in executed statement
# must be consistant with function implementation.
#
sub _objs_from_sth {
  my $self = shift;
  my $sth = shift;

  my ($trv_id, $tr_id, $vf_id, $cdna_start, $cdna_end, $tl_start, $tl_end,
      $pep_allele, $consequence_type);

  $sth->bind_columns(\$trv_id, \$tr_id, \$vf_id, \$cdna_start, \$cdna_end,
                     \$tl_start, \$tl_end, \$pep_allele, \$consequence_type);


  my %tr_hash;
  my %vf_hash;

  my @results;

  # construct all of the TranscriptVariation objects

  while($sth->fetch()) {

      my @consequences = split /,/,$consequence_type;
    
      my $trv = Bio::EnsEMBL::Variation::TranscriptVariation->new_fast
	  ( { 'dbID' => $trv_id,
	      'adaptor' => $self,
	      'cdna_start' => $cdna_start,
	      'cdna_end'   => $cdna_end,
	      'translation_start' => $tl_start,
	      'translation_end' => $tl_end,
	      'pep_allele_string' => $pep_allele,
	      'consequence_type' => \@consequences} );

    $trv->{'_vf_id'} = $vf_id; #add the variation feature
    $trv->{'_transcript_id'} = $tr_id; #add the transcript id
    push @results, $trv;
  }


  return \@results;
}

sub _tables {return ['transcript_variation','tv'];}

sub _columns {
    return qw (tv.transcript_variation_id tv.transcript_id
	       tv.variation_feature_id tv.cdna_start tv.cdna_end
	       tv.translation_start tv.translation_end
	       tv.peptide_allele_string tv.consequence_type
	       );
}
1;
