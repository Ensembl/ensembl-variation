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
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $trva = $reg->get_adaptor("human","variation","transcriptvariation");
  $vfa = $reg->get_adaptor("human","variation","variationfeature");
  $tra = $reg->get_adaptor("human","core","transcript");

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

use Bio::EnsEMBL::Variation::DBSQL::BaseVariationAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);
use Bio::EnsEMBL::Utils::TranscriptAlleles qw(type_variation);

use Bio::EnsEMBL::Variation::TranscriptVariation;

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseVariationAdaptor');

# size that we look either side of a variation for transcripts
our $UP_DOWN_SIZE = 5000;

sub _internal_fetch_all_by_Transcripts{
    my $self = shift;
    my $transcript_ref = shift;
    my $extra_constraint = shift;

    if (ref($transcript_ref) ne 'ARRAY' or ! $transcript_ref->[0]->isa('Bio::EnsEMBL::Transcript')){
      throw('Array Bio::EnsEMBL::Transcript expected');
    }
    
    my %tr_by_id;

    foreach my $tr (@{$transcript_ref}) {
      if (!$tr->isa('Bio::EnsEMBL::Transcript')){
		throw('Bio::EnsEMBL::Transcript is expected');
      }
      $tr_by_id{$tr->stable_id()} = $tr;
    }
	
	# hack if LRG transcript
	my @lrg_trs;
	foreach my $tr(values %tr_by_id) {
		if($tr->stable_id =~ /^LRG\_\d+/) {
			my $slice = $tr->feature_Slice;
			$slice = $slice->expand($UP_DOWN_SIZE, $UP_DOWN_SIZE);
			
			foreach my $vf(@{$slice->get_all_VariationFeatures}) {
				push @lrg_trs, @{$self->_calc_consequences($slice, $tr, $vf)};
			}
		}
	}

    my $instr = join (",", map{"'$_'"} keys( %tr_by_id));
    my $constraint = "tv.transcript_stable_id in ( $instr )";
    $constraint .= " AND $extra_constraint" if $extra_constraint; 
    my $transcript_variations = $self->generic_fetch($constraint);
    for my $tv (@{$transcript_variations}){
		#add to the TranscriptVariation object all the Transcripts
		$tv->{'transcript'} = $tr_by_id{$tv->{'_transcript_stable_id'}};
		delete $tv->{'_transcript_stable_id'}; #remove the transcript_id from the transcript_variation object
    }
	
	# add hacked LRG ones
	push @$transcript_variations, @lrg_trs;
	
    return $transcript_variations;
}

=head2 fetch_all_by_Transcripts

    Arg[1]      : listref of Bio::EnsEMBL::Transcript
    Example     : $tr = $ta->fetch_by_stable_id('ENST00000278995');
                  @tr_vars = @{$tr_var_adaptor->fetch_all_by_Transcripts([$tr])};
    Description : Retrieves all germline TranscriptVariation objects associated with
                  provided Ensembl Transcript. Attaches them to the TranscriptVariation.
    ReturnType  : ref to list of Bio::EnsEMBL::Variation::TranscriptVariations
    Exceptions  : throw on bad argument
    Caller      : general
    Status      : Stable

=cut

sub fetch_all_by_Transcripts {
    my $self = shift;
    my $transcripts = shift;
    my $constraint = 's.somatic = 0';
    return $self->_internal_fetch_all_by_Transcripts($transcripts, $constraint);
}

=head2 fetch_all_somatic_by_Transcripts

    Arg[1]      : listref of Bio::EnsEMBL::Transcript
    Example     : $tr = $ta->fetch_by_stable_id('ENST00000372348');
                  @tr_vars = @{$tr_var_adaptor->fetch_all_by_Transcripts([$tr])};
    Description : Retrieves all somatic TranscriptVariation objects associated with
                  provided Ensembl Transcript. Attaches them to the TranscriptVariation.
    ReturnType  : ref to list of Bio::EnsEMBL::Variation::TranscriptVariations
    Exceptions  : throw on bad argument
    Caller      : general
    Status      : Stable

=cut

sub fetch_all_somatic_by_Transcripts {
    my $self = shift;
    my $transcripts = shift;
    my $constraint = 's.somatic = 1';
    return $self->_internal_fetch_all_by_Transcripts($transcripts, $constraint);
}

=head2 fetch_all_by_VariationFeatures

  Arg [1]    : Bio::EnsEMBL::Variation::VariationFeature $vf
  Arg [2]    : (optional) Bio::EnsEMBL::Transcript $trs
  Example    :
     $vf = $vf_adaptor->fetch_by_dbID(1234);
     @tr_vars = @{$tr_var_adaptor->fetch_all_by_VariationFeatures([$vf])});
  Description: Retrieves all TranscriptVariation objects associated with
               provided Ensembl variation features. Attaches them to the given variation
               features. If a ref to an array of transcripts is specified as a second
	       argument, returns only the TranscriptVariations associated with those
	       transcripts.
  Returntype : ref to list of Bio::EnsEMBL::Variation::TranscriptVariations
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_VariationFeatures {
  my $self = shift;
  my $vf_ref = shift;
  my $tr_ref = shift;

  if(ref($vf_ref) ne 'ARRAY') {
    throw('Array Bio::EnsEMBL::Variation::VariationFeature expected');
  }
  if(defined($tr_ref) && ref($tr_ref) ne 'ARRAY') {
    throw('Array Bio::EnsEMBL::Transcript expected');
  }

  my %vf_by_id;
  my $tvs;
  
  my @have_dbID;
  my @no_dbID;
  my @lrg_vfs;
  
  foreach my $vf(@$vf_ref) {
	if($vf->seq_region_name =~ /^LRG\_\d+/) {
		push @lrg_vfs, $vf;
	}
	
	elsif(defined $vf->dbID) {
		push @have_dbID, $vf;
	}
	
	else {
		push @no_dbID, $vf;
	}
  }
  
  # try lrgs as no_dbIDs
  push @no_dbID, @lrg_vfs;
  
  # already existing VFs - get TVs from the database
  if(scalar @have_dbID) {

	%vf_by_id = map {$_->dbID(), $_ } @have_dbID;
	my $instr = join (",",keys( %vf_by_id));
	
	# If a list of transcripts was specified, get the stable_ids to restrict the query by
	my $tr_condition = "";
	if (defined($tr_ref)) {
	    my @stable_ids = map($_->stable_id(),@{$tr_ref});
	    $tr_condition = " AND tv.transcript_stable_id in ( '" . join("','",@stable_ids) . "' )";
	}
	
	$tvs = $self->generic_fetch( "tv.variation_feature_id in ( $instr )" . $tr_condition);
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
	
	my ($rf_adaptor, $ef_adaptor);
	
	# get a feature set adaptor to check whether there is regulatory build
	my $fs_adaptor = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(-species => $species, -type => "FeatureSet");
	if(defined $fs_adaptor && scalar @{$fs_adaptor->fetch_all_by_type('regulatory')}) {
		
		# get functional genomics adaptors
		$rf_adaptor = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(-species => $species, -type => "RegulatoryFeature");
		$ef_adaptor = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(-species => $species, -type => "ExternalFeature");
	}
	
	#unless(defined $rf_adaptor && defined $ef_adaptor) {
	#	warn("Must have functional genomics database attached to consider regulatory features");
	#}
	
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
		
		if(defined $rf_adaptor && defined $ef_adaptor && defined $gene_adaptor && defined $transcript_adaptor && $vf->seq_region_name !~ /^LRG\_\d+/) {
			
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
			#push @rf, @{$rf_adaptor->fetch_all_by_Slice($slice)};
			
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
		
		# iterate through them
		while(my $transcript = shift @transcripts) {
			
			# use _calc_consequences to get the TVs
			push @this_vf_tvs, @{$self->_calc_consequences($expanded_slice, $transcript, $vf)};
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
		
		# FIXME: Extremely ugly fix to limit the TranscriptVariations to a supplied list of transcripts. Just temporary to make the
		#Êfunctionality consistent when looking at de novo vfs. The filtering should be performed within the code above instead.
		if (!defined($tr_ref)) {
		    push @{$tvs}, @this_vf_tvs;
		}
		else {
		    while (my $tv = shift(@this_vf_tvs)) {
			next unless (grep($_ eq $tv->transcript_stable_id(),@{$tr_ref}));
			push(@{$tvs},$tv);
		    }
		}
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



=head2 _calc_consequences

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Arg [2]    : Bio::EnsEMBL::Transcript $transcript
  Arg [3]    : Bio::EnsEMBL::Variation::VariationFeature $vf
  Example    : my $tvs = $tva->_calc_consequences($slice, $transcript, $vf);
  Description: Internal method used by SNP Effect Predictor to calculate
               consqeuences on the fly. Not intended for external use.
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptVariation
  Exceptions : none
  Caller     : internal
  Status     : at risk

=cut

sub _calc_consequences {
	my $self = shift;
	my $expanded_slice = shift;
	my $transcript = shift;
	my $vf = shift;
	
	my @tvs;
	my $allele_string = $vf->allele_string;
	
	# HGMD
	if($allele_string eq 'HGMD_MUTATION') {
		my $trv = Bio::EnsEMBL::Variation::TranscriptVariation->new_fast( {
			'dbID' 				=> undef,
			'adaptor' 			=> $self,
			'consequence_type'	=> ['HGMD_MUTATION'],
			'transcript' 		=> $transcript,
			'variation_feature' => $vf,
			'_transcript_stable_id' => $transcript->stable_id,
			'_vf_id' => $vf->dbID,
		} );
		
		push @tvs, $trv;
		$vf->add_TranscriptVariation($trv);
		return \@tvs;
	}
	
	# expand the allele string
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
	my $start = $vf->seq_region_start - $expanded_slice->start + 1;
	my $end = $vf->seq_region_end - $expanded_slice->start + 1;
	
	# create a consequence type object
	my $consequence_type = Bio::EnsEMBL::Variation::ConsequenceType->new(
		$transcript->dbID,
		$vf->dbID || undef,
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
			'cdna_start'	    => $con->cdna_start,
			'cdna_end'			=> $con->cdna_end,
			'cds_start'         => $con->cds_start,
			'cds_end'           => $con->cds_end,
			'translation_start' => $con->aa_start,
			'translation_end'	=> $con->aa_end,
			'pep_allele_string' => join("/", @{$con->aa_alleles || []}),
			'consequence_type'	=> $con->type,
			'transcript' 		=> $transcript,
			'variation_feature' => $vf,
			'_transcript_stable_id' => $transcript->stable_id,
			'_vf_id'            => $vf->dbID,
		} );
		
		push @tvs, $trv;
		$vf->add_TranscriptVariation($trv);
	}
	
	return \@tvs;
}

#
# internal method responsible for constructing transcript variation objects
# from an executed statement handle.  Ordering of columns in executed statement
# must be consistant with function implementation.
#
sub _objs_from_sth {
  my $self = shift;
  my $sth = shift;

  my ($trv_id, $tr_stable_id, $vf_id, $cdna_start, $cdna_end, $cds_start, $cds_end, $tl_start, $tl_end,
      $pep_allele, $consequence_type);

  $sth->bind_columns(\$trv_id, \$tr_stable_id, \$vf_id, \$cdna_start, \$cdna_end, \$cds_start, \$cds_end,
                     \$tl_start, \$tl_end, \$pep_allele, \$consequence_type);


  my %tr_hash;
  my %vf_hash;

  my @results;

  # construct all of the TranscriptVariation objects

  while($sth->fetch()) {

      my @consequences = split /,/,$consequence_type;
	  
	  # hack for broken HGMD data in 59
	  if(defined($pep_allele)) {
		$pep_allele = undef if $pep_allele eq 'NULL';
	  }
    
      my $trv = Bio::EnsEMBL::Variation::TranscriptVariation->new_fast
	  ( { 'dbID' => $trv_id,
	      'adaptor' => $self,
	      'cdna_start' => $cdna_start,
	      'cdna_end'   => $cdna_end,
	      'cds_start'  => $cds_start,
		  'cds_end'    => $cds_end,
	      'translation_start' => $tl_start,
	      'translation_end' => $tl_end,
	      'pep_allele_string' => $pep_allele,
	      'consequence_type' => \@consequences} );

    $trv->{'_vf_id'} = $vf_id; #add the variation feature
    $trv->{'_transcript_stable_id'} = $tr_stable_id; #add the transcript stable id
    push @results, $trv;
  }


  return \@results;
}

sub _tables {
    return (
        ['transcript_variation','tv'],
        ['variation_feature', 'vf'],
        ['source', 's']
    );
}

sub _default_where_clause {
  return 'tv.variation_feature_id = vf.variation_feature_id AND vf.source_id = s.source_id';
}

sub _columns {
    return qw (tv.transcript_variation_id tv.transcript_stable_id
	       tv.variation_feature_id tv.cdna_start tv.cdna_end tv.cds_start tv.cds_end
	       tv.translation_start tv.translation_end
	       tv.peptide_allele_string tv.consequence_type
	       );
}
1;
