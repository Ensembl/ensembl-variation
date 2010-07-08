use strict;
use warnings;
use LWP::UserAgent;

package LRGImport;

our $dbCore;


# Insert entries into the analysis table if logic_name is not already in there
sub add_analysis {
  my $logic_name = shift;
  
  my $analysis_id = get_analysis_id($logic_name);
  if (!$analysis_id) {
    my $stmt = qq{
    INSERT INTO
      analysis (
        created,
        logic_name
      )
    VALUES
      (
        NOW(),
        '$logic_name'
      )
    };
    $dbCore->dbc->do($stmt);
    $analysis_id = $dbCore->dbc->db_handle->{'mysql_insertid'};
  }
  
  return $analysis_id;
}

# Adds a description for an analysis if it's not already present
sub add_analysis_description {
  my $analysis_id = shift;
  my $description = shift;
  my $display_label = shift;
  my $displayable = shift;
  my $web_data = shift;
  
  my $stmt = qq{
    INSERT IGNORE INTO
      analysis_description (
	analysis_id,
	description,
	display_label,
	displayable
      )
    VALUES (
      $analysis_id,
      '$description',
      '$display_label',
      $displayable
    )
  };
  $dbCore->dbc->do($stmt);
  
  if (defined($web_data)) {
    $stmt = qq{
      UPDATE TABLE
	analysis_description
      SET
	web_data = "$web_data"
      WHERE
	analysis_id = $analysis_id
    };
    $dbCore->dbc->do($stmt);
  }
}

#ÊAdd gene/transcript/translation annotation to a LRG
sub add_annotation {
  my $lrg = shift;
  my $lrg_name = shift;
  my $lrg_coord_system_name = shift;
  my $biotype = shift;
  my $analysis_logic_name = shift;
  
  # Get the coord_system_id, seq_region_id and analysis_id
  my $cs_id = LRGImport::get_coord_system_id($lrg_coord_system_name);
  my $seq_region_id = LRGImport::get_seq_region_id($lrg_name,$cs_id);
  my $analysis_id = LRGImport::get_analysis_id($analysis_logic_name);
  
  my $gene_id;
  my $sth;
  
  # Get the transcripts
  my $transnodes = $lrg->findNodeArray('fixed_annotation/transcript');

  foreach my $transnode (@{$transnodes}) {

    # Get the transcript start and end
    my $transcript_start = $transnode->data->{'start'};
    my $transcript_end = $transnode->data->{'end'};
    my $transcript_length = $transcript_end-$transcript_start+1;
    my $transcript_name = $transnode->data->{'name'};

    # Get the coding start and end coordinates
    my $cds_node = $transnode->findNode('coding_region');
    my $cds_start = $cds_node->data->{'start'};
    my $cds_end = $cds_node->data->{'end'};
  
    # Insert transcript entry into db
    my $transcript_id = LRGImport::add_transcript(undef,$gene_id,$analysis_id,$seq_region_id,$transcript_start,$transcript_end,1,$biotype,1);

    # Insert transcript_stable_id (we already know what it should be)
    my $transcript_stable_id = $lrg_name . '_' . $transcript_name;
    LRGImport::add_stable_id($transcript_id,$transcript_stable_id,'transcript');
    
    # If we are processing the first transcript of the gene, the gene needs to be added to the database
    if (!defined($gene_id)) {
      $gene_id = LRGImport::add_gene($biotype,$analysis_id,$seq_region_id,$transcript_start,$transcript_end,1,'LRG database',$transcript_id);
      # For genes, use the LRG identifier as gene name
      my $gene_stable_id = $lrg_name;
      # Insert gene stable id into gene_stable_id table
      LRGImport::add_stable_id($gene_id,$gene_stable_id,'gene');
      
      print "Gene:\t" . $gene_id . "\t" . $gene_stable_id . "\n";
    }
  
    print "\tTranscript:\t" . $transcript_id . "\t" . $transcript_stable_id . "\n";
  
    # Update transcript table to insert the gene id
    LRGImport::add_transcript($transcript_id,$gene_id);
      
    my $start_exon_id;
    my $end_exon_id;
    my $exon_id;
    my $exon_count = 0;
    
    my $end_phase = -1;
    my $phase;
    
    # Loop over the transcript nodes and pick out exon-intron pairs and insert the exons into the database
    my @transcript_nodes = @{$transnode->{'nodes'}};
    
    while (my $tr_node = shift(@transcript_nodes)) {
      if ($tr_node->{'name'} eq 'exon') {
	my $exon = $tr_node;
	
	#ÊThe phase of this exon will be the same as the end_phase of the last exon
	$phase = $end_phase;
	
	#ÊPeek at the next node, if it is an intron we extract the phase (end_phase of this exon). Otherwise, we put it back at the beginning of the array
	$tr_node = shift(@transcript_nodes);
	if (defined($tr_node) && $tr_node->{'name'} eq 'intron') {
	    $end_phase = $tr_node->{'data'}{'phase'};
	}
	else {
	    unshift(@transcript_nodes,$tr_node);
	}
	
	$exon_count++;
	my $lrg_coords = $exon->findNode('lrg_coords');
	my $exon_start = $lrg_coords->data->{'start'};
	my $exon_end = $lrg_coords->data->{'end'};
	my $exon_length = ($exon_end-$exon_start+1);
	$exon_id = LRGImport::add_exon($seq_region_id,$exon_start,$exon_end,1);
	
	# Get the next free exon stable id for this transcript
	my $exon_stable_id = LRGImport::get_next_stable_id('exon_stable_id',$transcript_stable_id . '_e');
	# Insert exon stable id into exon_stable_id table
	LRGImport::add_stable_id($exon_id,$exon_stable_id,'exon');
	
	print "\t\tExon:\t" . $exon_id . "\t" . $exon_stable_id;
      
	# If the coding start is within this exon, save this exon id as start_exon_id and calculate the coding start offset within the exon
	if (!defined($start_exon_id) && $exon_start <= $cds_start && $exon_end >= $cds_start) {
	  $start_exon_id = $exon_id;
	  $cds_start = $cds_start - $exon_start + 1;
	  print "\t[First coding]";
	}
	# If the coding end is within this exon, save this exon id as end exon id and calculate end offset within the exon
	if (!defined($end_exon_id) && $exon_start <= $cds_end && $exon_end >= $cds_end) {
	  $end_exon_id = $exon_id;
	  $cds_end = $cds_end - $exon_start + 1;
	  # The end phase will be -1 since translation stops within this exon
	  $end_phase = -1;
	  print "\t[Last coding]";
	}
	
	# Update the exon with the phases
	LRGImport::update_rows([qq{phase = '$phase'},qq{end_phase = '$end_phase'}],[qq{exon_id = '$exon_id'}],['exon']);
	#ÊLink the transcript and the exon
	LRGImport::add_exon_transcript($exon_id,$transcript_id,$exon_count);
	
	print "\t" . $phase . "\t" . $end_phase . "\n";
      }
    }
  
    # Insert translation into db (if available, the gene could be non-protein coding)
    if ($exon_count > 0) {
      my $translation_id = LRGImport::add_translation($transcript_id,$cds_start,$start_exon_id,$cds_end,$end_exon_id);
      
      # Get the next free translation stable_id
      my $translation_stable_id = LRGImport::get_next_stable_id('translation_stable_id',$transcript_stable_id . '_p');
      # Insert translation stable id into translation_stable_id table
      LRGImport::add_stable_id($translation_id,$translation_stable_id,'translation');
      print "\tTranslation:\t" . $translation_id . "\t" . $translation_stable_id . "\n";
      
      # Set the translation_id as the canonical_translation_id in transcript table
      LRGImport::update_rows([qq{canonical_translation_id = '$translation_id'}],[qq{transcript_id = '$transcript_id'}],['transcript']);
    }
  }
  
  #ÊIf necessary, update meta_coord entries for maximum lengths of tables
  foreach my $table (('gene','transcript','exon')) {
    my $max_length = LRGImport::get_max_length($seq_region_id,$table);
    if (defined($max_length)) {
      LRGImport::add_meta_coord($table,$cs_id,$max_length);
    }
  }
}

# Adds an assembly mapping
sub add_assembly_mapping {
    my $asm_sri = shift;
    my $cmp_sri = shift;
    my $asm_start = shift;
    my $asm_end = shift;
    my $cmp_start = shift;
    my $cmp_end = shift;
    my $ori = shift;
    
    my $stmt = qq{
	INSERT IGNORE INTO
	    assembly (
		asm_seq_region_id,
		cmp_seq_region_id,
		asm_start,
		asm_end,
		cmp_start,
		cmp_end,
		ori
	    )
	VALUES (
	    $asm_sri,
	    $cmp_sri,
	    $asm_start,
	    $asm_end,
	    $cmp_start,
	    $cmp_end,
	    $ori
	)
    };
    $dbCore->dbc->do($stmt);
}

#ÊAdd an attrib_type, if the code does not exist
sub add_attrib_type {
    my $code = shift;
    my $name = shift;
    my $description = shift;
    
    my $attrib_type_id = LRGImport::get_attrib_type_id($code);
    if (!defined($attrib_type_id)) {
	$attrib_type_id = LRGImport::insert_row({'code' => $code, 'name' => $name, 'description' => $description},'attrib_type');
    }
    
    return $attrib_type_id;
}

#ÊGets the coord_system_id of an existing coord_system or adds it as a new entry and returns the id
sub add_coord_system {
    my $name = shift;
    my $version = shift;
    my $rank = shift;
    my $species_id = shift;
    my $attrib = shift;
    
    # Set default values
    $species_id ||= 1;
    $attrib ||= 'default_version';
    
    # Check if coord_system already exists
    my $coord_system_id = get_coord_system_id($name,$version,$species_id);
    my $stmt;
    
    # If not add an entry for it
    if (!defined($coord_system_id)) {
	
	# If rank was not specified, get the highest rank in the table and add one
	if (!defined($rank)) {
	    $stmt = qq{
		SELECT
		    MAX(rank)
		FROM
		    coord_system
	    };
	    $rank = ($dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0]+1);
	}
	
	$stmt = qq{
	    INSERT INTO
		coord_system (
		    species_id,
		    name,
		    version,
		    rank,
		    attrib
		)
	    VALUES (
		$species_id,
		'$name',
		?,
		$rank,
		'$attrib'
	    )
	};
	my $sth = $dbCore->dbc->prepare($stmt);
	$sth->execute($version);
	$coord_system_id = $dbCore->dbc->db_handle->{'mysql_insertid'};
    }
    
    return $coord_system_id;
}

# Add DNA
sub add_dna {
    my $seq_region_id = shift;
    my $seq_region_seq = shift;
    
    my $stmt = qq{
	INSERT INTO
	    dna (
		seq_region_id,
		sequence
	    )
	VALUES (
	    $seq_region_id,
	    '$seq_region_seq'
	)
	ON DUPLICATE KEY UPDATE
	    sequence = '$seq_region_seq'
    };
    $dbCore->dbc->do($stmt);
}

# Add an exon
sub add_exon {
    my $seq_region_id = shift;
    my $seq_region_start = shift;
    my $seq_region_end = shift;
    my $seq_region_strand = shift;
    
    my $stmt = qq{
      INSERT INTO
	exon (
	  seq_region_id,
	  seq_region_start,
	  seq_region_end,
	  seq_region_strand
	)
      VALUES
	(
	  $seq_region_id,
	  $seq_region_start,
	  $seq_region_end,
	  $seq_region_strand
	)
    };
    $dbCore->dbc->do($stmt);
    my $exon_id = $dbCore->dbc->db_handle->{'mysql_insertid'};
    
    return $exon_id;
}

# Link an exon to a transcript
sub add_exon_transcript {
    my $exon_id = shift;
    my $transcript_id = shift;
    my $rank = shift;
    
    my $stmt = qq {
	INSERT INTO
	    exon_transcript (
		exon_id,
		transcript_id,
		rank
	    )
	VALUES (
	    $exon_id,
	    $transcript_id,
	    $rank
	)
    };
    $dbCore->dbc->do($stmt);
}

#ÊAdds an external db if it does not exist
sub add_external_db {
    # Required arguments
    my $dbname = shift;
    my $status = shift;
    my $priority = shift;
    
    # Optional arguments
    my %args = (
	'db_display_name' => shift,
	'db_release' => shift,
	'dbprimary_acc_linkable' => shift,
	'display_label_linkable' => shift,
	'type' => shift
    );
    
    my $external_db_ids = get_rows($dbname,'db_name','external_db','external_db_id');
    my $external_db_id;
    
    if (!defined($external_db_ids) || scalar(@{$external_db_ids}) == 0) {
	# external_db_id column is not auto_increment => get maximum id and add one
	my $stmt = qq{
	    SELECT
		MAX(external_db_id) + 1
	    FROM
		external_db
	};
	$external_db_id = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
	
	$stmt = qq{
	    INSERT INTO
		external_db (
		    external_db_id,
		    db_name,
		    status,
		    priority
		)
	    VALUES (
		$external_db_id,
		'$dbname',
		'$status',
		$priority
	    )
	};
	$dbCore->dbc->do($stmt);
	
	foreach my $field (keys %args) {
	    if (defined($args{$field})) {
		$stmt = qq{
		    UPDATE
			external_db
		    SET
			$field = '$args{$field}'
		    WHERE
			external_db_id = $external_db_id
		};
		$dbCore->dbc->do($stmt);
	    }
	}
    }
    else {
	$external_db_id = $external_db_ids->[0][0];
    }
    
    return $external_db_id;
}

# Add a gene
sub add_gene {
    my $biotype = shift;
    my $analysis_id = shift;
    my $seq_region_id = shift;
    my $seq_region_start = shift;
    my $seq_region_end = shift;
    my $seq_region_strand = shift;
    my $source = shift;
    my $canonical_transcript_id = shift;
    my $is_current = shift;
    my $status = shift;
    
    $is_current ||= 1;
    $status ||= 'KNOWN';
    
    my $stmt = qq{
	INSERT INTO
	    gene (
		biotype,
		analysis_id,
		seq_region_id,
		seq_region_start,
		seq_region_end,
		seq_region_strand,
		source,
		is_current,
		canonical_transcript_id,
		status
	    )
	VALUES (
	    '$biotype',
	    $analysis_id,
	    $seq_region_id,
	    $seq_region_start,
	    $seq_region_end,
	    $seq_region_strand,
	    '$source',
	    $is_current,
	    $canonical_transcript_id,
	    '$status'
	)
    };
    $dbCore->dbc->do($stmt);
    my $gene_id = $dbCore->dbc->db_handle->{'mysql_insertid'};
    
    return $gene_id;
}

# Add a gene attribute indicating that the specified Ensembl gene is overlapped by a LRG record
sub add_lrg_overlap {
    my $gene_stable_id = shift;
    my $lrg_name = shift;
    my $partial = shift;
    
    # Get the gene id for the stable id
    my $gene_id = get_object_id_by_stable_id('gene',$gene_stable_id);
    
    # Get the attrib type id to use
    my $attrib_type_id;
    if ($partial) {
	$attrib_type_id = add_attrib_type('GeneOverlapLRG','Gene overlaps LRG','This gene is partially overlapped by a LRG region (start or end outside LRG)');
    }
    else {
	$attrib_type_id = add_attrib_type('GeneInLRG','Gene in LRG','This gene is contained within an LRG region');
    }
    
    # As value, use the LRG name
    add_object_attrib($gene_id,$attrib_type_id,$lrg_name,'gene');
}

# Add all the necessary data to the db in order to project and transfer between the LRG and chromosome
sub add_mapping {
    my $lrg_name = shift;
    my $lrg_coord_system = shift;
    my $q_seq_length = shift;
    my $mapping = shift;
    my $assembly = shift;
    
    $assembly ||= get_assembly();
    
    my $pairs = $mapping->{'pairs'};
    my %additions;
  
    # Get the coord system id for this LRG or add a new entry if not present
    my $cs_id = add_coord_system($lrg_coord_system);

    my @meta_id;
    # Insert a general meta entry for 'LRG' if none is present (Eugene mentioned this, used to determine if a species has LRGS (?))
    push(@meta_id,(add_meta_key_value($lrg_coord_system,$lrg_coord_system)));

    # In order to project the slice, we need to add an entry to the meta table (if not present)
    my $meta_key = 'assembly.mapping';

    # Need the mapping to contig and also to chromosome via contig
    my $meta_value = $lrg_coord_system . '#contig';
    foreach my $append (('','#chromosome:' . $assembly)) {
	my $id = add_meta_key_value($meta_key,$meta_value . $append);
	push(@meta_id,$id);
    }
    
    # Add a seq_region for the LRG if it doesn't exist
    my @seq_region_ids;
    my $q_seq_region_id = add_seq_region($lrg_name,$cs_id,$q_seq_length);
    push(@seq_region_ids,$q_seq_region_id);
  
    #ÊAdd a seq_region_attrib indicating that it is toplevel
    my $attrib_type_id = get_attrib_type_id('toplevel') or warn("Could not get attrib_type_id for toplevel attribute!\n");
    LRGImport::insert_row({'seq_region_id' => $q_seq_region_id, 'attrib_type_id' => $attrib_type_id, 'value' => 1},'seq_region_attrib',1);
    
    #ÊAdd a seq_region_attrib indicating that it is non-reference
    $attrib_type_id = get_attrib_type_id('non_ref') or warn("Could not get attrib_type_id for non_ref attribute!\n");
    LRGImport::insert_row({'seq_region_id' => $q_seq_region_id, 'attrib_type_id' => $attrib_type_id, 'value' => 1},'seq_region_attrib',1);
    
    #ÊAdd a seq_region_attrib indicating that it is a LRG
    $attrib_type_id = get_attrib_type_id('LRG') or warn("Could not get attrib_type_id for LRG attribute!\n");
    LRGImport::insert_row({'seq_region_id' => $q_seq_region_id, 'attrib_type_id' => $attrib_type_id, 'value' => 1},'seq_region_attrib',1);
    
    # Get the seq_region_id for the target contig
    my $ctg_cs_id = get_coord_system_id('contig');
    
    $attrib_type_id = undef;
    my $contig_csi;
    my $sa = $dbCore->get_SliceAdaptor();
    
    #ÊLoop over the pairs array. For each mapping span, first get a chromosomal slice and project this one to contigs
    foreach  my $pair (sort {$a->[3]<=>$b->[3]} @{$pairs}) {
    
	# Each span is represented by a 'DNA' type
	if ($pair->[0] eq 'DNA') {
	    die("Distance between query and target is not the same, there is an indel which shouldn't be there. Check the code!") unless ($pair->[2]-$pair->[1] == $pair->[4]-$pair->[3]);
	    
	    # Get a chromosomal slice. Always in the forward orientation
	    my $chr_slice = $sa->fetch_by_region('chromosome',$mapping->{'chr_name'},$pair->[3],$pair->[4],1,$assembly);
	    
	    #ÊProject the chromosomal slice to contig
	    my $segments = $chr_slice->project('contig');
	    
	    # Die if the projection failed
	    die("Could not project " . $chr_slice->name() . " to contigs!") if (!defined($segments) || scalar(@{$segments}) == 0);
	    
	    # Loop over the projection segments and insert the corresponding mapping between the LRG and contig
	    foreach my $segment (@{$segments}) {
		
		#ÊThe projected slice on contig
		my $ctg_slice = $segment->[2];
		my $ctg_seq_region_id = $sa->get_seq_region_id($ctg_slice);
		my $ctg_start = $ctg_slice->start();
		my $ctg_end = $ctg_slice->end();
		# The orientation of the contig relative to the LRG is the orientation relative to chromosome multiplied by the orientation of chromosome releative to LRG
		my $ctg_strand = ($ctg_slice->strand() * $pair->[5]);
		
		# The offset of the chromosome mapping
		my $chr_offset = $segment->[0] - 1;
		my $chr_length = $segment->[1] - $segment->[0] + 1;
		
		my $lrg_start;
		my $lrg_end;
		
		# The LRG->contig mapping is straightforward if LRG and chromosome is the same orientation
		if ($pair->[5] > 0) {
		    $lrg_start = $pair->[1] + $chr_offset;
		    $lrg_end = $pair->[1] + $chr_offset + $chr_length - 1;
		}
		#ÊIf the LRG and chromosome are different orientation, the hits are backwards
		else {
		    $lrg_start = $pair->[2] - ($chr_offset + $chr_length) + 1;
		    $lrg_end = $pair->[2] - $chr_offset;    
		}
		add_assembly_mapping(
		    $q_seq_region_id,
		    $ctg_seq_region_id,
		    $lrg_start,
		    $lrg_end,
		    $ctg_start,
		    $ctg_end,
		    $ctg_strand
		);
	    }
	}
	# Else if this is a mismatch, put an rna_edit entry into the seq_region_attrib table
	elsif ($pair->[0] eq 'M') {
      
	    # Get the attrib_type_id for rna_edit
	    if (!defined($attrib_type_id)) {
		$attrib_type_id = get_attrib_type_id('_rna_edit');
	    }
	    my $value = $pair->[1] . ' ' . $pair->[2] . ' ' . $pair->[6];
	    add_object_attrib($q_seq_region_id,$attrib_type_id,$value,'seq_region');
	}
	# Else if this is a gap,
	elsif ($pair->[0] eq 'G') {
	    # If this is a deletion in the LRG, we just make a break in the assembly table. This means that we don't need to do anything at this point, a new 'DNA' span will be in the pairs array
	    # If this is an insertion in the LRG, we add the sequence as a contig and put a mapping between the LRG and the new contig
	    if ($pair->[2] >= $pair->[1]) {
		my $contig_seq = $pair->[6];
		my $contig_len = ($pair->[2] - $pair->[1] + 1);
		my $contig_name = $lrg_name . '_ins_' . $pair->[1] . '-' . $pair->[2];
		
		# Get the coord_system_id for contig
		if (!defined($contig_csi)) {
		    $contig_csi = get_coord_system_id('contig');
		}
	
		# Get or add a seq_region for this contig
		my $contig_sri = add_seq_region($contig_name,$contig_csi,$contig_len);
	    
		# Add DNA for this seq_region unless the database is not core. Ugly regexp to determine this
		add_dna($contig_sri,$contig_seq) if ($dbCore->dbc()->dbname() =~ m/\_core\_\d+\_\d+\s/);
		
		# Add a mapping between the LRG and the contig
		add_assembly_mapping(
		    $q_seq_region_id,
		    $contig_sri,
		    $pair->[1],
		    $pair->[2],
		    1,
		    $contig_len,
		    1
		);
	    }
	}	
    }
}

# Add a meta_coord entry
sub add_meta_coord {
    my $table = shift;
    my $cs_id = shift;
    my $max_length = shift;
    
    my $stmt = qq{
	INSERT INTO
	    meta_coord (
	      table_name,
	      coord_system_id,
	      max_length
	    )
	VALUES (
	    '$table',
	    $cs_id,
	    $max_length
	)
	ON DUPLICATE KEY UPDATE
	    max_length = $max_length
    };
    $dbCore->dbc->do($stmt);
}

# Add a meta key/value pair
sub add_meta_key_value {
    my $key = shift;
    my $value = shift;
    my $species_id = shift;
    
    #ÊSet default value
    $species_id ||= 1;
    
    # Check if a value already exists for this key/value pair
    my $stmt = qq{
	SELECT
	    meta_id
	FROM
	    meta
	WHERE
	    species_id = $species_id AND
	    meta_key = '$key' AND
	    meta_value = '$value'
	LIMIT 1
    };
    my $meta_id = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
    
    # If none exists, add it
    if (!defined($meta_id)) {
	$stmt = qq{
	    INSERT INTO
		meta (
		    species_id,
		    meta_key,
		    meta_value
		)
	    VALUES (
		$species_id,
		'$key',
		'$value'
	    )
		    
	};
	$dbCore->dbc->do($stmt);
	$meta_id = $dbCore->dbc->db_handle->{'mysql_insertid'};
    }
    
    return $meta_id;
}


#ÊAdd a seq_region or get the seq_region_id of an already existing one
sub add_seq_region {
    my $name = shift;
    my $coord_system_id = shift;
    my $length = shift;
    
    # Check if the seq_region already exists
    my $seq_region_id = get_seq_region_id($name,$coord_system_id);
    my $stmt;
    
    # If not, add it
    if (!defined($seq_region_id)) {
	$stmt = qq{
	    INSERT INTO
		seq_region (
		    name,
		    coord_system_id,
		    length
		)
	    VALUES (
		'$name',
		$coord_system_id,
		$length
	    )
	};
	$dbCore->dbc->do($stmt);
	$seq_region_id = $dbCore->dbc->db_handle->{'mysql_insertid'};
    }
    
    return $seq_region_id;
}

# Add a *_attrib to an object
sub add_object_attrib {
    my $object_id = shift;
    my $attrib_type_id = shift;
    my $value = shift;
    my $object_type = shift;
    
    my $table = $object_type . '_attrib';
    my $key = $object_type . '_id';
    
    my $stmt = qq{
	INSERT IGNORE INTO
	    $table (
		$key,
		attrib_type_id,
		value
	    )
	VALUES (
	    $object_id,
	    $attrib_type_id,
	    '$value'
	)
    };
    $dbCore->dbc->do($stmt);    
}

# Add stable id for a gene/transcript/translation. If one already exists, replace it
sub add_stable_id {
    my $object_id = shift;
    my $object_stable_id = shift;
    my $object_type = shift;
    my $version = shift;
    
    my $table = $object_type . "_stable_id";
    my $key = $object_type . "_id";
    
    $version ||= 1;
    
    # Insert object stable id into object_stable_id table
    my $stmt = qq{
      INSERT INTO
	$table (
	  $key,
	  stable_id,
	  version,
	  created_date,
	  modified_date
	)
      VALUES (
	$object_id,
	'$object_stable_id',
	$version,
	NOW(),
	NOW()
      )
      ON DUPLICATE KEY UPDATE
	stable_id = '$object_stable_id',
	version = $version,
	modified_date = NOW()
    };
    $dbCore->dbc->do($stmt);
}

# Add a transcript or update one with gene_id if a transcript_id is specified
sub add_transcript {
    my $transcript_id = shift;
    my $gene_id = shift;
    my $analysis_id = shift;
    my $seq_region_id = shift;
    my $seq_region_start = shift;
    my $seq_region_end = shift;
    my $seq_region_strand = shift;
    my $biotype = shift;
    my $is_current = shift;
    my $status = shift;
    
    $is_current ||= 1;
    $gene_id ||= 'NULL';
    $status ||= 'KNOWN';
    
    my $stmt;
    
    # If a transcript_id was specified, update that transcript with the gene_id
    if (defined($transcript_id)) {
	$stmt = qq{
	    UPDATE
		transcript
	    SET
		gene_id = $gene_id
	    WHERE
		transcript_id = $transcript_id
	};
	$dbCore->dbc->do($stmt);
    }
    # Else, insert transcript entry into db
    else {
	# Will this properly set gene_id to NULL if undef provided? Will not really matter though since it is updated later
	$stmt = qq{
	    INSERT INTO
	      transcript (
		gene_id,
		analysis_id,
		seq_region_id,
		seq_region_start,
		seq_region_end,
		seq_region_strand,
		biotype,
		is_current,
		status
	      )
	    VALUES (
	      $gene_id,
	      $analysis_id,
	      $seq_region_id,
	      $seq_region_start,
	      $seq_region_end,
	      $seq_region_strand,
	      '$biotype',
	      $is_current,
	      '$status'
	    )
	};
	$dbCore->dbc->do($stmt);
	$transcript_id = $dbCore->dbc->db_handle->{'mysql_insertid'};
    }
    
    return $transcript_id;
}

# Add a translation
sub add_translation {
    my $transcript_id = shift;
    my $cds_start = shift;
    my $start_exon_id = shift;
    my $cds_end = shift;
    my $end_exon_id = shift;
    
    my $stmt;
    my $translation_id = get_translation_id($transcript_id);
    # If no translation is stored for this transcript, add one
    if (!defined($translation_id)) {
	$stmt = qq{
	    INSERT INTO
		translation (
		    transcript_id,
		    seq_start,
		    start_exon_id,
		    seq_end,
		    end_exon_id
		)
	    VALUES (
		$transcript_id,
		$cds_start,
		$start_exon_id,
		$cds_end,
		$end_exon_id
	    ) 
	};
	$dbCore->dbc->do($stmt);
	$translation_id = $dbCore->dbc->db_handle->{'mysql_insertid'};
    }
    # Else, update the data for the one that is stored
    else {
	$stmt = qq{
	    UPDATE
		translation
	    SET
		seq_start = $cds_start,
		start_exon_id = $start_exon_id,
		seq_end = $cds_end,
		end_exon_id = $end_exon_id
	    WHERE
		translation_id = $translation_id
	};
	$dbCore->dbc->do($stmt);
    }
    
    return $translation_id;  
}

sub add_xref {
    my $external_db_name = shift;
    my $db_primary_acc = shift;
    my $display_label = shift;
    
    # Optional arguments
    my %args = (
	'version' => shift,
	'description' => shift,
	'info_type' => shift,
	'info_text' => shift
    );
    
    # Get external db id for database
    my $external_db_ids = get_rows($external_db_name,'db_name','external_db','external_db_id');
    if (!defined($external_db_ids) || scalar(@{$external_db_ids}) == 0) {
	warn("No external db entry found for database $external_db_name\n");
	return undef;
    }
    my $external_db_id = $external_db_ids->[0][0];
    
    # Check if xref entry already exists
    my $stmt = qq{
	SELECT
	    xref_id
	FROM
	    xref
	WHERE
	    external_db_id = $external_db_id AND
	    dbprimary_acc = '$db_primary_acc' AND
	    display_label = '$display_label'
	LIMIT 1
    };
    my $xref_id = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
    
    # Otherwise add an entry
    if (!$xref_id) {
	$stmt = qq{
	    INSERT INTO
		xref (
		    external_db_id,
		    dbprimary_acc,
		    display_label
		)
	    VALUES (
		$external_db_id,
		'$db_primary_acc',
		'$display_label'
	    )
	};
	$dbCore->dbc->do($stmt);
	$xref_id = $dbCore->dbc->db_handle->{'mysql_insertid'};
	
	foreach my $field (keys %args) {
	    if (defined($args{$field})) {
		$stmt = qq{
		    UPDATE
			xref
		    SET
			$field = '$args{$field}'
		    WHERE
			xref_id = $xref_id
		};
		$dbCore->dbc->do($stmt);
	    }
	}
    }
    
    return $xref_id;
}

sub add_object_xref {
    my $ensembl_id = shift;
    my $ensembl_object_type = shift;
    my $xref_id = shift;
    
    # Optional arguments
    my %args = (
	'linkage_annotation' => shift,
	'analysis_id' => shift
    );
    
    my $stmt = qq{
	INSERT IGNORE INTO
	    object_xref (
		ensembl_id,
		ensembl_object_type,
		xref_id
	    )
	VALUES (
	    $ensembl_id,
	    '$ensembl_object_type',
	    $xref_id
	)
    };
    
    $dbCore->dbc->do($stmt);
    my $object_xref_id = $dbCore->dbc->db_handle->{'mysql_insertid'};
    
    return $object_xref_id;
}

# Get a listing of available remote files from a supplied url
sub fetch_remote_lrg_ids {
  my $urls = shift;
  
  #ÊA hash holding the results
  my %result;
  
  # Create a user agent object
  my $ua = LWP::UserAgent->new;
  
  # Loop over the supplied URLs
  while (my $url = shift(@{$urls})) {
    # Create a request
    my $req = HTTP::Request->new(GET => $url);
    # Pass request to the user agent and get a response back
    my $res = $ua->request($req);
  
    # Store data on the request in the hash
    $result{$url} = {
      'success' => $res->is_success,
      'message' => $res->message
    };
  
    #ÊGet the listing if the request was successful
    if ($result{$url}{'success'}) {
      my @lrgs = ($res->content =~ m/(LRG\_\d+)\.xml/g);
      $result{$url}{'lrg_id'} = \@lrgs;
    }
  }
  
  return \%result;
}

#ÊAttempt to fetch the LRG XML file from the LRG server
sub fetch_remote_lrg {
  my $lrg_id = shift;
  my $urls = shift;
  my $destfile = shift;
  
  # This routine returns a hash with details of the request
  my %result;
  
  #ÊIf no destination file was specified, create one in the tmp directory
  $destfile ||= '/tmp/' . $lrg_id . '_' . time() . '.xml';
  #ÊRemove the destfile if it exists
  unlink($destfile);
  
  # Set the output file
  $result{'xmlfile'} = $destfile;
  
  # Create a user agent object
  my $ua = LWP::UserAgent->new;

  # Cycle through the supplied URLs to try until one is succsessful
  my $res;
  while (my $url = shift(@{$urls})) {
    # Create a request
    my $req = HTTP::Request->new(GET => "$url/$lrg_id\.xml");

    # Pass request to the user agent and get a response back
    $res = $ua->request($req);
    
    $result{$url} = {
	'message' => $res->message
    };
    $result{'success'} = $res->is_success;
    
    # Break if the request was successful
    last if ($result{'success'});
  }

  # If nothing could be found, return the result hash
  return \%result unless ($result{'success'});
  
  # Write the XML to a temporary file
  open(XML_OUT,'>',$destfile) or die("Could not open $destfile for writing");
  print XML_OUT $res->content;
  close(XML_OUT);
  
  # Return the result hash
  return \%result;
}

# Get the analysis_id for an analysis based on the logic_name
sub get_analysis_id {
  my $logic_name = shift;
  
  my $stmt = qq{
    SELECT
      analysis_id
    FROM
      analysis
    WHERE
      logic_name = '$logic_name'
    LIMIT 1
  };
  my $analysis_id = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
  
  return $analysis_id;
}

# Get the current assembly name
sub get_assembly {
    
    my $stmt = qq{
      SELECT
	meta_value
      FROM
	meta
      WHERE
	meta_key = 'assembly.name'
    };
    my $assembly = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
    
    return $assembly;
}

# Get the attrib_type_id for a name
sub get_attrib_type_id {
    my $code = shift;
    
    my $stmt = qq{
	SELECT
	    attrib_type_id
	FROM
	    attrib_type
	WHERE
	    code = '$code'
	LIMIT 1
    };
    my $attrib_type_id = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
    
    return $attrib_type_id;
}

# Get the coord_system_id for a species_id, name and version
sub get_coord_system_id {
    my $name = shift;
    my $version = shift;
    my $species_id = shift;
    
    # Set default value
    $species_id ||= 1;
    
    # Allow for version to be NULL
    my $condition;
    if (!defined($version) || $version eq 'NULL') {
	$condition = qq{version IS NULL AND };
    }
    else {
	$condition .= qq{version = '$version' AND };
    }
    $condition .= qq{
	species_id = $species_id AND
	name = '$name'
	LIMIT 1
    };
    
    my $stmt = qq{
        SELECT
	    coord_system_id
        FROM
	    coord_system
        WHERE
    };
    $stmt .= $condition;
    my $coord_system_id = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
    
    return $coord_system_id;
}

#ÊGet the maximum field value from specified tables and fields
sub get_max_key {
    my $table_fields = shift;
    
    # By default, list the tables that are affected by LRGs
    $table_fields ||= {
	'analysis' => 'analysis_id',
	'analysis_description' => 'analysis_id',
	'assembly' => 'asm_seq_region_id',
	'attrib_type' => 'attrib_type_id',
	'coord_system' => 'coord_system_id',
	'dna' => 'seq_region_id',
	'exon' => 'exon_id',
	'exon_stable_id' => 'exon_id',
	'exon_transcript' => 'exon_id',
	'external_db' => 'external_db_id',
	'gene' => 'gene_id',
	'gene_attrib' => 'gene_id',
	'gene_stable_id' => 'gene_id',
	'meta' => 'meta_id',
	'meta_coord' => 'coord_system_id',
	'object_xref' => 'object_xref_id',
	'seq_region' => 'seq_region_id',
	'seq_region_attrib' => 'seq_region_id',
	'transcript' => 'transcript_id',
	'transcript_stable_id' => 'transcript_id',
	'translation' => 'translation_id',
	'translation_stable_id' => 'translation_id',
	'xref' => 'xref_id'
    };
    
    my %max_values;
    my $stmt;
    while (my ($table,$field) = each(%{$table_fields})) {
	$stmt = qq{
	    SELECT
		MAX($field)
	    FROM
		$table
	};
	my $val = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
	$max_values{$table . '.' . $field} = $val;
    }
    
    return \%max_values;
}

# Get maximum length of a gene/transcript/translation on a seq_region
sub get_max_length {
    my $seq_region_id = shift;
    my $table = shift;
    
    my $stmt = qq{
      SELECT
        MAX(1 + t.seq_region_end - t.seq_region_start)
      FROM
	$table t
      WHERE
	t.seq_region_id = $seq_region_id
    };
    my $max_length = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
    
    return $max_length;    
}

# Get the next free stable id for a *_stable_id table.
sub get_next_stable_id {
  my $table = shift;
  my $prefix = shift;
  
  my $stable_id;
  my $suflength = 9;
  my $prelength = length($prefix);
  
  my $stmt = qq{
    SELECT
      MAX(CONVERT(SUBSTR(stable_id,$prelength+1),UNSIGNED INTEGER))
    FROM
      $table
    WHERE
      stable_id LIKE '$prefix%'
  };
  my $max = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
  $max ||= 0;

  $max++;
  $stable_id = $prefix . sprintf("%u",$max);
  
  return $stable_id;
};

sub get_object_ids_by_seq_region_id {
    my $object_type = shift;
    my $seq_region_id = shift;
    
    my $key = $object_type . '_id';
    my $table = $object_type;
    my @ids;
    
    my $stmt = qq{
	SELECT
	    $key
	FROM
	    $table
	WHERE
	    seq_region_id = $seq_region_id
    };
    my $sth = $dbCore->dbc->prepare($stmt);
    $sth->execute();
    
    while (my $row = $sth->fetchrow_arrayref()) {
	push(@ids,$row->[0]);
    }
    
    return \@ids;
}

# Get an object id by its stable_id
sub get_object_id_by_stable_id {
    my $object_type = shift;
    my $stable_id = shift;
    
    my $key = 'stable_id';
    my $table = $object_type . '_stable_id';
    my $id_field = $object_type . '_id';
    
    my $rows = get_rows($stable_id,$key,$table,$id_field);
    
    my $object_id;
    $object_id = $rows->[0][0] if (defined($rows) && scalar(@{$rows}) > 0);
    
    return $object_id;
}

# Get the seq_region_id for a seq_region name and coord_system_id
sub get_seq_region_id {
    my $name = shift;
    my $coord_system_id = shift;
    
    my $stmt = qq{
	SELECT
	    seq_region_id
	FROM
	    seq_region
	WHERE
	    name = '$name' AND
	    coord_system_id = $coord_system_id
	LIMIT 1
    };
    my $seq_region_id = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
    
    return $seq_region_id;
}
#ÊGet seq_region_ids using approximate name matching
sub get_seq_region_ids {
    my $name = shift;
    my $coord_system_id = shift;
    
    my $stmt = qq{
	SELECT
	    seq_region_id
	FROM
	    seq_region
	WHERE
	    name LIKE '\%$name\%' AND
	    coord_system_id = $coord_system_id
    };
    my @seq_region_ids;
    my $sth = $dbCore->dbc->prepare($stmt);
    $sth->execute();
    while (my $row = $sth->fetchrow_arrayref()) {
	push(@seq_region_ids,$row->[0]);
    }
    
    return \@seq_region_ids;
}

sub fetch_rows {
  my $fields = shift;
  my $tables = shift;
  my $conditions = shift;
  my $orders = shift;
  
  $orders ||= ["NULL"];
  
  my $field_string = join(', ',@{$fields});
  my $table_string = join(', ',@{$tables});
  # This will throw a MySQL error if no condition was given, which is good to prevent something nasty happening to the db
  my $condition_string = "(" . join(') AND (',@{$conditions}) . ")";
  my $order_string = join(", ",@{$orders});
  
  my $stmt = qq{
    SELECT
      $field_string
    FROM
      $table_string
    WHERE
      $condition_string
    ORDER BY
      $order_string
  };
  my $result = $dbCore->dbc->db_handle->selectall_arrayref($stmt);
  
  return $result;
}

sub get_rows {
    my $id = shift;
    my $key = shift;
    my $table = shift;
    my $fields = \@_;
    
    return fetch_rows($fields,[$table],["$key = '$id'"]);
}

# Get stable id for an object type and object id
sub get_stable_id {
    my $object_id = shift;
    my $object_type = shift;
    
    my $table = $object_type . "_stable_id";
    my $key = $object_type . "_id";
    
    my $stmt = qq{
	SELECT
	    stable_id
	FROM
	    $table
	WHERE
	    $key = $object_id
	LIMIT 1
    };
    my $stable_id = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
  
    return $stable_id;
}

# Get translation id for a transcript
sub get_translation_id {
    my $transcript_id = shift;
    
    my $stmt = qq{
	SELECT
	    translation_id
	FROM
	    translation
	WHERE
	    transcript_id = $transcript_id
	LIMIT 1
    };
    my $translation_id = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
    
    return $translation_id;
}

# Remove all traces of a LRG from database. If last argument is defined, will even remove coord_system, analysis and meta entries unless
# there still exists LRGs in the database
sub purge_db {
    my $lrg_name = shift;
    my $lrg_coord_system = shift;
    my $purge_all = shift;
    
    # Get the coord system id
    my $cs_id = get_coord_system_id($lrg_coord_system);
    if (!defined($cs_id)) {
      warn("Could not find coordinate system $lrg_coord_system in " . $dbCore->dbc()->dbname());
      $cs_id = -1;
    }
    my $ctg_coord_system_id = get_coord_system_id('contig');
    if (!defined($ctg_coord_system_id)) {
      warn("Could not find coordinate system contig in " . $dbCore->dbc()->dbname());
      $ctg_coord_system_id = -1;
    }
    
    # Get seq region ids on contigs for any seq regions associated with this LRG
    my @seq_region_ids = @{get_seq_region_ids($lrg_name . '_',$ctg_coord_system_id)};
    #ÊLastly, add the seq region id for this LRG itself
    push(@seq_region_ids,get_seq_region_id($lrg_name,$cs_id));
    
    # Delete entries from relevant tables
    foreach my $seq_region_id (@seq_region_ids) {
	next if (!defined($seq_region_id));
	# Delete from gene/transcript/exon/translation 
	foreach my $object_type (('gene','exon','transcript')) {
	    # Get the object_ids associated with this seq_region_id
	    my $object_ids = get_rows($seq_region_id,'seq_region_id',$object_type,$object_type . '_id');
	    while (my $object_id = shift(@{$object_ids})) {
	        my $oid = $object_id->[0];
		# Translations have to be removed by the corresponding transcript ids
		if ($object_type eq 'transcript') {
		    my $stmt = qq{
		    DELETE FROM
			t,
			tsi
		    USING
			translation t,
			translation_stable_id tsi
		    WHERE
			tsi.translation_id = t.translation_id AND
			t.transcript_id = $oid
		    };
		    $dbCore->dbc->do($stmt);
		    
		    my $translation_id = LRGImport::get_translation_id($oid);
		    # Remove from object_xref
		    remove_row([qq{ensembl_object_type = 'Translation'}, qq{ensembl_id = $translation_id}],'object_xref') if defined($translation_id);
		    #ÊRemove entries in the exon_transcript table
		    remove_row([qq{transcript_id = '$oid'}],'exon_transcript');
		}
		
		# Remove rows with this object id (also from associated stable_id table)
		remove_row([$object_type . '_id' . qq{ = $oid}],$object_type);
		remove_row([$object_type . '_id' . qq{ = $oid}],$object_type . '_stable_id');
		
		#ÊRemove xref objects linked to this object id. Does not remove the xref itself, just the row in object_xref
		remove_row([qq{ensembl_object_type = '} . ucfirst($object_type) . qq{'}, qq{ensembl_id = $oid}],'object_xref');
	    }
	}
	# Delete from assembly
	remove_row([qq{asm_seq_region_id = $seq_region_id}],'assembly');
	remove_row([qq{cmp_seq_region_id = $seq_region_id}],'assembly');
	# Delete from seq_region_attrib
	remove_row([qq{seq_region_id = $seq_region_id}],'seq_region_attrib');
	# Delete from dna
	remove_row([qq{seq_region_id = $seq_region_id}],'dna');
	# Delete from seq_region
	remove_row([qq{seq_region_id = $seq_region_id}],'seq_region');
    }
    
    #ÊGet any xrefs for this LRG
    my $xref_ids = fetch_rows([qq{x.xref_id}],[qq{xref x},qq{external_db edb}],[qq{x.external_db_id = edb.external_db_id},qq{edb.db_name = 'LRG' OR edb.db_name = 'ENS_LRG_transcript' OR edb.db_name = 'ENS_LRG_gene'},qq{x.display_label = '$lrg_name' OR x.display_label LIKE '$lrg_name\_t%' OR x.display_label LIKE '$lrg_name\_g%'}]);
    while (my $row = shift(@{$xref_ids})) {
      my $xref_id = $row->[0];
      # Remove from object_xref
      remove_row([qq{xref_id = '$xref_id'}],'object_xref');
      # Remove from xref
      remove_row([qq{xref_id = '$xref_id'}],'xref');
    }
    
    # Remove any gene attributes linked to this LRG
    remove_row([qq{value = '$lrg_name'}],'gene_attrib');
    
    # If specified, remove EVERYTHING LRG-related unless there still are LRG entries left
    if ($purge_all) {
      my $seq_regions = get_rows($cs_id,'coord_system_id','seq_region','seq_region_id');
      if (scalar(@{$seq_regions})) {
	warn("Seq regions belonging to LRG coordinate system still exist in database, will not clear LRG data!");
	return;
      }
      
      # Remove coord_system entry
      remove_row([qq{coord_system_id = $cs_id}],'coord_system');
      # Remove meta_coord entries
      remove_row([qq{coord_system_id = $cs_id}],'meta_coord');
      #ÊRemove meta table entries
      my $ids = fetch_rows(["meta_id"],["meta"],["meta_key = 'assembly.mapping' OR meta_key = 'lrg'","meta_value LIKE '%lrg%'"]);
      remove_row(["meta_id = '" . join("' OR meta_id = '",map($_->[0],@{$ids})) . "'"],"meta");
      # Remove analysis entries
      $ids = fetch_rows(["analysis_id"],["analysis"],["logic_name LIKE 'LRG_%'"]);
      remove_row(["analysis_id = '" . join("' OR analysis_id = '",map($_->[0],@{$ids})) . "'"],"analysis_description");
      remove_row(["analysis_id = '" . join("' OR analysis_id = '",map($_->[0],@{$ids})) . "'"],"analysis");
    }
}

# Insert row into table
sub insert_row {
    my $fields = shift;
    my $table = shift;
    my $ignore = shift;
    
    if ($ignore) {
	$ignore = 'IGNORE';
    }
    else {
	$ignore = '';
    }
    
    # keys() and values() are guaranteed to return in the same order so the code below is fine
    my $field_list = join(',',keys(%{$fields}));
    my $value_list = qq{'} . join(qq{','},values(%{$fields})) . qq{'};
    
    my $stmt = qq{
	INSERT $ignore INTO
	    $table (
		$field_list
	    )
	VALUES (
	    $value_list
	)
    };
    
    $dbCore->dbc->do($stmt);
    my $id = $dbCore->dbc->db_handle->{'mysql_insertid'};
    
    return $id;
}

# Remove rows from tables using specified conditions hash (joined by AND)
sub remove_row {
    my $conditions = shift;
    my $table = shift;
    
    # This will throw a MySQL error if no condition was given, which is good to prevent something nasty happening to the db
    my $condition_string = "(" . join(') AND (',@{$conditions}) . ")";
    
    my $stmt = qq{
	DELETE FROM
	    $table
	WHERE
	    $condition_string    
    };
    $dbCore->dbc->do($stmt);
}

#ÊUpdate a table where the specified conditions are fulfilled
sub update_rows {
    my $values = shift;
    my $conditions = shift;
    my $tables = shift;
    
    my $table = join(',',@{$tables});
    # This will throw a MySQL error if no condition was given, which is good to prevent something nasty happening to the db
    my $condition_string = "(" . join(') AND (',@{$conditions}) . ")";
    my $value_string = join(', ',@{$values});
    
    my $stmt = qq{
	UPDATE
	    $table
	SET
	    $value_string
	WHERE
	    $condition_string
    };
    $dbCore->dbc->do($stmt);
}
1;
