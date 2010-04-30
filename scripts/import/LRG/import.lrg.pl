#!perl

use strict;

use Getopt::Long;
use LRG;
use Bio::Seq;
use Bio::EnsEMBL::Registry;

# set default options
my $registry_file = "ensembl.registry";
my ($help);
my $clean;

usage() unless scalar @ARGV;

# get options from command line
GetOptions(
	   'registry_file=s' => \$registry_file,
	   'help' => \$help,
	   'clean=s' => \$clean,
);

usage() if $help;

# connect to core database
Bio::EnsEMBL::Registry->load_all( $registry_file );
#Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
our $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor('human','core');
die "ERROR: Could not connect to core database\n" unless defined $dbCore;

# Clean up data in the database if required
if ($clean) {
  my $stmt = qq{
    DELETE FROM
      meta
    WHERE
      meta_key = 'assembly.mapping' AND
      meta_value LIKE '$clean\#\%'
  };
 $dbCore->dbc->do($stmt);
  
  $stmt = qq{
    SELECT
      coord_system_id
    FROM
      coord_system
    WHERE
      name = '$clean'
  };
  my $result = $dbCore->dbc->db_handle->selectall_arrayref($stmt);
  foreach my $ref (@{$result}) {
    my $cs_id = $ref->[0];
    $stmt = qq{
      DELETE FROM
	coord_system
      WHERE
	coord_system_id = '$cs_id'
    };
    $dbCore->dbc->do($stmt);
    
    $stmt = qq{
      SELECT
	seq_region_id
      FROM
	seq_region
      WHERE
	coord_system_id = '$cs_id'
    };
    my $seq_region_id = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
    
    $stmt = qq{
      DELETE FROM
	seq_region
      WHERE
	seq_region_id = '$seq_region_id'
    };
    $dbCore->dbc->do($stmt);
    
    $stmt = qq{
      DELETE FROM
	assembly
      WHERE
	asm_seq_region_id = '$seq_region_id'
    };
    $dbCore->dbc->do($stmt);
    
    my $transcript_ids;
    
    foreach my $table (('exon','gene','transcript')) {
      my $key = $table . '_id';
      my $id_table = $table . '_stable_id';
      
      $stmt = qq{
	SELECT
	  $key
	FROM
	  $table
	WHERE
	  seq_region_id = '$seq_region_id'
      };
      my $ids = $dbCore->dbc->db_handle->selectall_arrayref($stmt);
      
      if ($table eq 'transcript') {
	$transcript_ids = \@{$ids};
      }
      
      foreach my $id (@{$ids}) {
	$stmt = qq{
	  DELETE FROM
	    $id_table
	  WHERE
	    $key = '$id->[0]'
	};
    	$dbCore->dbc->do($stmt);
	
	$stmt = qq{
	  DELETE FROM
	    $table
	  WHERE
	    $key = '$id->[0]'
	};
    	$dbCore->dbc->do($stmt);
      }
    }
    
    foreach my $id (@{$transcript_ids}) {
      $stmt = qq{
	SELECT
	  translation_id
	FROM
	  translation
	WHERE
	  transcript_id = '$id->[0]'
	LIMIT 1
      };
      my $translation_id = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
      
      $stmt = qq{
	DELETE FROM
	  translation_stable_id
	WHERE
	  translation_id = '$translation_id'
      };
      $dbCore->dbc->do($stmt);
      $stmt = qq{
	DELETE FROM
	  translation
	WHERE
	  translation_id = '$translation_id'
      };
      $dbCore->dbc->do($stmt);      
    }
    
  }
  exit(0);
}

# check for an input file
my $input_file = shift @ARGV;
die "ERROR: No input file specified\n" unless -e $input_file;

# create an LRG object from it
my $lrg = LRG::LRG::newFromFile($input_file, "\.$$\.temp");

# find the LRG ID
my $lrg_name = $lrg->findNode("fixed_annotation/id")->content();

die "ERROR: Problem with LRG identifier '$lrg_name'\n" unless $lrg_name =~ /^LRG\_[0-9]+$/;

# find mapping
my $updatable_annotation = $lrg->findNode("updatable_annotation");
my $annotation_sets = $updatable_annotation->findNodeArray("annotation_set");
my $mapping_node;
foreach my $annotation_set (@{$annotation_sets}) {
  $mapping_node = $annotation_set->findNode("mapping",{'most_recent' => '1'});
  last unless !$mapping_node;
}

# if we can't get GRCh37, get any mapping node
$mapping_node = $lrg->findNode("updatable_annotation/annotation_set/mapping") unless defined $mapping_node;

die "ERROR: Could not find mapping node in XML\n" unless defined $mapping_node;

my $assembly = $mapping_node->data->{'assembly'};
# Get the assembly from the database and compare it to the mapped assembly
my $stmt = qq{
  SELECT
    meta_value
  FROM
    meta
  WHERE
    meta_key = 'assembly.name'
};
my $sth_ref = $dbCore->dbc->db_handle->selectall_arrayref($stmt);
die "ERROR: Mapped assembly is not the same as core assembly (mapped is $assembly, core is " . $sth_ref->[0][0] . ")" unless $sth_ref->[0][0] eq $assembly;

# Statement to fetch the last auto_incremented key value
my $key_stmt = qq{
  SELECT LAST_INSERT_ID()
};

# get a slice adaptor
my $sa = $dbCore->get_SliceAdaptor();
my $slice = $sa->fetch_by_region("chromosome", $mapping_node->data->{'chr_name'});

# first we need to check if the coordinate system for this LRG exists
my $cs_id;
my $cs_name_ref = $dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT name from coord_system WHERE name="$lrg_name"});

# if it doesn't, insert into the coord_system table
if(! $cs_name_ref->[0][0]) {
  
  # get the max rank
  my $max_rank = $dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT MAX(rank) FROM coord_system;})->[0][0];
  $max_rank++;
  
  $dbCore->dbc->do(qq{INSERT INTO coord_system(species_id,name,rank,attrib) values(1,"$lrg_name",$max_rank,"default_version");});
}

# now get the ID
my $cs_id_ref = $dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT coord_system_id from coord_system WHERE name="$lrg_name";});
$cs_id = $cs_id_ref->[0][0];

print "Coord system:\t" . $cs_id . "\t" . $lrg_name . "\n";


# now we need to do the same for seq_region
my $q_seq_region_id;
my $q_seq_length = ($mapping_node->data->{'chr_end'} - $mapping_node->data->{'chr_start'}) + 1;

my $t_seq_region_id = $sa->get_seq_region_id($slice);
my $lrg_name_ref =$dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT name from seq_region WHERE name="$lrg_name"});
if (! $lrg_name_ref->[0][0]) {
  $dbCore->dbc->do(qq{INSERT INTO seq_region(name,coord_system_id,length)values("$lrg_name",$cs_id,$q_seq_length)});
}

# now get the ID
my $q_seqid_ref = $dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT seq_region_id from seq_region WHERE name="$lrg_name"});
$q_seq_region_id = $q_seqid_ref->[0][0];



# now we need to add entries to meta if needed
my $meta_ref = $dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT * FROM meta WHERE meta_key = "assembly.mapping" AND meta_value like "$lrg_name\#\%";});

if(! $meta_ref->[0][0]) {
  # get the assembly name
  my $ass_name = $dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT meta_value FROM meta WHERE meta_key = "assembly.name";})->[0][0];
  
  $dbCore->dbc->do(qq{INSERT INTO meta(species_id,meta_key,meta_value) VALUES(1,"assembly.mapping","$lrg_name\#chromosome:$ass_name");});
  $dbCore->dbc->do(qq{INSERT INTO meta(species_id,meta_key,meta_value) VALUES(1,"assembly.mapping","$lrg_name\#chromosome:$ass_name\#contig");});
  $dbCore->dbc->do(qq{INSERT INTO meta(species_id,meta_key,meta_value) VALUES(1,"assembly.mapping","$lrg_name\#chromosome:$ass_name\#supercontig");});
}

# Insert a general meta entry for 'LRG' if none is present
$stmt = qq{
  SELECT
    meta_id
  FROM
    meta
  WHERE
    meta_key = 'LRG'
};
if (!$dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0]) {
  $stmt = qq{
    INSERT INTO
      meta (
	species_id,
	meta_key,
	meta_value
      )
    VALUES (
      1,
      'LRG',
      'LRG'
    )
  };
  $dbCore->dbc->do($stmt);
}


# now we need to add entries to the assembly table
# and to the seq_region_attrib table for mismatches
foreach  my $span(@{$mapping_node->{'nodes'}}) {
  
  next unless $span->name eq 'mapping_span';
  
  # Do an extra check to make sure start and end coordinates are not flipped (start should ALWAYS be lower than end)
  my $chr_start = $span->data->{'start'};
  my $chr_end = $span->data->{'end'};
  if ($chr_end < $chr_start) {
    $chr_end = $chr_start;
    $chr_start = $span->data->{'end'};
  }
  
  # switched coords
  $dbCore->dbc->do(qq{INSERT IGNORE INTO assembly(asm_seq_region_id,cmp_seq_region_id,asm_start,asm_end,cmp_start,cmp_end,ori)values($q_seq_region_id,$t_seq_region_id,}.$span->data->{'lrg_start'}.",".$span->data->{'lrg_end'}.",".$chr_start.",".$chr_end.",".$span->data->{'strand'}."\)");
  
  # find mismatches etc
  foreach my $diff(@{$span->{'nodes'}}) {
	next unless $diff->name eq 'diff';
	
  	$dbCore->dbc->do(qq{INSERT IGNORE INTO seq_region_attrib(seq_region_id,attrib_type_id,value) VALUES($q_seq_region_id, 145, "}.$diff->data->{'lrg_start'}." ".$diff->data->{'lrg_end'}." ".$diff->data->{'lrg_sequence'}."\"\)");
  }
}


my $analysis_id;
my $analysis_description = q{
  Data from LRG database
};
my $analysis_display_label = q{
  LRG Genes
};
my $analysis_displayable = 1;
my $analysis_web_data = qq{
  {
    'colour_key' => 'rna_[status]',
    'caption' => 'LRG gene',
    'label_key' => '[text_label] [display_label]',
    'name' => 'LRG Genes',
    'default' => {
      'MultiTop' => 'gene_label',
      'contigviewbottom' => 'transcript_label',
      'MultiBottom' => 'collapsed_label',
      'contigviewtop' => 'gene_label',
      'alignsliceviewbottom' => 'as_collapsed_label',
      'cytoview' => 'gene_label'
    },
    'multi_caption' => 'LRG genes',
    'key' => 'ensembl'
  }
};
my $logic_name = 'LRG_download';
my $seq_region_id = $q_seq_region_id;
my $gene_id;
my $transcript_id;
my $translation_id;
my $exon_id;
my $biotype = 'LRG gene';

# Insert entries into the analysis and analysis_description tables if not already in there
$stmt = qq{
  SELECT
    analysis_id
  FROM
    analysis
  WHERE
    logic_name = '$logic_name'
  LIMIT 1
};
$analysis_id = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];

if (!$analysis_id) {
  $stmt = qq{
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
  $analysis_id = $dbCore->dbc->db_handle->selectall_arrayref($key_stmt)->[0][0];
  
  $stmt = qq{
    INSERT INTO
      analysis_description (
	analysis_id,
	description,
	display_label,
	displayable,
	web_data
      )
    VALUES (
      $analysis_id,
      '$analysis_description',
      '$analysis_display_label',
      $analysis_displayable,
      $analysis_web_data
    )
  };
  $dbCore->dbc->do($stmt);
}

print "Analysis:\t" . $analysis_id . "\t" . $logic_name . "\n";

my $max_exon_length = 0;
my $max_transcript_length = 0;
my $sth;

# Get the transcripts
my $transnodes = $lrg->findNodeArray('fixed_annotation/transcript');
foreach my $transnode (@{$transnodes}) {
  my $exons = $transnode->findNodeArray('exon');

# Get the transcript start and end
  my $transcript_start = $transnode->data->{'start'};
  my $transcript_end = $transnode->data->{'end'};
  my $transcript_length = $transcript_end-$transcript_start+1;
  $max_transcript_length = ($transcript_length > $max_transcript_length ? $transcript_length : $max_transcript_length);
  
# Get the start and end coordinates of the cDNA from the first and last exon (are the exons sorted in the array?)
  my $cdna_start = $exons->[0]->findNode('cdna_coords')->data->{'start'};
  my $cdna_end = $exons->[scalar(@{$exons})-1]->findNode('cdna_coords')->data->{'end'};

# Get the coding start and end coordinates
  my $cds_node = $transnode->findNode('coding_region');
  my $cds_start = $cds_node->data->{'start'};
  my $cds_end = $cds_node->data->{'end'};
  
# Insert transcript entry into db
  $stmt = qq{
    INSERT INTO
      transcript (
        analysis_id,
        seq_region_id,
        seq_region_start,
        seq_region_end,
        seq_region_strand,
        biotype,
        is_current
      )
    VALUES
      (
        $analysis_id,
        $seq_region_id,
        $transcript_start,
        $transcript_end,
        1,
        '$biotype',
        1
      )
  };
  $dbCore->dbc->do($stmt);
  $transcript_id = $dbCore->dbc->db_handle->selectall_arrayref($key_stmt)->[0][0];

  # Get the next free transcript stable_id
  my $transcript_stable_id = get_next_stable_id($dbCore->dbc,'transcript_stable_id',$lrg_name . '_T');
  # Insert transcript stable id into transcript_stable_id table
  $stmt = qq{
    INSERT INTO
      transcript_stable_id (
	transcript_id,
	stable_id,
	version,
	created_date,
	modified_date
      )
    VALUES (
      $transcript_id,
      '$transcript_stable_id',
      1,
      NOW(),
      NOW()
    ) 
  };
  $dbCore->dbc->do($stmt);
  
  if (!$gene_id) {
# Insert gene entry into db
    $stmt = qq{
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
	  canonical_transcript_id
	)
      VALUES
	(
	  '$biotype',
	  $analysis_id,
	  $seq_region_id,
	  $transcript_start,
	  $transcript_end,
	  1,
	  'LRG database',
	  1,
	  $transcript_id
	)
    };
    $dbCore->dbc->do($stmt);
    $gene_id = $dbCore->dbc->db_handle->selectall_arrayref($key_stmt)->[0][0];
 
# Get the next free gene stable_id
    my $gene_stable_id = get_next_stable_id($dbCore->dbc,'gene_stable_id',$lrg_name . '_G');
# Insert gene stable id into gene_stable_id table
    $stmt = qq{
      INSERT INTO
	gene_stable_id (
	  gene_id,
	  stable_id,
	  version,
	  created_date,
	  modified_date
	)
      VALUES (
	$gene_id,
	'$gene_stable_id',
	1,
	NOW(),
	NOW()
      ) 
    };
    $dbCore->dbc->do($stmt);
     
    print "Gene:\t" . $gene_id . "\t" . $gene_stable_id . "\n";
  }
  
  print "Transcript:\t" . $transcript_id . "\t" . $transcript_stable_id . "\n";
  
  # Update transcript table to insert the gene id
  $stmt = qq{
    UPDATE
      transcript
    SET
      gene_id = $gene_id
    WHERE
      transcript_id = $transcript_id
  };
  $dbCore->dbc->do($stmt);
  
  
  # Strand is always positive since we are in LRG-space
  $stmt = qq{
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
	?,
	?,
	1
      )
  };
  $sth = $dbCore->dbc->prepare($stmt);
  
  my $extr_stmt = qq {
    INSERT INTO
      exon_transcript (
	exon_id,
	transcript_id,
	rank
      )
    VALUES
      (
	?,
	?,
	?
      )
  };
  my $extr_sth = $dbCore->dbc->prepare($extr_stmt);
  
  my $start_exon_id;
  my $end_exon_id;
  my $exon_id;
  my $exon_count = 0;
  # Loop over the exons in the transcript and insert them into the database
  print "Exon:";
  
  foreach my $exon (@{$exons}) {
    $exon_count++;
    my $lrg_coords = $exon->findNode('lrg_coords');
    my $exon_start = $lrg_coords->data->{'start'};
    my $exon_end = $lrg_coords->data->{'end'};
    my $exon_length = ($exon_end-$exon_start+1);
    $max_exon_length = ($exon_length > $max_exon_length ? $exon_length : $max_exon_length);
    $sth->execute($exon_start,$exon_end);
    $exon_id = $dbCore->dbc->db_handle->selectall_arrayref($key_stmt)->[0][0];

    print "\t" . $exon_id;
    
# If the coding start is within this exon, save this exon id as start_exon_id and calculate the coding start offset within the exon
    if (!defined($start_exon_id) && $exon_start <= $cds_start && $exon_end >= $cds_start) {
      $start_exon_id = $exon_id;
      $cds_start = $cds_start - $exon_start + 1;
      print '[First coding]';
    }
# If the coding end is within this exon, save this exon id as end exon id and calculate end offset within the exon
    if (!defined($end_exon_id) && $exon_start <= $cds_end && $exon_end >= $cds_end) {
      $end_exon_id = $exon_id;
      $cds_end = $cds_end - $exon_start + 1;
      print '[Last coding]';
    }
    $extr_sth->execute($exon_id,$transcript_id,$exon_count);
  }
  print "\n";
  
  $stmt = qq{
    INSERT INTO
      translation (
	transcript_id,
	seq_start,
	start_exon_id,
	seq_end,
	end_exon_id
      )
    VALUES
      (
	$transcript_id,
	$cdna_start,
	$start_exon_id,
	$cdna_end,
	$end_exon_id
      )
  };
  $dbCore->dbc->do($stmt);
  $translation_id = $dbCore->dbc->db_handle->selectall_arrayref($key_stmt)->[0][0];
# Get the next free translation stable_id
# my $translation_stable_id = get_next_stable_id($dbCore->dbc,'translation_stable_id','ENSPLRG');
  my $translation_stable_id = get_next_stable_id($dbCore->dbc,'translation_stable_id',$lrg_name . '_P');
# Insert translation stable id into translation_stable_id table
  $stmt = qq{
    INSERT INTO
      translation_stable_id (
	translation_id,
	stable_id,
	version,
	created_date,
	modified_date
      )
    VALUES (
      $translation_id,
      '$translation_stable_id',
      1,
      NOW(),
      NOW()
    ) 
  };
  $dbCore->dbc->do($stmt);
  
  print "Translation:\t" . $translation_id . "\t" . $translation_stable_id . "\n";
}

my %h = (
  'gene' => $max_transcript_length,
  'transcript' => $max_transcript_length,
  'exon' => $max_exon_length,
);

while (my ($table,$length) = each(%h)) {
  $stmt = qq{
    SELECT
      max_length
    FROM
      meta_coord
    WHERE
      table_name = '$table' AND
      coord_system_id = '$cs_id'
    LIMIT 1
  };
  my $max_length = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
  
  if (!$max_length || $max_length < $length) {
    $stmt = qq{
      INSERT INTO
	meta_coord (
	  table_name,
	  coord_system_id,
    	  max_length
	)
      VALUES
	(
	  '$table',
	  $cs_id,
	  $length
	)
      ON DUPLICATE KEY UPDATE
        max_length = $length
    };
    $sth = $dbCore->dbc->do($stmt);
  }
}

# Get the next free stable id for a *_stable_id table.
sub get_next_stable_id {
  my $dbh = shift;
  my $table = shift;
  my $prefix = shift;
  
  my $stable_id;
  my $max = 0;
  my $suflength = 9;
  
 # $prefix .= 0;
  
  my $stmt = qq{
    SELECT
      MAX(stable_id)
    FROM
      $table
    WHERE
      stable_id LIKE '$prefix%'
  };
  my $max_stable_id = $dbh->db_handle->selectall_arrayref($stmt)->[0][0];
 # $prefix = substr($prefix,0,-1);
  if ($max_stable_id) {
    my ($suffix) = $max_stable_id =~ m/^$prefix([0-9]+)$/;
    $suflength = length($suffix);
 #   ($max) = $suffix =~ m/^0*([1-9].*)$/;
    $max = int($suffix);
  }
  $max++;
 # $stable_id = $prefix . sprintf("%0" . $suflength . "d",$max);
  $stable_id = $prefix . sprintf("%u",$max);
  
  return $stable_id;
};

sub usage() {
	
	print
	"Usage:
	
	perl import.lrg.pl [options] LRG.xml
	
	-r || --registry_file	registry file
	-h || --help		print this message\n";
	
	die;
}