=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     https://www.apache.org/licenses/LICENSE-2.0

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


=head1 BasePhenotypeAnnotation

This module contains the main variables and methods for the Phenotype Annotation hive pipeline.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation;

use strict;
use warnings;

use DBI qw(:sql_types);
use String::Approx qw(amatch adist);
use Algorithm::Diff qw(diff);

use Bio::EnsEMBL::Variation::Utils::SpecialChar qw(replace_char);

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');


=head2 new

  Arg [-debug] : boolean - run in debug mode, default: 0 - off
  Arg [-skip_synonyms] : boolean - default: 0
  Arg [-skip_sets] : boolean - default: 0
  Example    : $basePheno = Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation->new
                    (-debug   => 1,
                     -skip_synonyms => 1);

  Description: Constructor. Instantiates a new BasePhenotypeAnnotation object.
  Returntype : Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation
  Exceptions : none

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = bless {
    "debug"         => 0,
    "skip_synonyms" => 1,
    "skip_sets"     => 1,
  }, $class;
  $self->{pubmed_prefix} = "PMID:";

  return $self;
}

=head2 workdir

  Arg [1]    : string $workdir (optional)
               The working directory path.
  Example    : $wkdir = $obj->workdir()
  Description: Get/set the working directory path
  Returntype : string
  Exceptions : none

=cut

sub workdir {
  my ($self, $wkdir) = @_;
  $self->{workdir} = $wkdir if defined $wkdir;
  return $self->{workdir};
}


=head2 skip_synonyms

  Arg [1]    : boolean $skip_synonyms (optional)
               The new skip_synonyms flag status
  Example    : $skip_syn = $obj->skip_synonyms()
  Description: Get/set the skip_synonyms flag status
  Returntype : boolean
  Exceptions : none

=cut

sub skip_synonyms {
  my ($self, $skip_synonyms) = @_;
  $self->{skip_synonyms} = $skip_synonyms if defined $skip_synonyms;
  return $self->{skip_synonyms};
}


=head2 get_pubmed_prefix

  Example    : $pmid_prefix = $obj->get_pubmed_prefix()
  Description: Get the pubmed_prefix string ("PMID:")
  Returntype : string
  Exceptions : none

=cut

sub get_pubmed_prefix {
  my $self = shift;
  return $self->{pubmed_prefix};
}

=head2 core_db_adaptor

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor $db_adaptor (optional)
               The new core_db_adaptor
  Example    : $core_dba = $obj->core_db_adaptor()
  Description: Get/set the core_db_adaptor
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none

=cut

sub core_db_adaptor {
  my ($self, $db_adaptor) = @_;
  $self->{core_dba} = $db_adaptor if defined $db_adaptor;
  return $self->{core_dba};
}

=head2 variation_db_adaptor

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor $db_adaptor (optional)
               The new core_db_adaptor
  Example    : $variation_dba = $obj->variation_db_adaptor()
  Description: Get/set the variation_db_adaptor
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none

=cut

sub variation_db_adaptor {
  my ($self, $db_adaptor) = @_;
  $self->{variation_dba} = $db_adaptor if defined $db_adaptor;
  return $self->{variation_dba};
}


=head2 ontology_db_adaptor

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor $db_adaptor (optional)
               The new ontology_db_adaptor
  Example    : $ontology_dba = $obj->ontology_db_adaptor()
  Description: Get/set the ontology_db_adaptor
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none

=cut

sub ontology_db_adaptor {
  my ($self, $db_adaptor) = @_;
  $self->{ontology_dba} = $db_adaptor if defined $db_adaptor;
  return $self->{ontology_dba};
}


=head2 debug

  Arg [1]    : boolean $debug (optional)
               The new debug mode flag
  Example    : $debug = $obj->debug()
  Description: Get/set debug mode flag
  Returntype : boolean
  Exceptions : none

=cut

sub debug {
  my ($self, $debug) = @_;
  $self->{debug} = $debug if defined $debug;
  return $self->{debug};
}


=head2 skip_sets

Arg [1]    : boolean $skip_sets (optional)
             The new skip_sets flag
Example    : $skip_sets = $obj->skip_sets()
Description: Get/set skip_sets flag
Returntype : boolean
Exceptions : none

=cut

sub skip_sets {
  my ($self, $skip_set) = @_;
  $self->{skip_sets} = $skip_set if defined $skip_set;
  return $self->{skip_sets};
}


=head2 logFH

  Arg [1]    : FileHandle $logFH (optional)
               The new logFH
  Example    : $logF = $obj->logFH()
  Description: Get/Set the logging file handle
  Returntype : FileHandle
  Exceptions : none

=cut

sub logFH {
  my ($self, $logFH) = @_;
  $self->{logFH} = $logFH if defined $logFH;
  return $self->{logFH};
}

=head2 errFH

  Arg [1]    : FileHandle $errFH (optional)
               The new errFH
  Example    : $errF = $obj->errFH()
  Description: Get/Set for the error logging file handle
  Returntype : FileHandle
  Exceptions : none

=cut

sub errFH {
  my ($self, $errFH) = @_;
  $self->{errFH} = $errFH if defined $errFH;
  return $self->{errFH};
}


=head2 pipelogFH

  Arg [1]    : FileHandle $pipelogFH (optional)
               The new logFH
  Example    : $logF = $obj->pipelogFH()
  Description: Get/Set for the pipeline logging file handle
  Returntype : FileHandle
  Exceptions : none

=cut

sub pipelogFH {
  my ($self, $pipelogFH) = @_;
  $self->{pipelogFH} = $pipelogFH if defined $pipelogFH;
  return $self->{pipelogFH};
}


=head2 print_logFH

  Example    : $obj->print_logFH("new log message")
  Description: Write message to logging file handle
  Returntype : none
  Exceptions : none

=cut

sub print_logFH {
  my ($self, $msg) = @_;

  my $logFH = $self->logFH;

  if ($logFH) {
    print $logFH $msg;
  }
  else {
    die("Missing logFH variable\n");
  }
}


=head2 print_errFH

  Example    : $obj->print_errFH("new log message")
  Description: Write message to error logging file handle
  Returntype : none
  Exceptions : none

=cut

sub print_errFH {
  my ($self, $msg) = @_;

  my $errFH = $self->errFH;

  if ($errFH) {
    print $errFH $msg;
  }
  else {
    die("Missing errFH variable\n");
  }
}


=head2 print_pipelogFH

  Example    : $obj->print_pipelogFH("new log message")
  Description: Write message to pipeline logging file handle
  Returntype : none
  Exceptions : none

=cut

sub print_pipelogFH {
  my $self = shift;
  my $msg = shift;
  my $pipelogFH = $self->pipelogFH;

  if ($pipelogFH) {
    print $pipelogFH $msg;
  }
  else {
    die("Missing pipelogFH variable\n");
  }
}



#----------------------------
# PUBLIC METHODS

=head2 get_seq_region_ids

  Example    : $obj->get_seq_region_ids()
  Description: Fetch seq_region_id(s) from varaition db
  Returntype : Hash of seq_region_id identifiers
  Exceptions : none


=cut

sub get_seq_region_ids {
  my $self = shift;

  my $sth = $self->variation_db_adaptor->dbc->prepare(qq{
    SELECT seq_region_id, name
    FROM seq_region
  });
  $sth->execute;

  my (%seq_region_ids, $id, $name);
  $sth->bind_columns(\$id, \$name);
  $seq_region_ids{$name} = $id while $sth->fetch();
  $sth->finish;

  return \%seq_region_ids;
}


=head2 get_or_add_source

  Arg [1]    : hashref $source_info
  Example    : $obj->get_or_add_source($source_info)
  Description: Inserts or Updates the source in the variation database and returns the source_id.
  Returntype : Integer source_id
  Exceptions : none

=cut

sub get_or_add_source {
  my ($self, $source_info) = @_;

  my $db_adaptor = $self->variation_db_adaptor;

  my $stmt = qq{
    SELECT
      source_id
    FROM
      source
    WHERE
      name = '$source_info->{source_name}'
    LIMIT 1
  };
  my $sth = $db_adaptor->dbc->prepare($stmt);
  $sth->execute();
  my $source_id;
  $sth->bind_columns(\$source_id);
  $sth->fetch();

  if (!defined($source_id)) {
    $stmt = qq{
      INSERT INTO
        source (name, description, url, somatic_status, version )
      VALUES (
        '$source_info->{source_name}',
        '$source_info->{source_description}',
        '$source_info->{source_url}',
        '$source_info->{source_status}',
        $source_info->{source_version}
      )
    };
    $db_adaptor->dbc->do($stmt);
    $source_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};

    $self->print_logFH("Added source for $source_info->{source_name} (source_id = $source_id)\n") if ($self->{debug});
  }
  else {
    $stmt = qq{
      UPDATE
        source
      SET name=?,
        description=?,
        url=?,
        version=?
      WHERE
        source_id=?
    };
    my $update_source_sth =$db_adaptor->dbc->prepare($stmt);
    $update_source_sth->execute($source_info->{source_name},
                                $source_info->{source_description},
                                $source_info->{source_url},
                                $source_info->{source_version},
                                $source_id);
  }

  return $source_id;
}


=head2 save_phenotypes

  Arg [1]    : hashref $source_info
  Arg [2]    : hashref $input_data
  Example    : $self->save_phenotypes(\%source_info, $input_data);
  Description: Save the phenotype data ($input_data) for the given source ($source_info).
  Returntype : none
  Exceptions : none

=cut

sub save_phenotypes {
  my ($self, $source_info, $input_data) = @_;

  my $core_dba = $self->core_db_adaptor;
  my $phenotype_dba =  $self->variation_db_adaptor->get_PhenotypeAdaptor;

  my $set = defined($source_info->{set}) ? $source_info->{set} : undef;
  $source_info->{source_status} = 'germline' unless defined($source_info->{source_status});

  my $initial_phenotype_count = $phenotype_dba->generic_count;
  my @attrib_types = @{$self->_get_attrib_types()};
  my %phenotype_cache = map {$_->description() => $_->dbID()} @{$phenotype_dba->fetch_all};
  $self->{phenotype_cache} = \%phenotype_cache;

  my @ids;
  my %synonym;
  my @phenotypes;
  if (exists($input_data->{'synonyms'})) {
    %synonym = %{$input_data->{'synonyms'}};
    # To get all the ids of the source (Uniprot)
    @ids = keys(%synonym);
  }
  if (exists($input_data->{'phenotypes'})) {
    @phenotypes = @{$input_data->{'phenotypes'}};
  }
  $self->print_logFH( "Got ".(scalar @phenotypes)." phenotype objects\n") if ($self->{debug});

  # Get internal variation ids for the rsIds
  $self->print_logFH( "Retrieving internal variation IDs\n") if ($self->{debug});
  if (scalar @ids == 0) {
    @ids = map {$_->{'id'}} @phenotypes;
  }
  my $variation_ids = $source_info->{object_type} =~ /Variation/ ? $self->_get_dbIDs(\@ids) : {};

  # Get coordinates of objects
  my $coords;

  $self->print_logFH( "Retrieving object coordinates\n") if ($self->{debug});

  # might be able to copy them from data (QTLs and Genes)
  if(defined($phenotypes[0]->{seq_region_id})) {
    foreach my $p(@phenotypes) {
      my $coord = {};
      foreach my $key(qw(id start end strand)) {
        $coord->{'seq_region_'.$key} = $p->{'seq_region_'.$key};
      }
      push @{$coords->{$p->{id}}}, $coord;
    }
  }

  $coords ||= $self->_get_coords(\@ids,$variation_ids, $source_info->{object_type}, $source_info->{object_type} eq 'Gene' ? $self->core_db_adaptor : $self->variation_db_adaptor);

  # uniquify coords
  foreach my $id (keys %$coords) {
    my %tmp = map {
      $_->{seq_region_id}."_".
      $_->{seq_region_start}."_".
      $_->{seq_region_end}."_".
      $_->{seq_region_strand}."_" => $_
    } @{$coords->{$id}};

    $coords->{$id} = [values %tmp];
  }

  # Get or add a source
  my $source_id = $self->get_or_add_source($source_info);
  $self->print_logFH( "$source_info->{source_name} source_id is $source_id\n");

  # Add the synonyms if required
  unless ($self->{skip_synonyms}) {
    $self->print_logFH( "Adding synonyms\n") if ($self->{debug});
    $self->_add_synonyms(\%synonym,$variation_ids,$source_id);
  }

  # Now, insert phenotypes
  die("ERROR: No phenotypes or objects retrieved from input\n") unless scalar @phenotypes;

  $self->print_logFH( "Adding phenotypes\n") if ($self->{debug});
  $source_info->{source_id}=$source_id;
  my %data = (phenotypes => \@phenotypes, coords => $coords,
              variation_ids => $variation_ids, attrib_types=> \@attrib_types);
  $self->_add_phenotypes(\%data,$source_info);

  $self->print_logFH( "$initial_phenotype_count initial phenotypes\n") if ($self->{debug});
  my $added_phenotypes = $phenotype_dba->generic_count - $initial_phenotype_count;
  $self->print_logFH("$added_phenotypes new phenotypes added\n") if ($self->{debug});


  # Add the variation sets if required
  unless ($self->{skip_sets}) {
    if (%synonym && !$self->{skip_synonyms}) {
      $self->print_logFH("Adding variation sets for synonyms\n") if ($self->{debug});
      $self->_add_set($set,$source_id,'synonym');
    }
    elsif (@phenotypes) {
      $self->print_logFH("Adding variation sets for phenotypes\n") if ($self->{debug});
      $self->_add_set($set,$source_id,'phenotype');
    }
  }
}


=head2 dump_phenotypes

  Arg [1]    : String $source_name
  Arg [2]    : boolean $clean (optional) - delete phenotype data from db (default: 0)
  Example    : $self->dump_phenotypes($source_id, 1);
  Description: Dump the existing phenotype_features, phenotype_features_attribs
               for the particular source and removes them if clean option selected.
               $clean option removes the phenotype feature data including phenotypes
               and phenotype_ontology_accessions that is not attached to any phenotype_feature.
  Returntype : none
  Exceptions : none

=cut

sub dump_phenotypes {
  my ($self, $source_name, $clean) = @_;

  die ("source_name needs to be specified!\n") unless defined $source_name && defined $self->workdir;
  $clean ||= 0;

  # Prepared statements
  my $pfa_select_stmt = qq{
    SELECT pfa.*
    FROM phenotype_feature_attrib pfa, phenotype_feature pf, source s
    WHERE pfa.phenotype_feature_id = pf.phenotype_feature_id
        AND pf.source_id = s.source_id
        AND s.name = \"$source_name\"
  };
  my $pfa_delete_stmt = qq{
    DELETE pfa.*
    FROM phenotype_feature_attrib pfa, phenotype_feature pf, source s
    WHERE pfa.phenotype_feature_id = pf.phenotype_feature_id
      AND pf.source_id = s.source_id
      AND s.name = \"$source_name\"
  };

  my $pf_select_stmt = qq{
    SELECT pf.*
    FROM phenotype_feature pf, source s
    WHERE pf.source_id = s.source_id
    AND s.name = \"$source_name\"
  };
  my $pf_delete_stmt = qq{
    DELETE pf.*
    FROM phenotype_feature pf, source s
    WHERE pf.source_id = s.source_id
    AND s.name = \"$source_name\"
  };

  my $p_extra_select_stmt = qq{
    SELECT *
    FROM phenotype p LEFT JOIN phenotype_feature pf ON pf.phenotype_id = p.phenotype_id
    WHERE pf.phenotype_id IS null;
  };
  my $p_extra_delete_stmt = qq{
    DELETE p.*
    FROM phenotype p LEFT JOIN phenotype_feature pf ON pf.phenotype_id = p.phenotype_id
    WHERE pf.phenotype_id IS null;
  };

  my $poa_extra_select_stmt = qq{
    SELECT *
    FROM phenotype_ontology_accession poa LEFT JOIN phenotype_feature pf ON pf.phenotype_id = poa.phenotype_id
    WHERE pf.phenotype_id IS null;
  };
  my $poa_extra_delete_stmt = qq{
    DELETE poa.*
    FROM phenotype_ontology_accession poa LEFT JOIN phenotype_feature pf ON pf.phenotype_id = poa.phenotype_id
    WHERE pf.phenotype_id IS null;
  };

  my $db_adaptor    = $self->variation_db_adaptor;

  _sql_to_file($pfa_select_stmt, $db_adaptor, $self->workdir."/"."pfa_".$source_name.".txt");
  _sql_to_file($pf_select_stmt, $db_adaptor, $self->workdir."/"."pf_".$source_name.".txt");
  _sql_to_file($p_extra_select_stmt, $db_adaptor, $self->workdir."/"."p_extra_".$source_name.".txt");
  _sql_to_file($poa_extra_select_stmt, $db_adaptor, $self->workdir."/"."poa_extra_".$source_name.".txt");

  if ($clean) {
    my $sth = $db_adaptor->dbc->prepare($pfa_delete_stmt);
    $sth->execute();

    $sth = $db_adaptor->dbc->prepare($pf_delete_stmt);
    $sth->execute();

    $sth = $db_adaptor->dbc->prepare($p_extra_delete_stmt);
    $sth->execute();

    $sth = $db_adaptor->dbc->prepare($poa_extra_delete_stmt);
    $sth->execute();
  }
}

#----------------------------
# PRIVATE METHODS

# getter for the internal phenotype_cache shared between methods
sub _phenotype_cache {
  my $self = shift;
  return $self->{phenotype_cache};
}

# getter for the internal sql statments handles shared between methods
sub _sql_statements {
  my $self = shift;
  return $self->{sql_statements};
}

# run the sql (sql_stmt) on the given db (db_adaptor) and save results to file (outFile)
sub _sql_to_file{
  my ($sql_stmt, $db_adaptor, $outFile) = @_;

  my $sth = $db_adaptor->dbc->prepare($sql_stmt);
  $sth->execute();

  open (my $fhDump, ">", $outFile) || die ("Failed to open file".$outFile.": $!\n");
  while (my $row = $sth->fetchrow_arrayref()) {
    for(@$row) { $_ = "NULL" if !defined $_; } #some values are undef e.g. study_id, this is expected
    print $fhDump join( "\t", map {qq($_)} @$row), "\n";
  }
  close($fhDump);
}

# insert the phenotype data (including phenotype feature attribs) into the db
sub _add_phenotypes {
  my ($self, $data, $source_info) = @_;

  my $db_adaptor    = $self->variation_db_adaptor;

  my $phenotypes    = $data->{phenotypes};
  my $coords        = $data->{coords};
  my @attrib_types  = @{$data->{attrib_types}};
  my $variation_ids = $data->{variation_ids};

  my $object_type = $source_info->{object_type};
  my $source_id   = $source_info->{source_id};
  my $threshold   = $source_info->{threshold};

  # Prepared statements
  my $st_ins_stmt = qq{
    INSERT INTO
      study ( source_id, external_reference, study_type, description
    )
    VALUES (
      $source_id, ?, ?, ?
    )
  };

  my $extra_cond = '';
  my $left_join = '';

  # add special joins for checking certain sources
  if ($source_info->{source_name_short} =~ m/GWAS/i) {
    $left_join = qq{
      LEFT JOIN
      (
        phenotype_feature_attrib pfa
        JOIN attrib_type at
        ON pfa.attrib_type_id = at.attrib_type_id
      )
      ON pf.phenotype_feature_id = pfa.phenotype_feature_id
    };

    $extra_cond = qq{ AND at.code = "p_value" AND pfa.value = ? };
  }

  my $pf_check_stmt = qq{
    SELECT
      pf.phenotype_feature_id
    FROM
      phenotype_feature pf
    $left_join
    WHERE
      pf.object_id = ? AND
      pf.type = ? AND
      pf.phenotype_id = ? AND
      pf.source_id = ? AND
      (pf.study_id = ? OR pf.study_id IS NULL)
      $extra_cond
    LIMIT 1
  };

  my $pf_ins_stmt = qq{
    INSERT INTO phenotype_feature (
      phenotype_id, source_id, study_id,
      type, object_id, is_significant,
      seq_region_id, seq_region_start, 
      seq_region_end, seq_region_strand
    )
    VALUES (
      ?, ?, ?,
      ?, ?, ?,
      ?, ?, 
      ?, ?
    )
  };

  my $attrib_ins_stmt = qq{
    INSERT INTO phenotype_feature_attrib (
      phenotype_feature_id,
      attrib_type_id,
      value
    )
    SELECT ?, at.attrib_type_id, ?
    FROM attrib_type at
    WHERE at.code = ?
  };

  my $attrib_ins_cast_stmt = qq{
    INSERT INTO phenotype_feature_attrib (
      phenotype_feature_id,
      attrib_type_id,
      value
    )
    SELECT ?, at.attrib_type_id, CAST(? AS CHAR)
    FROM attrib_type at
    WHERE at.code = ?
  };

  my $attrib_id_ext_stmt = qq{
    SELECT attrib_id, value from attrib, attrib_type
    WHERE attrib_type.code = 'ontology_mapping'
    AND attrib.attrib_type_id = attrib_type.attrib_type_id
  };

  my $ontology_accession_ins_stmt = qq{
    INSERT IGNORE INTO phenotype_ontology_accession
    (phenotype_id, accession, mapped_by_attrib, mapping_type)
    values (?,?,?,?)
   };

  my $st_ins_sth   = $db_adaptor->dbc->prepare($st_ins_stmt);
  my $pf_check_sth = $db_adaptor->dbc->prepare($pf_check_stmt);
  my $pf_ins_sth   = $db_adaptor->dbc->prepare($pf_ins_stmt);
  my $attrib_ins_sth = $db_adaptor->dbc->prepare($attrib_ins_stmt);
  my $attrib_ins_cast_sth = $db_adaptor->dbc->prepare($attrib_ins_cast_stmt);
  my $ontology_accession_ins_sth = $db_adaptor->dbc->prepare($ontology_accession_ins_stmt);

  $self->{sql_statements}{st_ins_sth} = $st_ins_sth;

  ## get the attrib id for the type of description to ontology term linking
  my $attrib_id_ext_sth = $db_adaptor->dbc->prepare($attrib_id_ext_stmt);
  $attrib_id_ext_sth->execute();
  my $ont_attrib_type = $attrib_id_ext_sth->fetchall_hashref("value");

  if ($source_info->{source_name} eq 'RGD'){ #required as db source:name is 'RGD' and attrib_type_id 509 for ontology_mapping is 'Rat Genome Database'
    $source_info->{source_attrib_type} = 'Rat Genome Database';
  }
  my $mapped_by = 'Data source';
  if (exists $source_info->{source_mapped_attrib_type} && exists $ont_attrib_type->{$source_info->{source_mapped_attrib_type}}){
    $mapped_by = $source_info->{source_mapped_attrib_type};
  }

  # First, sort the array according to the phenotype description
  my @sorted = sort {($a->{description} || $a->{name}) cmp ($b->{description} || $b->{name})} @{$phenotypes};
  $self->{study_count} = 0;
  my $phenotype_feature_count = 0;

  my $total = scalar @sorted;
  my $i = 0;

  while (my $phenotype = shift(@sorted)) {
    $self->_progress($i++, $total);

    $object_type = $phenotype->{type} if defined($phenotype->{type});

    my $object_id = $phenotype->{"id"};
    # If the rs could not be mapped to a variation id, skip it
    if ($object_type =~ /Variation/ && (!defined($variation_ids->{$object_id}))){
      $self->print_errFH( "WARN: skipping rsID: could not be mapped to a variation id : $object_id \n");
      next;
    }

    $object_id = $variation_ids->{$object_id}[1] if ($object_type eq 'Variation');

    my $study_id = $self->_get_study_id($phenotype,$source_id);

    # get phenotype ID
    my $phenotype_id = $self->_get_phenotype_id($phenotype);

    foreach my $acc (@{$phenotype->{accessions}}){
      $acc =~ s/\s+//g;
      $acc = $self->_iri2acc($acc) if $acc =~ /^http/;
      $ontology_accession_ins_sth->execute( $phenotype_id, $acc,  $ont_attrib_type->{$mapped_by}->{'attrib_id'}, $phenotype->{ontology_mapping_type} ) || die ("Failed to import phenotype accession\n");
    }
    if ($phenotype->{"associated_gene"}) {
      $phenotype->{"associated_gene"} =~ s/\s//g;
      $phenotype->{"associated_gene"} =~ s/;/,/g;
    }

    $phenotype->{"p_value"} = $self->_convert_p_value($phenotype->{"p_value"}) if (defined($phenotype->{"p_value"}));

    # Check if this phenotype_feature already exists for this variation and source, in that case we probably want to skip it
    my $pf_id;

    $pf_check_sth->bind_param(1,$object_id,SQL_VARCHAR);
    $pf_check_sth->bind_param(2,$object_type,SQL_VARCHAR);
    $pf_check_sth->bind_param(3,$phenotype_id,SQL_INTEGER);
    $pf_check_sth->bind_param(4,$source_id,SQL_INTEGER);
    $pf_check_sth->bind_param(5,$study_id,SQL_INTEGER);

    # For nhgri-ebi gwas data
    if ($source_info->{source_name_short} =~ m/GWAS/i) {
      $pf_check_sth->bind_param(6,$phenotype->{"p_value"},SQL_VARCHAR);
    }

    $pf_check_sth->execute();
    $pf_check_sth->bind_columns(\$pf_id);
    $pf_check_sth->fetch();
    if (defined($pf_id)) {
      $self->print_errFH( "WARN: skipping phenotype_feature ($pf_id), phenotype($phenotype_id), source ($source_id) as it already exists for this variation and source.\n");
      next;
    }

    my $is_significant = defined($threshold) ? ($phenotype->{"p_value"} && $phenotype->{"p_value"} < $threshold ? 1 : 0) : 1;

    # Else, insert this variation annotation
    foreach my $coord(@{$coords->{$object_id}}) {
      $pf_ins_sth->bind_param(1,$phenotype_id,SQL_INTEGER);
      $pf_ins_sth->bind_param(2,$source_id,SQL_INTEGER);
      $pf_ins_sth->bind_param(3,$study_id,SQL_INTEGER);
      $pf_ins_sth->bind_param(4,$object_type,SQL_VARCHAR);
      $pf_ins_sth->bind_param(5,$object_id,SQL_VARCHAR);
      $pf_ins_sth->bind_param(6,$is_significant,SQL_INTEGER);
      $pf_ins_sth->bind_param(7,$coord->{seq_region_id},SQL_INTEGER);
      $pf_ins_sth->bind_param(8,$coord->{seq_region_start},SQL_INTEGER);
      $pf_ins_sth->bind_param(9,$coord->{seq_region_end},SQL_INTEGER);
      $pf_ins_sth->bind_param(10,$coord->{seq_region_strand},SQL_INTEGER);
      $pf_ins_sth->execute();
      $phenotype_feature_count++;

      # get inserted ID
      $pf_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};

      # add attribs
      foreach my $attrib_type(grep {defined($phenotype->{$_}) && $phenotype->{$_} ne ''} @attrib_types) {
        my $value = $phenotype->{$attrib_type};
        my $sth = $value =~ m/^\d+(\.\d+)?$/ ? $attrib_ins_cast_sth : $attrib_ins_sth;
        $sth->bind_param(1,$pf_id,SQL_INTEGER);
        $sth->bind_param(2,$value,SQL_VARCHAR);
        $sth->bind_param(3,$attrib_type,SQL_VARCHAR);
        $sth->execute();
      }
    }
  }

  $self->_end_progress();
  $self->print_logFH( "$self->{study_count} new studies added\n") if ($self->{debug});
  $self->print_logFH( "$phenotype_feature_count new phenotype_features added\n") if ($self->{debug});

}

# insert the synonym data into variation_synonym table
sub _add_synonyms {
  my ($self, $synonyms, $variation_ids, $source_id) = @_;

  my $db_adaptor = $self->variation_db_adaptor;

  # If we actually didn't get any synonyms, just return
  return if (!defined($synonyms) || !scalar(keys(%{$synonyms})));

  # Some prepeared statements needed for inserting the synonyms into database
  my $ins_stmt = qq{
    INSERT IGNORE INTO
      variation_synonym (
      variation_id,
      source_id,
      name
      )
    VALUES (
      ?,
      $source_id,
      ?
    )
  };
  my $ins_sth = $db_adaptor->dbc->prepare($ins_stmt);

  my $alt_count = 0;
  my $variation_count = 0;

  foreach my $rs_id (keys %{$variation_ids}) {

    my $var_id = $variation_ids->{$rs_id}[0];

    # If we have a variation id, we can proceed
    if (defined($var_id)) {

      $variation_count++;

      $ins_sth->bind_param(1,$var_id,SQL_INTEGER);

      # Handle all synonym ids for this rs_id
      while (my $alt_id = shift(@{$synonyms->{$rs_id}})) {

        # Add the id as synonym, if it is already present, it will just be ignored
        $ins_sth->bind_param(2,$alt_id,SQL_VARCHAR);
        $ins_sth->execute();
        $alt_count++;
      }
    } else {
      $self->print_logFH( "failed to add synonyms as variation was missing for $rs_id \n");
    }
  }

  $self->print_logFH( "Added $alt_count synonyms for $variation_count rs-ids\n");
}

# update variation_set_variation
sub _add_set {
  my ($self, $set, $source_id, $set_from) = @_;

  my $db_adaptor = $self->variation_db_adaptor;
  return if (!defined($set));
  return if (!defined($set_from) && ($set_from ne 'phenotype' || $set_from ne 'synonym'));

  my $variation_set_id;

  # Get variation_set_id
  my $select_set_stmt = qq{
  SELECT v.variation_set_id
  FROM variation_set v, attrib a
  WHERE v.short_name_attrib_id=a.attrib_id 
  AND a.value = ?
  };
  my $sth1 = $db_adaptor->dbc->prepare($select_set_stmt);
  $sth1->bind_param(1,$set,SQL_VARCHAR);
  $sth1->execute();
  $sth1->bind_columns(\$variation_set_id);
  $sth1->fetch();
  return if (!defined($variation_set_id));

  # Insert into variation_set_variation
  my $insert_set_stmt = qq{
    INSERT IGNORE INTO variation_set_variation (variation_id,variation_set_id)
    SELECT distinct v.variation_id, ? 
  };
  if ($set_from eq 'phenotype') {
    $insert_set_stmt .= qq{ 
      FROM phenotype_feature pf, variation v WHERE 
        v.name=pf.object_id AND
        pf.type='Variation' AND
        pf.source_id=?
    };
  }
  elsif ($set_from eq 'synonym') {
    $insert_set_stmt .= qq{ 
      FROM variation_synonym vs, variation v WHERE 
        v.variation_id=vs.variation_id AND
        vs.source_id=?
    };
  }

  my $sth2 = $db_adaptor->dbc->prepare($insert_set_stmt);
  $sth2->bind_param(1,$variation_set_id,SQL_INTEGER);
  $sth2->bind_param(2,$source_id,SQL_INTEGER);
  $sth2->execute();
}

# fetch coordinates for variation, gene, or structural_variation
sub _get_coords {
  my ($self, $ids, $variation_ids, $type, $db_adaptor) = @_;

  my $tables;
  my $where_clause;

  my @object_ids = ($type eq 'Variation') ? map { $variation_ids->{$_}[1] } keys(%$variation_ids): @$ids;

  if($type eq 'Variation') {
    $tables = 'variation_feature f, variation v';
    $where_clause = 'v.variation_id = f.variation_id AND v.name = ?';
  }
  elsif($type =~ /Structural/) {
    $tables = 'structural_variation_feature f, structural_variation v';
    $where_clause = 'v.structural_variation_id = f.structural_variation_id AND v.variation_name = ?';
  }
  elsif($type eq 'Gene') {
    $tables = 'gene f';
    $where_clause = 'f.stable_id = ?';
  }

  my $sth = $db_adaptor->dbc->prepare(qq{
    SELECT
      f.seq_region_id, f.seq_region_start, f.seq_region_end, f.seq_region_strand
    FROM
      $tables
    WHERE
      $where_clause
  });

  my $coords = {};
  my ($sr_id, $start, $end, $strand);

  foreach my $id (@object_ids) {
    $sth->bind_param(1,$id,SQL_VARCHAR);
    $sth->execute();
    $sth->bind_columns(\$sr_id, \$start, \$end, \$strand);

    push @{$coords->{$id}}, {
      seq_region_id => $sr_id,
      seq_region_start => $start,
      seq_region_end => $end,
      seq_region_strand => $strand
    } while $sth->fetch();
  }

  $sth->finish();

  return $coords;
}

# look up variation ids via rsIDs and synonyms
sub _get_dbIDs {
  my ($self, $rs_ids) = @_;

  my $id_stmt = qq{
    SELECT DISTINCT
      v.variation_id, v.name
    FROM
      variation v,
      variation_feature vf
    WHERE
      v.name = ? AND
      v.variation_id=vf.variation_id
    LIMIT 1
  };
  my $syn_stmt = qq{
    SELECT DISTINCT
      v.variation_id, v.name
    FROM
      variation_feature vf,
      variation_synonym vs JOIN
      variation v ON vs.variation_id = v.variation_id
    WHERE
      vs.name = ? AND
      v.variation_id=vf.variation_id
    LIMIT 1
  };
  my $id_sth = $self->variation_db_adaptor->dbc->prepare($id_stmt);
  my $syn_sth = $self->variation_db_adaptor->dbc->prepare($syn_stmt);

  my %mapping;

  foreach my $rs_id (@{$rs_ids}) {
    $id_sth->bind_param(1,$rs_id,SQL_VARCHAR);
    $id_sth->execute();
    my ($var_id,$var_name);
    $id_sth->bind_columns(\$var_id,\$var_name);
    $id_sth->fetch();

    # If we couldn't find the rs_id, look in the synonym table
    if (!defined($var_id)) {
      $syn_sth->bind_param(1,$rs_id,SQL_VARCHAR);
      $syn_sth->execute();
      $syn_sth->bind_columns(\$var_id,\$var_name);
      $syn_sth->fetch();
    }
    $self->print_errFH("$rs_id - no mapping found in variation db\n") unless $var_id ;
    $mapping{$rs_id} = [$var_id,$var_name] if $var_id && $var_name;
  }

  return \%mapping;
}

# fetch attrib types from db
sub _get_attrib_types {
  my $self = shift;

  my $sth = $self->variation_db_adaptor->dbc->prepare(qq{
    SELECT code
    FROM attrib_type
  });
  $sth->execute();

  my ($attrib_type, @tmp_types);
  $sth->bind_columns(\$attrib_type);
  while ($sth->fetch()) {
    push @tmp_types, $attrib_type if ($attrib_type ne 'description');
  }
  $sth->finish;

  return \@tmp_types;
}

# clean up + search for phenotype in db, if not found it gets inserted
sub _get_phenotype_id {
  my ($self, $phenotype) = @_;

  my %phenotype_cache = %{$self->_phenotype_cache};

  my ($name, $description);
  $name = $phenotype->{name};
  $description = $phenotype->{description};

  # Clean up
  $description =~ s/^\s+|\s+$//g; # Remove spaces at the beginning and the end of the description
  $description =~ s/\n//g; # Remove 'new line' characters
  $description =~ s/[\(\)]//g; # Remove characters ( )

  # Replace special characters in the phenotype description
  $description = replace_char($description);

  # Check phenotype description in the format "description; name"
  if (!defined($name) || $name eq '') {
    my ($p_desc,$p_name) = split(";",$description);
    if ($p_name) {
      $p_name =~ s/ //g;
      if ($p_name =~ /^\w+$/) {
        $description = $p_desc;
        $name = $p_name;
      }
    }
  }

  if(scalar keys %phenotype_cache) {

    # check cache first
    return $phenotype_cache{$description} if defined $phenotype_cache{$description};

    my @tmp = keys %phenotype_cache;

    # lc everything
    my $description_bak = $description;
    $description = lc($description);

    # store mapped
    my %mapped;
    $mapped{lc($_)} = $_ for @tmp;
    @tmp = keys %mapped;

    # check if it matches lc
    if(defined($mapped{$description})) {
      $phenotype_cache{$description_bak} = $phenotype_cache{$mapped{$description}};
      return $phenotype_cache{$description_bak};
    }

    # try a fuzzy match using String::Approx
    my @matches = amatch($description, [10], @tmp);

    if(@matches) {

      # we only want the best match
      my $best = scalar @matches == 1 ? $matches[0] : (sort {abs(adist($description, $a)) <=> abs(adist($description, $b))} @matches)[0];
      $self->print_errFH( "\nPHENOTYPE:\n\tINPUT: $description\n\tBEST:  $best\n\tDIST: ".adist($description, $best)."\n");

      my $skip = 0;

      # Assuming a perfect match check has been done in the previous lines, e.g. if ($phenotype_cache{$description})
      $skip = 1 if (adist($description, $best) == 0);

      # find characters that differ
      my $diff = diff([split(//,$description)], [split(//,$best)]);

      # skip if mismatch is anything word-like
      my $diff_string = '';
      my $previous_diff_pos;
      foreach(map {@$_} @$diff) {
        my $diff_pos  = $_->[1];
        my $diff_char = $_->[2];

        $skip = 1 if $diff_char =~ /\w/i;

        if (!$previous_diff_pos) {
          $diff_string .= "[$diff_char";
        }
        elsif ($diff_pos == $previous_diff_pos) {
          $diff_string .= "/$diff_char]";
        }
        elsif ($diff_pos == ($previous_diff_pos+1)) {
          $diff_string .= $diff_char;
        }
        else {
          $diff_string .= "] [$diff_char";
        }
        $previous_diff_pos = $diff_pos;
      }
      if ( $diff_string ne '') {
        $diff_string .= ']' if ($diff_string !~ /\]$/);
        $self->print_errFH( "\tDIFFERENCE(S): $diff_string\n") if ( $diff_string ne '');
      }

      # cache this match so we don't have to fuzz again
      $phenotype_cache{$description_bak} = $phenotype_cache{$mapped{$best}};
      unless ($skip) {
        $self->print_errFH( "\tUSED (with diff): $best\n");
        return $phenotype_cache{$mapped{$best}};
      }
    }

    # restore from backup before inserting
    $description = $description_bak;
  }

  # finally if no match, do an insert
  my $sth = $self->variation_db_adaptor->dbc->prepare(qq{
    INSERT IGNORE INTO phenotype ( name, description ) VALUES ( ?,? )
  });

  $sth->bind_param(1,$name,SQL_VARCHAR);
  $sth->bind_param(2,$description,SQL_VARCHAR);
  $sth->execute();
  my $phenotype_id = $self->variation_db_adaptor->dbc->db_handle->{'mysql_insertid'};

  # update cache
  $phenotype_cache{$description} = $phenotype_id;
  $self->{phenotype_cache} = \%phenotype_cache;

  return $phenotype_id;
}

# clean up (truncated pubmed ids) + search for study in db, if not found it gets inserted
sub _get_study_id {
  my ($self, $phenotype, $source_id) = @_;

  my $st_ins_sth = $self->_sql_statements->{st_ins_sth};

  my $study_id;
  if(defined($phenotype->{study}) || defined($phenotype->{study_description}) || defined($phenotype->{study_type})) {

    my $sql_study = '= ?';
    my $sql_type  = '= ?';

    # To avoid duplication of study entries
    if (!defined $phenotype->{"study"}) {$sql_study = 'IS NULL'; }
    if (!defined $phenotype->{"study_type"}) {$sql_type = 'IS NULL'; }

    my $st_check_stmt = qq{
      SELECT study_id
      FROM study
      WHERE
      source_id = $source_id AND
      external_reference $sql_study AND
      study_type $sql_type
      LIMIT 1
    };

    my $st_check_sth = $self->variation_db_adaptor->dbc->prepare($st_check_stmt);
    my $param_num = 1;

    if (defined $phenotype->{"study"}) {
      if (length($phenotype->{"study"}) > 255) {
        $self->print_errFH( "WARNING: study external_references truncated search FROM:>".$phenotype->{"study"}. "<\n");
        $phenotype->{"study"} = substr($phenotype->{"study"}, 0, 254);
        $phenotype->{"study"} = substr($phenotype->{"study"}, 0,rindex($phenotype->{"study"}, ",PMID"));
        $self->print_errFH( "WARNING: study external_references truncated search TO  :>".$phenotype->{"study"}. "<\n");
      }
      $st_check_sth->bind_param($param_num,$phenotype->{"study"},SQL_VARCHAR);
      $param_num++;
    }
    if (defined $phenotype->{"study_type"}) {
      $st_check_sth->bind_param($param_num,$phenotype->{"study_type"},SQL_VARCHAR);
      $param_num++;
    }
    $st_check_sth->execute();
    $st_check_sth->bind_columns(\$study_id);
    $st_check_sth->fetch();

    if (!defined($study_id)) {
      if (length($phenotype->{"study"}) > 255) {
        $self->print_errFH( "WARNING: study external_references truncated FROM:>".$phenotype->{"study"}. "<\n");
        $phenotype->{"study"} = substr($phenotype->{"study"}, 0, 254);
        $phenotype->{"study"} = substr($phenotype->{"study"}, 0,rindex($phenotype->{"study"}, ",PMID"));
        $self->print_errFH( "WARNING: study external_references truncated TO  :>".$phenotype->{"study"}. "<\n");
      }
      $st_ins_sth->bind_param(1,$phenotype->{"study"},SQL_VARCHAR);
      $st_ins_sth->bind_param(2,$phenotype->{"study_type"},SQL_VARCHAR);
      $st_ins_sth->bind_param(3,$phenotype->{"study_description"},SQL_VARCHAR);
      $st_ins_sth->execute();
      $study_id = $self->variation_db_adaptor->dbc->db_handle->{'mysql_insertid'};
      $self->{study_count}++;
    }
  }
  return $study_id;
}


# update or initiate progress bar
sub _progress {
  my ($self, $i, $total) = @_;

  my $width = 60;
  my $percent = int(($i/$total) * 100);
  my $numblobs = int((($i/$total) * $width) - 2);

  # this ensures we're not writing to the terminal too much
  return if defined($self->{prev_prog}) && $numblobs.'-'.$percent eq $self->{prev_prog};
  $self->{prev_prog} = $numblobs.'-'.$percent;

  $self->print_pipelogFH( sprintf ("\r% -${width}s% 1s% 10s", '['.('=' x $numblobs).($numblobs == $width - 2 ? '=' : '>'), ']', "[ " . $percent . "% ]"));
}

# end progress bar
sub _end_progress {
  my $self = shift;
  $self->_progress(1,1);
  $self->print_pipelogFH("\n");
}

# standardize p_values
sub _convert_p_value {
  my ($self, $pval) = @_;

  my $sci_pval = '';
  # If a scientific format is not found, then ...
  if ($pval !~ /^\d+.*e.+$/i) {
    # If a range format is found (e.g. 10^-2 > p > 10^-3)
    if ($pval =~ /^\d+\^(-\d+)/) {
      if (length("$1")==1) { $1 = "0$1"; }
      $sci_pval = "1.00e$1"; # e.g 10^-2 > p > 10^-3 => 1.00e-2
    }
    # If a decimal format is found (e.g. 0.0023)
    elsif ($pval =~ /(\d+.*)/){
      $sci_pval = $1;
    #$sci_pval = sprintf("%.2e",$pval); # e.g. 0.002 => 2,30e-3
    }
    elsif ($pval =~ /^\w+/) {
      $sci_pval = "NULL";
    }
  }
  else {
    $pval =~ tr/E/e/;
    if ($pval =~ /^(\d+)(e-?\d+)/) {
      $pval="$1.00$2";
    }
    if ($pval =~ /^(\d+\.\d{1})(e-?\d+)/) {
      $pval="$1"."0$2";
    }
    $sci_pval = $pval;
  }
  return $sci_pval;
}

## format ontology accessions
sub _iri2acc{
  my ($self, $iri) = @_;

  my @a = split/\//, $iri;
  my $acc = pop @a;
  $acc =~ s/\_/\:/;

  return $acc;
}

1;
