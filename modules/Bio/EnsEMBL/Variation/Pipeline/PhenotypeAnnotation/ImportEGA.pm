=head1 LICENSE
#TODO: check correct licence

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2018] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportEGA;

use warnings;
use strict;

use File::Path qw(make_path);
use LWP::Simple;
use POSIX 'strftime';
use DBI qw(:sql_types);
use Data::Dumper; #TODO: remove if not needed
use Bio::EnsEMBL::Registry;
#use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants qw(GWAS NONE);
#use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation qw($variation_dba);

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

my %source_info;
my $workdir;
my $core_dba;
my $variation_dba;
my $phenotype_dba;

my $pubmed_prefix = 'PMID:';
my $dbh;

my $debug;

sub fetch_input {
    my $self = shift;

    my $pipeline_dir = $self->required_param('pipeline_dir');
    my $species      = $self->required_param('species');
    my $conf_file    = $self->required_param('ega_database_conf');


    $core_dba    = $self->get_species_adaptor('core');
    $variation_dba  = $self->get_species_adaptor('variation'); #TODO: why did the init -> Base class -> this not work?
    $phenotype_dba  = $variation_dba->get_PhenotypeAdaptor; 

    $debug        = $self->param('debug_mode');

    #TODO: 
    # original call: perl ensembl-variation/scripts/import/import_phenotype_data.pl -host ${host} -user ${user} -pass ${pass}  -dbname ${dbname} \
    #-source nhgri -infile ${datadir}/gwascatalog.txt -verbose -version 20121031

    my $dateStr = strftime "%Y%m", localtime;

    %source_info = (source_description => 'Variants imported from the European Genome-phenome Archive with phenotype association',
                    source_url => 'https://www.ebi.ac.uk/ega/',
                    object_type => 'Variation',

                    source_name => 'EGA', #TODO: figure out where this is used
                    source_version => $dateStr, # it is current month
                    source_status => 'germline',
                    source => 'ega',

                    #  object_type => 'QTL', #default type, import code will switch to Gene for the Gene-type ones
                  #  source_mapped_attrib_type => 'Rat Genome Database', #for ontology mapping (attr_type_id 509) entry in phenotype_ontology_accession (attr_id 588)
                    
                  #  threshold => $self->required_param('threshold_qtl'),
                    );

    $workdir = $pipeline_dir."/".$source_info{source_name}."/".$species;
    make_path($workdir);
    my $file_ega = "ega.studies.csv";

    #parse database connection details:
    my %database_conf;
    open(CONF,'<',$conf_file) or die ("Could not open $conf_file for reading");
    while (<CONF>) {
        chomp;                  # no newline
        s/#.*//;                # no comments
        s/^\s+//;               # no leading white
        s/\s+$//;               # no trailing white
        next unless length;     # anything left?
        my ($var, $value) = split(/\s*=\s*/, $_, 2);
        $database_conf{$var} = $value;
    }

    #get input file EGA:
    local (*STDOUT, *STDERR);
    open STDOUT, ">", $workdir."/".'log_import_out_'.$file_ega; #TODO: what is best error/out log naming convention?
    open STDERR, ">", $workdir."/".'log_import_err_'.$file_ega;
    my $dsn = "dbi:mysql:$database_conf{DATABASE}:$database_conf{HOST}:$database_conf{PORT}";
    $dbh = DBI->connect($dsn,$database_conf{USER},$database_conf{PASS});
    get_ega_file($workdir."/".$file_ega) unless -e $workdir."/".$file_ega;

    print "Found files (".$workdir."/".$file_ega."), will skip new fetch\n" if -e $workdir."/".$file_ega;
    $self->param('ega_file', $file_ega);
}

sub run {
  my $self = shift;

  my $file_ega = $self->required_param('ega_file');

  local (*STDOUT, *STDERR);
  open STDOUT, ">>", $workdir."/".'log_import_out_'.$file_ega; #TODO: what is best error/out log naming convention?
  open STDERR, ">>", $workdir."/".'log_import_err_'.$file_ega;

  #get source id
  my $source_id = $self->get_or_add_source(\%source_info,$variation_dba);
  print STDOUT "$source_info{source} source_id is $source_id\n" if ($debug);

  # get phenotype data + save it (all in one method)
  my $results = parse_ega($workdir."/".$file_ega, $source_id);
  print "Got ".(scalar @{$results->{'studies'}})." new studies \n" if $debug ;

  my %param_source = (source_name => $source_info{source_name},
                      type => ['Variation']);
  $self->param('output_ids', { source => \%param_source,
                               species => $self->required_param('species')
                             });
}

sub write_output {
  my $self = shift;

  $self->dataflow_output_id($self->param('output_ids'), 1);
  print "Passing EGA import for checks\n" if $self->param('debug_mode');
}


# EGA specific data fetch methods
sub get_ega_file {
  my $outfile = shift;

  my $studies = get_all_study_stable_ids();
  open (FILE, ">".$outfile);

  foreach my $stable_id (sort {$a cmp $b } @$studies) {
    my $study_id = get_study_id($stable_id);
    my $pmid = get_publication($study_id);
    if (!$pmid) {
      warn "Cannot find a publication for $stable_id\n";
      next;
    }
    #TODO: can I replace it with source_url:http://www.ebi.ac.uk/ega/ OR does it matter the https/ http?
  #  my $url = "https://www.ebi.ac.uk/ega/studies/".$stable_id;

    print FILE "$stable_id,$pmid,$source_info{source_url}studies/$stable_id\n";
  }

  close (FILE);
}

sub get_all_study_stable_ids {
  my @studies;

  my $sth = $dbh->prepare(qq[SELECT stable_id FROM study]);
  $sth->execute();

  while (my $id = $sth->fetchrow() ) {
    push @studies, $id;
  }

  $sth->finish();

  return \@studies;
}

sub get_study_id {
  my ($stable_id) = @_;

  my $sth = $dbh->prepare(qq[SELECT study_id FROM study WHERE stable_id = ?]);
  $sth->execute($stable_id);

  my $study_id = $sth->fetchrow();
  $sth->finish();

  return $study_id;
}

sub get_publication {
    my ($study_id) = @_;

    my $sth = $dbh->prepare(qq[SELECT publication_id FROM study_publication WHERE study_id = ?]);
    $sth->execute($study_id);

    my $publication_id = $sth->fetchrow();
    $sth->finish();

    my $pmid;
    if ($publication_id) {
      $sth = $dbh->prepare(qq[SELECT pmid from publication WHERE publication_id = ?]);
      $sth->execute($publication_id);

      $pmid = $sth->fetchrow();
      $sth->finish();
    }

    return $pmid;
}



# EGA specific phenotype parsing method
sub parse_ega {
  my $infile = shift;
  my $source_id = shift;
  
  my $study_check_stmt = qq{
    SELECT
      study_id
    FROM
      study
    WHERE
      name=? AND source_id=$source_id
    LIMIT 1
  };
  my $nhgri_check_stmt = qq{
    SELECT
      st.study_id,st.study_type
    FROM
      study st, source s
    WHERE
      external_reference=? AND s.name like '%nhgri%'
      AND s.source_id=st.source_id
    LIMIT 1
  };
  my $study_ins_stmt = qq{
    INSERT INTO
      study (
      name,
      source_id,
      external_reference,
      url,
      study_type
      )
    VALUES (
      ?,
      $source_id,
      ?,
      ?,
      ?
    )
  };
  # NHGRI and EGA associated studies
  my $asso_study_check_stmt = qq{
    SELECT
      study1_id
    FROM
      associate_study
    WHERE
      study1_id = ? AND study2_id = ?
    LIMIT 1
  };
  my $asso_study_ins_stmt = qq{
    INSERT INTO
      associate_study (study1_id,study2_id)
    VALUES (?,?)
  };
  
  my $nhgri_check_sth = $variation_dba->dbc->prepare($nhgri_check_stmt);
  my $study_check_sth = $variation_dba->dbc->prepare($study_check_stmt);
  my $study_ins_sth   = $variation_dba->dbc->prepare($study_ins_stmt);
  my $asso_study_check_sth = $variation_dba->dbc->prepare($asso_study_check_stmt);
  my $asso_study_ins_sth   = $variation_dba->dbc->prepare($asso_study_ins_stmt);
  
  # Open the input file for reading
  open(IN,'<',$infile) or die ("Could not open $infile for reading");
  
  # Read through the file and parse out the desired fields
  my @new_studies;
  while (<IN>) {
    chomp $_;
    my @attributes = split(",",$_);
    next if ($attributes[1] eq '');
    my $name = $attributes[0];
    my $pubmed = $pubmed_prefix.$attributes[1];
    my $url = $attributes[2];
    
    # NHGRI study
    my $nhgri_study_id;
    my $study_type;
    $nhgri_check_sth->bind_param(1,$pubmed,SQL_VARCHAR);
    $nhgri_check_sth->execute();
    $nhgri_check_sth->bind_columns(\$nhgri_study_id,\$study_type);
    $nhgri_check_sth->fetch();
    
    if (!defined($nhgri_study_id)) {
      warn "No NHGRI-EBI study found for the EGA $name | $pubmed !\n";
      next;
    }
    
    # EGA study
    my $study_id;
    $study_check_sth->bind_param(1,$name,SQL_VARCHAR);
    $study_check_sth->execute();
    $study_check_sth->bind_columns(\$study_id);
    $study_check_sth->fetch();
    if (!defined($study_id)) {
      $study_ins_sth->bind_param(1,$name,SQL_VARCHAR);
      $study_ins_sth->bind_param(2,$pubmed,SQL_VARCHAR);
      $study_ins_sth->bind_param(3,$url,SQL_VARCHAR);
      $study_ins_sth->bind_param(4,$study_type,SQL_VARCHAR);
      $study_ins_sth->execute();
      $study_id = $variation_dba->dbc->db_handle->{'mysql_insertid'};
      push (@new_studies, $name)
    }
    
    my $is_associated;
    $asso_study_check_sth->bind_param(1,$nhgri_study_id,SQL_INTEGER);
    $asso_study_check_sth->bind_param(2,$study_id,SQL_INTEGER);
    $asso_study_check_sth->execute();
    $asso_study_check_sth->bind_columns(\$is_associated);
    $asso_study_check_sth->fetch();
    
    if (!defined($is_associated)) {
      $asso_study_ins_sth->bind_param(1,$nhgri_study_id,SQL_INTEGER);
      $asso_study_ins_sth->bind_param(2,$study_id,SQL_INTEGER);
      $asso_study_ins_sth->execute();
    }
  }
  close(IN);

  my %result = ('studies' => \@new_studies);
  return \%result;
}

1;
