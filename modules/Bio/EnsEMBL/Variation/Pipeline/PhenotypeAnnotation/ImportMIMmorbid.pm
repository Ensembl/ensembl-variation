=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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


=head1 ImportMIMmorbid

This module imports MIM morbid data. The module fetches the data from
core xrefs.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportMIMmorbid;

use warnings;
use strict;

use File::Path qw(make_path);
use File::stat;
use POSIX qw(strftime);

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

my %source_info;

sub fetch_input {
  my $self = shift;

  my $pipeline_dir = $self->required_param('pipeline_dir');
  my $species      = $self->required_param('species');

  $self->debug($self->param('debug_mode'));

  my $dateStr = strftime("%Y%m%d", localtime);

  %source_info = (source_description => 'Online Mendelian Inheritance in Man (OMIM) database',
                  source_url => 'https://www.omim.org/',
                  object_type => 'Gene',
                  source_status => 'germline',
                  source_name => 'MIM morbid',      #source name in the variation db
                  source_name_short => 'MIMmorbid', #source identifier in the pipeline
                  source_version => $dateStr,
                  data_types => 'phenotype_feature',
                  );

  my $workdir = $pipeline_dir."/".$source_info{source_name_short}."/".$species;
  make_path($workdir) or die "Failed to create $workdir $!\n";
  $self->workdir($workdir);

  open(my $logFH, ">", $workdir."/".'log_import_out_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $errFH, ">", $workdir."/".'log_import_err_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $pipelogFH, ">", $workdir."/".'log_import_debug_pipe_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  $self->logFH($logFH);
  $self->errFH($errFH);
  $self->pipelogFH($pipelogFH);

  # get input file for MIM import
  my $file_mim = "mim_dump.txt";
  if ( -e $workdir."/".$file_mim ){
    print $logFH "Found files (".$workdir."/".$file_mim."), will skip new fetch\n";
    my $fileTime = strftime("%Y%m%d", localtime(stat($workdir."/".$file_mim)->mtime)); #get file date
    print $errFH "WARNING: File $file_mim to be imported has a different date than today!: $fileTime \n" if $fileTime ne $dateStr;
  } else {
    my $st_getdata = qq{
      SELECT g.stable_id, g.seq_region_id,
            g.seq_region_start, g.seq_region_end, g.seq_region_strand,
            x.dbprimary_acc, x.description
      FROM external_db e, xref x, object_xref o, gene g
      WHERE e.external_db_id = x.external_db_id AND
            x.xref_id = o.xref_id AND
            o.ensembl_id = g.gene_id AND e.db_name = 'MIM_MORBID'
    };
    my $sth = $self->core_db_adaptor->dbc->prepare($st_getdata);
    $sth->execute();
    open(OUT, ">$workdir/$file_mim") || die ("ERROR: Unable to write to file $workdir/$file_mim\n");
    print OUT join("\t", @{$sth->{NAME}})."\n";
    while(my @row = $sth->fetchrow_array()) {
      print OUT join("\t", @row)."\n";
    }
    close OUT;
  }

  $self->param('mim_file', $file_mim);
}

sub run {
  my $self = shift;

  my $file_mim = $self->required_param('mim_file');

  # dump and clean pre-existing phenotype features
  $self->dump_phenotypes($source_info{source_name}, 1);

  #get source id
  my $source_id = $self->get_or_add_source(\%source_info);
  $self->print_pipelogFH("$source_info{source_name} source_id is $source_id\n") if ($self->debug);

  # get phenotype data
  my $results = parse_input_file($self->workdir."/".$file_mim);
  $self->print_pipelogFH("Got ".(scalar @{$results->{'phenotypes'}})." new phenotypes \n") if ($self->debug);

  # save phenotypes
  $self->save_phenotypes(\%source_info, $results);

  my %param_source = (source_name => $source_info{source_name_short},
                      type => $source_info{object_type});
  $self->param('output_ids', { source => \%param_source,
                               species => $self->required_param('species'),
                               run_type => $self->required_param('run_type'),
                             });
}

sub write_output {
  my $self = shift;

  $self->print_pipelogFH("Passing $source_info{source_name_short} import (".$self->required_param('species').") for checks (check_phenotypes)\n") if ($self->debug);
  close($self->logFH) if defined $self->logFH ;
  close($self->errFH) if defined $self->errFH ;
  close($self->pipelogFH) if defined $self->pipelogFH ;

  $self->dataflow_output_id($self->param('output_ids'), 2);
}


=head2 parse_input_file

  Arg [1]    : string $infile
               The input file path
  Example    : $results = parse_input_file($infile)
  Description: Parse phenotypes from MIM morbid input file
  Returntype : hashref with results (key 'phenotypes')
  Exceptions : none

=cut

sub parse_input_file {
  my $infile = shift;

  open(IN, ($infile =~ /(z|gz)$/i ? "zcat $infile | " : $infile)) || die ("Could not open $infile for reading\n");

  my @phenotypes;
  # first record is empty
  <IN>;

  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;

    my @content = split(/\t/,$_);

    my $desc = $content[6];
    next if (!$desc || $desc eq '');
    $desc =~ s/^\s+//;
    $desc = (split(';', $desc))[0]; # Only keep the first term

    my $gene = $content[0];
    next if (!$gene || $gene eq '');

    push @phenotypes, {
      'id'                => $gene,
      'description'       => $desc,
      'external_id'       => $content[5],
      'seq_region_id'     => $content[1],
      'seq_region_start'  => $content[2],
      'seq_region_end'    => $content[3],
      'seq_region_strand' => $content[4],
      'type'              => 'Gene',
    };
  }
  close IN;

  my %result = ('phenotypes' => \@phenotypes);
  return (\%result);
}

1;
