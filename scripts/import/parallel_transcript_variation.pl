# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

use warnings;

use Getopt::Long;
use Fcntl ':flock';
use DBI;
use DBH;
use FindBin qw( $Bin );

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);
use Bio::EnsEMBL::Variation::ConsequenceType;
use Bio::EnsEMBL::Utils::TranscriptAlleles qw(type_variation);  #function to calculate the consequence type of a Variation in a transcript
use ImportUtils qw(dumpSQL debug load create_and_load);

my ($TMP_DIR, $TMP_FILE, $LIMIT,$status_file);


my ($species,$limit,$num_processes,$registry_file);
my $run_as_last;

GetOptions('species=s'   => \$species,
           'tmpdir=s'  => \$ImportUtils::TMP_DIR,
           'tmpfile=s' => \$ImportUtils::TMP_FILE,
           'limit=s'   => \$limit,
           'num_processes=i' => \$num_processes,
           'status_file=s' => \$status_file,
           'registry_file=s' => \$registry_file,
           'run_as_last' => \$run_as_last);

$num_processes ||= 1;

$LIMIT = ($limit) ? " $limit " : ''; #the limit will refer to the slices the process must used (1,4-5,7,.....)

usage('-num_processes must at least be 1') if ($num_processes == 0);
usage('-status_file argument is required') if (!$status_file);

$registry_file ||= $Bin . "/ensembl.registry";

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $fdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'funcgen') if $species =~ /hum|hom|mouse|mus/i;
my $species_name = $cdba->species();#used by funcg table

my $dbVar = $vdba->dbc;
my $dbCore = $cdba;
my $dbFunc = $fdba;

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

if ($run_as_last)
{
  last_process($dbCore,$dbVar, $run_as_last);
}
else
{
  transcript_variation($dbCore, $dbVar, $dbFunc);

  open STATUS, ">>$TMP_DIR/$status_file"
    or throw("Could not open tmp file: $TMP_DIR/$status_file\n"); 
  flock(STATUS,LOCK_EX);
  seek(STATUS, 0, 2); #move to the end of the file
  print STATUS "process finished\n";
  flock(STATUS,LOCK_UN);
  close STATUS;

  #check if it is the last process
  my $processes = `cat $TMP_DIR/$status_file | wc -l`;
  if ($processes == $num_processes){
    #if is the last process, finish it
    debug("Last process: ready to import data");
    last_process($dbCore,$dbVar);
  }
}

#
# Loads the transcript variation table.  Retrieves every transcript in the
# the database and types all of the snps in the vicinity of the transcript.
# The amino acid changes for coding snps is also set.
#
#

sub transcript_variation {
  my $dbCore = shift;
  my $dbVar  = shift;
  my $dbFunc = shift;

  my $UPSTREAM = 5000;
  my $DNSTREAM = 5000;
  my $MAX_FEATURE_LENGTH  = 1000;
 
  my $sth = $dbVar->prepare
    (qq{SELECT vf.variation_feature_id, vf.seq_region_start, vf.seq_region_end,
               vf.seq_region_strand, vf.allele_string
        FROM   variation_feature vf
        WHERE  vf.seq_region_id = ?
        AND    vf.seq_region_end >= ?
        AND    vf.seq_region_start <= ?
        AND    vf.seq_region_start >= ?
	});

  my $sa = $dbCore->get_SliceAdaptor();

  my $dbname = $dbVar->dbname(); #get the name of the database to create the file
  my $host = `hostname`;
  chop $host;
  open FH, ">$TMP_DIR/$dbname.transcript_variation_$host\:$$\.txt";
  #open FH, ">/tmp/$dbname.transcript_variation_$host\_$$\.txt" or die "Could not open file: $!\n";
  my $inc_non_ref = 1;

  my $slices = $sa->fetch_all('toplevel', undef, $inc_non_ref);
  #order the slices
  my @slices_ordered = sort {$a->seq_region_name cmp $b->seq_region_name }  @{$slices};
  #
  #my @slice_names = map {$_->seq_region_name}, @slices_ordered;
  #my $name_string = join " ", @slice_names;
  #debug("Slices are\n$name_string\n");
  
  my ($offset,$length) = split /,/,$LIMIT; #get the offset and length of the elements we have to get from the slice
  # assumes that variation features have already been pushed to toplevel


  my %done;
  foreach my $slice (splice (@slices_ordered,$offset,$length)) {
    #next unless $slice->dbID =~ /27509|27515|27524|27527/;
    debug("Processing transcript variations for ",
	  $slice->seq_region_name(), "\n");
    #if it is human, mouse $dbFunc is defined, therefore go through below to calculate regulatory feature
    if ($dbFunc) {
      my $efa = $dbFunc->get_ExternalFeatureAdaptor();
      my $rfa = $dbFunc->get_RegulatoryFeatureAdaptor();
      my $gene_adaptor = $dbCore->get_GeneAdaptor();
      my $transcript_adaptor = $dbCore->get_TranscriptAdaptor();
      #first consider the regulatory region that overlap with SNPs
      my (@rf);
      foreach my $f  (@{$efa->fetch_all_by_Slice($slice)}) {
        if ($f->feature_set->name =~ /miRanda/ or $f->feature_set->name =~ /VISTA\s+enhancer\s+set/i or $f->feature_set->name =~ /cisRED\s+motifs/i) {
	  #print "feature_name is ",$f->feature_set->name,"\n";
          push @rf, $f;
	}
      }

      #print "slice name is ",$slice->name,"\n";
      my @rf_fg;# = @{$rfa->fetch_all_by_Slice($slice)};
      push @rf, @rf_fg if @rf_fg >0;
      foreach my $rf (@rf) { 
	# request all variations which lie in the region of a regulate feature
	$sth->execute($slice->get_seq_region_id(),
		      $rf->seq_region_start(),
		      $rf->seq_region_end(),
		      $rf->seq_region_start());
	my $rows = $sth->fetchall_arrayref();
	my ($start,$end, $rf_start, $rf_end); #start, end of the variation feature in the slice

	foreach my $row (@$rows) {
	  $start = $row->[1];
	  $end = $row->[2];
	  $rf_start = $rf->seq_region_start();
	  $rf_end = $rf->seq_region_end();

	  if ($end >= $rf_start and $start <= $rf_end) {
	    #print "start is $start and rf_start is $rf_start and rf_end is $rf_end\n";
	    foreach my $dbEntry (@{$rf->get_all_DBEntries("$species_name\_core_Gene")}) {
	      #print "gene is ",$dbEntry->primary_id,"\n";
              my $g = $gene_adaptor->fetch_by_stable_id($dbEntry->primary_id); #get the gene for the stable_id
	      if(!defined $g) {
		#some of the genes do not seem to be in the core database due to core gene update
		#print "Gene : ",$dbEntry->primary_id," is not in coredb\n";
		next;
	      }
	      foreach my $tr (@{$g->get_all_Transcripts()}) {
		my @arr = ($tr->stable_id,$row->[0],'\N','\N','\N','\N','\N','\N','\N','REGULATORY_REGION');
		if (! $done{$row->[0]}{$tr->dbID}) {#get rid of duplicated transcript entries
		  print FH join("\t", @arr), "\n";
		  $done{$row->[0]}{$tr->dbID}=1;
		}
	      }
	    }
	    foreach my $dbEntry (@{$rf->get_all_DBEntries("$species_name\_core_Transcript")}) {	
	      my $tr = $transcript_adaptor->fetch_by_stable_id($dbEntry->primary_id); #get transcript for stable_id
	      next if(!defined $tr); #some of the transcripts do not seem to be in the core database

	      my @arr = ($tr->stable_id,$row->[0],'\N','\N','\N','\N','\N','\N','\N','REGULATORY_REGION');
	      if (! $done{$row->[0]}{$tr->dbID}) {#get rid of duplicated transcript entries
		print FH join("\t", @arr), "\n";
		$done{$row->[0]}{$tr->dbID}=1;
	      }
	    }
	  }
	}
      }
    }

    # then compute the effect of the variation on each of the transcripts
    # of the gene
    #next if $slice->seq_region_name() ne '17'; 
    my $load_transcripts = 1;

    #we want all genes, not just genes with protein_coding biotype
    my $transcripts = $slice->get_all_Transcripts();

    # request all variations which lie in the region of a gene
    # in future, could use get_all_Transcript instead
    foreach my $tr (@$transcripts) {
      #debug("Time to do transcript ",$tr->stable_id," : ",scalar(localtime(time)));


      #debug( "seq_region_start is ",$tr->seq_region_start()," seq_region_id is ",
      #$slice->get_seq_region_id()," seq_region_end is ",$tr->seq_region_end());
      $sth->execute($slice->get_seq_region_id(),
                    $tr->seq_region_start() - $UPSTREAM,
                    $tr->seq_region_end()   + $DNSTREAM,
                    $tr->seq_region_start() - $MAX_FEATURE_LENGTH - $UPSTREAM);
      my $rows = $sth->fetchall_arrayref();


      my $tr_name = $tr->display_id;

      my ($start,$end, $strand); #start, end and strand of the variation feature in the slice

      foreach my $row (@$rows) {
		# HGMD MUTATION
		if($row->[4] eq 'HGMD_MUTATION') {
		  print FH $tr->stable_id, "\t", $row->[0], (('\N'."\t") x 8), 'HGMD_MUTATION', "\n";
		  next;
		}
		
	next if ($row->[4] =~ /LARGE|INS|DEL|CNV|MUTATION/);#for LARGEINSERTION and LARGEDELETION alleles we don't calculate transcripts
	# put variation in slice coordinates
	$start = $row->[1] - $slice->start() + 1;
	$end = $row->[2] - $slice->start() + 1;
	$strand = $row->[3];
	expand(\$row->[4]);	#expand the alleles
	my @alleles = split('/',$row->[4]);

	if ($strand != $tr->strand()) {
	  # flip feature onto same strand as transcript
	  for (my $i = 0; $i < @alleles; $i++) {
	    reverse_comp(\$alleles[$i]) unless $alleles[$i] =~ /INS|DEL/;
	  }
	  $strand = $tr->strand();
	}
	shift @alleles;		#remove reference allele
	
	my $consequence_type;

	$consequence_type = Bio::EnsEMBL::Variation::ConsequenceType->new($tr->dbID,$row->[0],$start,$end,$strand,\@alleles);

	# compute the effect of the variation on each of the transcripts
	# of the gene


	my ($consequences);
	if ($row->[4] !~ /\+/) {
	  #print "Getting consequences for ", (join " ", @$row), " in transcript ", $tr->stable_id, " ", $tr->dbID;
          #print "\n";
	  $consequences = type_variation($tr, "", $consequence_type);
	}
	foreach my $ct (@$consequences) {
	  my $final_ct = join(",",@{$ct->type});
	  my @arr = (#$ct->transcript_id,
		     $tr->stable_id,
		     $ct->variation_feature_id,
		     join("/", @{$ct->aa_alleles||[]}),
		     $ct->aa_start,
		     $ct->aa_end,
		     $ct->cdna_start,
		     $ct->cdna_end,
		     $ct->cds_start,
		     $ct->cds_end,
		     $final_ct); 
	  @arr = map {($_) ? ($_) = $_ =~ /^(.*)[,]*$/ : '\N'} @arr;
	  print FH join("\t", @arr), "\n";
	}
      }
    }
  }

  close FH;
  return;
}

#will know if it is the last process running counting the lines of the status_file.If so, load all data
sub last_process{
    my $dbCore = shift;
    my $dbVar = shift;
    my $run_as_last = shift;



    my $dbname = $dbVar->dbname(); #get the name of the database to create the file
    if (!$run_as_last)
    {
      my $call = "cat $TMP_DIR/$dbname.transcript_variation*.txt > $TMP_DIR/$TMP_FILE";
      system($call);

      unlink(<$TMP_DIR/$dbname.transcript_variation*.txt>);
    }

    debug("Deleting existing transcript variation");

    $dbVar->do("DELETE FROM transcript_variation");

    #load transcript_variation data
    debug("Reimporting processed transcript variation");
    load($dbVar, qw(transcript_variation
 		    transcript_stable_id variation_feature_id peptide_allele_string
 		    translation_start translation_end cdna_start cdna_end cds_start cds_end consequence_type));



    debug("Preparing to update consequence type in variation feature table");
    # create a file which contains the var_feat_id and the max consequence type
    dumpSQL($dbVar, qq{SELECT variation_feature_id, consequence_type FROM transcript_variation});
 
    # sort the file
    system("sort $TMP_DIR/$TMP_FILE |uniq > $TMP_DIR/$TMP_FILE.su");
    open IN, "$TMP_DIR/$TMP_FILE\.su" or die "can't open tmp_file $!";
    open(FH, ">" . $TMP_DIR . "/" . $TMP_FILE);

    my $previous_variation_feature_id = 0;
    my %consequence_types = %Bio::EnsEMBL::Variation::ConsequenceType::CONSEQUENCE_TYPES;
    #my %splice_sites =  %Bio::EnsEMBL::Variation::ConsequenceType::SPLICE_SITES;
    my ($variation_feature_id, $consequence_type);
    my $highest_priority_type = 'INTERGENIC'; #by default, has this type
    my $highest_splice_site = '';
    my $regulatory_region;
    my $type;
    my $splice_site;
    my @types;
    my $final_type;
    
    while (<IN>) {
      ($variation_feature_id, $consequence_type) = split;
      if (($variation_feature_id != $previous_variation_feature_id) && ($previous_variation_feature_id != 0)){
	$final_type = "$highest_splice_site," if $highest_splice_site;
	$final_type .= "$regulatory_region," if $regulatory_region;
	$final_type .= "$highest_priority_type";
	print FH join("\t",$previous_variation_feature_id,$final_type),"\n";
      
	$regulatory_region = '';
	$splice_site = '';
	$highest_splice_site = '';
	$highest_priority_type = 'INTERGENIC';
	$final_type = '';
        $type = '';
      }
      $previous_variation_feature_id = $variation_feature_id;

      if (defined $consequence_type) {#get types, there might be a splice_site and a normal one and regulatory_region
	@types = split(',',$consequence_type);
	foreach my $ct (@types) {
	  if ($ct =~ /regulatory/i) {
	    $regulatory_region = $ct;
	  }
	  elsif ($ct =~ /splice/i) {
	    $splice_site = $ct;
	  }
	  else {
	    $type = $ct;
	  }
	}

	#find the highest consequence type
	if ($type and $consequence_types{$type}  and $consequence_types{$type} < $consequence_types{$highest_priority_type}){
	  $highest_priority_type = $type;
	}
	#and the highest splice site
	if ($splice_site and $highest_splice_site eq ''){
	  $highest_splice_site = $splice_site;
	}
	if ($splice_site and $highest_splice_site and $consequence_types{$splice_site} < $consequence_types{$highest_splice_site}){
	  $highest_splice_site = $splice_site;
	}
      }
  }
    #and print the last variation

    if (defined $consequence_type) {
      @types = split(',',$consequence_type);
      #print "type is @types\n";
      foreach my $ct (@types) {
	if ($ct =~ /regulatory/i) {
	  $regulatory_region = $ct;
	}
	elsif ($ct =~ /splice/i) {
	  $splice_site = $ct;
	}
	else {
	  $type = $ct;
	}
      }

      #find the highest consequence type
      if ($type and $consequence_types{$type} and $consequence_types{$type} < $consequence_types{$highest_priority_type}){
	$highest_priority_type = $type;
      }
      #and the highest splice site
      if (defined $splice_site and $splice_site ne '' and $highest_splice_site eq ''){
	$highest_splice_site = $splice_site;
      }
      if (defined $splice_site and $splice_site ne '' and $highest_splice_site ne '' and $consequence_types{$splice_site} < $consequence_types{$highest_splice_site}){
	$highest_splice_site = $splice_site;
      }
    }

    $final_type = "$highest_splice_site," if defined $highest_splice_site;
    $final_type .= "$regulatory_region," if defined $regulatory_region;
    $final_type .= "$highest_priority_type";
    print FH join("\t",$previous_variation_feature_id,$final_type),"\n";

    close(FH);
    

    # upload into tmp table
    debug("creating temporary table with consequence type");
    create_and_load($dbVar,"tmp_consequence_type", "variation_feature_id *", "consequence_type");
    # do a multi table update with that one.
    $dbVar->do(qq{UPDATE variation_feature vf, tmp_consequence_type tct 
		      SET vf.consequence_type = tct.consequence_type
		      WHERE vf.variation_feature_id = tct.variation_feature_id
                      AND vf.consequence_type != 'NO_CONSEQUENCE'
                  });
    $dbVar->do(qq{DROP TABLE tmp_consequence_type});
    
    #set consequence_type to INTERGENIC if previous transcript has been removed from current release
    $dbVar->do(qq{UPDATE variation_feature vf 
	          SET vf.consequence_type = 'INTERGENIC'
		  WHERE NOT EXISTS
		  (SELECT * FROM transcript_variation tv
		  WHERE vf.variation_feature_id = tv.variation_feature_id)
		  AND vf.consequence_type != 'INTERGENIC'
		});
    
    unlink("$TMP_DIR/$status_file");
    unlink("$TMP_DIR/$TMP_FILE.su");
    update_meta_coord($dbCore, $dbVar, 'transcript_variation');
}

#
# updates the meta coord table
#
sub update_meta_coord {
  my $dbCore = shift;
  my $dbVar  = shift;
  my $table_name = shift;
  my $csname = shift || 'chromosome';

  my $csa = $dbCore->get_CoordSystemAdaptor();

  my $cs = $csa->fetch_by_name($csname);

# Added IGNORE to statement
  my $sth = $dbVar->prepare
    ('INSERT IGNORE INTO meta_coord set table_name = ?, coord_system_id = ?');

  $sth->execute($table_name, $cs->dbID());

  $sth->finish();

  return;
}

sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: perl parallel_transcript_variation.pl <options>

options:
    -chost <hostname>    hostname of core Ensembl MySQL database (default = ecs2)
    -cuser <user>        username of core Ensembl MySQL database (default = ensro)
    -cpass <pass>        password of core Ensembl MySQL database
    -cport <port>        TCP port of core Ensembl MySQL database (default = 3364)
    -cdbname <dbname>    dbname of core Ensembl MySQL database
    -vhost <hostname>    hostname of variation MySQL database to write to
    -vuser <user>        username of variation MySQL database to write to (default = ensadmin)
    -vpass <pass>        password of variation MySQL database to write to
    -vport <port>        TCP port of variation MySQL database to write to (default = 3306)
    -vdbname <dbname>    dbname of variation MySQL database to write to
    -limit <number>      limit the number of rows for testing
    -tmpdir <dir>        temp directory to use (with lots of space!)
    -tmpfile <filename>   name of temp file to use
    -num_processes <number> number of processes that are running (default = 1)
    -status_file <filename> name of a temp file where all the processes write when they finish
EOF

  die("\n$msg\n\n");
}
