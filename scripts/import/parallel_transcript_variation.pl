use strict;
use warnings;

use Getopt::Long;
use Fcntl ':flock';
use DBI;
use DBH;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::VariationFeature; #to get the consequence_types priorities
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);
use Bio::EnsEMBL::Utils::TranscriptAlleles qw(type_variation);  #function to calculate the consequence type of a Variation in a transcript
use ImportUtils qw(debug load create_and_load);

my ($TMP_DIR, $TMP_FILE, $LIMIT,$status_file);

{
  my ($vhost, $vport, $vdbname, $vuser, $vpass,
      $chost, $cport, $cdbname, $cuser, $cpass,
      $limit, $num_processes);

  GetOptions('chost=s'   => \$chost,
             'cuser=s'   => \$cuser,
             'cpass=s'   => \$cpass,
             'cport=i'   => \$cport,
             'cdbname=s' => \$cdbname,
             'vhost=s'   => \$vhost,
             'vuser=s'   => \$vuser,
             'vpass=s'   => \$vpass,
             'vport=i'   => \$vport,
             'vdbname=s' => \$vdbname,
             'tmpdir=s'  => \$ImportUtils::TMP_DIR,
             'tmpfile=s' => \$ImportUtils::TMP_FILE,
             'limit=s'   => \$limit,
	     'num_processes=i' => \$num_processes,
	     'status_file=s' => \$status_file);

  #added default options
  $chost    ||= 'ecs2';
  $cuser    ||= 'ensro';
#  $cport    ||= 3364;

  $vport    ||= 3306;
  $vuser    ||= 'ensadmin';
  
  $num_processes ||= 1;

  $LIMIT = ($limit) ? " $limit " : ''; #the limit will refer to the slices the process must used (1,4-5,7,.....)

  usage('-vdbname argument is required') if(!$vdbname);
  usage('-cdbname argument is required') if(!$cdbname);

  usage('-num_processes must at least be 1') if ($num_processes == 0);
  usage('-status_file argument is required') if (!$status_file);

  my $dbCore = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (-host   => $chost,
     -user   => $cuser,
     -pass   => $cpass,
     -port   => $cport,
     -dbname => $cdbname);

  my $dbVar = DBH->connect
    ("DBI:mysql:host=$vhost;dbname=$vdbname;port=$vport",$vuser, $vpass );
  die("Could not connect to variation database: $!") if(!$dbVar);


  $TMP_DIR  = $ImportUtils::TMP_DIR;
  $TMP_FILE = $ImportUtils::TMP_FILE;

  transcript_variation($dbCore, $dbVar);

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

  my $inc_non_ref = 1;
  my $slices = $sa->fetch_all('toplevel', undef, $inc_non_ref);
  #order the slices
  my @slices_ordered = sort {$a->seq_region_name cmp $b->seq_region_name }  @{$slices};
  my ($offset,$length) = split /,/,$LIMIT; #get the offset and length of the elements we have to get from the slice
  # assumes that variation features have already been pushed to toplevel
  foreach my $slice (splice (@slices_ordered,$offset,$length)) {
    debug("Processing transcript variations for ",
          $slice->seq_region_name(), "\n");
    my $genes = $slice->get_all_Genes();
    # request all variations which lie in the region of a gene
    foreach my $g (@$genes) {

      $sth->execute($slice->get_seq_region_id(),
                    $g->seq_region_start() - $UPSTREAM,
                    $g->seq_region_end()   + $DNSTREAM,
                    $g->seq_region_start() - $MAX_FEATURE_LENGTH - $UPSTREAM);
      my $rows = $sth->fetchall_arrayref();

      foreach my $tr (@{$g->get_all_Transcripts()}) {

        next if(!$tr->translation()); # skip pseudogenes
	my $utr3 = $tr->three_prime_utr();
        my $utr5 = $tr->five_prime_utr();


# compute the effect of the variation on each of the transcripts
        # of the gene
	my $consequence_type;
	my ($start,$end, $strand); #start, end and strand of the variation feature in the slice
        foreach my $row (@$rows) {
	  next if ($row->[4] =~ /LARGE/); #for LARGEINSERTION and LARGEDELETION alleles we don't calculate transcripts
          # put variation in slice coordinates
          $start = $row->[1] - $slice->start() + 1;
          $end = $row->[2] - $slice->start() + 1;
	  $strand = $row->[3];
	  expand(\$row->[4]);#expand the alleles
          my @alleles = split('/',$row->[4]);

          if($strand != $tr->strand()) {
            # flip feature onto same strand as transcript
            for(my $i = 0; $i < @alleles; $i++) {
              reverse_comp(\$alleles[$i]);
            }
            $strand = $tr->strand();
          }
	  shift @alleles; #remove reference allele
	  $consequence_type = Bio::EnsEMBL::Variation::ConsequenceType->new($tr->dbID,$row->[0],$start,$end,$strand,\@alleles);
	  my $consequences;
	  if ($row->[4] !~ /\+/){
	      $consequences = type_variation($tr, $consequence_type);
	  }
          foreach my $ct (@$consequences) {
	      my $type; 
	      $type = join(',',$ct->splice_site,$ct->type) if ($ct->splice_site != ''); #when there is splice site
	      $type = $ct->type if ($ct->splice_site == ''); #when there is no splice site;
            my @arr = ($ct->transcript_id,
                       $ct->variation_feature_id,
                       join("/", @{$ct->aa_alleles||[]}),
                       $ct->aa_start,
                       $ct->aa_end,
                       $ct->cdna_start,
                       $ct->cdna_end,
		       $type);
            @arr = map {($_) ? $_ : '\N'} @arr;
            print FH join("\t", @arr), "\n";
          }
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
    
    debug("Reimporting processed transcript variation");
    my $dbname = $dbVar->dbname(); #get the name of the database to create the file
    my $call = "cat $TMP_DIR/$dbname.transcript_variation*.txt > $TMP_DIR/$TMP_FILE";
    system($call);
    
    unlink(<$TMP_DIR/$dbname.transcript_variation*.txt>);
    
    load($dbVar, qw(transcript_variation
		    transcript_id variation_feature_id peptide_allele_string
		    translation_start translation_end cdna_start cdna_end consequence_type));
    
    
    update_meta_coord($dbCore, $dbVar, 'transcript_variation');
    debug("Preparing to update consequence type in variation feature table");
    #and delete the status file
    my $sth = $dbVar->prepare( "SELECT STRAIGHT_JOIN vf.variation_feature_id, tv.consequence_type ".
			       "FROM variation_feature vf LEFT JOIN transcript_variation tv ON ".
			       "vf.variation_feature_id = tv.variation_feature_id ".
			       "ORDER BY vf.variation_feature_id" );
    $sth->{mysql_use_result} = 1;
    # create a file which contains the var_feat_id and the max consequence type
    $sth->execute();
    my ($variation_feature_id, $consequence_type);
    $sth->bind_columns(\$variation_feature_id,\$consequence_type);
    my $previous_variation_feature_id = 0;

    my %consequence_types = %Bio::EnsEMBL::Variation::ConsequenceType::CONSEQUENCE_TYPES;
    my %splice_sites =  %Bio::EnsEMBL::Variation::ConsequenceType::SPLICE_SITES;

    my $highest_priority_type = 'INTERGENIC'; #by default, has this type
    my $highest_splice_site = '';
    my $type;
    my $splice_site;
    my @types;

    open(FH, ">" . $TMP_DIR . "/" . $TMP_FILE);
    while($sth->fetch()){
	if (($variation_feature_id != $previous_variation_feature_id) && ($previous_variation_feature_id != 0)){
	    print FH join("\t",$previous_variation_feature_id,join(',',$highest_splice_site,$highest_priority_type)), "\n" if ($highest_splice_site ne '');
	    print FH join("\t",$previous_variation_feature_id,$highest_priority_type), "\n" if ($highest_splice_site eq '');
	    $highest_splice_site = '';
	    $highest_priority_type = 'INTERGENIC';
	}
 	$previous_variation_feature_id = $variation_feature_id;

	@types = split(',',$consequence_type) if (defined $consequence_type); #get types, there might be a splice_site and a normal one
	if (@types > 1){
	    $splice_site = shift @types; #and the splice_site
	    $type = shift @types; #get the consequence type

	}
	elsif (@types == 1){
	    $type = shift @types; #get the consequence type
	    $splice_site = '';
	}
	else{
	    $highest_priority_type ||= 'INTERGENIC';
	    $splice_site = '';
	}
	#find the highest consequence type
	if (defined $consequence_type and ($consequence_types{$type} < $consequence_types{$highest_priority_type})){
	    $highest_priority_type = $type;
	}
	#and the highest splice site
	if ($splice_site ne '' and $highest_splice_site eq ''){
	    $highest_splice_site = $splice_site;
	}
	if ($splice_site ne '' and $highest_splice_site ne '' and $splice_sites{$splice_site} < $splice_sites{$highest_splice_site}){
	    $highest_splice_site = $splice_site;
	}
    }
    #and print the last variation

    @types = split(',',$consequence_type) if (defined $consequence_type); #get types, there might be a splice_site and a normal one
    if (@types > 1){
	$type = shift @types; #get the consequence type
	$splice_site = shift @types; #and the splice_site
    }
    elsif (@types == 1){
	$type = shift @types; #get the consequence type
	$splice_site ||= '';
    }
    else{
	$type ||= 'INTERGENIC';
	$splice_site ||= '';
    }
    #find the highest consequence type
    if (defined $consequence_type and ($consequence_types{$type} < $consequence_types{$highest_priority_type})){
	$highest_priority_type = $type;
    }
    #and the highest splice site
    if ($splice_site ne '' and $highest_splice_site eq ''){
	$highest_splice_site = $splice_site;
    }
    if ($splice_site ne '' and $highest_splice_site ne '' and $splice_sites{$splice_site} < $splice_sites{$highest_splice_site}){
	$highest_splice_site = $splice_site;
    }
    print FH join("\t",$previous_variation_feature_id,join(',',$highest_splice_site,$highest_priority_type)), "\n" if ($highest_splice_site ne '');
    print FH join("\t",$previous_variation_feature_id,$highest_priority_type), "\n" if ($highest_splice_site eq '');
    
    close(FH);
    # upload into tmp table
    debug("creating temporary table with consequence type");
    create_and_load($dbVar,"tmp_consequence_type", "variation_feature_id *", "consequence_type");
    # do a multi table update with that one.
    $dbVar->do(qq{UPDATE variation_feature vf, tmp_consequence_type tct 
		      SET vf.consequence_type = tct.consequence_type
		      WHERE vf.variation_feature_id = tct.variation_feature_id
		  });
    $dbVar->do(qq{DROP TABLE tmp_consequence_type});
    unlink("$TMP_DIR/$status_file");
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

  my $sth = $dbVar->prepare
    ('INSERT INTO meta_coord set table_name = ?, coord_system_id = ?');

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
