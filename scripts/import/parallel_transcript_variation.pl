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
use ImportUtils qw(debug load create_and_load);
use Data::Dumper;

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

        foreach my $row (@$rows) {
          my %var;
          $var{'vf_id'}  = $row->[0];
          # put variation in slice coordinates
          $var{'start'}  = $row->[1] - $slice->start() + 1;
          $var{'end'}    = $row->[2] - $slice->start() + 1;
          $var{'strand'} = $row->[3];
          $var{'tr_id'}  = $tr->dbID();
	  next if ($row->[4] =~ /LARGE/); #for LARGEINSERTION and LARGEDELETION alleles we don't calculate transcripts
	  expand(\$row->[4]);#expand the alleles
          my @alleles = split('/',$row->[4]);

          if($var{'strand'} != $tr->strand()) {
            # flip feature onto same strand as transcript
            for(my $i = 0; $i < @alleles; $i++) {
              reverse_comp(\$alleles[$i]);
            }
            $var{'strand'} = $tr->strand();
          }
          $var{'alleles'} = \@alleles;
	  my $vars;
	  if ($row->[4] !~ /\+/){
	      $vars = type_variation($tr, \%var);
	  }
          foreach my $v (@$vars) {
            my @arr = ($v->{'tr_id'},
                       $v->{'vf_id'},
                       join("/", @{$v->{'aa_alleles'}||[]}),
                       $v->{'aa_start'},
                       $v->{'aa_end'},
                       $v->{'cdna_start'},
                       $v->{'cdna_end'},
                       $v->{'type'});
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

#
# Classifies a variation which is in the vicinity of a transcript
#
sub type_variation {
  my $tr  = shift;
  my $var = shift;

  my $tm = $tr->get_TranscriptMapper();

  my @coords = $tm->genomic2cdna($var->{'start'},
                                 $var->{'end'},
                                 $var->{'strand'});

  # Handle simple cases where the variation is not split into parts.
  # Call method recursively with component parts in complicated case.
  # E.g. a single multi-base variation may be both intronic and coding


  if(@coords > 1) {
    my @out;

    foreach my $c (@coords) {
      my %new_var = %{$var};
      $new_var{'end'} = $var->{'start'} + $c->length() - 1;
      $var->{'start'} = $new_var{'end'} + 1;
      push @out, @{type_variation($tr, \%new_var)};
    }

    return \@out;
  }

  my $c = $coords[0];

  if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {

    # check if the variation is completely outside the transcript:

    if($var->{'end'} < $tr->start()) {
      $var->{'type'} = ($tr->strand() == 1) ? 'UPSTREAM' : 'DOWNSTREAM';
      return [$var];
    }
    if($var->{'start'} > $tr->end()) {
      $var->{'type'} = ($tr->strand() == 1) ? 'DOWNSTREAM' : 'UPSTREAM';
      return [$var];
    }

    # variation must be intronic since mapped to cdna gap, but is within
    # transcript
    $var->{'type'} = 'INTRONIC';
    return [$var];
  }

  $var->{'cdna_start'} = $c->start();
  $var->{'cdna_end'}   = $c->end();

  @coords = $tm->genomic2cds($var->{'start'}, $var->{'end'},$var->{'strand'});

  if(@coords > 1) {
    my @out;

    foreach my $c (@coords) {
      my %new_var = %{$var};
      $new_var{'end'} = $var->{'start'} + $c->length() - 1;
      $var->{'start'} = $new_var{'end'} + 1;
      push @out, @{type_variation($tr, \%new_var)};
    }
    return \@out;
  }

  $c = $coords[0];

  if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {
    # mapped successfully to CDNA but not to CDS, must be UTR

    if($var->{'end'} < $tr->coding_region_start()) {
      $var->{'type'} = ($tr->strand() == 1) ? '5PRIME_UTR' : '3PRIME_UTR';
    }
    elsif($var->{'start'} > $tr->coding_region_end()) {
      $var->{'type'} = ($tr->strand() == 1) ? '3PRIME_UTR' : '5PRIME_UTR';
    }
    else {
      throw('Unexpected: CDNA variation which is not in CDS is not in UTR');
    }
    return [$var];
  }

  $var->{'cds_start'} = $c->start();
  $var->{'cds_end'}   = $c->end();

  @coords = $tm->genomic2pep($var->{'start'}, $var->{'end'}, $var->{'strand'});

  if(@coords != 1 || $coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
    throw("Unexpected: Could map to CDS but not to peptide coordinates.");
  }

  $c = $coords[0];

  $var->{'aa_start'} = $c->start();
  $var->{'aa_end'}   = $c->end();

  apply_aa_change($tr, $var);

  return [$var];
}



#
# Determines the effect of a coding variation on the peptide sequence
#

sub apply_aa_change {
  my $tr = shift;
  my $var = shift;

  my $peptide = $tr->translate->seq();

  my ($attrib) = @{$tr->slice()->get_all_Attributes('codon_table')}; #for mithocondrial dna it is necessary to change the table

  my $codon_table;
  $codon_table = $attrib->value() if($attrib);
  $codon_table ||= 1; # default vertebrate codon table 

  my $len = $var->{'aa_end'} - $var->{'aa_start'} + 1;
  my $old_aa = substr($peptide, $var->{'aa_start'} -1 , $len);

  my $codon_cds_start = $var->{'aa_start'} * 3 - 2;
  my $codon_cds_end   = $var->{'aa_end'}   * 3;
  my $codon_len = $codon_cds_end - $codon_cds_start + 1;

  my @alleles = @{$var->{'alleles'}};
  
  shift(@alleles); # ignore reference allele

  my $var_len = $var->{'cds_end'} - $var->{'cds_start'} + 1;

  my @aa_alleles = ($old_aa);

  foreach my $a (@alleles) {
    $a =~ s/\-//;
    my $cds = $tr->translateable_seq();
    
    if($var_len != length($a)) {
      if(abs(length($a) - $var_len) % 3) {
        # frameshifting variation, do not set peptide_allele string
        # since too complicated and could be very long
        $var->{'type'} = 'FRAMESHIFT_CODING';
        return;
      }

      if($var_len == 0) { # insertion
        $aa_alleles[0] = '-';
        $old_aa    = '-';
      }
    }

    my $new_aa;

    if(length($a)) {
      substr($cds, $var->{'cds_start'}-1, $var_len) = $a;
      my $codon_str = substr($cds, $codon_cds_start-1, $codon_len);

      my $codon_seq = Bio::Seq->new(-seq      => $codon_str,
                                    -moltype  => 'dna',
                                    -alphabet => 'dna');


      $new_aa = $codon_seq->translate(undef,undef,undef,$codon_table)->seq();
    } else {
      $new_aa = '-'; # deletion
    }

    if(uc($new_aa) ne uc($old_aa)) {
      push @aa_alleles, $new_aa;
    }
  }

  if(@aa_alleles > 1) {
    $var->{'type'} = 'NON_SYNONYMOUS_CODING';
  } else {
    $var->{'type'} = 'SYNONYMOUS_CODING';
  }

  $var->{'aa_alleles'} = \@aa_alleles;
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
    my %consequence_types = %Bio::EnsEMBL::Variation::VariationFeature::CONSEQUENCE_TYPES;
    my @consequence_types_ordered = sort {$consequence_types{$a} <=> $consequence_types{$b}} keys %consequence_types;
    my $var_consequence_type = {}; #will contain all the consequence types for the variation_feature
    my $highest_priority;
    open(FH, ">" . $TMP_DIR . "/" . $TMP_FILE);
    while($sth->fetch()){
	if (!$consequence_type){$consequence_type = 'INTERGENIC'}
	#when we have seen all consequence_types for the variation_feature, find the highest
	if (($variation_feature_id != $previous_variation_feature_id) && ($previous_variation_feature_id != 0)){
	    foreach my $type (@consequence_types_ordered){
		if ($var_consequence_type->{$type}){
		    $highest_priority = $type;
		    last;
		}
	    }
	    print FH join("\t",$previous_variation_feature_id,$highest_priority), "\n";
	    $var_consequence_type= {}; #initialize the consequence type for the next variation_feature
	}
	$previous_variation_feature_id = $variation_feature_id;
	$var_consequence_type->{$consequence_type}++;
    }
    #and print the last variation
    foreach my $type (@consequence_types_ordered){
	if ($var_consequence_type->{$type}){
	    $highest_priority = $type;
	    last;
	}
    }
    print FH join("\t",$previous_variation_feature_id,$highest_priority), "\n";

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
