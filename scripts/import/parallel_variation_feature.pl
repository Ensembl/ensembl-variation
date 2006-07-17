use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use Fcntl ':flock';
use DBI;
use DBH;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use Bio::EnsEMBL::Utils::Sequence qw(expand);
use ImportUtils qw(debug load);

use constant MITHOCONDRIAL_CONTIG => 'MT'; #the human mithocondrial contigs are already in top_level

my ($TMP_DIR, $TMP_FILE, $LIMIT,$status_file);

{
  my ($vhost, $vport, $vdbname, $vuser, $vpass,
      $chost, $cport, $cdbname, $cuser, $cpass,
      $limit, $num_processes, $top_level);

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
	     'toplevel' => \$top_level,
	     'status_file=s' => \$status_file);

  #added default options
  $chost    ||= 'ecs2';
  $cuser    ||= 'ensro';
  $cport    ||= 3365;

  $vport    ||= 3306;
  $vuser    ||= 'ensadmin';
  
  $num_processes ||= 1;

  $LIMIT = ($limit) ? " $limit " : '';

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

  load_asm_cache($dbCore) unless defined $top_level;
  variation_feature($dbCore, $dbVar, $top_level);

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
      #if is the last process, delete the variation table and upload with the new coordinate information from the different tables
      last_process($dbCore,$dbVar);
  }
}




#
# preloads the mapper cache
#
sub load_asm_cache {
  my $dbCore = shift;

  debug("Building assembly mapping data cache");

  my $slice_adaptor = $dbCore->get_SliceAdaptor();
  my $asma = $dbCore->get_AssemblyMapperAdaptor();
  my $csa  = $dbCore->get_CoordSystemAdaptor();

  my $top_cs  = $csa->fetch_by_name('chromosome');
#  my $sctg_cs = $csa->fetch_by_name('supercontig');
  my $sctg_cs = $csa->fetch_by_rank(2); #replaced to support other species without supercontig
  my $seq_cs  = $csa->fetch_by_name('seqlevel');

  my $mapper1 = $asma->fetch_by_CoordSystems($top_cs, $sctg_cs);
  my $mapper2 = $asma->fetch_by_CoordSystems($top_cs, $seq_cs);

  debug("Registering all chromosome/contig assembly information");

  $mapper2->max_pair_count(10e6);
  $mapper2->register_all();

  debug("Registering all superctg/chromosome assembly information");

  $mapper1->max_pair_count(10e6);
  $mapper1->register_all();


  debug("Registration DONE");

  return;
}


#
# moves variation features to toplevel,
# sets map_weight, sets allele_string
#
sub variation_feature {
  my $dbCore = shift;
  my $dbVar  = shift;
  my $top_level = shift;
  my $mithocondrial = 0;
  my $slice_adaptor = $dbCore->get_SliceAdaptor();
  my $asma = $dbCore->get_AssemblyMapperAdaptor();
  my $csa  = $dbCore->get_CoordSystemAdaptor();

  my $top_cs  = $csa->fetch_by_name('toplevel');
  my $sctg_cs;
  if ($dbCore->dbc()->dbname() =~ /anopheles|danio/i){
      #there are no supercontig in these species, but scaffold
      $sctg_cs = $csa->fetch_by_name('scaffold');
  }
  elsif ($dbCore->dbc()->dbname() =~ /canis/i){
      $sctg_cs = $csa->fetch_by_rank(5);
  }
  else{
      $sctg_cs = $csa->fetch_by_name('supercontig');
  }

#  my $sctg_cs = $csa->fetch_by_rank(2); #replaced to support other species without supercontig

  my $mapper = $asma->fetch_by_CoordSystems($top_cs, $sctg_cs) if (! $top_level);
  
  debug("Processing variation features");

  my $sth = $dbVar->prepare
    (qq{SELECT STRAIGHT_JOIN vf.variation_feature_id, vf.seq_region_id,
               vf.seq_region_start, vf.seq_region_end, vf.seq_region_strand,
               vf.variation_id, a.allele, tmw.count, v.name, vf.flags, vf.source_id, vf.validation_status, vf.consequence_type
        FROM   variation_feature vf FORCE INDEX(PRIMARY), allele_string a, tmp_map_weight tmw,
               variation v
        WHERE  a.variation_id = vf.variation_id 
        AND    vf.variation_id = tmw.variation_id
        AND    vf.variation_id = v.variation_id 
	$LIMIT 
        ORDER BY vf.seq_region_id, vf.seq_region_start 
	},{mysql_use_result=>1});

  $sth->execute();

  my ($vf_id, $sr_id, $sr_start, $sr_end, $sr_strand,
      $v_id, $allele, $map_weight, $v_name, $vf_flags, $vf_source_id, $vf_validation_status, $vf_consequence_type);

  $sth->bind_columns(\$vf_id, \$sr_id, \$sr_start, \$sr_end, \$sr_strand,
                     \$v_id, \$allele, \$map_weight, \$v_name, \$vf_flags, \$vf_source_id, \$vf_validation_status, \$vf_consequence_type);

  my ($cur_vf_id, $cur_map_weight, $cur_v_id, $cur_v_name,
     $top_coord, $top_sr_id, $top_sr_start, $top_sr_end, $top_sr_strand, 
      $ref_allele,$cur_vf_flags, $cur_source_id, $cur_validation_status, $cur_consequence_type);
  my %alleles;
  my %alleles_expanded; #same hash as before, but with the expanded alleles: $alleles_expanded{AGAGAG} = (AG)3

  my $dbname = $dbVar->dbname(); #get the name of the database to create the file
  my $host = `hostname`;
  chop $host;
  #open FH, ">$TMP_DIR/$dbname.variation_feature_$host\:$$\.txt"
  open FH, ">/tmp/$dbname.variation_feature_$host\_$$\.txt"
    or throw("Could not open tmp file: $TMP_DIR/variation_feature_$$\n");

  while($sth->fetch()) {
    next if $map_weight >3; #excluding SNPs with map_weight > 3
    if(!defined($cur_vf_id) || $cur_vf_id != $vf_id) {
      if($top_coord) {
        my $allele_str;
        # construct an allele string. Remember to check the expanded version

        if($alleles_expanded{$ref_allele}) {
          # make sure the reference allele is first. Remember to convert ref_allele (in expanded version) to a compressed allele
          delete $alleles{$alleles_expanded{$ref_allele}};
          $allele_str = join('/', ($alleles_expanded{$ref_allele}, keys %alleles));
        } else {
	  $allele_str = undef;
          warn("Reference allele $ref_allele for $cur_v_name not found in alleles: " .
               join("/", keys %alleles), " discarding feature");
        }
	
	if($allele_str) {
	  if ($top_level) {
	    print FH join("\t", $cur_vf_id, $top_sr_id, $top_sr_start, $top_sr_end, $top_sr_strand,
			  $cur_v_id, $allele_str, $cur_v_name,
			  $cur_map_weight,$cur_vf_flags,$cur_source_id, $cur_validation_status,$cur_consequence_type), "\n";
	  }
	  else {
	    print FH join("\t", $cur_vf_id, $top_sr_id, $top_coord->start(),
			  $top_coord->end(), $top_coord->strand(),
			  $cur_v_id, $allele_str, $cur_v_name,
			  $cur_map_weight,$cur_vf_flags,$cur_source_id,$cur_validation_status,$cur_consequence_type), "\n";
	  }
	}
      }
      
      %alleles = ();
      %alleles_expanded = ();
      $ref_allele = undef;
      $top_sr_id = undef;
      $top_coord = undef;
      $cur_vf_id = $vf_id;
      $cur_map_weight = $map_weight;
      $cur_v_id  = $v_id;
      $cur_v_name = $v_name;
      $cur_vf_flags = 1 if (defined $vf_flags);
      $cur_vf_flags = '\N' if (! defined $vf_flags);
      $cur_source_id = $vf_source_id;
      $cur_validation_status = $vf_validation_status if (defined $vf_validation_status);
      $cur_validation_status = '\N' if (! defined $vf_validation_status);
      $cur_consequence_type = $vf_consequence_type;
      $top_sr_start = $sr_start;
      $top_sr_end = $sr_end;
      $top_sr_strand = $sr_strand;
      $allele = uc $allele;    #convert allele to uppercase format

      # map the variation coordinates to toplevel

      my $slice = $slice_adaptor->fetch_by_seq_region_id($sr_id);
      
      if(!$slice) {
        warning("Could not locate seq_region with id=$sr_id");
        next;
      }
      if ($slice->seq_region_name() eq MITHOCONDRIAL_CONTIG()){$mithocondrial = 1} #the human mithocondrial conigs are already top_level
      else{$mithocondrial = 0} #for the res
      ###if it is a toplevel coordinates already, we need find ref_allele for it
      if ($top_level || $mithocondrial) {
	$top_coord = $slice->sub_Slice($sr_start, $sr_end,$sr_strand);
	if (!$top_coord) {
	  ###if start = end+1, this indicate a indel
	  if ($sr_start == $sr_end+1) {
	    $ref_allele = '-';
	    $top_sr_id = $sr_id;
	    $top_coord=1;
	    #warning("Could not locate $sr_id, $sr_start, $sr_end, $sr_strand, maybe indels");
	  }
	  else {
	    next;
	  }
	}
	else {
	  $ref_allele = $top_coord->seq();
	  $ref_allele = '-' if(!$ref_allele);
	  $ref_allele = uc $ref_allele;   #convert reference allele to uppercase
	  $top_sr_id = $sr_id;
	  debug ("ref_allele is '-' for $sr_id, $sr_start, $sr_end, $sr_strand") if ($ref_allele eq '-');
	}
      }
      else {
	#debug("seq_region_name is ",$slice->seq_region_name(),"$sr_start, $sr_end, $sr_strand\n");
	my @coords = $mapper->map($slice->seq_region_name(), $sr_start, $sr_end,
				  $sr_strand, $sctg_cs);
	
	if (@coords != 1 || $coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
	  $top_coord = undef;
	  warn ("top_coord is not found for $vf_id\n");
	} else {
	  
	  # obtain the seq_region_id of the seq_region we mapped to
	  # and the reference allele from the genome sequence
	  ($top_coord) = @coords;
	  $slice = $slice_adaptor->fetch_by_region
	    ($top_coord->coord_system()->name(),$top_coord->id(),
	     $top_coord->start(),$top_coord->end(), $top_coord->strand(),
	     $top_coord->coord_system()->version());
	  
	  $ref_allele = $slice->seq();
	  $ref_allele = '-' if(!$ref_allele);
	  $ref_allele = uc $ref_allele;   #convert reference allele to uppercase

	  $top_sr_id = $slice->get_seq_region_id();
	}
      }
    }

    $alleles{$allele} = 1;
    my $allele_copy = $allele;
    #make a copy of the allele, but in the expanded version
    &expand(\$allele) if ($allele !~ /LARGE/); #only expand alleles with the (AG)5 format
    $alleles_expanded{$allele} = $allele_copy; 
  }

  $sth->finish();

  # print the last row, excluding SNPs with map_weight > 3
  
  if($top_coord and $cur_map_weight <=3) {

    my $allele_str;

    if($alleles{$ref_allele}) {
      # make sure the reference allele is first
      delete $alleles{$ref_allele};
      $allele_str = join('/', ($ref_allele, keys %alleles));
    } else {
      $allele_str = undef;
      warn("Reference allele $ref_allele for $cur_v_name not found in alleles: " .
           join("/", keys %alleles), " discarding feature\n");
    }
    
    if($allele_str) {
      if ($top_level) {
	print FH join("\t", $cur_vf_id, $top_sr_id, $top_sr_start, $top_sr_end, $top_sr_strand,
		      $cur_v_id, $allele_str, $cur_v_name,
		      $cur_map_weight,$cur_vf_flags, $cur_source_id,$cur_validation_status,$cur_consequence_type), "\n";
      }
      else {
	print FH join("\t", $cur_vf_id, $top_sr_id, $top_coord->start(),
		      $top_coord->end(), $top_coord->strand(),
		      $cur_v_id, $allele_str, $cur_v_name, $cur_map_weight, $cur_vf_flags, $cur_source_id, $cur_validation_status,$cur_consequence_type), "\n";
      }
    }
  }

  close FH;
  my $call = "lsrcp /tmp/$dbname.variation_feature_$host\_$$\.txt $TMP_DIR/$dbname.variation_feature_$host\_$$\.txt";
  system($call);
  unlink("/tmp/$dbname.variation_feature_$host\_$$\.txt");
}

#will know if it is the last process running counting the lines of the status_file.If so, load all data
sub last_process{
    my $dbCore = shift;
    my $dbVar = shift;
    
    debug("Deleting existing variation features");
    
    $dbVar->do("DELETE FROM variation_feature");
    
    debug("Reimporting processed variation features");
    my $dbname = $dbVar->dbname(); #get the name of the database to create the file    
    #group all the fragments in 1 file
    my $call = "cat $TMP_DIR/$dbname.variation_feature*.txt > $TMP_DIR/$TMP_FILE";
    system($call);

    #and finally, load the information
    load($dbVar, qw(variation_feature variation_feature_id seq_region_id
		    seq_region_start seq_region_end seq_region_strand variation_id
		    allele_string variation_name map_weight flags source_id validation_status consequence_type));

    #unlink(<$TMP_DIR/$dbname.variation_feature*>);
    $dbVar->do("DROP TABLE tmp_map_weight");

    update_meta_coord($dbCore, $dbVar, 'variation_feature');
    #and delete the status file

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
    ('INSERT INTO meta_coord set table_name = ?, coord_system_id = ?, max_length =500');

  $sth->execute($table_name, $cs->dbID());

  $sth->finish();

  return;
}

sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: perl parallel_variation_feature.pl <options>

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
