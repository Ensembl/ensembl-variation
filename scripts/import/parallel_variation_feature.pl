use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use DBI;
use DBH;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use Bio::EnsEMBL::Utils::Sequence qw(expand);
use ImportUtils qw(debug load);


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
#  $cport    ||= 3364;

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
  print STATUS "process finished\n";
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
  my $sctg_cs = $csa->fetch_by_name('supercontig');
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

  my $slice_adaptor = $dbCore->get_SliceAdaptor();
  my $asma = $dbCore->get_AssemblyMapperAdaptor();
  my $csa  = $dbCore->get_CoordSystemAdaptor();

  my $top_cs  = $csa->fetch_by_name('toplevel');
  my $sctg_cs = $csa->fetch_by_name('supercontig');

  my $mapper = $asma->fetch_by_CoordSystems($top_cs, $sctg_cs) if (! $top_level);
  
  debug("Processing variation features");

  my $sth = $dbVar->prepare
    (qq{SELECT vf.variation_feature_id, vf.seq_region_id,
               vf.seq_region_start, vf.seq_region_end, vf.seq_region_strand,
               vf.variation_id, a.allele, tmw.count, v.name
        FROM   variation_feature vf, allele a, tmp_map_weight tmw,
               variation v
        WHERE  a.variation_id = vf.variation_id
        AND    vf.variation_id = tmw.variation_id
        AND    vf.variation_id = v.variation_id
	$LIMIT
        ORDER BY vf.seq_region_id, vf.seq_region_start, vf.variation_feature_id},{mysql_use_result => 1});

  $sth->execute();

  my ($vf_id, $sr_id, $sr_start, $sr_end, $sr_strand,
      $v_id, $allele, $map_weight, $v_name);

  $sth->bind_columns(\$vf_id, \$sr_id, \$sr_start, \$sr_end, \$sr_strand,
                     \$v_id, \$allele, \$map_weight, \$v_name);

  my ($cur_vf_id, $cur_map_weight, $cur_v_id, $cur_v_name,
     $top_coord, $top_sr_id, $top_sr_start, $top_sr_end, $top_sr_strand, $ref_allele);
  my %alleles;
  my %alleles_expanded; #same hash as before, but with the expanded alleles: $alleles_expanded{AGAGAG} = (AG)3

  open FH, ">$TMP_DIR/variation_feature_$$\.txt"
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
          warn("Reference allele $ref_allele not found in alleles: " .
               join("/", keys %alleles), " discarding feature");
        }
	
        if($allele_str) {
	  if ($top_level) {
	    print FH join("\t", $cur_vf_id, $top_sr_id, $top_sr_start, $top_sr_end, $top_sr_strand,
			  $cur_v_id, $allele_str, $cur_v_name,
			  $cur_map_weight), "\n";
	  }
	  else {
	    print FH join("\t", $cur_vf_id, $top_sr_id, $top_coord->start(),
			  $top_coord->end(), $top_coord->strand(),
			  $cur_v_id, $allele_str, $cur_v_name,
			  $cur_map_weight), "\n";
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
      
      $top_sr_start = $sr_start;
      $top_sr_end = $sr_end;
      $top_sr_strand = $sr_strand;

      # map the variation coordinates to toplevel

      my $slice = $slice_adaptor->fetch_by_seq_region_id($sr_id);

      if(!$slice) {
        warning("Could not locate seq_region with id=$sr_id");
        next;
      }

      ###if it is a toplevel coordinates already, we need find ref_allele for it
      if ($top_level) {
	$top_coord = $slice->sub_Slice($sr_start, $sr_end,$sr_strand);
	if (!$top_coord) {
	  ###if start = end+1, this indicate a indel
	  if ($sr_start == $sr_end+1) {
	    $ref_allele = '-';
	    $top_sr_id = $sr_id;
	    $top_coord=1;
	    warning("Could not locate $sr_id, $sr_start, $sr_end, $sr_strand, maybe indels");
	  }
	  else {
	    next;
	  }
	}
	else {
	  $ref_allele = $top_coord->seq();
	  $ref_allele = '-' if(!$ref_allele);
	  $top_sr_id = $sr_id;
	  debug ("ref_allele is '-' for $sr_id, $sr_start, $sr_end, $sr_strand") if ($ref_allele eq '-');
	}
      }
      else {
	my @coords = $mapper->map($slice->seq_region_name(), $sr_start, $sr_end,
				  $sr_strand, $sctg_cs);
	
	if (@coords != 1 || $coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
	  $top_coord = undef;
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

  # print the last row
  if($top_coord) {
    my $allele_str;

    if($alleles{$ref_allele}) {
      # make sure the reference allele is first
      delete $alleles{$ref_allele};
      $allele_str = join('/', ($ref_allele, keys %alleles));
    } else {
      $allele_str = undef;
      warn("Reference allele $ref_allele not found in alleles: " .
           join("/", keys %alleles), " discarding feature\n");
    }

    if($allele_str) {
      if ($top_level) {
	print FH join("\t", $cur_vf_id, $top_sr_id, $top_sr_start, $top_sr_end, $top_sr_strand,
		      $cur_v_id, $allele_str, $cur_v_name,
		      $map_weight), "\n";
      }
      else {
	print FH join("\t", $vf_id, $top_sr_id, $top_coord->start(),
		      $top_coord->end(), $top_coord->strand(),
		      $cur_v_id, $allele_str, $v_name, $map_weight), "\n";
      }
    }
  }

  close FH;

}

#will know if it is the last process running counting the lines of the status_file.If so, load all data
sub last_process{
    my $dbCore = shift;
    my $dbVar = shift;
    
    debug("Deleting existing variation features");
    
    $dbVar->do("DELETE FROM variation_feature");
    
    debug("Reimporting processed variation features");
    
    #group all the fragments in 1 file
    my $call = "cat $TMP_DIR/variation_feature*.txt > $TMP_DIR/$TMP_FILE";
    system($call);

    unlink(<$TMP_DIR/variation_feature*.txt>);
    #and finally, load the information
    load($dbVar, qw(variation_feature variation_feature_id seq_region_id
			     seq_region_start seq_region_end seq_region_strand variation_id
			     allele_string variation_name map_weight));

    $dbVar->do("DROP TABLE tmp_map_weight");

    update_meta_coord($dbCore, $dbVar, 'variation_feature');
    #and delete the status file

#    unlink("$TMP_DIR/$status_file");
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
