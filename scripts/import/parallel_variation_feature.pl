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
#use Data::Dumper;
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

  my $dbVar_write = DBH->connect
    ("DBI:mysql:host=$vhost;dbname=$vdbname;port=$vport",$vuser, $vpass );
  die("Could not connect to variation database: $!") if(!$dbVar_write);

  $TMP_DIR  = $ImportUtils::TMP_DIR;
  $TMP_FILE = $ImportUtils::TMP_FILE;

  load_asm_cache($dbCore) unless defined $top_level;
  variation_feature($dbCore, $dbVar, $dbVar_write, $top_level);

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
       #if is the last process, delete the variation table and upload with the new coordinate information from the different tables and the failed_Variation information
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
  my $dbVar_write = shift;
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

  #my $sctg_cs = $csa->fetch_by_rank(2); #replaced to support other species without supercontig

  my $mapper = $asma->fetch_by_CoordSystems($top_cs, $sctg_cs) if (! $top_level);
  
  debug("Processing variation features");

  my $sth = $dbVar->prepare
    (qq{SELECT STRAIGHT_JOIN vf.variation_feature_id, vf.seq_region_id,
               vf.seq_region_start, vf.seq_region_end, vf.seq_region_strand,
               vf.variation_id, a.allele, tmw.count, v.name, vf.flags, vf.source_id, vf.validation_status, vf.consequence_type
        FROM   variation_feature vf FORCE INDEX(PRIMARY), variation v, tmp_map_weight tmw,
               allele_string a
        WHERE  a.variation_id = vf.variation_id 
        AND    vf.variation_id = tmw.variation_id
        AND    vf.variation_id = v.variation_id 
        $LIMIT
        ORDER BY vf.seq_region_id, vf.seq_region_start 
        #limit 10
	},{mysql_use_result=>1});

  $sth->execute();

  my ($vf_id, $sr_id, $sr_start, $sr_end, $sr_strand,
      $v_id, $allele, $map_weight, $v_name, $vf_flags, $vf_source_id, $vf_validation_status, $vf_consequence_type);

  $sth->bind_columns(\$vf_id, \$sr_id, \$sr_start, \$sr_end, \$sr_strand,
                     \$v_id, \$allele, \$map_weight, \$v_name, \$vf_flags, \$vf_source_id, \$vf_validation_status, \$vf_consequence_type);

  my ($cur_vf_id, $cur_map_weight, $cur_v_id, $cur_v_name, $cur_allele,
     $top_coord, $top_sr_id, $top_sr_start, $top_sr_end, $top_sr_strand, 
      $ref_allele,$cur_vf_flags, $cur_source_id, $cur_validation_status, $cur_consequence_type, $special);
  my %alleles;
  my %alleles_expanded; #same hash as before, but with the expanded alleles: $alleles_expanded{AGAGAG} = (AG)3

  my $dbname = $dbVar->dbname(); #get the name of the database to create the file
  my $host = `hostname`;
  chop $host;
  open FH, ">$TMP_DIR/$dbname.variation_feature_$host\:$$\.txt"
    or throw("Could not open tmp file: $TMP_DIR/variation_feature_$$\n");

  while($sth->fetch()) {
	
    #excluding SNPs with map_weight>3
    if ($map_weight>3) {	
      #needs to be written to the failed_variation table
      $dbVar_write->do(qq{INSERT IGNORE INTO failed_variation (variation_id,failed_description_id) VALUES ($v_id,1)});
      #next; # don't skip now we are flagging not failing
    }
    if(!defined($cur_vf_id) || $cur_vf_id != $vf_id) {
      
	  # check top level coord status
	  if($top_coord) {
        
		my $allele_str;
		
		# check coord length
		my $coord_length = $top_level ? $top_sr_end - $top_sr_start + 1 : $top_coord->end - $top_coord->start + 1;
		my $tmp_allele_str = join("/", keys %alleles);
		my $ok = 0;
		
		foreach my $allele(values %alleles_expanded) {
		  my $tmp_allele = $allele;
		  $tmp_allele =~ s/\-//g;
		  $ok = 1 if length($tmp_allele) == $coord_length;
		  
		  # special case alleles
		  $ok = 1 if $allele !~ /^[ACGT-]+$/;
		  
		  last if $ok;
		}
		
		if(!$ok) {
		  warn("None of the alleles $tmp_allele_str for $cur_v_name match the given coordinate length of $coord_length");
      
		  #needs to be written to the failed_variation table
		  $dbVar_write->do(qq{INSERT IGNORE INTO failed_variation (variation_id,failed_description_id) VALUES ($cur_v_id,15)}) if ($cur_map_weight==1);
		}
		
        # construct an allele string. Remember to check the expanded version
        if($alleles_expanded{$ref_allele}) {
          # make sure the reference allele is first. Remember to convert ref_allele (in expanded version) to a compressed allele
          delete $alleles{$alleles_expanded{$ref_allele}};
          $allele_str = join('/', ($alleles_expanded{$ref_allele}, keys %alleles));
        }
		elsif (defined($special) && scalar keys %alleles <= 2) {		  
		  $allele_str = $special;
		}
		else {
		  #$allele_str = undef;  # don't delete allele string now that we're flagging not failing
		  $allele_str = join("/", keys %alleles);
			  warn("Reference allele $ref_allele for $cur_v_name not found in alleles: " .
			  $allele_str, " flagging feature cur_allele was $cur_allele");
		  $dbVar_write->do(qq{INSERT IGNORE INTO failed_variation (variation_id,failed_description_id) VALUES ($cur_v_id,2)}) if ($cur_map_weight ==1) && $ok; #only put into failed_variation when map_weight==1
		}
		
		#if($allele_str) {
		  $allele_str ||= '\N';
		
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
		#}
      }
      else {
		print "vf_id is $vf_id is no top_coord\n";
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
      $cur_allele = $allele;
      
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
		}
		
		else {
		  
		  # obtain the seq_region_id of the seq_region we mapped to
		  # and the reference allele from the genome sequence
		  ($top_coord) = @coords;
	
		  $slice = $slice_adaptor->fetch_by_seq_region_id
			($top_coord->id(),
			 $top_coord->start(),$top_coord->end(), $top_coord->strand(),
			 $top_coord->coord_system()->version());
		  if (! $slice) {
			print "$cur_v_id ",join " ", $top_coord->coord_system()->name(),$top_coord->id(),
			  $top_coord->start(),$top_coord->end(), $top_coord->strand(),
			$top_coord->coord_system()->version(),"\n";
			$ref_allele = "";
		  }
		  else {
			$ref_allele = $slice->seq();
			$ref_allele = '-' if(!$ref_allele);
			$ref_allele = uc $ref_allele;   #convert reference allele to uppercase
			$top_sr_id = $top_coord->id();
		    #$top_sr_id = $slice->get_seq_region_id();
		  }
		}
      }
	  
	  $special = undef;
    }

    $alleles{$allele} = 1;
    my $allele_copy = $allele;
	
	$special = $allele if $allele =~ /LARGE|INS|DEL|CNV|MUTA/;
	
    #make a copy of the allele, but in the expanded version
    &expand(\$allele) if ($allele !~ /LARGE|INS|DEL|CNV|MUTA/); #only expand alleles with the (AG)5 format
    $alleles_expanded{$allele} = $allele_copy; 
  }

  # check for problems which may have terminated the fetch early
  die $sth->errstr if $sth->err;


  # print the last row
  if($top_coord ) {
    my $allele_str;
	
	# check coord length
	my $coord_length = $top_level ? $top_sr_end - $top_sr_start + 1 : $top_coord->end - $top_coord->start + 1;
	my $tmp_allele_str = join("/", keys %alleles);
	my $ok = 0;
	foreach my $allele(values %alleles_expanded) {
	  my $tmp_allele = $allele;
	  $tmp_allele =~ s/\-//g;
	  $ok = 1 if length($tmp_allele) == $coord_length;
	  
	  # special case alleles don't need to match length
	  $ok = 1 if $allele !~ /^[ACGT-]+$/;
	  
	  last if $ok;
	}
	
	if(!$ok) {
	  warn("None of the alleles $tmp_allele_str for $cur_v_name match the given coordinate length of $coord_length");
      
	  #needs to be written to the failed_variation table
      $dbVar_write->do(qq{INSERT IGNORE INTO failed_variation (variation_id,failed_description_id) VALUES ($cur_v_id,15)}) if ($cur_map_weight==1);
	}

    if($alleles{$ref_allele}) {
      # make sure the reference allele is first
      delete $alleles{$ref_allele};
      $allele_str = join('/', ($ref_allele, keys %alleles));
    }
	elsif (defined($special) && scalar keys %alleles <= 2) {		  
	  $allele_str = $special;
	}
    else {
      #$allele_str = undef; # don't delete allele string now that we're flagging not failing
	  $allele_str = join("/", keys %alleles);
      warn("Reference allele $ref_allele for $cur_v_name not found in alleles: $allele_str flagging feature\n");
      
	  #needs to be written to the failed_variation table
      $dbVar_write->do(qq{INSERT IGNORE INTO failed_variation (variation_id,failed_description_id) VALUES ($cur_v_id,2)}) if ($cur_map_weight==1) && $ok;
    }
    
    #if($allele_str) {
	  $allele_str ||= '\N';
	
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
    #}
  }
  else {
	print "again vf_id is $vf_id is no top_coord\n";
  }

  close FH;
}

#will know if it is the last process running counting the lines of the status_file.If so, load all data
sub last_process{
    my $dbCore = shift;
    my $dbVar = shift;
    
    debug("Deleting existing variation features");
    $dbVar->do("TRUNCATE TABLE variation_feature");

    debug("Reimporting processed variation features");
    my $dbname = $dbVar->dbname(); #get the name of the database to create the file    
    #group all the fragments in 1 file
    my $call = "cat $TMP_DIR/$dbname.variation_feature*.txt > $TMP_DIR/$TMP_FILE";
    system($call);

    #and finally, load the information
    load($dbVar, qw(variation_feature variation_feature_id seq_region_id
		    seq_region_start seq_region_end seq_region_strand variation_id
		    allele_string variation_name map_weight flags source_id validation_status consequence_type));

    unlink(<$TMP_DIR/$dbname.variation_feature*>);
    $dbVar->do("DROP TABLE tmp_map_weight");

    #and delete the status file

    unlink("$TMP_DIR/$status_file");
    update_meta_coord($dbCore, $dbVar, 'variation_feature');
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
    ('INSERT IGNORE INTO meta_coord set table_name = ?, coord_system_id = ?, max_length =500');

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
