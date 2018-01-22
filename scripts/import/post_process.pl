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
use Data::Dumper;
use Benchmark;
use DBI;
use DBH;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use ImportUtils qw(dumpSQL debug create_and_load load);


my ($TMP_DIR, $TMP_FILE, $LIMIT);

{
  my ($vhost, $vport, $vdbname, $vuser, $vpass,
      $chost, $cport, $cdbname, $cuser, $cpass,
      $limit);

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
             'limit=i'   => \$limit);

  #added default options
  $chost    ||= 'ecs2';
  $cuser    ||= 'ensro';
#  $cport    ||= 3364;

  $vport    ||= 3306;
  $vuser    ||= 'ensadmin';

  $LIMIT = ($limit) ? " LIMIT $limit " : '';

  usage('-vdbname argument is required') if(!$vdbname);
  usage('-cdbname argument is required') if(!$cdbname);

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

#  load_asm_cache($dbCore);
#  variation_feature($dbCore, $dbVar);
  flanking_sequence($dbCore, $dbVar);
#  transcript_variation($dbCore, $dbVar);
#  ld_populations($dbCore,$dbVar);

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

  $dbVar->do(qq{CREATE TEMPORARY TABLE tmp_map_weight
                SELECT variation_id, count(*) as count
                FROM   variation_feature
                GROUP BY variation_id});

  $dbVar->do(qq{ALTER TABLE tmp_map_weight 
                ADD INDEX variation_idx(variation_id)});

  my $slice_adaptor = $dbCore->get_SliceAdaptor();
  my $asma = $dbCore->get_AssemblyMapperAdaptor();
  my $csa  = $dbCore->get_CoordSystemAdaptor();

  my $top_cs  = $csa->fetch_by_name('toplevel');
  my $sctg_cs = $csa->fetch_by_name('supercontig');

  my $mapper = $asma->fetch_by_CoordSystems($top_cs, $sctg_cs);

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
        GROUP BY vf.variation_feature_id, a.allele
        ORDER BY vf.seq_region_id, vf.seq_region_start,
                 variation_feature_id});


  $sth->execute();

  my ($vf_id, $sr_id, $sr_start, $sr_end, $sr_strand,
      $v_id, $allele, $map_weight, $v_name);

  $sth->bind_columns(\$vf_id, \$sr_id, \$sr_start, \$sr_end, \$sr_strand,
                     \$v_id, \$allele, \$map_weight, \$v_name);

  my ($cur_vf_id, $cur_map_weight, $cur_v_id, $cur_v_name,
     $top_coord, $top_sr_id, $ref_allele);
  my %alleles;

  open FH, ">$TMP_DIR/$TMP_FILE"
    or throw("Could not open tmp file: $TMP_DIR/$TMP_FILE\n");

  while($sth->fetch()) {
    if(!defined($cur_vf_id) || $cur_vf_id != $vf_id) {
      if($top_coord) {
        my $allele_str;

        # construct an allele string
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
          print FH join("\t", $cur_vf_id, $top_sr_id, $top_coord->start(),
                        $top_coord->end(), $top_coord->strand(),
                        $cur_v_id, $allele_str, $cur_v_name,
                        $cur_map_weight), "\n";
        }
      }

      %alleles = ();
      $ref_allele = undef;
      $top_sr_id = undef;
      $cur_vf_id = $vf_id;
      $cur_map_weight = $map_weight;
      $cur_v_id  = $v_id;
      $cur_v_name = $v_name;

      # map the variation coordinates to toplevel

      my $slice = $slice_adaptor->fetch_by_seq_region_id($sr_id);

      if(!$slice) {
        warning("Could not locate seq_region with id=$sr_id");
        next;
      }

      my @coords = $mapper->map($slice->seq_region_name(), $sr_start, $sr_end,
                                $sr_strand, $sctg_cs);

      if(@coords != 1 || $coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
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

    $alleles{$allele} = 1;
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
      print FH join("\t", $vf_id, $top_sr_id, $top_coord->start(),
                        $top_coord->end(), $top_coord->strand(),
                        $cur_v_id, $allele_str, $v_name, $map_weight), "\n";
    }
  }

  close FH;

   debug("Deleting existing variation features");

   $dbVar->do("DELETE FROM variation_feature");

   debug("Reimporting processed variation features");

   load($dbVar, qw(variation_feature variation_feature_id seq_region_id
           seq_region_start seq_region_end seq_region_strand variation_id
           allele_string variation_name map_weight));

  $dbVar->do("DROP TABLE tmp_map_weight");

  update_meta_coord($dbCore, $dbVar, 'variation_feature');
}



#
# Compresses flanking sequence storage by only storing genomic coordinates
# when the genomic sequence exactly matches the flanking sequence.
#
sub flanking_sequence {
  my $dbCore = shift;
  my $dbVar  = shift;

  debug("Compressing storage of flanking sequence");

  my $slice_adaptor = $dbCore->get_SliceAdaptor();

  my $update_sth = $dbVar->prepare
    (qq{UPDATE flanking_sequence
        SET up_seq   = ?,
            down_seq = ?,
            up_seq_region_start = ?,
            up_seq_region_end = ?,
            down_seq_region_start = ?,
            down_seq_region_end = ?,
            seq_region_id = ?,
            seq_region_strand = ?});

  my $sth = $dbVar->prepare(qq{SELECT fs.variation_id, fs.up_seq, fs.down_seq,
                                      vf.seq_region_id, vf.seq_region_start,
                                      vf.seq_region_end, vf.seq_region_strand
                               FROM variation_feature vf FORCE INDEX(variation_idx), flanking_sequence fs
                               WHERE vf.variation_id = fs.variation_id
                               GROUP BY vf.variation_id
#                               ORDER BY vf.seq_region_id, vf.seq_region_start
                               $LIMIT},
                            );


  $sth->execute();

  my ($var_id, $up_seq, $dn_seq, $sr_id, $sr_start, $sr_end, $sr_strand);
  $sth->bind_columns(\$var_id, \$up_seq, \$dn_seq, \$sr_id, \$sr_start,
                     \$sr_end, \$sr_strand);

  while($sth->fetch()) {
    my $up_len = length($up_seq);
    my $dn_len = length($dn_seq);

    # figure out the coordinates of the flanking sequence

    my ($up_sr_start, $up_sr_end, $dn_sr_start, $dn_sr_end);

    if($sr_strand == 1) {
      $up_sr_start = $sr_start - $up_len;
      $up_sr_end   = $sr_start - 1;
      $dn_sr_start = $sr_end + 1;
      $dn_sr_end   = $sr_end + $dn_len;
    } else {
      $up_sr_start = $sr_end + 1;
      $up_sr_end   = $sr_end + $up_len;
      $dn_sr_start = $sr_start - $dn_len;
      $dn_sr_end   = $sr_start - 1;
    }

    my $slice = $slice_adaptor->fetch_by_seq_region_id($sr_id);

    if(!$slice) {
      warning("Could not obtain slice for seq_region_id $sr_id\n");
      next;
    }

    # compare sequence in database to flanking sequence
    # if it matches store only coordinates, otherwise still store flanking

    # there are sometimes off by ones in dbSNP, try to take this into account
    # with a 'wobble' of one base on either side
    my @wobble = (0, 1, -1);
    my $i = 0;

    while(defined($up_seq) && defined($dn_seq) && $i < 3) {
      my $w = $wobble[$i++];

      if(defined($up_seq)) {
        my $up = $slice->subseq($up_sr_start+$w, $up_sr_end+$w, $sr_strand);
        $up_seq = undef if(uc($up) eq uc($up_seq));
        print STDERR "*"  if($w && !defined($up_seq));
      }

      if(defined($dn_seq)) {
        my $dn = $slice->subseq($dn_sr_start+$w, $dn_sr_end+$w, $sr_strand);
        $dn_seq = undef if(uc($dn) eq uc($dn_seq));
        print STDERR "*" if($w && !defined($dn_seq));
      }
    }

    $up_sr_start = $up_sr_end = undef if(defined($up_seq));
    $dn_sr_start = $dn_sr_end = undef if(defined($dn_seq));

    # if we changed something update the row in the database
    if(!defined($dn_seq) || !defined($up_seq)) {
      $update_sth->execute($up_seq, $dn_seq, $up_sr_start, $up_sr_end,
                           $dn_sr_start, $dn_sr_end, $sr_id, $sr_strand);
    }
  }

  $sth->finish();
  $update_sth->finish();

  return;
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

  my $sth = $dbVar->prepare
    (qq{SELECT vf.variation_feature_id, vf.seq_region_start, vf.seq_region_end,
               vf.seq_region_strand, vf.allele_string
        FROM   variation_feature vf
        WHERE  vf.seq_region_id = ?
        AND    vf.seq_region_end >= ?
        AND    vf.seq_region_start <= ?});

  my $sa = $dbCore->get_SliceAdaptor();

  open FH, ">$TMP_DIR/$TMP_FILE";

  my $inc_non_ref = 1;
  my $slices = $sa->fetch_all('toplevel', undef, $inc_non_ref);

  # assumes that variation features have already been pushed to toplevel
  foreach my $slice (@$slices) {
    debug("Processing transcript variations for ",
          $slice->seq_region_name(), "\n");
    my $genes = $slice->get_all_Genes();

    # request all variations which lie in the region of a gene

    foreach my $g (@$genes) {
      $sth->execute($slice->get_seq_region_id(),
                    $g->seq_region_start() - $UPSTREAM,
                    $g->seq_region_end()   + $DNSTREAM);

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

          my @alleles = split('/', $row->[4]);
          if($var{'strand'} != $tr->strand()) {
            # flip feature onto same strand as transcript
            for(my $i = 0; $i < @alleles; $i++) {
              reverse_comp(\$alleles[$i]);
            }
            $var{'strand'} = $tr->strand();
          }
          $var{'alleles'} = \@alleles;

          my $vars = type_variation($tr, \%var);

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

  load($dbVar, qw(transcript_variation
                  transcript_id variation_feature_id peptide_allele_string
                  translation_start translation_end cdna_start cdna_end type));


  


  return;
}


#for a list of populations, creates the pairwise_ld information
sub ld_populations{
    my $dbCore = shift;
    my $dbVar = shift;
    my $population_id;

    my $sth = $dbVar->prepare(qq{SELECT distinct(i.population_id)
				     FROM individual i, individual_genotype ig 
				     WHERE ig.individual_id = i.individual_id
				     $LIMIT
    });
    $sth->execute();
    $sth->bind_columns(\$population_id);
    #for all the populations, we calculate the pairwise_ld informaiton
    while ($sth->fetch()){
	my $t0 = new Benchmark;
	&pairwise_ld($dbCore,$dbVar,$population_id);
	my $t1 = new Benchmark;
	my $td = timediff ($t1,$t0);
	print "the code took:",timestr($td),"\n";
    }
    $sth->finish;
}
#
# populates the table pairwise_ld with the information for the genotypes for a certain population
#
sub pairwise_ld{
    my $dbCore = shift;
    my $dbVar = shift;
    my $population_id = shift;

    my %seq_region; #hash containing the mapping between seq_region_id->name region
    my %alleles_variation = (); #will contain a record of the alleles in the variation. A will be the major, and a the minor. When more than 2 alleles
    my %genotype_information; #will contain all the genotype information to write in the file, if necessary
    #, the genotypes for that variation will be discarded
    my $previous_variation_id = ''; #to know if it is a new variation and we can get the new alleles
    my $genotype; #genotype in the AA, Aa or aa format to write to the file
    debug("Loading pairwise_ld table for population $population_id\n");

    #necessary the order to know when we change variation
    my $sth = $dbVar->prepare
	(qq{SELECT ig.variation_id, vf.variation_feature_id, vf.seq_region_id, vf.seq_region_start, ig.individual_id, ig.allele_1, ig.allele_2, vf.seq_region_end
		FROM   variation_feature vf, individual_genotype ig, individual i
		WHERE  ig.variation_id = vf.variation_id
		AND i.individual_id = ig.individual_id
		AND i.population_id = ?
		ORDER BY ig.variation_id,vf.seq_region_id	       
		});
    $sth->execute($population_id);
    open ( FH, ">$TMP_DIR/temp_genotype.txt");

    #retrieve all the genotypes and format them into the calculation_genotype required format:
    # snp_id chr position individual_id genotype
    my ($variation_id, $variation_feature_id, $seq_region_id, $seq_region_start, $individual_id, $allele_1,$allele_2,$seq_region_end);
    my $seq_region_name;
    $sth->bind_columns(\$variation_id, \$variation_feature_id, \$seq_region_id, \$seq_region_start, \$individual_id, \$allele_1,\$allele_2,\$seq_region_end);
    my $count_not_discarded = 0;
    my $count_discarded = 0;
    while ($sth->fetch()){
	if ($previous_variation_id eq ''){
	    $previous_variation_id = $variation_id;
	}
	#if it is a new variation, write to the file (if necessary) and empty the hash
	if ($previous_variation_id ne $variation_id){
	    #if the variation has 2 or 1 alleles, print all the genotypes to the file
	    if (keys %alleles_variation <= 2){		
		$count_not_discarded++;
		&convert_genotype(\%alleles_variation,\%genotype_information);
		foreach my $individual_id (keys %genotype_information){
		    print FH join("\t",$previous_variation_id,$genotype_information{$individual_id}{seq_region_start},$individual_id,$genotype_information{$individual_id}{genotype},$genotype_information{$individual_id}{variation_feature_id},$genotype_information{$individual_id}{seq_region_id},$genotype_information{$individual_id}{seq_region_end},$population_id),"\n";
		}
	    }
	    else{
		$count_discarded++;
	    }
	    $previous_variation_id = $variation_id;
	    %alleles_variation = (); #new variation, flush the hash
	    %genotype_information = (); #new variation, flush the hash
	}
	#we store the genotype information for the variation
	$genotype_information{$individual_id}{variation_feature_id} = $variation_feature_id;
	$genotype_information{$individual_id}{seq_region_start} = $seq_region_start;
	$genotype_information{$individual_id}{allele_1} = $allele_1;
	$genotype_information{$individual_id}{allele_2} = $allele_2;
	$genotype_information{$individual_id}{seq_region_end} = $seq_region_end;
	$genotype_information{$individual_id}{seq_region_id} = $seq_region_id;

	#and the alleles
	$alleles_variation{$allele_1}++;
	$alleles_variation{$allele_2}++;
    }
    print "Number of variations discarded for population $population_id are: $count_discarded and used $count_not_discarded\n\n";
    $sth->finish();
    close FH;
    debug("calculating ld distance");
    #and, finally, create the table with the information, calling the calc_genotypes.pl script
    system("perl calc_genotypes.pl $TMP_DIR/temp_genotype.txt $TMP_DIR/temp_genotype_out.txt");
    #delete the original file
    unlink("$TMP_DIR/temp_genotype.txt");
    #rename the formated file
    rename("$TMP_DIR/temp_genotype_out.txt", "$TMP_DIR/$TMP_FILE");
    debug("Importing pairwise data");
    #and import the data in the database
    load($dbVar, qw(pairwise_ld variation_feature_id_1 variation_feature_id_2 population_id seq_region_id seq_region_start seq_region_end snp_distance_count r2 Dprime));
    unlink("$TMP_DIR/$TMP_FILE");
    return;
}
#
# Converts the genotype into the required format for the calculation of the pairwise_ld value: AA, Aa or aa
# From the Allele table, will select the alleles and compare to the alleles in the genotype
#
sub convert_genotype{
    my $alleles_variation = shift; #reference to the hash containing the alleles for the variation present in the genotypes
    my $genotype_information = shift; #reference to a hash containing the values to be written to the file
    my @alleles_ordered; #the array will contain the alleles ordered by apparitions in the genotypes (only 2 values possible)
    
    @alleles_ordered = sort({$alleles_variation->{$b} <=> $alleles_variation->{$a}} keys %{$alleles_variation});

    #let's convert the allele_1 allele_2 to a genotype in the AA, Aa or aa format, where A corresponds to the major allele and a to the minor
    foreach my $individual_id (keys %{$genotype_information}){
	#if both alleles are different, this is the Aa genotype
	if ($genotype_information->{$individual_id}{allele_1} ne $genotype_information->{$individual_id}{allele_2}){
	    $genotype_information->{$individual_id}{genotype} = 'Aa';
	}
	#when they are the same, must find out which is the major
	else{	    
	    if ($alleles_ordered[0] eq $genotype_information->{$individual_id}{allele_1}){
		#it is the major allele
		$genotype_information->{$individual_id}{genotype} = 'AA';
	    }
	    else{
		$genotype_information->{$individual_id}{genotype} = 'aa';
	    }
	    
	}
    }
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
      substr($cds, $var->{'cds_start'}, $var_len) = $a;
      my $codon_str = substr($cds, $codon_cds_start-1, $codon_len);

      my $codon_seq = Bio::Seq->new(-seq      => $codon_str,
                                    -moltype  => 'dna',
                                    -alphabet => 'dna');

      $new_aa = $codon_seq->translate()->seq();
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



sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: perl post_process.pl <options>

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
EOF

  die("\n$msg\n\n");
}
