
use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);

use DBI;


my $LIMIT = ' LIMIT 10000 ';

{
  my $vhost   = 'ecs4';
  my $vport   = 3352;
  my $vdbname = 'mcvicker_variation';
  my $vuser   = 'ensadmin';
  my $vpass   = 'ensembl';

  my $dbCore = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (-host   => 'ecs2',
     -user   => 'ensro',
     -port   => 3364,
     -dbname => 'mus_musculus_core_22_32b');

  my $dbVar = DBI->connect
    ("DBI:mysql:host=$vhost;dbname=$vdbname;port=$vport",$vuser, $vpass );
  die("Could not connect to variation database: $!") if(!$dbVar);

  variation_feature($dbCore, $dbVar);
  flanking_sequence($dbCore, $dbVar);
  variation_group_feature($dbCore, $dbVar);
  transcript_variation($dbCore, $dbVar);

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
        ORDER BY variation_feature_id});


  $sth->execute();

  my ($vf_id, $sr_id, $sr_start, $sr_end, $sr_strand,
      $v_id, $allele, $map_weight, $v_name);

  $sth->bind_columns(\$vf_id, \$sr_id, \$sr_start, \$sr_end, \$sr_strand,
                     \$v_id, \$allele, \$map_weight, \$v_name);

  my ($cur_vf_id, $cur_map_weight, $cur_v_id, $cur_v_name,
     $top_coord, $top_sr_id, $ref_allele);
  my %alleles;

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
          print STDERR join("\t", $cur_vf_id, $top_sr_id, $top_coord->start(),
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
      print STDERR join("\t", $vf_id, $top_sr_id, $top_coord->start(),
                        $top_coord->end(), $top_coord->strand(),
                        $cur_v_id, $allele_str, $v_name, $map_weight), "\n";
    }
  }

  $dbVar->do("DROP TABLE tmp_map_weight");

}



#
# Compresses flanking sequence storage by only storing genomic coordinates
# when the genomic sequence exactly matches the flanking sequence.
#
sub flanking_sequence {
  my $dbCore = shift;
  my $dbVar  = shift;

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
                               FROM variation_feature vf, flanking_sequence fs
                               WHERE vf.variation_id = fs.variation_id
                               GROUP BY fs.variation_id
                               ORDER BY vf.seq_region_id, vf.seq_region_start
                               $LIMIT},
                            {'mysql_use_result' => 1});


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
      #      $update_sth->execute($up_seq, $dn_seq, $up_sr_start, $up_sr_end,
      #                          $dn_sr_start, $dn_sr_end, $sr_id, $sr_strand);
      print STDERR "updated\n";
    } else {
      print STDERR "not updated\n";
    }
  }

  $sth->finish();
  $update_sth->finish();
}



sub variation_group_feature {

}


sub transcript_variation {

}



sub debug {
  print STDERR @_,"\n";
}
