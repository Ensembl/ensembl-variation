#! /usr/local/ensembl/bin/perl
#

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Data::Dumper;
use Bio::SeqIO;
use Getopt::Long;

# try to use named options here and write a sub usage() function
# eg -host -user -pass -port -snp_dbname -core_dbname etc
# optional chromosome name or genomic sequence file
# optional more than one genomic sequence file
# optional a directory or sequence files (for unknown placing)

our ($cuser,$cpass,$chost,$cport,$cdbname,$vuser,$vpass,$vhost,$vport,$vdbname,$chr_name);

GetOptions('chost=s'          => \$chost,
	   'cuser=s'          => \$cuser,
	   'cpass=s'          => \$cpass,
	   'cport=i'          => \$cport,
	   'cdbname=s'        => \$cdbname,
	   'vhost=s'          => \$vhost,
	   'vuser=s'          => \$vuser,
	   'vpass=s'          => \$vpass,
	   'vport=i'          => \$vport,
	   'vdbname=s'        => \$vdbname,
	   'chr_name=s'       => \$chr_name,
	  );

#added default options
$chost    ||='ecs2';
$cport    ||= 3365;
$cuser    ||='ensro';

$vhost    ||='ecs2';
$vport    ||= 3365;
$vuser    ||= 'ensro';

usage('-cdbname argument is required') if(!$cdbname);
usage('-vdbname argument is required') if(!$vdbname);


my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host =>$chost,
					      -user  =>$cuser,
					      -port  => $cport,
					      -dbname =>$cdbname);


my $vdb = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new (
							  -host =>$vhost,
							  -user  =>$vuser,
							  -port  =>$vport,
							  -dbname =>$vdbname,
							 );
print Dumper "cdb is $cdb\n";
print Dumper "vdb is $vdb\n";


generate_input_seq($cdb, $vdb, $chr_name);

sub generate_input_seq {

  my $cdb = shift;
  my $vdb = shift;
  my $chr_name = shift;

  my $var_adaptor = $vdb->get_VariationAdaptor();

  my ($sthc, $sthv, @seq_region_ids,%variation_ids);

  if ($chr_name) {
    $sthc = $cdb->dbc->prepare(qq{select seq_region_id from seq_region where name = ?}) if $chr_name;
    $sthc->execute($chr_name);

    while (my $seq_region_ids = $sthc->fetchrow_array()) {
      push @seq_region_ids, $seq_region_ids;
    }

    my $id_str = "IN (".join(',',map{"'$_'"} @seq_region_ids).");";

    $sthv = $vdb->dbc->prepare(qq{select variation_name,variation_id,seq_region_id from variation_feature where seq_region_id $id_str});
  }
  else {
    $sthv = $vdb->dbc->prepare(qq{select v.name,f.variation_id,f.seq_region_id,f.up_seq,f.down_seq from variation v, flanking_sequence f where v.variation_id=f.variation_id;});
  }
  
  $sthv->{"mysql_use_result"} = 1; # Get results as soon as ready; V.fast!
  $sthv->execute();

  my $file_size = 100000;
  my $file_count=1;
  my $i = 0;

  while (my ($variation_name,$variation_id,
             $seq_region_id,$up_seq,$down_seq) = $sthv->fetchrow_array()) {
    #print "$variation_name,$variation_id,$seq_region_id,$up_seq,$down_seq\n";
    $variation_ids{$variation_id}{'var_name'}=$variation_name;
    $variation_ids{$variation_id}{'up_seq'}=$up_seq;
    $variation_ids{$variation_id}{'down_seq'}=$down_seq;
    $variation_ids{$variation_id}{'seq_region_id'}=$seq_region_id;
    $i++;
    if( $i == $file_size ){
      my @ids = keys %variation_ids;
      print_seqs ($var_adaptor, $file_count, \@ids,\%variation_ids);
      @ids = ();
      %variation_ids = ();
      $file_count++;
      $i = 0;
    }
  }

  ###print the last batch of seq
  my @ids = keys %variation_ids;
  print_seqs ($var_adaptor, $file_count, \@ids,\%variation_ids);

}

sub print_seqs {

  my ($var_adaptor, $file_count, $ids, $variation_ids) = @_;
  my @ids = @$ids;
  my %variation_ids = %$variation_ids;

  open OUT, ">$file_count\_query_seq" or die "can't open query_seq file : $!";

  foreach my $var_id (@ids) {
    my ($flanking_sequence,$down_seq,$up_seq);
    if ($variation_ids{$var_id}{'seq_region_id'}) {
      $flanking_sequence = $var_adaptor->get_flanking_sequence($var_id);
      ($down_seq,$up_seq) = @{$flanking_sequence} if $flanking_sequence;
    }
    else {
      $up_seq = $variation_ids{$var_id}{'up_seq'};
      $down_seq = $variation_ids{$var_id}{'down_seq'};
    }

    my $seq = lc($up_seq)."W".lc($down_seq);
    print OUT ">$variation_ids{$var_id}{'var_name'}\n$seq\n";
  }

  close OUT;
}

sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: generate_input_seq.pl  <options>

options:
    -chost <hostname>    hostname of core Ensembl MySQL database (default = ecs2)
    -cuser <user>        username of core Ensembl MySQL database (default = ensro)
    -cpass <pass>        password of core Ensembl MySQL database
    -cport <port>        TCP port of core Ensembl MySQL database (default = 3365)
    -cdbname <dbname>    dbname of core Ensembl MySQL database
    -vhost <hostname>    hostname of variation MySQL database to write to
    -vuser <user>        username of variation MySQL database to write to (default = ensro)
    -vpass <pass>        password of variation MySQL database to write to
    -vport <port>        TCP port of variation MySQL database to write to (default = 3365)
    -vdbname <dbname>    dbname of variation MySQL database to write to
    -chr_name <chromosomename> chromosome name for which flanking sequences are mapping to

EOF

  die("\n$msg\n\n");
}

###example:
###./generate_input_seq.pl -cdbname homo_sapiens_core_27_35a -vdbname homo_sapiens_variation_27_35a -chr_name 21
