use strict;
use warnings;

#generic object for the dbSNP data. Contains the general methods to dump the data into the new Variation database. Any change in the methods
# will need to overload the correspondent method in the subclass for the specie

package SARA::GetSARA;

use Bio::Index::Fastq;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use ImportUtils qw(dumpSQL debug create_and_load load);
use Time::HiRes qw(tv_interval gettimeofday);

#creates the object and assign the attributes to it (connections, basically)
sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($dbCore, $dbVar, $dbSara, $tmp_dir, $tmp_file, $species) =
        rearrange([qw(DBCORE DBVAR DBSARA TMPDIR TMPFILE SPECIES)],@_);

  return bless {'dbSara' => $dbSara,
		'dbCore' => $dbCore,
		'dbVar' => $dbVar, ##this is a dbconnection
		'tmpdir' => $tmp_dir,
		'tmpfile' => $tmp_file,
		'species' => $species,
	       }, $class;
}

#main and only function in the object that dumps all dbSNP data 
sub get_sara{ 

  my $self = shift;
  my $var_dbname = ($self->{'dbVar'}) ? $self->{'dbVar'}->dbname : "var_dbname";

  my $sth = $self->{'dbSara'}->prepare(qq{SELECT distinct seq_region_id
                                          FROM $var_dbname.variation_feature
                              });
   my ($seq_region_id,%rec_seq_region_id);
   $sth->execute();
   $sth->bind_columns(\$seq_region_id);

   while ($sth->fetch()) {
     $rec_seq_region_id{$seq_region_id}=1;
   }

  my $queue_hugemem = "-q hugemem -R'select[mem>2000] rusage[mem=2000]'";
  my $queue_linux64 = "-q normal -R'select[type==LINUX64]'";
  
  #$self->update_feature_table();
  #my $t0 = [gettimeofday];
  #debug("Time to start get_query_snp_pos $t0");
  #$self->get_query_snp_pos(\%rec_seq_region_id,$var_dbname);
  #my $t1 = [gettimeofday];
  #debug("Time to start get_flanking_seq_qual $t1");
  #$self->get_flanking_seq_qual(\%rec_seq_region_id,$var_dbname,$queue_linux64);
  #print "Time to run GetSARA ", tv_interval($t0,$t1),"\n";
  #$self->make_genotype_allele_table(\%rec_seq_region_id,$var_dbname);
  $self->insert_allele_gtype_to_vardb(\%rec_seq_region_id,$var_dbname);
  #$self->read_coverage($var_dbname);
}

sub update_feature_table {
   my $self = shift;

   #cigar line has query_start and query_end reversed when query_strand=-1, but we don't want this, we reverse them back so we can get right base and right quality values (we can alse reverse sequence but we can't reverse quality data)

   $self->{'dbSara'}->do(qq{UPDATE ssahaSNP_feature sf, query_match_length ql 
                           SET sf.query_start_new = ql.length-sf.query_end,
                               sf.query_end_new = ql.length-sf.query_start
                           WHERE sf.query_strand =-1
                           AND sf.query_name=ql.query_name
                          });

}

sub get_query_snp_pos {
  my $self = shift;
  my $rec_seq_region_id = shift;
  my $var_dbname = shift;
  my $tmp_dir = $self->{'tmpdir'};
  my $tmp_file = $self->{'tmpfile'};


  my $species = $self->{'species'};

  foreach my $seq_region_id (keys %$rec_seq_region_id) {
    my $call = "bsub -q normal -J $var_dbname\_ssaha_feature_job_$seq_region_id -o $tmp_dir/ssaha_feature_out\_$seq_region_id /usr/local/ensembl/bin/perl parallel_sara_feature.pl -species $species -seq_region_id $seq_region_id -job snp_pos -tmpdir $tmp_dir -tmpfile $tmp_file";
    print "call is $call\n";
    system($call);
  }

  my $call1 = "bsub -q normal -K -w 'done($var_dbname\_ssaha_feature_job*)' -J waiting_process sleep 1"; #waits until all variation features have finished to continue
  system($call1);

}

sub get_flanking_seq_qual {

  my $self = shift;
  my $rec_seq_region_id = shift;
  my $var_dbname = shift;
  my $queue = shift;

  my $tmp_dir = $self->{'tmpdir'};
  my $tmp_file = $self->{'tmpfile'};

  my $species = $self->{'species'};
  my $count;
  foreach my $seq_region_id (keys %$rec_seq_region_id) {
    my $call = "bsub $queue -J $var_dbname\_ssaha_flank_job_$seq_region_id -o $tmp_dir/ssaha_flank_out\_$seq_region_id /usr/local/ensembl64/bin/perl parallel_sara_feature.pl -species $species -seq_region_id $seq_region_id -job flank -tmpdir $tmp_dir -tmpfile $tmp_file";
    $count++;
    print "call is $call $count\n";
    system($call);
    sleep(1200);
  }

  my $call1 = "bsub -q normal -K -w 'done($var_dbname\_ssaha_flank_job*)' -J waiting_process sleep 1"; #waits until all variation features have finished to continue
  system($call1);
}

sub make_genotype_allele_table {

  my $self = shift;
  my $rec_seq_region_id = shift;
  my $var_dbname = shift;
  my $tmp_dir = $self->{'tmpdir'};
  my $tmp_file = $self->{'tmpfile'};

  my $species = $self->{'species'};

  foreach my $seq_region_id (keys %$rec_seq_region_id) {
    my $call = "bsub -q bigmem -R'select[mem>5000] rusage[mem=5000]' -J $var_dbname\_ssaha_gtype_job_$seq_region_id -o $tmp_dir/ssaha_gtype_out\_$seq_region_id /usr/local/ensembl/bin/perl parallel_sara_feature.pl -species $species -seq_region_id $seq_region_id -job gtype_allele -tmpdir $tmp_dir -tmpfile $tmp_file";
    print "call is $call\n";
    system($call);
  }

  my $call1 = "bsub -q normal -K -w 'done($var_dbname\_ssaha_gtype_job*)' -J waiting_process sleep 1"; #waits until all variation features have finished to continue
  system($call1);
}

sub insert_allele_gtype_to_vardb {

  my ($self, $rec_seq_region_id, $var_dbname) = @_;

  #$self->{'dbVar'}->do(qq{ALTER TABLE allele add unique index unique_allele_idx(var+iation_id,allele,sample_id)});
  
  foreach my $seq_region_id (keys %$rec_seq_region_id) {
    #insert individual allele first
    $self->{'dbSara'}->do(qq{INSERT IGNORE  INTO $var_dbname.allele (variation_id,allele,sample_id) select ga.variation_id,ga.allele,ip.population_sample_id from gtype_allele_$seq_region_id ga, $var_dbname.sample s, $var_dbname.individual_population ip where ga.sample_name = s.name and s.sample_id = ip.individual_sample_id});
    #then insert reference allele
    $self->{'dbSara'}->do(qq{INSERT IGNORE INTO $var_dbname.allele (variation_id,allele,sample_id) select ga.variation_id,ga.allele,s.sample_id from gtype_allele_$seq_region_id ga, $var_dbname.sample s where ga.sample_name = s.name and s.name like "ref%"}) ;
    $self->{'dbVar'}->do(qq{DROP INDEX unique_allele_idx ON allele});
    $self->{'dbSara'}->do(qq{INSERT INTO $var_dbname.tmp_individual_genotype_single_bp
    (variation_id,allele_1,allele_2,sample_id) select ig.variation_id,ig.allele_1,ig.allele_2,s.sample_id from tmp_individual_genotype_single_bp_$seq_region_id ig, $var_dbname.sample s where ig.sample_name = s.name});
  }

  #delete entries from variation, variation_feature and flanking_sequence tables that variation_id is not in allele and tmp_genotype_single_bp table

  $self->{'dbSara'}->do(qq{CREATE TABLE variation_old SELECT * FROM $var_dbname.variation});
  $self->{'dbSara'}->do(qq{CREATE TABLE variation_feature_old SELECT * FROM $var_dbname.variation_feature});
  $self->{'dbSara'}->do(qq{CREATE TABLE flanking_sequence_old SELECT * FROM $var_dbname.flanking_sequence});

  $self->{'dbVar'}->do(qq{DELETE FROM v USING variation v LEFT JOIN tmp_individual_genotype_single_bp t ON v.variation_id = t.variation_id
                           WHERE t.variation_id IS NULL
                          });
  $self->{'dbVar'}->do(qq{DELETE FROM vf USING variation_feature vf LEFT JOIN tmp_individual_genotype_single_bp t ON vf.variation_id = t.variation_id
                           WHERE t.variation_id IS NULL
                          });
  $self->{'dbVar'}->do(qq{DELETE FROM f USING flanking_sequence f LEFT JOIN tmp_individual_genotype_single_bp t ON f.variation_id = t.variation_id
                           WHERE t.variation_id IS NULL
                          });
}

sub read_coverage {

  my $self = shift;
  my $var_dbname = shift;
  my $tmp_dir = $self->{'tmpdir'};
  my $tmp_file = $self->{'tmpfile'};

  my $species = $self->{'species'};


  my $alignment_file ="/turing/mouse129_extra/yuan/human_celera/analysis/chimp/ALIGNMENT_FILE";
  #my $alignment_file ="/ecs2/scratch4/yuan/hum/CELERA/read_coverage/chimp/TEST";
  my $sth = $self->{'dbCore'}->prepare(qq{SELECT sr.seq_region_id,sr.name
                                 FROM   seq_region_attrib sra, attrib_type at, seq_region sr
                                 WHERE sra.attrib_type_id=at.attrib_type_id 
                                 AND at.code="toplevel"
                                 AND sr.seq_region_id = sra.seq_region_id 
                                });
  my ($seq_region_id,$seq_region_name,%rec_seq_region);
  $sth->execute();
  $sth->bind_columns(\$seq_region_id,\$seq_region_name);

  while ($sth->fetch()) {
    $rec_seq_region{$seq_region_id}=$seq_region_name;
  }

  foreach my $seq_region_id (keys %rec_seq_region) {
    my $seq_region_name = $rec_seq_region{$seq_region_id};
    my $call = "bsub -q hugemem -R'select[mem>3000] rusage[mem=3000]' -J $var_dbname\_read_coverage_job_$seq_region_id -o /$tmp_dir/read_coverage_out\_$seq_region_id /usr/local/bin/perl /nfs/acari/yuan/ensembl/src/ensembl-variation/scripts/import/parallel_sara_feature.pl -species $species -seq_region_name $seq_region_name -job read_coverage -alignment_file $alignment_file -tmpdir $tmp_dir -tmpfile $tmp_file";
    print "call is $call\n";
    system($call);
  }
  my $call1 = "bsub -q normal -K -w 'done($var_dbname\_read_coverage_job*)' -J waiting_process sleep 1"; #waits until all variation features have finished to continue
  system($call1);
}

1;
