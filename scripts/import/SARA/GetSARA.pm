=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

use strict;


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

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
use Getopt::Long;

# configure these before running!!!
our $ssahabuild;
our $ssaha2;
our $getseqreads;
our $ssahaSNP_cons;
our $search_read;
our $fastq_dir;
our $output_dir;
our $target_dir;

die("Can't run without variables set - please edit script") unless defined($ssahaSNP_cons) and defined($getseqreads) and defined($search_read) and defined($fastq_dir) and defined($output_dir) and defined($target_dir);

#creates the object and assign the attributes to it (connections, basically)
sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($dbCore, $dbVar, $dbSara, $tmp_dir, $tmp_file, $species, $source_id) =
        rearrange([qw(DBCORE DBVAR DBSARA TMPDIR TMPFILE SPECIES SOURCE_ID)],@_);

  return bless {'dbSara' => $dbSara,
		'dbCore' => $dbCore,
		'dbVar' => $dbVar, ##this is a dbconnection
		'tmpdir' => $tmp_dir,
		'tmpfile' => $tmp_file,
		'species' => $species,
		'source_id' => $source_id,
	       }, $class;
}

#main and only function in the object that dumps all dbSNP data 
sub get_sara{ 

  my $self = shift;
  my $var_dbname = ($self->{'dbVar'}) ? $self->{'dbVar'}->dbname : "var_dbname";
  my %rec_seq_region_id;

  my $sth = $self->{'dbSara'}->prepare(qq{SELECT distinct seq_region_id
                                          FROM $var_dbname.variation_feature
                                    });
  my ($seq_region_id);
  $sth->execute();
  $sth->bind_columns(\$seq_region_id);

  while ($sth->fetch()) {
    $rec_seq_region_id{$seq_region_id}=1;
  }

  #can't built index_file for human in lustre file system, ask Tim? (need 2-3 days), so need to use turing for flanking_qual part of the comput
  my $queue_hugemem = "-q hugemem -R'select[mem>5000] rusage[mem=5000]'";
  my $queue_linux64 = "-q normal -R'select[type == LINUX64 && myens_genomics2 < 200] rusage[myens_genomics1=10]'";
  my $queue_long = "-q long -M5000000 -R'select[mem>5000] rusage[mem=5000]'";
  my $queue;
  $queue = $queue_hugemem if $self->{'tmpdir'} =~ /turing/;
  $queue = $queue_long if $self->{'tmpdir'} !~ /turing/;
  #$queue = $queue_linux64 if $self->{'tmpdir'} !~ /turing/;
  
  #my $t0 = [gettimeofday];
  #debug("Time to start get_query_snp_pos $t0");
  #$self->get_query_snp_pos(\%rec_seq_region_id,$var_dbname);
  #my $t1 = [gettimeofday];
  #debug("Time to start get_flanking_seq_qual $t1");
  #$self->make_reads_file(\%rec_seq_region_id,$queue);
  #$self->make_pileup_reads_file();
  #$self->get_flanking_seq_qual(\%rec_seq_region_id,$var_dbname,$queue);
  #print "Time to run GetSARA ", tv_interval($t0,$t1),"\n";
  #$self->make_genotype_allele_table(\%rec_seq_region_id,$var_dbname,$queue);
  #$self->merge_tables();
  #$self->insert_allele_gtype_to_vardb(\%rec_seq_region_id,$var_dbname);
  #$self->remove_faulse_variation_in_multiple_strains();
  #$self->read_coverage($var_dbname);##this is been calculated ealier
  #$self->remove_empty_tables();
  $self->remove_tables();
  #$self->merge_table_with_same_seq_region_id();
}

sub get_query_snp_pos {
  my $self = shift;
  my $rec_seq_region_id = shift;
  my $var_dbname = shift;
  my $source_id = $self->{'source_id'};
  my $tmp_dir = $self->{'tmpdir'};
  my $tmp_file = $self->{'tmpfile'};
  my $count;

  my $species = $self->{'species'};
  
  foreach my $seq_region_id (keys %$rec_seq_region_id) {
  #foreach my $seq_region_id (226032,226034) {
    my $call = "bsub -q long -R'select[myens_genomics1 < 200] rusage[myens_genomics1=10]' -J $var_dbname\_ssaha_feature_job_$seq_region_id -o $tmp_dir/ssaha_feature_out\_$seq_region_id perl parallel_sara_zm.pl -species $species -source_id $source_id -seq_region_id $seq_region_id -job snp_pos -tmpdir $tmp_dir -tmpfile $tmp_file";
    $count++;
    
    print "call is $call $count\n";
    system($call);
  }

  my $call1 = "bsub -q normal -K -w 'done($var_dbname\_ssaha_feature_job*)' -J waiting_process sleep 1"; #waits until all variation features have finished to continue
  system($call1);

}

sub get_flanking_seq_qual {
  ##use farm run smaller jobs or use turing run all jobs
  my $self = shift;
  my $rec_seq_region_id = shift;
  my $var_dbname = shift;
  my $queue = shift;
  my $species = $self->{'species'};
  my $tmp_dir = $self->{'tmpdir'};
  my $tmp_file = $self->{'tmpfile'};
  my $count;
  my $defined_table_row = 100000;
  my $reads_dir;
  $reads_dir = "[reads_dir]/reads_out";
  #job 1 somehow failed with memory 7 GB, so sent to turing, but only maxmem 3GB

  my %rec_id_name;
  #note the read_file is seq_region_name rather than seq_region_id,so use this
  #my $sth = $self->{'dbSara'}->prepare(qq{SELECT seq_region_id,seq_region_name
  #                           FROM $var_dbname.tmp_seq_region 
  #                           });  
  #my ($seq_region_id,$seq_region_name);  
  #$sth->execute();  $sth->bind_columns(\$seq_region_id,\$seq_region_name);  
  #while ($sth->fetch()) {    
  #  $rec_id_name{$seq_region_id}=$seq_region_name;  
  #}
 
  my $call; 
  foreach my $seq_region_id (keys %$rec_seq_region_id) {
  #foreach my $seq_region_id (226032,226034) {
    my $row_count_ref = $self->{'dbSara'}->db_handle->selectall_arrayref(qq{SELECT COUNT(*) FROM snp_pos_$seq_region_id});
    my $row_count = $row_count_ref->[0][0];
    if ($row_count > $defined_table_row) {
      #my $reads_file = "$reads_dir/reads_out_$rec_id_name{$seq_region_id}";
      my $reads_file = "$reads_dir/reads_out_$seq_region_id";
      $call = "bsub $queue -J $var_dbname\_ssaha_flank_job_$seq_region_id -o $tmp_dir/ssaha_flank_out\_$seq_region_id perl parallel_sara_zm.pl -species $species -seq_region_id $seq_region_id -job flank -reads_file $reads_file -tmpdir $tmp_dir -tmpfile $tmp_file";
      #$call = "bsub $queue -J $var_dbname\_ssaha_flank_job_$seq_region_id -o $tmp_dir/ssaha_flank_out\_$seq_region_id perl parallel_sara_feature.pl -species $species -seq_region_id $seq_region_id -job flank -index_file $index_file -tmpdir $tmp_dir -tmpfile $tmp_file";
    }
    else {
      my $reads_file = "$reads_dir/reads_out_small";
      $call = "bsub $queue -J $var_dbname\_ssaha_flank_job_$seq_region_id -o $tmp_dir/ssaha_flank_out\_$seq_region_id perl parallel_sara_zm.pl -species $species -seq_region_id $seq_region_id -job flank -reads_file $reads_file -tmpdir $tmp_dir -tmpfile $tmp_file";
    }
    $count++;      
    print "call is $call $count\n";      
    system($call);      
    sleep(5);  
  }

  my $call1 = "bsub -q normal -K -w 'done($var_dbname\_ssaha_flank_job*)' -J waiting_process sleep 1"; #waits until all variation features have finished to continue
  system($call1);
}

sub make_reads_file {
  ##use turing run all seq_region_ids
  my $self = shift;
  my $rec_seq_region_id = shift;
  my $queue = shift;
  my $tmp_dir = $self->{'tmpdir'};
  my $tmp_file = $self->{'tmpfile'};
  my $reads_dir = "[reads_dir]";

  my $defined_table_row = 100000;
  my (@big_seq_region_id,@small_seq_region_id);

  foreach my $seq_region_id (keys %$rec_seq_region_id) {
  #foreach my $seq_region_id (226052,226062,226046) {
    my $row_count_ref = $self->{'dbSara'}->db_handle->selectall_arrayref(qq{SELECT COUNT(*) FROM snp_pos_$seq_region_id});
    my $row_count = $row_count_ref->[0][0];
    if ($row_count > $defined_table_row) {
      push @big_seq_region_id,$seq_region_id;
    }
    else {
      push @small_seq_region_id, $seq_region_id;
    }
  }
  if ($tmp_dir !~ /turing/) {
    die "search_read job has to run on turing\n";
  }

  foreach my $seq_region_id (@big_seq_region_id) {
  #foreach my $seq_region_id (226053) {
    debug("Dumping reads names for $seq_region_id...");

    dumpSQL($self->{'dbSara'},qq{SELECT DISTINCT query_name FROM snp_pos_$seq_region_id});

    system("sort $tmp_dir/$tmp_file |uniq >$reads_dir/reads_dir/reads_name_$seq_region_id");

    #my $seq_region_id="missed";
    ###needs 50GB memeory to run search_read (for any size of reads_name)
    system("bsub -q hugemem -R'select[mem>50000] rusage[mem=50000]' -J reads_file_$seq_region_id -o $reads_dir/reads_dir/out_reads_file_$seq_region_id $getseqreads $reads_dir/reads_dir/reads_name_$seq_region_id $reads_dir/genome-reads.fastq $reads_dir/reads_dir/reads_out_$seq_region_id");
    if ($? == -1) {
      print "failed to execute: $!\n";
    }
    elsif ($? & 127) {
      printf "child died with signal %d, %s coredump\n",
	($? & 127),  ($? & 128) ? 'with' : 'without';
    }
    else {
      printf "child exited with value %d\n", $? >> 8;
    }
  }

  foreach my $seq_region_id (@small_seq_region_id) {
    debug("Dumping reads names for $seq_region_id...");    

    dumpSQL($self->{'dbSara'},qq{SELECT query_name FROM snp_pos_$seq_region_id});    
    system("sort $tmp_dir/$tmp_file |uniq >>$reads_dir/reads_dir/reads_name_small");

  }
  ###needs 50GB memeory to run search_read (for any size of reads_name)

  system("bsub -q hugemem -R'select[mem>50000] rusage[mem=50000]' -J reads_file_small -o $reads_dir/reads_dir/out_reads_file_small $getseqreads $reads_dir/reads_dir/reads_name_small $reads_dir/genome-reads.fastq $reads_dir/reads_dir/reads_out_small");
  if ($? == -1) {
    print "failed to execute: $!\n";
  }
  elsif ($? & 127) {
    printf "child died with signal %d, %s coredump\n",
    ($? & 127),  ($? & 128) ? 'with' : 'without';
  }
  else {
    printf "child exited with value %d\n", $? >> 8;
  }

  my $call1 = "bsub -q hugemem -R'select[mem>100] rusage[mem=100]' -K -w 'done(reads_file*)' -J waiting_process sleep 1"; #waits until all variation features have finished to continue
  system($call1);
  system("scp $reads_dir/reads_dir/reads_out_* [reads_out_dir]");
  #return ("$reads_dir/reads_dir/reads_out_$seq_region_id$t\_000.fastq");

}

sub make_pileup_reads_file {
  ##use turing run all chromosomes
  my $self = shift;
  my $tmp_dir = $self->{'tmpdir'};
  my $tmp_file = $self->{'tmpfile'};
  
  opendir DIR, "$output_dir" or die "Failed to open dir : $!";
  my @reads_dirs = grep /dir$/, readdir(DIR);
  print "files are @reads_dirs\n";
  
  #foreach my $read_dir (@reads_dirs) {
  foreach my $read_dir("1_dir") { 
    my ($chr) = $read_dir =~ /(\S+)\_dir/; 
    debug("chr is $chr Get reads fastq for $read_dir...");

    ###needs 50GB memeory to run search_read (for any size of reads_name)
    system("bsub -q hugemem -R'select[mem>60000] rusage[mem=60000]' -J reads_file_$read_dir -o $output_dir/$read_dir/out_reads_file_$chr $search_read -fastq 1 $output_dir/$read_dir/$chr\_read_name $output_dir/$read_dir/reads_out_$chr $fastq_dir/readname.tag $fastq_dir/*fastq");

    if ($? == -1) {
      print "failed to execute: $!\n";
    }
    elsif ($? & 127) {
      printf "child died with signal %d, %s coredump\n",
	($? & 127),  ($? & 128) ? 'with' : 'without';
    }
    else {
      printf "child exited with value %d\n", $? >> 8;
    }

    debug("Running pileup SNP...");
    if (! -e "$output_dir/$read_dir/reads_out_$chr\.fastq") {
      system("cat $output_dir/$read_dir/reads_out_$chr\_*.fastq >$output_dir/$read_dir/reads_out_$chr\.fastq");
    }
    system("bsub -q hugemem -M20000000 -R'select[mem>20000] rusage[mem=20000]' -o $output_dir/$read_dir/out_$chr\_SNP $ssahaSNP_cons $output_dir/$read_dir/$chr\_align $output_dir/$read_dir/$chr\_cigar $target_dir/$chr\.fa $output_dir/$read_dir/reads_out_$chr\.fastq");
  
  }
}

sub make_genotype_allele_table {
  ##use farm run these jobs  
  my $self = shift;
  my $rec_seq_region_id = shift;
  my $var_dbname = shift;
  my $queue = shift;
  my $tmp_dir = $self->{'tmpdir'};
  my $tmp_file = $self->{'tmpfile'};
  my $count;
  my $species = $self->{'species'};
  
  foreach my $seq_region_id (keys %$rec_seq_region_id) {
  #foreach my $seq_region_id (226055) { 
    my $call = "bsub -q normal -R'select[myens_genomics1 < 200] rusage[myens_genomics1=10]' -J $var_dbname\_ssaha_gtype_job_$seq_region_id -o $tmp_dir/ssaha_gtype_out\_$seq_region_id perl parallel_sara_zm.pl -species $species -seq_region_id $seq_region_id -job gtype_allele -tmpdir $tmp_dir -tmpfile $tmp_file";
    $count++;
    print "call is $call $count\n";
    system($call);
  }

  my $call1 = "bsub -q normal -K -w 'done($var_dbname\_ssaha_gtype_job*)' -J waiting_process sleep 1"; #waits until all variation features have finished to continue
  system($call1);
}

sub insert_allele_gtype_to_vardb {

  my ($self, $rec_seq_region_id, $var_dbname) = @_;

  $self->{'dbVar'}->do(qq{ALTER TABLE allele add unique index unique_allele_idx(variation_id,allele,sample_id)});
 
  $self->{'dbVar'}->do(qq{CREATE TABLE tmp_individual_genotype_single_bp (
                          variation_id int not null,allele_1 varchar(255),allele_2 varchar(255),sample_id int,
                          key variation_idx(variation_id),
                          key sample_idx(sample_id)
                          ) MAX_ROWS = 100000000}
   );
 
  $self->{'dbVar'}->do(qq{CREATE UNIQUE INDEX ind_genotype_idx ON tmp_individual_genotype_single_bp(variation_id,sample_id,allele_1,allele_2)});

  debug("Insert individual allele first...");

  $self->{'dbSara'}->do(qq{INSERT IGNORE  INTO $var_dbname.allele (variation_id,allele,sample_id) select ga.variation_id,ga.allele,ip.population_sample_id from gtype_allele ga, $var_dbname.sample s, $var_dbname.individual_population ip where ga.sample_name = s.name and s.sample_id = ip.individual_sample_id});
  
  debug("Then insert reference allele...");
  $self->{'dbSara'}->do(qq{INSERT IGNORE INTO $var_dbname.allele (variation_id,allele,sample_id) select ga.variation_id,ga.allele,s.sample_id from gtype_allele ga, $var_dbname.sample s where ga.sample_name = s.name and s.name like "ENS%"}) ;

  debug("insert into tmp_genotype table...");
  $self->{'dbSara'}->do(qq{INSERT IGNORE INTO $var_dbname.tmp_individual_genotype_single_bp
    (variation_id,allele_1,allele_2,sample_id) select ig.variation_id,ig.allele_1,ig.allele_2,s.sample_id from tmp_individual_genotype_single_bp ig, $var_dbname.sample s where ig.sample_name = s.name});


  $self->{'dbVar'}->do(qq{DROP INDEX unique_allele_idx ON allele});
  $self->{'dbVar'}->do(qq{DROP INDEX ind_genotype_idx ON tmp_individual_genotype_single_bp});


  #delete entries from variation, variation_feature and flanking_sequence tables that variation_id is not in allele and tmp_genotype_single_bp table

  #$self->{'dbSara'}->do(qq{CREATE TABLE variation_old SELECT * FROM $var_dbname.variation});
  #$self->{'dbSara'}->do(qq{CREATE TABLE variation_feature_old SELECT * FROM $var_dbname.variation_feature});
  #$self->{'dbSara'}->do(qq{CREATE TABLE flanking_sequence_old SELECT * FROM $var_dbname.flanking_sequence});

  $self->{'dbVar'}->do(qq{CREATE TABLE uniq_var_id_tmp_gtype SELECT DISTINCT variation_id FROM tmp_individual_genotype_single_bp});
  $self->{'dbVar'}->do(qq{ALTER TABLE uniq_var_id_tmp_gtype ADD INDEX variation_idx(variation_id)});

  debug("delete from variation table...");
  $self->{'dbVar'}->do(qq{DELETE FROM v USING variation v LEFT JOIN uniq_var_id_tmp_gtype u ON v.variation_id = u.variation_id
                            WHERE u.variation_id IS NULL
                           });
  debug("delete from variation_feature table...");
  $self->{'dbVar'}->do(qq{DELETE FROM vf USING variation_feature vf LEFT JOIN uniq_var_id_tmp_gtype u ON vf.variation_id = u.variation_id
                            WHERE u.variation_id IS NULL
                           });
  debug("delete from flanking_sequence table...");
  $self->{'dbVar'}->do(qq{DELETE FROM f USING flanking_sequence f LEFT JOIN uniq_var_id_tmp_gtype u ON f.variation_id = u.variation_id
                            WHERE u.variation_id IS NULL
                           });
}

sub remove_faulse_variation_in_multiple_strains {

  my $self = shift;
  $self->{'dbVar'}->do(qq{CREATE TABLE varid_remove 
                          SELECT t.variation_id,group_concat(t.allele_1) as allele_1,group_concat(t.allele_2) as allele_2 
                          FROM tmp_individual_genotype_single_bp t 
                          GROUP BY variation_id
                          });
  $self->{'dbVar'}->do(qq{CREATE TABLE varid_remove1
                          SELECT * FROM varid_remove
                          WHERE allele_1=allele_2
                          });
  $self->{'dbVar'}->do(qq{CREATE TABLE varid_same_as_ref 
                         SELECT DISTINCT v.variation_id from varid_remove1 v, variation_feature vf 
                         WHERE substring(v.allele_1,1,1) like substring(vf.allele_string,1,1) 
                         AND substring(v.allele_1,3,1) like substring(vf.allele_string,1,1) 
                         AND substring(v.allele_1,5,1) like substring(vf.allele_string,1,1) 
                         AND v.variation_id=vf.variation_id 
                         AND length(v.allele_1)=5
                         });
  $self->{'dbVar'}->do(qq{INSERT INTO varid_same_as_ref
                          SELECT DISTINCT v.variation_id from varid_remove1 v, variation_feature vf
                          WHERE substring(v.allele_1,1,1) like substring(vf.allele_string,1,1)
                          AND substring(v.allele_1,3,1) like substring(vf.allele_string,1,1)
                          AND v.variation_id=vf.variation_id
                          AND length(v.allele_1)=3
                          });
  $self->{'dbVar'}->do(qq{INSERT INTO varid_same_as_ref
                         SELECT DISTINCT v.variation_id from varid_remove1 v, variation_feature vf
                         WHERE substring(v.allele_1,1,1) like substring(vf.allele_string,1,1)
                         AND v.variation_id=vf.variation_id
                         AND length(v.allele_1)=1
                         });
  $self->{'dbVar'}->do(qq{ALTER TABLE varid_same_as_ref ADD INDEX variation_id(variation_id)});
  #remove these variation_id from variation/variation_feature/flanking_sequence/allele/tmp_gtype tables
  debug("delete from variation table...");
  $self->{'dbVar'}->do(qq{DELETE FROM v USING variation v, varid_same_as_ref u 
                          WHERE v.variation_id = u.variation_id
                          });
  debug("delete from variation_feature table...");
  $self->{'dbVar'}->do(qq{DELETE FROM vf USING variation_feature vf, varid_same_as_ref u 
                          WHERE vf.variation_id = u.variation_id
                          });
  debug("delete from flanking_sequence table...");
  $self->{'dbVar'}->do(qq{DELETE FROM f USING flanking_sequence f, varid_same_as_ref u 
                          WHERE f.variation_id = u.variation_id
                          });
  debug("delete from allele table...");
  $self->{'dbVar'}->do(qq{DELETE FROM v USING allele a, varid_same_as_ref u 
                          WHERE a.variation_id = u.variation_id
                          });
  debug("delete from tmp_individual_genotype_single_bp table...");
  $self->{'dbVar'}->do(qq{DELETE FROM v USING tmp_individual_genotype_single_bp t, varid_same_as_ref u 
                          WHERE t.variation_id = u.variation_id
                          });
  debug("Drop tables varid_remove varid_remove1 varid_same_as_ref");
  #$self->{'dbVar'}->do(qq{DROP TABLES varid_remove varid_remove1 varid_same_as_ref});

}
sub read_coverage {
  #not used here
  my $self = shift;
  my $var_dbname = shift;
  my $tmp_dir = $self->{'tmpdir'};
  my $tmp_file = $self->{'tmpfile'};

  my $species = $self->{'species'};


  my $alignment_file ="[alignment_file]";
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
    my $call = "bsub -q bigmem -R'select[mem>3000] rusage[mem=3000]' -J $var_dbname\_read_coverage_job_$seq_region_id -o /$tmp_dir/read_coverage_out\_$seq_region_id perl parallel_sara_feature.pl -species $species -seq_region_name $seq_region_name -job read_coverage -alignment_file $alignment_file -tmpdir $tmp_dir -tmpfile $tmp_file";
    print "call is $call\n";
    system($call);
  }
  my $call1 = "bsub -q normal -K -w 'done($var_dbname\_read_coverage_job*)' -J waiting_process sleep 1"; #waits until all variation features have finished to continue
  system($call1);
}

sub merge_tables {

  my ($self) = @_;

  #foreach my $table ("combine_feature","snp_pos","flanking_qual","gtype_allele","failed_flanking_qual","failed_gtype","tmp_individual_genotype_single_bp") {
  foreach my $table ("gtype_allele","failed_gtype","tmp_individual_genotype_single_bp") {
    debug("Merge table $table...");
    my $single_table_ref = $self->{'dbSara'}->db_handle->selectall_arrayref(qq{SHOW tables like "$table\_%"});
    my @tables = map {$_->[0] } @$single_table_ref;

    foreach my $t(@tables) {
      if ($t eq $table) {
	$self->{'dbSara'}->do(qq{DROP TABLE $table});
      }
    }
    $self->{'dbSara'}->do(qq{CREATE TABLE  $table like $tables[0]});

    my $table_names = join ',',@tables;
    #print "table_names is $table_names\n";
    $self->{'dbSara'}->do(qq{ALTER TABLE $table engine = merge union($table_names) insert_method=last});
  }
}

sub remove_tables {

  my ($self) = @_;

  #foreach my $table ("combine_feature","snp_pos","failed_flanking_qual","flanking_qual","gtype_allele","failed_gtype","tmp_individual_genotype_single_bp") {
  foreach my $table ("gtype_allele","failed_gtype","tmp_individual_genotype_single_bp") { 
    my $single_table_ref = $self->{'dbSara'}->db_handle->selectall_arrayref(qq{SHOW tables like "$table%"});
    my @tables = map {$_->[0] } @$single_table_ref;

    my $table_names = join ',',@tables;
    print "table_names is $table_names\n";
    $self->{'dbSara'}->do(qq{DROP TABLES $table_names});
  }
}

sub remove_empty_tables {

  my ($self) = @_;

  #foreach my $table ("combine_feature","snp_pos","failed_flanking_qual","flanking_qual","gtype_allele","failed_gtype","tmp_individual_genotype_single_bp") {
  foreach my $table ("flanking_qual","failed_flanking_qual") { 
    my $single_table_ref = $self->{'dbSara'}->db_handle->selectall_arrayref(qq{SHOW tables like "$table%"});
    my @tables = map {$_->[0] } @$single_table_ref;
    my @empty_tables;
    foreach my $table (@tables) {
      my $table_row_ref = $self->{'dbSara'}->db_handle->selectall_arrayref(qq{SELECT COUNT(*) FROM $table});
      #print "table is $table\n";
      my $table_row = $table_row_ref->[0][0];
      if (! $table_row) {
	push @empty_tables, $table;
      }
    }
    my $table_names = join ',',@empty_tables;
    print "empty_table_names is $table_names\n";
    $self->{'dbSara'}->do(qq{DROP TABLES $table_names});
  }
}

sub merge_table_with_same_seq_region_id {

  my ($self) = @_;

  #foreach my $table ("combine_feature","snp_pos","flanking_qual","gtype_allele","failed_qual","failed_gtype","tmp_individual_genotype_single_bp") {
  foreach my $table ("flanking_qual","failed_qual") {
    my $single_table_ref = $self->{'dbSara'}->db_handle->selectall_arrayref(qq{SHOW tables like "$table\_%"});
    my @tables = map {$_->[0] } @$single_table_ref;
    my $main_tab_name;
    $self->{'dbSara'}->do(qq{CREATE TABLE IF NOT EXISTS $main_tab_name like $tables[0]});
    foreach my $tab (@tables) {
      if ($tab =~ /(.*\_.*).*\_(\d+)(\_\d+)$/) {
	$main_tab_name = $1.$3;
	$self->{'dbSara'}->do(qq{CREATE TABLE IF NOT EXISTS $main_tab_name like $tab});
	print "table is $tab and main table is $main_tab_name\n";
	#$self->{'dbSara'}->do(qq{INSERT INTO $main_tab_name select * from $tab});
	#$self->{'dbSara'}->do(qq{DROP TABLE $tab});
      }
    }
  }
}

1;
