#! /usr/local/bin/perl

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Data::Dumper;
use Bio::SeqIO;
use FindBin qw( $Bin );
use Getopt::Long;
use ImportUtils qw(dumpSQL debug create_and_load load);

our ($species, $input_file, $output_file, $align_file, $TMP_DIR, $TMP_FILE, $individual);

GetOptions('species=s'    => \$species,
           'input_file=s' => \$input_file,
           'output_file=s' => \$output_file,
           'align_file=s' => \$align_file, #for read_coverage
           'tmpdir=s'     => \$ImportUtils::TMP_DIR,
           'tmpfile=s'    => \$ImportUtils::TMP_FILE,
           'individual=s' => \$individual
          );
#my $registry_file ||= $Bin . "/ensembl.registry2";

#usage('-species argument is required') if(!$species);

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

my $registry_file;
$registry_file ||= $Bin . "/ensembl.registry";
#$registry_file ||= $Bin . "/ensembl.registry" if ($TMP_FILE =~ /venter/);
#$registry_file ||= $Bin . "/ensembl.registry2" if ($TMP_FILE =~ /watson/);

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbCore = $cdb->dbc;
my $dbVar = $vdb->dbc;
my $slice_adaptor = $cdb->get_SliceAdaptor;
my $buffer = {};

#make_pileup_reads_file();
parse_pileup_snp_file();
#create_vdb();
#PAR_regions();

sub make_pileup_reads_file {
  ##use turing run all chromosomes
  ##before run, run this first: /nfs/team71/psg/zn1/bin/tag.pl try.fastq > readname.tag
  my $self = shift;
  my $tmp_dir = $self->{'tmpdir'};
  my $tmp_file = $self->{'tmpfile'};
  my $fastq_dir = "/turing/mouse129_extra/yuan/watson/fastq";
  my $output_dir = "/turing/mouse129_extra/yuan/watson/output_dir";
  my $target_dir = "/turing/mouse129_extra/yuan_watson/target_dir";
  
  opendir DIR, "$output_dir" or die "Failed to open dir : $!";
  my @reads_dirs = grep /dir$/, readdir(DIR);
  print "files are @reads_dirs\n";
  
  #foreach my $read_dir (@reads_dirs) {
  foreach my $read_dir("17_dir") { 
    my ($chr) = $read_dir =~ /(\S+)\_dir/; 
    print "chr is $chr Get reads fastq for $read_dir...\n";
=head
    ###needs 50GB memeory to run search_read (for any size of reads_name)
    system("bsub -q hugemem -R'select[mem>60000] rusage[mem=60000]' -J reads_file_$read_dir -o $output_dir/$read_dir/out_reads_file_$chr /nfs/team71/psg/zn1/bin/search_read -fastq 1 $output_dir/$read_dir/$chr\_read_name $output_dir/$read_dir/reads_out_$chr $fastq_dir/readname.tag $fastq_dir/*fastq");

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
=cut
    print "Running pileup SNP...\n";
    if (! -e "$output_dir/$read_dir/reads_out_$chr\.fastq") {
      system("cat $output_dir/$read_dir/reads_out_$chr\_*.fastq >$output_dir/$read_dir/reads_out_$chr\.fastq");
    }
    chdir("$output_dir/$read_dir");
    #system("bsub -q hugemem -M20000000 -R'select[mem>20000] rusage[mem=20000]' -o $output_dir/$read_dir/out_$chr\_SNP /nfs/team71/psg/zn1/src/ssahaSNP2/parse_SNP/view_ssahaSNP/cons/ABI/ssahaSNP_cons $chr\_align $chr\_cigar ../../target_dir/$chr\.fa reads_out_$chr\.fastq");
    system("bsub -q hugemem -M20000000 -R'select[mem>20000] rusage[mem=20000]' -o $output_dir/$read_dir/out_$chr\_pileup /nfs/team71/psg/zn1/src/ssahaSNP2/parse_SNP/view_ssahaSNP/cons/cons/ssahaSNP_cons $chr\_align $chr\_cigar ../../target_dir/$chr\.fa reads_out_$chr\.fastq");
  }
}

sub parse_pileup_snp_file {

  my $input_file = "/turing/mouse129_extra/yuan/watson/output_dir/SNP_file";
  #my $input_file = "/lustre/work1/ensembl/yuan/SARA/watson/output_dir/TEST_SNP";
  my $individual_name = "watson";
  my $snp_count;

  open IN, "$input_file" or die "can't open snp_file : $!";

  while (<IN>) {
    my (%base_big,%base_small,%all_letter);
    my ($null,$target_name,$pos,$num_reads,$ref_base,$snp_base,$reads_bases,$num_same_ref,$num_diff_ref,$num_same_ref_cap,$num_diff_ref_cap) = split;
    $num_same_ref_cap ||='\N';
    $num_diff_ref_cap ||='\N';
    my $allele_string = "$ref_base/$snp_base";
    my ($chr,$start,$end) = split /\-/, $target_name;
    my $seq_region_start = $pos + $start -1; #most chr start =1, but haplotype chr start >1
    #print_buffered($buffer, "$TMP_DIR/pileup_file",join "\t", $target_name,$chr,$seq_region_start,$num_reads,$ref_base,$snp_base,$reads_bases,$num_same_ref,$num_diff_ref,$num_same_ref_cap,$num_diff_ref_cap,"\n");

    $reads_bases =~ s/\-+//g;
    my @bases = split "", $reads_bases;
    map {/[A-Z]/ ? $base_big{$_}++ : $base_small{$_}++} @bases;

    #if without this if{}, i.e keep the single capital base, it will reduce 25746 entries
    if (scalar keys %base_big >2) {
      #possible tri-alleleic alleles
      foreach my $B (keys %base_big) {
	if ($base_big{$B}==1) {
	  delete $base_big{$B};
	}
      }
    }

    if (scalar keys %base_big >2) {
      #get rid of tri-alleleic alleles
      print_buffered($buffer,"$TMP_DIR/failed_file", join ("\t", $chr,$seq_region_start,$allele_string,join (",",keys %base_big))."\n");
      next;
    }
    elsif (scalar keys %base_big == 2) {
      my @keys = keys %base_big;
      if ($keys[0] ne $snp_base and $keys[1] ne $snp_base) {
	#get rid of tri-alleleic alleles
	print_buffered($buffer, "$TMP_DIR/failed_file",join ("\t", $chr,$seq_region_start,$allele_string,join (",",keys %base_big))."\n");
	next;
      }
      my $num1 = $base_big{$keys[0]};
      my $num2 = $base_big{$keys[1]};
      ($num1,$num2) = ($num2,$num1) if $num2 < $num1 ;
      my $ratio = $num1/$num2;
      my ($allele_1,$allele_2) = @keys;
      #if change 0.25 to 0.2, give 48745 more entries which is 1.3% 48745/3818462
      if ($ratio >=0.25 and $allele_string =~ /$allele_1/i and $allele_string =~ /$allele_2/i) {
	#if have both alleles capital, they need to satisfy this ratio
	$snp_count++;
	print_buffered($buffer,"$TMP_DIR/output_file", join ("\t", $snp_count,$chr,$seq_region_start,$allele_string,$allele_1,$allele_2,$individual_name)."\n");
      }
      else {
	#if have two alleles, but not satisfy above ratio, get the allele with big number as genotypes
	#also if only one allele, then this allele has to equal snp_base, i.e different from ref_base to quality for a snp
	#177515 difference if use allele_string rather then snp_base which is 4.6% 177515/3818462
	#add $ratio < 0.34 will reduce 5k entries which is 0.01% 5002/3818462
	#add $ratio <= 0.25 will reduce another 2k entries it's safer to use <=0.25, i.e at least 1 to 4 to ignore 1
	if ($ratio <= 0.25 and $base_big{$allele_1}==$num2 and $snp_base =~ /$allele_1/i) {
	  $snp_count++;
	  print_buffered($buffer, "$TMP_DIR/output_file",join ("\t", $snp_count,$chr,$seq_region_start,$allele_string,$allele_1,$allele_1,$individual_name)."\n");
	}
	elsif ($ratio <= 0.25 and $base_big{$allele_2}==$num2 and $snp_base =~ /$allele_2/i) {
	  $snp_count++;
	  print_buffered($buffer, "$TMP_DIR/output_file",join ("\t", $snp_count,$chr,$seq_region_start,$allele_string,$allele_2,$allele_2,$individual_name)."\n");
	}
	else {
	  #case like G       T       AAAAAAAT
	  print_buffered($buffer, "$TMP_DIR/failed_file",join ("\t", $chr,$seq_region_start,$allele_string,join (",",keys %base_big))."\n");
	  next;
	}
      }
    }
    elsif (scalar keys %base_big < 2) {
      my @cap_letters = keys %base_big;
      my @small_letters = keys %base_small;
      #if only one big letter and it's different from ref_base, i.e is a snp (use allele_string instead of snp_base affect 1895 entries, which is 0.05% 1895/3818462
      #if use >2 instead of >=2, lose 252031 entries which is 6.6% 252031/3818462
      if ($cap_letters[0] and $base_big{$cap_letters[0]} >=2 and $snp_base =~ /$cap_letters[0]/i) {
	#if only have one capital letter, but have more than 2, use it
	$snp_count++;
	print_buffered($buffer,"$TMP_DIR/output_file", join "\t", $snp_count,$chr,$seq_region_start,$allele_string,$cap_letters[0],$cap_letters[0],$individual_name,"\n");
	next;
      }
      $reads_bases =~ tr/[A-Z]/[a-z]/;
      my @all_letters = split "", $reads_bases;
      map {$all_letter{$_}++,1}  @all_letters;

      if (scalar keys %all_letter >2) {
      #possible tri-alleleic alleles, delete if only one exist
	foreach my $B (keys %all_letter) {
	  if ($all_letter{$B}==1) {
	    delete $all_letter{$B};
	  }
	}
      }
      if (keys %all_letter >2) {
	#ignore it if have more than 2 alleles
	print_buffered($buffer,"$TMP_DIR/failed_file", join "\t", $chr,$seq_region_start,$allele_string,join(",",keys %all_letter),"\n");
	next;
      }
      my ($key1,$key2) = keys %all_letter;
      my $n1 = $all_letter{$key1};
      my $n2 = $key2 ? $all_letter{$key2} : 0;

      #make n1 is always smaller than n2
      ($n1,$n2) = ($n2,$n1) if $n2 < $n1 ;
      my $ratio = $n1/$n2 ;
      if ($n1 >2 and $ratio >=0.25 and $allele_string =~ /$key1/i and $allele_string =~ /$key2/i) {
	#if small number of alleles >2 and have correct ratio, take it
	$snp_count++;
	print_buffered($buffer,"$TMP_DIR/output_file", join "\t", $snp_count,$chr,$seq_region_start,$allele_string,uc($key1),uc($key2),$individual_name,"\n");
      }
      else {
	if ($n2 > 3) {#current watson has $n2>2
	  #if use n2>2, there are 32125 failed, if use n2>3, there are 46680 failed out of total 3818462
	  #if only one letter OR have less than three first letter, and have at least 3 of more letters, take it
	  $snp_count++;
	  #this allele must be different from ref_base to quality as a snp
	  if ($all_letter{$key1}==$n2 and $snp_base =~ /$key1/i) {
	    print_buffered($buffer,"$TMP_DIR/output_file", join "\t", $snp_count,$chr,$seq_region_start,$allele_string,uc($key1),uc($key1),$individual_name,"\n");
	  }
	  elsif ($all_letter{$key2}==$n2 and $snp_base =~ /$key2/) {
	    print_buffered($buffer,"$TMP_DIR/output_file", join "\t", $snp_count,$chr,$seq_region_start,$allele_string,uc($key2),uc($key2),$individual_name,"\n");
	  }
	}
	else {
	  #these are failed ones
	  print_buffered($buffer,"$TMP_DIR/failed_file", join "\t", $chr,$seq_region_start,$allele_string,join(",",keys %all_letter),"\n");
	}
      }
    }
  }

  print_buffered($buffer);
=head
  system("mv $TMP_DIR/failed_file $TMP_DIR/$TMP_FILE");
  create_and_load($dbVar,"tri_allelic_alleles","chr *","seq_region_start i*","alleles");
  system("mv $TMP_DIR/output_file $TMP_DIR/$TMP_FILE");
  create_and_load($dbVar,"parse_snp","snp_count i*","chr *","seq_region_start i*","allele_string","allele_1","allele_2","individual_name");
  system("mv $TMP_DIR/pileup_file $TMP_DIR/$TMP_FILE");
  create_and_load($dbVar,"pileup_snp_out","target_name","chr *","seq_region_start i*","num_reads i","ref_base","snp_base","reads_bases","num_same_ref i","num_diff_ref i","num_same_ref_cap i","num_diff_ref_cap i");
=cut
}

sub create_vdb {

  my %rec_strain = ("watson" => 1);
  my $individual_name = "watson";
  my $individual_type_id = 1;
  my $pop_size = keys %rec_strain;
  $pop_size ||=1;
  my $ind_pop_name = "ENSEMBL:Watson";
  my $ind_sample_pop_desc = "Population for $pop_size individual(s)";
  my $ind_sample_desc = "Individual within population $ind_pop_name";

  debug("Inserting into population, individual and sample tables");

  $dbVar->do(qq{INSERT INTO sample (name,description) values ("$ind_pop_name","$ind_sample_pop_desc")});
  my $population_ind_sample_id = $dbVar->db_handle->{'mysql_insertid'};
  $dbVar->do(qq{INSERT INTO population (sample_id) values ($population_ind_sample_id)});

  foreach my $strain_name (keys %rec_strain) {
    $dbVar->do(qq{INSERT INTO sample (name,size,description) values ("$strain_name",NULL,"$ind_sample_desc")});
    my $individual_sample_id = $dbVar->db_handle->{'mysql_insertid'};
    $dbVar->do(qq{INSERT INTO individual (sample_id,individual_type_id) values ($individual_sample_id, $individual_type_id)});
    $dbVar->do(qq{INSERT INTO individual_population (individual_sample_id,population_sample_id) values ($individual_sample_id,$population_ind_sample_id)});
  }

  debug("Creating database now...");

  my $source_name = "ENSEMBL"; #needs change everytime
  my $var_pre_name = "ENSEMBL";
  my $variation_name = "ENSSNP";###needs change for different species
  my $length_name = length($variation_name);

    $dbVar->do(qq{INSERT INTO source (name) values ("$source_name")});
    my $source_id = $dbVar->{'mysql_insertid'};
    debug("Insert into variation table...");

    $dbVar->do(qq{ALTER TABLE variation add column internal_name varchar(50)});
    $dbVar->do(qq{INSERT INTO variation (source_id,name,internal_name) select s.source_id as source_id,concat("$variation_name",snp.snp_count) as name,concat("$variation_name\_",snp.chr,"-",snp.seq_region_start) as internal_name from source s, parse_snp snp});


    dumpSQL($dbCore,qq{SELECT sr.seq_region_id, sr.name
                                FROM   seq_region_attrib sra, attrib_type at, seq_region sr
                                WHERE sra.attrib_type_id=at.attrib_type_id 
                                AND at.code="toplevel" 
                                AND sr.seq_region_id = sra.seq_region_id 
                     });
   create_and_load($dbVar,"tmp_seq_region","seq_region_id i*","seq_region_name *");

 debug("Insert into variation_feature table... NOTE ABOUT Y CHROMOSOME");
 $dbVar->do(qq{INSERT INTO variation_feature (seq_region_id,seq_region_start,seq_region_end,seq_region_strand,variation_id,allele_string,variation_name,map_weight,flags,source_id,validation_status,consequence_type)
             SELECT ts.seq_region_id,snp.seq_region_start,snp.seq_region_start,1,v.variation_id,snp.allele_string,v.name,1,"genotyped",s.source_id,NULL,"INTERGENIC"
             FROM tmp_seq_region ts,variation v,parse_snp snp, source s
             WHERE s.name = "$source_name"
             AND snp.chr = ts.seq_region_name
             AND substring(v.name,$length_name+1) = snp.snp_count
             });

  debug("Insert into variation_feature table... ABOUT Y CHROMOSOME");
  $dbVar->do(qq{create table vf_y_top select * from variation_feature where seq_region_id=226054 and seq_region_start>=1 and seq_region_start<=2709520});
  $dbVar->do(qq{insert into vf_y_top select * from variation_feature where seq_region_id=226054 and seq_region_start>=57443438 and seq_region_start<=57772954});
  $dbVar->do(qq{update vf_y_top set seq_region_id=226031});
  $dbVar->do(qq{insert into variation_feature (seq_region_id,seq_region_start,seq_region_end,seq_region_strand,variation_id,allele_string,variation_name,map_weight,flags,source_id,validation_status,consequence_type) select seq_region_id,seq_region_start,seq_region_end,seq_region_strand,variation_id,allele_string,variation_name,map_weight,flags,source_id,validation_status,consequence_type from vf_y_top});
  #Query OK, 5020 rows affected (0.10 sec)
  #Records: 5020  Duplicates: 0  Warnings: 0

  debug("Insert into flanking_sequence table...");
  $dbVar->do(qq{INSERT IGNORE INTO flanking_sequence (variation_id,up_seq,down_seq,up_seq_region_start,up_seq_region_end,down_seq_region_start,down_seq_region_end,seq_region_id,seq_region_strand)
              SELECT vf.variation_id,NULL,NULL,vf.seq_region_start-101,vf.seq_region_start-1,vf.seq_region_end+1,vf.seq_region_end+101,vf.seq_region_id,vf.seq_region_strand
              FROM variation_feature vf
               });

  debug("Insert into allele table...");
  $dbVar->do(qq{CREATE UNIQUE INDEX unique_allele_idx ON allele (variation_id,allele(2),frequency,sample_id)});
  $dbVar->do(qq{INSERT IGNORE INTO allele (variation_id,allele,sample_id)
                SELECT v.variation_id,substring(snp.allele_string,1,1) as allele,null as sample_id
                FROM variation v, parse_snp snp
                WHERE substring(v.name,$length_name+1) = snp.snp_count
                });

  $dbVar->do(qq{INSERT IGNORE INTO allele (variation_id,allele,sample_id)
                SELECT v.variation_id,substring(snp.allele_string,3,1) as allele,null as sample_id
                FROM variation v, parse_snp snp
                WHERE substring(v.name,$length_name+1) = snp.snp_count
                });
  #$dbVar->do("DROP INDEX unique_allele_idx ON allele"); need it in the following par regions

  $dbVar->do(qq{CREATE TABLE tmp_individual_genotype_single_bp (
                            variation_id int not null,allele_1 varchar(255),allele_2 varchar(255),sample_id int,
                            key variation_idx(variation_id),
                            key sample_idx(sample_id)
                            ) MAX_ROWS = 100000000
               });

  debug("Insert into tmp_individual_genotype_single_bp table...");
  $dbVar->do(qq{INSERT INTO tmp_individual_genotype_single_bp (variation_id,allele_1,allele_2,sample_id)
                SELECT  v.variation_id,snp.allele_1,snp.allele_2,s.sample_id
                FROM variation v, parse_snp snp, sample s
                WHERE substring(v.name,$length_name+1) = snp.snp_count
                AND snp.individual_name = s.name
                });
}

#PAR regions are identical regions between chr Y and X. We need to copy variations from X to Y
sub PAR_regions{

  my $ae_adaptor =  $cdb->get_AssemblyExceptionFeatureAdaptor( $species, 'Core', 'AssemblyExceptionFeature' );
  my $slice = $slice_adaptor->fetch_by_region( 'Chromosome', 'Y' );
  my @exceptions = @{ $ae_adaptor->fetch_all_by_Slice($slice) }; #get the PAR regions between chr Y->X

  my $last_variation_id = get_last_table_id("variation");
  my $last_variation_feature_id = get_last_table_id("variation_feature");
  my $last_variation_name = get_last_variation_name();

  get_variations(\@exceptions,$last_variation_id,$last_variation_feature_id,$last_variation_name ); #method to get all variations in PAR region

  #import the 3 table
#   my $sql = qq{LOAD DATA LOCAL INFILE "$TMP_DIR/variation.txt" INTO TABLE variation (variation_id,source_id,name)};
#   $dbVar->do($sql);
#   #unlink("$TMP_DIR/variation.txt");
#   $sql = qq{LOAD DATA LOCAL INFILE "$TMP_DIR/variation_feature.txt" INTO TABLE variation_feature (variation_feature_id,seq_region_id,seq_region_start,seq_region_end,seq_region_strand,variation_id,allele_string,variation_name,map_weight,flags, source_id)};
#   $dbVar->do($sql);
#   #unlink("$TMP_DIR/variation_feature.txt");
#   $sql = qq{LOAD DATA LOCAL INFILE "$TMP_DIR/flanking_sequence.txt" INTO TABLE flanking_sequence (variation_id,up_seq_region_start,up_seq_region_end,down_seq_region_start,down_seq_region_end,seq_region_id,seq_region_strand)};
#   $dbVar->do($sql);
#   #unlink("$TMP_DIR/flanking_sequence.txt");
#   get_allele_genptype();
}

sub get_variations{
  my $exceptions = shift;
  my $last_variation_id = shift;
  my $last_variation_feature_id = shift;
  my $last_variation_name = shift;

  $last_variation_id = 17310793;
  $dbVar->do(qq{CREATE TABLE IF NOT EXISTS tmp_varid (
				       variation_id int not null,
				       variation_id_new int,
				       primary key variation_idx(variation_id))});

  foreach my $exception (@$exceptions) {
    my $target_seq_region_id = $exception->slice->get_seq_region_id();
    my $target_start = $exception->start(); #Y coordinates
    my $target_end = $exception->end();
    my $seq_region_id = $exception->alternate_slice->get_seq_region_id();
    my $start = $exception->alternate_slice()->start(); #X coordinates
    my $end = $exception->alternate_slice()->end();

    my ($variation_id,$seq_region_start,$seq_region_strand,$allele_string);
    my $snp_pos; #position of the SNP relative to the beginning to the region
    my $sth = $dbVar->prepare(qq{SELECT variation_id,seq_region_start,seq_region_strand,allele_string FROM variation_feature WHERE seq_region_id = ? and seq_region_start >= ? and seq_region_end <= ? and seq_region_start <= ?});
    $sth->bind_param(1,$seq_region_id);
    $sth->bind_param(2,$start);
    $sth->bind_param(3,$end);
    $sth->bind_param(4,$end);
    $sth->execute();
    $sth->bind_columns(\$variation_id,\$seq_region_start,\$seq_region_strand,\$allele_string);

    while ($sth->fetch){
        $snp_pos = $seq_region_start - $start; #the position of the snp relative to the beginning of the PAR region
        write_file("variation.txt",$last_variation_id+1,1,"ENSSNP".($last_variation_name+1));
        write_file("variation_feature.txt",$last_variation_feature_id+1,$target_seq_region_id,$snp_pos+$target_start,$snp_pos+$target_start,$seq_region_strand,$last_variation_id+1,$allele_string,"ENSSNP".($last_variation_name+1),1,"genotyped",1);
        write_file("flanking_sequence.txt",$last_variation_id+1,$snp_pos+$target_start - 1 - 100,$snp_pos+$target_start - 1,$snp_pos+$target_start  +1, $snp_pos+$target_start + 1 + 100,$target_seq_region_id,$seq_region_strand);
        $last_variation_id++;
        $last_variation_feature_id++;
        $last_variation_name++;
	$dbVar->do(qq{INSERT INTO tmp_varid (variation_id,variation_id_new) values($variation_id,$last_variation_id)});
    }
  }
}

sub get_allele_genptype {


  $dbVar->do(qq{INSERT INTO allele (variation_id,allele,sample_id)
                SELECT t.variation_id_new as variation_id,a.allele,a.sample_id
                FROM tmp_varid t, allele a
                WHERE t.variation_id=a.variation_id
               });

  $dbVar->do(qq{INSERT INTO tmp_individual_genotype_single_bp (variation_id,allele_1,allele_2,sample_id)
                SELECT tv.variation_id_new as variation_id,tg.allele_1,tg.allele_2,tg.sample_id
                FROM tmp_varid tv, tmp_individual_genotype_single_bp tg
                WHERE tv.variation_id=tg.variation_id
               });
}

#function to return the last id used in the table tablename (the id must be called "tablename_id")
sub get_last_table_id{

    my $tablename = shift;

    my $max_id;
    my $sth = $dbVar->prepare(qq{SELECT MAX($tablename\_id) from $tablename});
    $sth->execute();
    $sth->bind_columns(\$max_id);
    $sth->fetch;
    $sth->finish();

    return $max_id if (defined $max_id);
    return 0 if (!defined $max_id);
}

sub write_file{
    my $filename = shift;
    my @values = @_;

    open FH, ">>$TMP_DIR/$filename" || die "Could not open file with information: $!\n
";
    my @a = map {defined($_) ? $_ : '\N'} @values; #to replace undefined values by \N in the file
    print FH join("\t", @a), "\n";
    close FH || die "Could not close file with information: $!\n";
    
}


#function to return the last id used in the table tablename (the id must be called "tablename_id")
sub get_last_variation_name{
    #my $dbSanger = shift;


    my $max_id;
    my $sth = $dbVar->prepare(qq{SELECT max(round(substring(name,7))) from variation});
    $sth->execute();
    $sth->bind_columns(\$max_id);
    $sth->fetch;
    $sth->finish();

    return $max_id if (defined $max_id);
    return 0 if (!defined $max_id);
}

sub read_coverage {

  debug("reading read coverage data...");
  open IN, "$align_file" or die "can't open alignment file\n";
  my $buffer={}; #buffer to store the lines to be written
  while (<IN>) {
    #ALIGNMENT 764 gnl|ti|900540082 10-1-110718848 5 794 39881618 39882407 F 790 99.75 795 SD
    my @all = split;
    my ($chr,$null2,$null3,$null4);
    if (@all == 13) {
      my ($null1,$score,$read_name,$target_name,$read_start,$read_end,$target_start,$target_end,$dir,$match_length,$identity,$read_length,$strain_name) = split;
      ($chr) = split /\-/, $target_name if $target_name =~ /\-/;
      ($null2,$null3,$chr,$null4) = split /\:/, $target_name if $target_name =~ /\:/;

      ($target_end, $target_start) = ($target_start, $target_end) if ($target_end < $target_start);
      my $my_score = $score/$read_length if $read_length != 0;
      if ($read_length == 0) {
        print "read_length=0 $_\n";
      }
      if ($my_score > 0.5) {
        my $file;
        if ($dbVar->dbname =~ /platypus/) {
          if ($chr =~ /Contig/) {
            $file = "Contig\.mapped";
          }
          elsif ($chr =~ /Ultra/) {
            $file = "Ultra\.mapped";
          }
          else {
            $file = "Chr\.mapped";
          }
        }
        else {
          if ($chr =~ /^NT/) {
            $file = "NT.mapped";
          }
          else {
            $file = "$chr.mapped";
          }
        }
        print_buffered($buffer,"$TMP_DIR/$file",join("\t",$strain_name,$target_start,$target_end)."\n");
      }
    }
  }
  print_buffered($buffer); #flush the buffer
}

sub print_buffered {
    my $buffer = shift;
    my $filename = shift;
    my $text = shift;

    local *FH;

    if( ! $filename ) {
        # flush the buffer
        foreach my $file (keys %{$buffer}){
            open( FH, ">>$file" ) or die;
            print FH $buffer->{ $file };
            close FH;
        }

    } else {
        $buffer->{ $filename } .= $text;
        if( length( $buffer->{ $filename } ) > 10_000 ) {
            open( FH, ">>$filename" ) or die;
            print FH $buffer->{ $filename };
            close FH;
            $buffer->{ $filename } = '';
        }
    }
}
