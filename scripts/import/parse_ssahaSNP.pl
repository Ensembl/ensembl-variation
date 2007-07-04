#! /usr/local/bin/perl
#
#bsub -q hugemem -R'select[mem>10000] rusage[mem=10000]' -o output_parse_ssahaSNP ./parse_ssahaSNP.pl -species rat -input_file /ecs4/scratch3/yuan/rat/CELERA/out/celera_rat_map1-copy20-idt92-cover100-match80_20060703.out -align_file /ecs4/scratch3/yuan/rat/CELERA/out/ALIGNMENT_file -tmpdir /ecs4/scratch3/yuan/rat/CELERA/out

use strict;
#use DBH;
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

# try to use named options here and write a sub usage() function
# eg -host -user -pass -port -snp_dbname -core_dbname etc
# optional chromosome name or genomic sequence file
# optional more than one genomic sequence file
# optional a directory or sequence files (for unknown placing)


our ($species, $input_file, $output_file, $align_file, $TMP_DIR, $TMP_FILE);

GetOptions('species=s'    => \$species,
	   'input_file=s' => \$input_file,
	   'output_file=s' => \$output_file,
	   'align_file=s' => \$align_file, #for read_coverage
	   'tmpdir=s'     => \$ImportUtils::TMP_DIR,
	   'tmpfile=s'    => \$ImportUtils::TMP_FILE,
	  );
my $registry_file ||= $Bin . "/ensembl.registry";

usage('-species argument is required') if(!$species);

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

Bio::EnsEMBL::Registry->load_all( $registry_file );
#my  $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species,'core','slice');
my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbCore = $cdb->dbc;
my $dbVar = $vdb->dbc;
my $slice_adaptor = $cdb->get_SliceAdaptor;

my %strain_name = ('GK/Ox'  =>1,
                   'SS/Jr'  =>1,
                   'SHRSP/mdc' =>1,
                   'WKY/mdc' =>1,
 		   'BN/Crl' =>1,
                  );

#&check_tri_allelic();
&create_vdb();
#&read_coverage();


sub check_tri_allelic {
  my (%allele,%rec_line,%rec_strain,%rec_var_uniq);

  open IN, "$input_file" or die "can't open input_file:$!!";
  open OUT, ">$output_file" or die "can't open output_file:$!!";
  open OUT1, ">$output_file\_all" or die "can't open output_file:$!!";
  open TMP, ">$TMP_DIR/$TMP_FILE" or die "can't open tmp_file:$!!";

  debug("reading parse_SNP file...");
  
  while (<IN>) {
    next if /Warning/i;
    if (/^SNP|^\s+\d+\s+\d+\s+\d+\s+/) {
      s/^SNP\:\s+|^\s+|\[|\]//g;
      my $line = $_;
      #SNP: 0 0 1247 1 1-1-16827091 G/A G41P615599RD4.T0 320 [ WIBR ] 4 1 2
      #SNP: 0 1 17289 0 220498186 1 1-1-267910886 A/T SHRSPa130g10.r1 881 [ SHRSP ] 1 1
      #SNP: 17 0 76064 2 1-1-267910886 C/A gnl|ti|912259493 -557 [ SD ] 6 2 6
      #     17 0 76064 2 1-1-267910886 C/A gnl|ti|928066838 308 [ SD ] 6 2 6

      #Read SHRSPa130g10.r1 is from strain SHRSP. Total read coverage at position
      #220498186 is 1 and coverage of srain SHRSP is also 1.

      my ($null1,$null2,$ref_pos,$null3,$subj_name,$allele_string,$read_name,$read_pos,$strain_name,$tot_cover_num,$num_strain_snp,$strain_tot_coverage) = split;
      next if $tot_cover_num ==0;
      my ($ref_allele,$read_allele) = split /\//, $allele_string;
      $strain_name =~ s/\[|\s+\]//g;
      my ($chr,$chr_start,$chr_end) = split /\-/, $subj_name;
      my $uniq_pos = "$chr\-$ref_pos";

      $allele{$uniq_pos}{$strain_name}{$ref_allele}++;
      $allele{$uniq_pos}{$strain_name}{$read_allele}++;
      $rec_line{$uniq_pos}{$strain_name} = $line;
      $rec_strain{$strain_name}=1;
    }
  }
   my $ref_strain = "refstrain";

   my $strain_sample_pop_desc = "population including 2 individuals from wgs and nisc";
   my $strain_sample_pop_name = "ENSEMBL:Platypus";
   my $strain_sample_pop_size = 2;
   my $individual_type_id = 2;
   debug("Insert into population and sample tables");

   $dbVar->do(qq{INSERT INTO sample (name,description) values ("$ref_strain","reference population sample")});
   my $population_ref_sample_id = $dbVar->db_handle->{'mysql_insertid'};
   $dbVar->do(qq{INSERT INTO population (sample_id) values ($population_ref_sample_id)});

   $dbVar->do(qq{INSERT INTO sample (name,size,description) values ("$strain_sample_pop_name",$strain_sample_pop_size,"$strain_sample_pop_desc")});
   my $population_strain_sample_id = $dbVar->db_handle->{'mysql_insertid'};
   $dbVar->do(qq{INSERT INTO population (sample_id) values ($population_strain_sample_id)});

   foreach my $strain_name (keys %rec_strain) {
     $dbVar->do(qq{INSERT INTO sample (name,description) values ("$strain_name","individual strain")});
     my $sample_id = $dbVar->db_handle->{'mysql_insertid'};
    
     $dbVar->do(qq{INSERT INTO individual (sample_id,gender,individual_type_id) values ($sample_id,"Unknown", $individual_type_id)});
     $dbVar->do(qq{INSERT INTO individual_population (individual_sample_id,population_sample_id) values ($sample_id,$population_strain_sample_id)}) if ($strain_name ne $ref_strain) ;
   }

  debug("Calculate tri-alleleic alleles");
  my $count=0;
  my $count1=0;
  foreach my $uniq_pos (keys %{allele}) {
    $count1++;
    my ($chr,$chr_start) = split /\-/, $uniq_pos;
    foreach my $strain_name (keys %{$allele{$uniq_pos}}) {
      my $line = $rec_line{$uniq_pos}{$strain_name};
      my ($null1,$null2,$ref_pos,$null3,$subj_name,$allele_string,$read_name,$read_pos,$strain_name,$tot_cover_num,$num_strain_snp,$strain_tot_coverage) = split /\s+/,$line;
      my @alleles = keys %{$allele{$uniq_pos}{$strain_name}};
      if (@alleles >2 ) {
	print TMP "$uniq_pos\t$strain_name\t@alleles\n";
      }
      else {
	print OUT1 "$count1\t$uniq_pos\t$chr\t$chr_start\t$strain_name\t$allele_string\t$read_name\t$read_pos\t$tot_cover_num\t$num_strain_snp\t$strain_tot_coverage\n";
          if (! $rec_var_uniq{$uniq_pos}) {
  	  $rec_var_uniq{$uniq_pos}=++$count;
  	  print OUT "$count\t$uniq_pos\t$chr\t$chr_start\t$strain_name\t$allele_string\t$read_name\t$read_pos\t$tot_cover_num\t$num_strain_snp\t$strain_tot_coverage\n";
  	}
      }
    }
  }
  create_and_load($dbVar,"tri_allelic_alleles","uniq_pos *","strain_name","alleles");
  system("mv $output_file $TMP_DIR/$TMP_FILE");
  create_and_load($dbVar,"ssaha_snp_out","count i*","uniq_pos *","chr *","target_start i*","strain_name","allele_string","read_name *","read_pos i*","tot_cover_num i","strain_snp_num i","strain_tot_coverage i");
  system("mv $output_file\_all $TMP_DIR/$TMP_FILE");
  create_and_load($dbVar,"ssaha_snp_out_all","count i*","uniq_pos *","chr *","target_start i*","strain_name","allele_string","read_name *","read_pos i*","tot_cover_num i","strain_snp_num i","strain_tot_coverage i");
}

sub create_vdb {
 
  debug("Creating database now...");
  my (%rec_var,%rec_strain,%sample_id,%allele);
  my $source_name = "ENSEMBL";
  my $var_pre_name = "ENSEMBL";
  my $length_pre_name = length($var_pre_name);
  my $variation_name = "ENSOANSNP";


    $dbVar->do(qq{INSERT INTO source (name) values ("$source_name")});
    my $source_id = $dbVar->{'mysql_insertid'};
    debug("Insert into variation table...");

    $dbVar->do(qq{ALTER TABLE variation add column internal_name varchar(50)});
    $dbVar->do(qq{INSERT INTO variation (source_id,name,internal_name) select s.source_id as source_id,concat("$variation_name",snp.count),concat("$var_pre_name\_",snp.uniq_pos) as name from source s, ssaha_snp_out snp});


   dumpSQL($dbCore,qq{SELECT sr.seq_region_id,sr.name
                     FROM   seq_region_attrib sra, attrib_type at, seq_region sr
                     WHERE sra.attrib_type_id=at.attrib_type_id 
                     AND at.code="seqlevel"
                     AND sr.seq_region_id = sra.seq_region_id 
        });

#    dumpSQL($dbCore,qq{SELECT sr.seq_region_id,sr.name
#                       FROM seq_region sr
#                       WHERE coord_system_id = 4
#                     });
   create_and_load($dbVar,"tmp_seq_region","seq_region_id i*","seq_region_name *");

  debug("Insert into variation_feature table...");
  $dbVar->do(qq{INSERT INTO variation_feature (seq_region_id,seq_region_start,seq_region_end,seq_region_strand,variation_id,allele_string,variation_name,map_weight,flags,source_id,validation_status,consequence_type)
              SELECT ts.seq_region_id,snp.target_start,snp.target_start,1,v.variation_id,snp.allele_string,v.name,1,"genotyped",s.source_id,NULL,"INTERGENIC"
              FROM tmp_seq_region ts,variation v,ssaha_snp_out snp, source s
              WHERE s.name = "$source_name"
              AND snp.chr = ts.seq_region_name
              AND substring(v.internal_name,$length_pre_name+2) = snp.uniq_pos
              });

  debug("Insert into flanking_sequence table...");
  $dbVar->do(qq{INSERT INTO flanking_sequence (variation_id,up_seq,down_seq,up_seq_region_start,up_seq_region_end,down_seq_region_start,down_seq_region_end,seq_region_id,seq_region_strand)
              SELECT vf.variation_id,NULL,NULL,vf.seq_region_start-101,vf.seq_region_start-1,vf.seq_region_end+1,vf.seq_region_end+101,vf.seq_region_id,vf.seq_region_strand
              FROM variation_feature vf
               });

}

sub read_coverage {

  debug("reading read coverage data...");
  open IN, "$align_file" or die "can't open alignment file\n";
  while (<IN>) {
    #ALIGNMENT 764 gnl|ti|900540082 10-1-110718848 5 794 39881618 39882407 F 790 99.75 795 SD
    my @all = split;
    if (@all == 13) {
      my ($null1,$score,$read_name,$target_name,$read_start,$read_end,$target_start,$target_end,$dir,$match_length,$identity,$read_length,$strain_name) = split;
      my ($chr) = split /\-/, $target_name;
      ($target_end, $target_start) = ($target_start, $target_end) if ($target_end < $target_start);
      my $my_score = $score/$read_length if $read_length != 0;
      if ($read_length == 0) {
	print "read_length=0 $_\n";
      }
      if ($my_score > 0.5) {
        if ($dbVar->dbname =~ /platypus/) {
          if ($chr =~ /Contig/) {
	    open MAP, ">>$TMP_DIR/Contig\.mapped" or die "can't open mapping file:$!";
          }
          elsif ($chr =~ /Ultra/) {
            open MAP, ">>$TMP_DIR/Ultra\.mapped" or die "can't open mapping file:$!";
          }
          else {
            open MAP, ">>$TMP_DIR/Chr\.mapped" or die "can't open mapping file:$!";
          }
          print MAP "$strain_name\t$chr\t$target_start\t$target_end\n";
        }
        else {
          open MAP, ">>$TMP_DIR/$chr\.mapped" or die "can't open mapping file:$!";
          print MAP "$strain_name\t$target_start\t$target_end\n";
        }
      }
    }
  }
}


