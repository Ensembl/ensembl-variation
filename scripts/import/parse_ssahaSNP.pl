#!/usr/bin/env perl
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



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


#perl parse_ssahaSNP_new.pl -species human -input_file [ssaha_ooutput_dir]/output_parse.txt -output_file [ssaha_ooutput_dir]/watson_output.txt -align_file [ssaha_ooutput_dir]/ALIGNMENT_file -tmpdir [tmpdir] -tmpfile temp_Watson.txt -individual 'Venter'

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

usage('-species argument is required') if(!$species);

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

#&check_tri_allelic();
#&create_vdb();
&read_coverage();

#method to read all allele data and detect situations when more than 2 different alleles for same variation
sub check_tri_allelic{

    #the input file should already be sorted by position
    open IN, "$input_file" or die "can't open input_file:$!!";
    open OUT, ">$output_file" or die "can't open output_file:$!!";
    #loop the file with the SNP information
    #SNP: 76 0 121186939 3 1-1-247249719 C/A gnl|ti|1759736436 -32 [  ] 0 3 0
    #     76 0 121186939 3 1-1-247249719 C/A gnl|ti|1759736680 220 [  ] 0 3 0
    #     76 0 121186939 3 1-1-247249719 C/A gnl|ti|1759732524 202 [  ] 0 3 0
    #SNP: 77 0 141615092 1 1-1-247249719 G/A gnl|ti|1759735512 56 [  ] 0 1 0
    #SNP: 78 0 141765872 1 1-1-247249719 A/T gnl|ti|1759732978 -148 [  ] 0 1 0
    my $previous_snp_pos=0; #variable to determinate wether we are in the same snp but different read.
    my %alleles; #hash to store the different alleles for this variation
    my %rec_strain;
    my $previous_line; #storing the previous line, the one we need to print to the file
    my $line;
    my ($null1,$null2,$ref_pos,$null3,$subj_name,$allele_string,$read_name,$read_pos,$strain_name,$tot_cover_num,$num_strain_snp,$strain_tot_coverage); #values in the line
    my ($ref_allele,$read_allele);
    my ($chr,$chr_start,$chr_end);
    my $uniq_pos;
    my $count=0;
    my $buffer={}; #buffer to store the lines to be written
    while(<IN>){
	next if /Warning/i;
	if (/^SNP|^\s+\d+\s+\d+\s+\d+\s+/) {
	    s/^SNP\:\s+|^\s+|\[|\]//g;
	    ($null1,$null2,$ref_pos,$null3,$subj_name,$allele_string,$read_name,$read_pos,$strain_name,$tot_cover_num,$num_strain_snp,$strain_tot_coverage) = split;
	    next if ($tot_cover_num == 0); #not sure what it means
	    #there is a SNP, record it
	    $strain_name =~ s/\[|\s+\]//g if (!defined $individual);
	    $strain_name = $individual if (defined $individual);
	    ($ref_allele,$read_allele) = split /\//, $allele_string;
	    ($chr,$chr_start,$chr_end) = split /\-/, $subj_name;
	    $uniq_pos = "$chr\-$ref_pos";	    
	    if ($previous_snp_pos ne $uniq_pos && $previous_snp_pos ne 0){
		#different SNP, time to print to file the information
		($null1,$null2,$ref_pos,$null3,$subj_name,$allele_string,$read_name,$read_pos,$strain_name,$tot_cover_num,$num_strain_snp,$strain_tot_coverage) = split /\s+/,$previous_line;
		foreach my $str_name (keys %alleles) {
		  if (keys %{$alleles{$str_name}}>2){
		    #this is tri-allelic
		    print_buffered($buffer,"$TMP_DIR/$TMP_FILE",join("\t",$previous_snp_pos,$str_name,keys %alleles) . "\n");
		  }
		  else{
		    $count++; #get the number of unique variation in the reads
		    ($chr,$chr_start) = split /\-/, $previous_snp_pos;
		    print_buffered($buffer,"$output_file",join("\t",$count,$previous_snp_pos,$chr,$chr_start,$str_name,$allele_string,$read_name,$read_pos,$tot_cover_num,$num_strain_snp,$strain_tot_coverage)."\n");
		    last;# we only need one record for each variation
		  }
		}
		#flush hash
		%alleles=();
	    }
	    $previous_snp_pos = $uniq_pos;
	    $alleles{$strain_name}{$ref_allele}++;
	    $alleles{$strain_name}{$read_allele}++;
	    $rec_strain{$strain_name}=1;
	    $previous_line = $_;
	}
    }
    #print the last variation
    $count++; #get the number of unique variation in the reads
    #different SNP, time to print to file the information

    foreach my $str_name (keys %alleles) {
      if (keys %alleles{$str_name}>2){
	#this is tri-allelic
	print_buffered($buffer,"$TMP_DIR/$TMP_FILE",join("\t",$previous_snp_pos,$str_name,keys %alleles) . "\n");
      }
      else{
	$count++;
        ($chr,$chr_start,$chr_end) = split /\-/, $subj_name;
	print_buffered($buffer,"$output_file",join("\t",$count,$previous_snp_pos,$chr,$chr_start,$str_name,$allele_string,$read_name,$read_pos,$tot_cover_num,$num_strain_snp,$strain_tot_coverage)."\n");
	last;# we only need one record for each variation
      }
    }

    close IN || die "Could not close SNP file:$!\n";

    print_buffered($buffer);#flush the buffer;

    create_and_load($dbVar,"tri_allelic_alleles","uniq_pos *","strain_name","alleles");
    system("mv $output_file $TMP_DIR/$TMP_FILE");
    create_and_load($dbVar,"ssaha_snp_out","count i*","uniq_pos *","chr *","target_start i*","strain_name","allele_string","read_name *","read_pos i*","tot_cover_num i","strain_snp_num i","strain_tot_coverage i");

    my $individual_type_id = 1;
    my $individual_size = keys %rec_strain;
    my $ref_strain = "refstrain";
    my $ind_pop_name = "ENSEMBL:strain";
    my $ind_sample_pop_desc = "Population for $individual_size individual(s)";
    MY $ind_sample_desc = "Individual within population $ind_pop_name";
    debug("Inserting into population, individual and sample tables");

    $dbVar->do(qq{INSERT INTO sample (name,description) values ("$ref_strain","reference population sample")});
    my $population_ref_sample_id = $dbVar->db_handle->{'mysql_insertid'};
    $dbVar->do(qq{INSERT INTO population (sample_id) values ($population_ref_sample_id)});
    $dbVar->do(qq{INSERT INTO sample (name,description) values ("$ind_pop_name","$ind_sample_pop_desc")});
    my $population_ind_sample_id = $dbVar->db_handle->{'mysql_insertid'};
    $dbVar->do(qq{INSERT INTO population (sample_id) values ($population_ind_sample_id)});


    foreach my $strain_name (keys %rec_strain) {
      $dbVar->do(qq{INSERT INTO sample (name,size,description) values ("$strain_name",NULL,"$ind_sample_desc")});
      my $individual_sample_id = $dbVar->db_handle->{'mysql_insertid'};
      $dbVar->do(qq{INSERT INTO individual (sample_id,individual_type_id) values ($individual_sample_id, $individual_type_id)});
      $dbVar->do(qq{INSERT INTO individual_population (individual_sample_id,population_sample_id) values ($individual_sample_id,$population_ind_sample_id)});
    }
}


sub create_vdb {
 
  debug("Creating database now...");

  my $source_name = "ENSEMBL"; #needs change everytime
  my $var_pre_name = "ENSEMBL";
  my $length_pre_name = length($var_pre_name);
  my $variation_name = "ENSSNP";###needs change for different species


    $dbVar->do(qq{INSERT INTO source (name) values ("$source_name")});
    my $source_id = $dbVar->{'mysql_insertid'};
    debug("Insert into variation table...");

    $dbVar->do(qq{ALTER TABLE variation add column internal_name varchar(50)});
    $dbVar->do(qq{INSERT INTO variation (source_id,name,internal_name) select s.source_id as source_id,concat("$variation_name",snp.count) as name,concat("$var_pre_name\_",snp.uniq_pos) as internal_name from source s, ssaha_snp_out snp});


    dumpSQL($dbCore,qq{SELECT sr.seq_region_id, sr.name
  		                FROM   seq_region_attrib sra, attrib_type at, seq_region sr
 		                WHERE sra.attrib_type_id=at.attrib_type_id 
	                        AND at.code="toplevel" 
                                AND sr.seq_region_id = sra.seq_region_id 
                     });
   create_and_load($dbVar,"tmp_seq_region","seq_region_id i*","seq_region_name *");

  debug("Insert into variation_feature table... NOTE ABOUT Y CHROMOSOME");
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
  my $buffer={}; #buffer to store the lines to be written
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
