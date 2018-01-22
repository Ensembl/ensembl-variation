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


my %rec_file;
my %rec_name;

use strict;
use Getopt::Long;

our ($input_dir, $output_dir, $individual_name, $file_dir);

GetOptions('input_dir=s'    => \$input_dir,
           'output_dir=s'    => \$output_dir,
           'individual_name=s'   => \$individual_name,
		   'filedir=s'			=> \$file_dir,
          );

#($input_dir and $output_dir) || die "We need input_dir, output_dir defined\n";
my %strain_name = ('GKO_Ox'    =>  'GK/Ox',
                   'DAHLS'  =>  'SS/Jr',
                   'SHRSP'  =>  'SHRSP/mdc',
                   'WKY'    =>  'WKY/mdc',
                   'SS_jr'     =>  'SS/Jr',
		   'ref_rat' => 'BN/Crl',
		   'SD'    => 'SD',
                  );


#get_human_individual_name();
#get_tetraodon_individual_name();
#get_tetraodon_pair_end();
make_alignment_file($input_dir,$output_dir,$individual_name);
#change_to_id_pair($input_dir,$output_dir);


sub get_human_individual_name {

  LINE : while (<$file_dir/cra-homo_sapiens-wgs-*\.TRACEINFO.xml.gz>) {
    #print "my file is #$_#\n";
    $_ =~ /cra-homo_sapiens-wgs-(\d+)\.TRACEINFO.xml.gz/;
    my $file_number = $1;
    #print "my file_num is #$file_number#\n";
    my $file = $_;
    open IN, "gunzip -c $file |" or die "can't open $file : $!";
    while (<IN>) {
      #<seq_lib_id>HuBB2kb</seq_lib_id>
      if (/\<seq_lib_id\>(.*)\<\/seq_lib_id\>/) {
	my ($individual) = $1 =~ /(Hu[ABCDF]+)\d+/;
	#print "individual is $individual\n";
	$rec_file{$file_number}=$individual;
	next LINE;
      }
    }
  }

  #foreach my $file_number (keys %rec_file) {
  #  print "$file_number\t$rec_file{$file_number}\n";
  #}
}

sub get_tetraodon_individual_name {

  while (<$file_dir/*gsc*\.TRACEINFO.xml.gz>) {
    print "my file is #$_#\n";
    my $file = $_;
    my ($read_name,$insert_size,$lib_id);
    open IN, "gunzip -c $file |" or die "can't open $file : $!";
    while (<IN>) {
      #<trace_name>G41P605616FB4.T0</trace_name>
      if (/\<trace_name\>(.*)\<\/trace_name\>/) {
	($read_name) = $1 =~ /(.*)\.?.*/;
      	#print "read_name is $read_name\n";
      }
      #<insert_size>F</insert_size>
      elsif (/\<insert_size\>(.*)\<\/insert_size\>/) {
	($insert_size) = $1;
 	#print "$read_name\t$lib_id\t$insert_size\n";
	$rec_name{$read_name}{'insert_size'} = $insert_size;
      }
    }
  }
}
sub get_tetraodon_pair_end {

  LINE : while (<$file_dir/*\.TRACEINFO.xml.gz>) {
    print "my file is #$_#\n";
    my $file = $_;
    my ($read_name,$trace_end);
    open IN, "gunzip -c $file |" or die "can't open $file : $!";
    while (<IN>) {
      #<trace_name>G41P605616FB4.T0</trace_name>
      if (/\<trace_name\>(.*)\<\/trace_name\>/) {
	($read_name) = $1 =~ /(.*)\.?.*/;
      	#print "read_name is $read_name\n";
      }
      #<trace_end>F</trace_end>
      elsif (/\<trace_end\>(.*)\<\/trace_end\>/) {
	($trace_end) = $1;
	if ($trace_end eq "F") {
	  $trace_end = ".p1k";
	}
	else {
	  $trace_end = ".q1k";
	}
	$rec_name{$read_name}{'pair_end'}=$read_name.$trace_end;
	#print "$read_name and $read_name.$trace_end\n";
      }
    }
  }
}


sub make_alignment_file {
  my ($input_dir,$output_dir,$individual_name) = @_;

  open SNP, ">$output_dir/SNP_file" or die "can't open SNP file:$!";
  open OUT, ">$output_dir/ALIGNMENT_file" or die "can't open output file:$!";
  open CIG, ">$output_dir/cigar_file" or die "can't open cigar file:$!";

  my $score;

  while (<$input_dir/ssaha_out\*>) {
    #my ($file_number) = $_ =~ /cra-homo_sapiens-wgs-(\d+)\.fasta/;
    my ($file_name) = $_ =~ /ssaha_out_(\S+)\.fastq\_\d+.*$/;
    $individual_name = $strain_name{$file_name} if ($strain_name{$file_name});
    open IN, "$_" or die "can't open input file:$!";
    while (<IN>) {
      if (/^Score/) {
	$score=1;
      }
      elsif (/^$|^ProcessSNP_start|\=+/gi) {
	$score=0;
      }
      elsif ($score==1 and /^\d+\s+\S+\s+\S+\s+\d+\s+\d+\s+\d+\s+\d+\s+\S+\s+\d+\s+\S+\s+\d+$/) {
	chomp;
	my @all = split;
# 	if ($rec_name{$all[1]}{'insert_size'} >=130000) {
# 	  $individual_name = "BAC";
# 	}
# 	else {
# 	  $individual_name = "plasmid";
# 	}
# 	$all[1] = $rec_name{$all[1]}{'pair_end'};
	my $line = join " ", @all;
	my $line = "ALIGNMENT ".$line." $individual_name" if $_;
	print OUT "$line\n";
      }
      elsif ($score==1 and /^\S+\s+\d+\s+\d+\s+\S+\s+\S+\s+\d+\s+\d+\s+\S+\s+\d+/) {
        print CIG "$_";
      }
      elsif (/^ssaha:SNP/) {
	my @all=split;
#	$all[3] = $rec_name{$all[3]}{'pair_end'};
	my $line = join " ",@all;
	print SNP "$line\n";
      }
    }
  }
}

sub change_to_id_pair {

  my %rec_name;

  my ($input_dir,$output_dir) = @_;
  open IN, "$input_dir/id_to_pair_insert_lib_fixed";
  #open IN, "$input_dir/id_pair_test";

  while (<IN>) {
    chomp;
    #17000092642284 19000081944854.p1k 40139 HuAA50.1b
    my ($read_name,$template) = split;
    $rec_name{$read_name}=$template;
  }

  open ALI, "$input_dir/ALIGNMENT_file" or die "can't open file align_pair_input :$!";
  open SNP, "$input_dir/SNP_file" or die "can't open file snp_pair_input :$!";;
  open OUTA,">$input_dir/ALIGNMENT_file_pair" or die "can't open file align_pair :$!";
  open OUTS,">$input_dir/SNP_file_pair" or die "can't open file snp_pair :$!";;

  while (<ALI>) {
    my @all = split;
    if (@all == 13) {
      $all[2] = $rec_name{$all[2]};
      print OUTA join " ", @all,"\n";
    }
  }

  while (<SNP>) {
    my @all = split;
    if (@all == 15) {
      $all[3] = $rec_name{$all[3]};
      print OUTS join " ", @all,"\n";
    }
  }
}
