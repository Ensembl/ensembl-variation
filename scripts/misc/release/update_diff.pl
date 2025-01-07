=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2025] EMBL-European Bioinformatics Institute

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


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

#!/usr/bin/perl -w

#This script compares the annotation between two Ensembl GFF3 files 
#and reports transcripts that have been created, updated or deleted,
#based on transcript, exon and CDS coordinate differences.
#Both GFF3 files are linked uning the gene and transcript stable ids.
#This script does not aim to detect gene or transcript stable id 
#changes between releases.

use strict;
use Getopt::Long;


my $old_file;
my $new_file;
my $outfile;
&GetOptions(
            'old=s' => \$old_file,
            'new=s' => \$new_file,
            'out=s' => \$outfile
           );

my $usage =  <<_USAGE_;

This script compares the annotation between two Ensembl GFF3 files 
and reports transcripts that have been created, updated or deleted,
based on transcript, exon and CDS coordinate differences.
The input files can be compressed (.gz).
The output file contains one line per transcript and includes 
a header with the field names.

perl diff_ensembl_gff3_files.pl -old OLD_GFF3_FILE -new NEW_GFF3_FILE -out OUTPUT_FILE
 
_USAGE_

die $usage unless ($old_file and $new_file and $outfile);


my %old_annot = parse_gff3_file($old_file);
my %new_annot = parse_gff3_file($new_file);

open (OUT, ">$outfile") or die "Can't open file $outfile: $!";
print OUT join("\t", "transcript_id", "status", "gene_id", "version", "chr", "start", "end", "strand", "biotype", "exon_starts", "exon_ends", "CDS_start", "CDS_end")."\n";

foreach my $gene_id (keys %new_annot){
  #New gene
  if (!exists($old_annot{$gene_id})){
    foreach my $transcript_id (keys %{$new_annot{$gene_id}{"transcripts"}}){
      print OUT join("\t", $transcript_id, "new", @{transcript_info(undef, $new_annot{$gene_id}{"transcripts"}{$transcript_id}, "new")})."\n"; 
    }
  }
  else{
    foreach my $transcript_id (keys %{$new_annot{$gene_id}{"transcripts"}}){
	  #Novel transcript
      if (!exists($old_annot{$gene_id}{"transcripts"}{$transcript_id})){
        print OUT join("\t", $transcript_id, "new", @{transcript_info(undef, $new_annot{$gene_id}{"transcripts"}{$transcript_id}, "new")})."\n"; 
      }
      else{
		#Transcript comparison - check for updated transcripts
		my $transcript_status = "unchanged";
		my @key_fields = ("chr", "start", "end", "strand", "exon_starts", "exon_ends");		
		foreach my $key (@key_fields){
		  if ($new_annot{$gene_id}{"transcripts"}{$transcript_id}{$key} ne $old_annot{$gene_id}{"transcripts"}{$transcript_id}{$key}){
            $transcript_status = "updated";
          }
        }
        @key_fields = ("CDS_start", "CDS_end");	
        foreach my $key (@key_fields){
          if ((exists($new_annot{$gene_id}{"transcripts"}{$transcript_id}{$key}) and !exists($old_annot{$gene_id}{"transcripts"}{$transcript_id}{$key})) or
              (exists($old_annot{$gene_id}{"transcripts"}{$transcript_id}{$key}) and !exists($new_annot{$gene_id}{"transcripts"}{$transcript_id}{$key})) or
              (exists($old_annot{$gene_id}{"transcripts"}{$transcript_id}{$key}) and exists($new_annot{$gene_id}{"transcripts"}{$transcript_id}{$key}) and
               $new_annot{$gene_id}{"transcripts"}{$transcript_id}{$key} ne $old_annot{$gene_id}{"transcripts"}{$transcript_id}{$key})
             ){
            $transcript_status = "updated";
          }
        }       
        if ($transcript_status eq "updated"){
          print OUT join("\t", $transcript_id, "updated", @{transcript_info($old_annot{$gene_id}{"transcripts"}{$transcript_id}, $new_annot{$gene_id}{"transcripts"}{$transcript_id}, "updated")})."\n";
        }
      }
    }
    #Deleted transcript
    foreach my $transcript_id (keys %{$old_annot{$gene_id}{"transcripts"}}){
      if (!exists($new_annot{$gene_id}{"transcripts"}{$transcript_id})){
        print OUT join("\t", $transcript_id, "deleted", @{transcript_info($old_annot{$gene_id}{"transcripts"}{$transcript_id}, undef, "deleted")})."\n";
      } 
    }
  }
}
#Deleted gene
foreach my $gene_id (keys %old_annot){
  if (!exists($new_annot{$gene_id})){
    foreach my $transcript_id (keys %{$old_annot{$gene_id}{"transcripts"}}){
      print OUT join("\t", $transcript_id, "deleted", @{transcript_info($old_annot{$gene_id}{"transcripts"}{$transcript_id}, undef, "deleted")})."\n";
    }
  }
}

close (OUT);


########

sub parse_gff3_file {
  my $file = shift;
  my %annot;
  my %tid2gid;
  
  if ($file =~ /\.gz$/){
    open (GFF3, "zcat $file | ") or die "Can't open file $file: $!";
  }
  else{
    open (GFF3, $file) or die "Can't open file $file: $!";
  }
  while (<GFF3>){
	chomp;
    next if /^#/;
    my ($chr, $source, $feat, $start, $end, $score, $strand, $frame, $attributes) = split(/\t/);
    my %attribs;
    foreach my $attrib (split(/;/, $attributes)){
      my ($key, $value) = split(/=/, $attrib);
      $attribs{$key} = $value;
    }
    if ($_ =~ /ID=gene/){
	  my $gene_id = $attribs{"gene_id"};
	  $annot{$gene_id}{"chr"} = $chr;
	  $annot{$gene_id}{"start"} = $start;
	  $annot{$gene_id}{"end"} = $end;
	  $annot{$gene_id}{"strand"} = $strand;
	  $annot{$gene_id}{"source"} = $source;
	  foreach my $key (keys %attribs){
        $annot{$gene_id}{$key} = $attribs{$key};
      }
    }
    elsif ($_ =~ /ID=transcript/){
	  if ($attribs{"Parent"} =~ /gene:(ENS.+)/){
        my $parent_gene_id = $1;
        my $transcript_id = $attribs{"transcript_id"};
	    $annot{$parent_gene_id}{"transcripts"}{$transcript_id}{"chr"} = $chr;
	    $annot{$parent_gene_id}{"transcripts"}{$transcript_id}{"start"} = $start;
	    $annot{$parent_gene_id}{"transcripts"}{$transcript_id}{"end"} = $end;
	    $annot{$parent_gene_id}{"transcripts"}{$transcript_id}{"strand"} = $strand;
	    $annot{$parent_gene_id}{"transcripts"}{$transcript_id}{"source"} = $source;      
	    foreach my $key (keys %attribs){
          $annot{$parent_gene_id}{"transcripts"}{$transcript_id}{$key} = $attribs{$key};
        }
        $tid2gid{$transcript_id} = $parent_gene_id;
      }
    }
    elsif ($feat eq "exon"){
      if ($attribs{"Parent"} =~ /transcript:(ENS.+)/){
		my $parent_transcript_id = $1; 
		my $parent_gene_id = $tid2gid{$parent_transcript_id};
        my $exon_id = $attribs{"exon_id"};
        $annot{$parent_gene_id}{"transcripts"}{$parent_transcript_id}{"exons"}{$exon_id}{"start"} = $start;
        $annot{$parent_gene_id}{"transcripts"}{$parent_transcript_id}{"exons"}{$exon_id}{"end"} = $end;        
        foreach my $key (keys %attribs){
          $annot{$parent_gene_id}{"transcripts"}{$parent_transcript_id}{"exons"}{$exon_id}{$key} = $attribs{$key};
        }
        my @starts = ();
        my @ends = ();
        foreach my $exon_id (keys %{$annot{$parent_gene_id}{"transcripts"}{$parent_transcript_id}{"exons"}}){
          push(@starts, $annot{$parent_gene_id}{"transcripts"}{$parent_transcript_id}{"exons"}{$exon_id}{"start"});
          push(@ends, $annot{$parent_gene_id}{"transcripts"}{$parent_transcript_id}{"exons"}{$exon_id}{"end"});
        }
        $annot{$parent_gene_id}{"transcripts"}{$parent_transcript_id}{"exon_starts"} = join(",", sort {$a<=>$b} @starts);
        $annot{$parent_gene_id}{"transcripts"}{$parent_transcript_id}{"exon_ends"} = join(",", sort {$a<=>$b} @ends);
      } 
    }
    elsif ($feat eq "CDS"){
      if ($attribs{"Parent"} =~ /transcript:(ENS.+)/){
		my $parent_transcript_id = $1; 
		my $parent_gene_id = $tid2gid{$parent_transcript_id};
        $annot{$parent_gene_id}{"transcripts"}{$parent_transcript_id}{"protein_id"} = $attribs{"protein_id"};
        unless (defined($annot{$parent_gene_id}{"transcripts"}{$parent_transcript_id}{"CDS_start"}) and $annot{$parent_gene_id}{"transcripts"}{$parent_transcript_id}{"CDS_start"} <= $start){
          $annot{$parent_gene_id}{"transcripts"}{$parent_transcript_id}{"CDS_start"} = $start;
        }
        unless (defined($annot{$parent_gene_id}{"transcripts"}{$parent_transcript_id}{"CDS_end"}) and $annot{$parent_gene_id}{"transcripts"}{$parent_transcript_id}{"CDS_end"} >= $end){
          $annot{$parent_gene_id}{"transcripts"}{$parent_transcript_id}{"CDS_end"} = $end;
        }
      } 
    }    
  }
  close (GFF3);
  
  return %annot;
}


sub transcript_info {
  my ($old_transcript_data, $new_transcript_data, $transcript_status) = @_;
  my @fields = ("Parent", "version", "chr", "start", "end", "strand", "biotype", "exon_starts", "exon_ends", "CDS_start", "CDS_end");
  my @values;
  foreach my $field (@fields){
    my $value;   
    my $old_value = $old_transcript_data->{$field} || "null";
    my $new_value = $new_transcript_data->{$field} || "null";
    if ($transcript_status eq "updated"){
      if ($old_value ne $new_value){
        $value = $old_value."->".$new_value;
      }
      else{
        $value = $old_value;
      }
    }
    elsif ($transcript_status eq "new"){
      $value = $new_value;
    }
    elsif ($transcript_status eq "deleted"){
      $value = $old_value;
    }    
    if ($field eq "Parent"){
      $value =~ s/gene://;
    }    
	push(@values, $value);	
  }
  return \@values;
}
