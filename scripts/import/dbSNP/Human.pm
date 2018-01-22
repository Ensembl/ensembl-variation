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
#object that contains the specific methods to dump data when the specie is a HUMAN (adds HGVbase and TSC information). 
package dbSNP::Human;

#use dbSNP::GenericContig;
use dbSNP::GenericChromosome;

use vars qw(@ISA);
use ImportUtils qw(debug load dumpSQL create_and_load);

#@ISA = ('dbSNP::GenericContig');
@ISA = ('dbSNP::GenericChromosome');

sub dump_dbSNP{
    my $self = shift;
    #first, dump all dbSNP data as usual
    $self->SUPER::dump_dbSNP();
    #then, get HGVbase IDs from Yuans file
    #$self->dump_HGVbaseIDs();
    #and finally, get TSC data from dbSNP
    #$self->dump_TSCIDs();
    $self->dump_AFFYIDs();
    #get mitochondrial SNPs provided by Yuan in .tb file formats---DON'T RUN THIS ANYMORE
    #$self->dump_mitocondrialSNPs();
}

#specific function to get the HGVbase IDs from a file provided by Yuan and add them to the variation_synonym table
sub dump_HGVbaseIDs{
    my $self = shift;
    
    #copy the file with the rs-> HGVbaseID information to the temp folder
    system "gunzip -c dbSNP/rs_hgvbase.txt.gz > " . $self->{'tmpdir'} . "/" . $self->{'tmpfile'};
    debug("Loading HGVbase data");
    create_and_load($self->{'dbVar'},"tmp_rs_hgvbase","rsID *","HGVbaseID");
    #add a new source to the Source table
    $self->{'dbVar'}->do(qq{INSERT INTO source (name,version) values ('HGVbase',15)
				  });
    debug("Adding HGVbaseIDs to synonym table");
    my $source_id = $self->{'dbVar'}->db_handle->{'mysql_insertid'}; #get the last autoinc id from the database (the one from the HGVbase source)
    #add the HGVbaseIDs to the Variation_Synonym table
    $self->{'dbVar'}->do(qq{INSERT INTO variation_synonym (variation_id,source_id,name)
				      SELECT v.variation_id, $source_id, trh.HGVbaseID
				      FROM variation v, tmp_rs_hgvbase trh
				      WHERE v.name = trh.rsID
				  });
    #and finally, remove the temporary table
    $self->{'dbVar'}->do(qq{DROP TABLE tmp_rs_hgvbase
				  });
}

#specific function to get the TSC IDs from dbSNP and add them to the variation_synonym table
sub dump_TSCIDs{
    my $self = shift;
    
         #add the TSC source to the table
     $self->{'dbVar'}->do(qq{INSERT INTO source (name,version) values ('TSC',1)
     });
     my $source_id = $self->{'dbVar'}->db_handle->{'mysql_insertid'}; #get the last autoinc id in the database (the one from the TSC source)
     #and finally add the TSC ids to the synonyms table
     debug("Dumping TSC information from dbSNP");
     
     my $stmt = qq{
	SELECT
	    'rs'+s.snp_id,
	    $source_id,
	    s.loc_snp_id
	FROM
	    SubSNP s
	WHERE
	    s.loc_snp_id LIKE 'TSC%'
     };
     dumpSQL($self->{'dbSNP'},$stmt);
     debug("Loading TSC ids into temporary table");
     create_and_load($self->{'dbVar'},"tmp_rs_TSC","rsID *","source_id","TSCid");
     $self->{'dbVar'}->do(qq{ INSERT IGNORE INTO variation_synonym (variation_id, source_id, name)
                                      SELECT v.variation_id, trt.source_id, trt.TSCid
                                      FROM variation v, tmp_rs_TSC trt
                                      WHERE v.name = trt.rsID
                                 }
                              );
     #and finally, remove the temporary table
     $self->{'dbVar'}->do(qq{DROP TABLE tmp_rs_TSC
                                   });
 }
 
 sub dump_AFFYIDs{
 
   my $self = shift;
   my ($source_name,$set_name);
 
   debug("Dumping AFFY information from dbSNP");
=head
    my $stmt = qq{
	SELECT
	    'rs'+s.snp_id,
	    s.loc_snp_id,
	    loc_batch_id
	FROM
	    SubSNP s,
	    Batch b
	WHERE
	    s.batch_id = b.batch_id AND
	    b.handle = 'AFFY'
    };
   dumpSQL($self->{'dbSNP'},$stmt);
   debug("Loading  ids into temporary table");
   create_and_load($self->{'dbVar'},"tmp_rs_AFFY","rsID *","AFFYid", "affy_name");
=cut
   foreach my $table ("yuan_aff_100k_var_46","yuan_aff_500k_var_46","yuan_aff_genome6_var_47") {
     if ($table =~ /100k/i) {
       $source_name = "Affy GeneChip 100K Array";
       $set_name = "Mapping50K";
     }
     elsif ($table =~ /500k/i) {
       $source_name = "Affy GeneChip 500K Array";
       $set_name = "Mapping250K";
     }
     elsif ($table =~ /genome6/i) {
       $source_name = "Affy GenomeWideSNP_6.0";
       $set_name = "6.0";
     }
 
     debug("Creating name_pair table with source_name $source_name...");
 
 #    $self->{'dbVar'}->do(qq{CREATE TABLE $table\_name_pair like $table.name_pair});
 #    $self->{'dbVar'}->do(qq{insert into $table\_name_pair select * from $table.name_pair});
 #    $self->{'dbVar'}->do(qq{insert ignore into $table\_name_pair
 #                            select a.AFFYid as affy_name,a.rsID as rs_name
 #                            from tmp_rs_AFFY a, $table.name_pair c
 #                            where a.AFFYid=c.affy_name});
 
 
     my $source_id_ref = $self->{'dbVar'}->db_handle->selectall_arrayref(qq{
                      SELECT source_id from source where name = "$source_name"});
     my $source_id = $source_id_ref->[0][0];
 
     if (!$source_id) {
       $self->{'dbVar'}->do(qq{insert into source (name) values("$source_name")});
       $source_id = $self->{'dbVar'}->db_handle->{'mysql_insertid'};
     }
 
     debug("Inserting in variation_synonym table from $table\_name_pair...");
     $self->{'dbVar'}->do(qq{ INSERT IGNORE INTO variation_synonym_test (variation_id, source_id, name)
                                      SELECT v.variation_id, $source_id as source_id, t.affy_name as name
                                      FROM variation v, $table\_name_pair t
                                      WHERE v.name = t.rs_name
                                 }
                         );
 
     #update rs_ID to rsCurrent from rsHigh
     #$self->{'dbVar'}->do(qq{update tmp_rs_AFFY_test t, rsHist h set t.rsID=h.rsCurrent where t.rsID=h.rsHigh});
     debug("Inserting in variation_synonym table from table  tmp_rs_AFFY with $source_name...");
     $self->{'dbVar'}->do(qq{ INSERT IGNORE INTO variation_synonym_test (variation_id, source_id, name)
                                      SELECT v.variation_id, $source_id as source_id, t.AFFYid as name
                                      FROM variation v, tmp_rs_AFFY_test t
                                      WHERE v.name = t.rsID
                                      AND t.affy_name like "%$set_name%"
                                 }
                         );
 
    #and finally, remove the temporary table
    #$self->{'dbVar'}->do(qq{DROP TABLE $table\_name_pair});
 
   }
 }
 
#will get from the RefSNP.tb and ContigHit files the information about the Simon mapped mitochondrial SNPs, and add the information to the relevant
# tables: Variation, Allele, Source, Flanking_sequence and Variation_Feature
sub dump_mitocondrialSNPs{
    my $self = shift;
    my %mitoSNPs; #hash containing all the information in the files RefSNP.tb and ContigHit.tb provided by Yuan at:
    my $variation_id; #internal id of the variation added to the database
    my $region; #hash that will contain, for a certain region, the seq_region_id in the database
    my $slice_adaptor = $self->{'dbCore'}->get_SliceAdaptor();
    my $slice;
    my $seq_region_id; #region for the mitocontig MT_NC_001807, extracted from the core database
    my $status;
    #/ecs2/scratch4/yuan/hum/MT_35
    #first of all, add the new source of information
     $self->{'dbVar'}->do(qq{INSERT INTO source (name) values ('mitomap.com')
     });
     my $source_id = $self->{'dbVar'}->dbh()->{'mysql_insertid'}; #get the last autoinc id in the database (the one from the mitomap.com source)
     
    #reads and loads into a hash table all the information in the RefSNP file
    $self->read_RefSNP(\%mitoSNPs,'/ecs2/scratch4/yuan/hum/MT_35/RefSNP.tb');
    #reads and loads into a hash table all the information in the ContigHit table referent to the location
    $self->read_ContigHit(\%mitoSNPs,'/ecs2/scratch4/yuan/hum/MT_35/ContigHit.tb');
    #and finally, add the information to the database
    foreach my $snp (keys %mitoSNPs){
	if ($mitoSNPs{$snp}{'status'} eq 'by-other-pop'){
	    $status = 4;
	}
	#insert in the Variation table
         $self->{'dbVar'}->do(qq{INSERT INTO variation (source_id,name,validation_status) VALUES ($source_id, "$mitoSNPs{$snp}{'name'}", $status);
                                   });
         $variation_id = $self->{'dbVar'}->dbh()->{'mysql_insertid'}; #get the last autoinc id in the database (the in the variation table)
 
	if (!exists $region->{$mitoSNPs{$snp}{'region'}}){
	    $slice = $slice_adaptor->fetch_by_region('toplevel',$mitoSNPs{$snp}{'region'}); #will get the slice for the region where the SNP is present
	    $region->{$mitoSNPs{$snp}{'region'}} = $slice_adaptor->get_seq_region_id($slice); #get the seq_region_id and store it in a hash
	}

	$seq_region_id = $region->{$mitoSNPs{$snp}{'region'}};
	#insert in the Flanking_sequence table
         $self->{'dbVar'}->do(qq{INSERT INTO flanking_sequence (variation_id,seq_region_id,seq_region_strand,up_seq,down_seq)
                                           VALUES ($variation_id,$seq_region_id,$mitoSNPs{$snp}{'strand'},"$mitoSNPs{$snp}{'up_seq'}",
                                                   "$mitoSNPs{$snp}{'down_seq'}")
                                       });
         #insert all the alleles
         foreach my $allele (split /\//,$mitoSNPs{$snp}{'alleles'}){
             $self->{'dbVar'}->do(qq{INSERT INTO allele (variation_id, allele) VALUES ($variation_id,"$allele")
                                           });
         }
         #and finally, insert the variation_feature table
         $self->{'dbVar'}->do(qq{INSERT INTO variation_feature (variation_id, seq_region_id,	
								     seq_region_start, seq_region_end, seq_region_strand, variation_name,source_id,validation_status)
					  VALUES ($variation_id, $seq_region_id, $mitoSNPs{$snp}{'start'}, $mitoSNPs{$snp}{'end'}, 
						  $mitoSNPs{$snp}{'strand'}, "$mitoSNPs{$snp}{'name'}",$source_id, $status)
				      });
    }
	
}

#will read from Yuans directory /ecs2/scratch4/yuan/hum/MT_35 the RefSNP.tb file. Important that in each release the format of the data and the location
#of the file might change
sub read_RefSNP{
    my $self = shift;
    my $snps = shift;
    my $file_location = shift;

    my ($snp_id, $snp_name, $alleles, $up_seq, $down_seq, $status); #values we want to get from the file
    my @line;
    open IN,$file_location || die "file with mitochondrial SNP information doesn't exist at $file_location\n";
    while (<IN>){
	chomp; #remove the last character
	@line = split /\t/;
	$snp_id = $line[0]; #internal id of the SNP
	$snp_name = $line[1]; #name of the snp
	$alleles = $line[5]; #alleles in C/T format
	$up_seq = $line[6]; #up_seq
	$down_seq = $line[7]; #down_seq
	$status = $line[13]; #validation_status
	$snps->{$snp_id}->{'name'} = $snp_name;
	$snps->{$snp_id}->{'alleles'} = $alleles;
	$snps->{$snp_id}->{'up_seq'} = $up_seq;
	$snps->{$snp_id}->{'down_seq'} = $down_seq;
	$snps->{$snp_id}->{'status'} = $status;
	
    }
    close IN;
}

#will read from Yuans directory /ecs2/scratch4/yuan/hum/MT_35 the ContigHit.tb file. Important that in each release the format of the data and the 
#location of the file might change
sub read_ContigHit{
    my $self = shift;
    my $snps = shift;
    my $file_location = shift;

    my ($snp_id,$region, $strand, $start, $end);
    my @line;
    open IN,$file_location || die "Could not open file with mitochondrial SNP information at $file_location\n";
    while (<IN>){
	chomp; #remove last character
	@line = split /\t/;
	$snp_id = $line[0]; #snp_id
	$region = $line[7]; #name of the region MT_NC....
	$strand = $line[9]; #strand of the SNP in the sequence
	$start = $line[10]; #physmapstart of the SNP
	$end = $line[11];  #physmapend of the SNP
	#get region name withou th MT
	$region =~ /MT_(NC_\d+)/;
	$snps->{$snp_id}->{'region'} = $1;
	$snps->{$snp_id}->{'strand'} = $strand;
	$snps->{$snp_id}->{'start'} = $start;
	$snps->{$snp_id}->{'end'} = $end;
	
    }
    close IN;
}

1;
