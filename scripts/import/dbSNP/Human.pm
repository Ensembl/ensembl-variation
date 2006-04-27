use strict;
use warnings;
#object that contains the specific methods to dump data when the specie is a HUMAN (adds HGVbase and TSC information). 
package dbSNP::Human;

use dbSNP::GenericContig;
use vars qw(@ISA);
use ImportUtils qw(debug load dumpSQL create_and_load);

@ISA = ('dbSNP::GenericContig');

sub dump_dbSNP{
    my $self = shift;
    #first, dump all dbSNP data as usual
    $self->SUPER::dump_dbSNP();
    #then, get HGVbase IDs from Yuans file
    #$self->dump_HGVbaseIDs();
    #and finally, get TSC data from dbSNP
    #$self->dump_TSCIDs();    
    #get mitochondrial SNPs provided by Yuan in .tb file formats---DON'T RUN THIS ANYMORE
    #$self->dump_mitocondrialSNPs();
}

#specific function to get the HGVbase IDs from a file provided by Yuan and add them to the variation_synonym table
sub dump_HGVbaseIDs{
    my $self = shift;
    
    #copy the file with the rs-> HGVbaseID information to the temp folder
    system "gunzip -c dbSNP/rs_hgvbase.txt.gz > " . $self->{'tmpdir'} . "/" . $self->{'tmpfile'};
    debug("Loading HGVbase data");
    create_and_load($self->{'dbVariation'},"tmp_rs_hgvbase","rsID *","HGVbaseID");
    #add a new source to the Source table
    $self->{'dbVariation'}->do(qq{INSERT INTO source (name,version) values ('HGVbase',15)
				  });
    debug("Adding HGVbaseIDs to synonym table");
    my $source_id = $self->{'dbVariation'}->{'mysql_insertid'}; #get the last autoinc id from the database (the one from the HGVbase source)
    #add the HGVbaseIDs to the Variation_Synonym table
    $self->{'dbVariation'}->do(qq{INSERT INTO variation_synonym (variation_id,source_id,name)
				      SELECT v.variation_id, $source_id, trh.HGVbaseID
				      FROM variation v, tmp_rs_hgvbase trh
				      WHERE v.name = trh.rsID
				  });
    #and finally, remove the temporary table
    $self->{'dbVariation'}->do(qq{DROP TABLE tmp_rs_hgvbase
				  });
}

#specific function to get the TSC IDs from dbSNP and add them to the variation_synonym table
sub dump_TSCIDs{
    my $self = shift;
    
    #add the TSC source to the table
    $self->{'dbVariation'}->do(qq{INSERT INTO source (name,version) values ('TSC',1)
    });
    my $source_id = $self->{'dbVariation'}->{'mysql_insertid'}; #get the last autoinc id in the database (the one from the TSC source)
    #and finally add the TSC ids to the synonyms table
    debug("Dumping TSC information from dbSNP");
    dumpSQL($self->{'dbSNP'}, qq{SELECT concat('rs',ss.snp_id), $source_id, s.loc_snp_id 
				     FROM SubSNP s, SNPSubSNPLink ss 
				     WHERE ss.subsnp_id = s.subsnp_id 
				     AND s.loc_snp_id like 'TSC%'
				 }
	    );
    debug("Loading TSC ids into temporary table");
    create_and_load($self->{'dbVariation'},"tmp_rs_TSC","rsID *","source_id","TSCid");
    $self->{'dbVariation'}->do(qq{ INSERT IGNORE INTO variation_synonym (variation_id, source_id, name)
				     SELECT v.variation_id, trt.source_id, trt.TSCid 
				     FROM variation v, tmp_rs_TSC trt
				     WHERE v.name = trt.rsID 
				}
			     );
    #and finally, remove the temporary table
    $self->{'dbVariation'}->do(qq{DROP TABLE tmp_rs_TSC
				  });
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
    $self->{'dbVariation'}->do(qq{INSERT INTO source (name) values ('mitomap.com')
    });
    my $source_id = $self->{'dbVariation'}->dbh()->{'mysql_insertid'}; #get the last autoinc id in the database (the one from the mitomap.com source)
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
	$self->{'dbVariation'}->do(qq{INSERT INTO variation (source_id,name,validation_status) VALUES ($source_id, "$mitoSNPs{$snp}{'name'}", $status);
				  });
	$variation_id = $self->{'dbVariation'}->dbh()->{'mysql_insertid'}; #get the last autoinc id in the database (the in the variation table)

	if (!exists $region->{$mitoSNPs{$snp}{'region'}}){
	    $slice = $slice_adaptor->fetch_by_region('toplevel',$mitoSNPs{$snp}{'region'}); #will get the slice for the region where the SNP is present
	    $region->{$mitoSNPs{$snp}{'region'}} = $slice_adaptor->get_seq_region_id($slice); #get the seq_region_id and store it in a hash
	}

	$seq_region_id = $region->{$mitoSNPs{$snp}{'region'}};
	#insert in the Flanking_sequence table
	$self->{'dbVariation'}->do(qq{INSERT INTO flanking_sequence (variation_id,seq_region_id,seq_region_strand,up_seq,down_seq) 
					  VALUES ($variation_id,$seq_region_id,$mitoSNPs{$snp}{'strand'},"$mitoSNPs{$snp}{'up_seq'}",
						  "$mitoSNPs{$snp}{'down_seq'}")
				      });
	#insert all the alleles
	foreach my $allele (split /\//,$mitoSNPs{$snp}{'alleles'}){
	    $self->{'dbVariation'}->do(qq{INSERT INTO allele (variation_id, allele) VALUES ($variation_id,"$allele")
					  });
	}
	#and finally, insert the variation_feature table
	$self->{'dbVariation'}->do(qq{INSERT INTO variation_feature (variation_id, seq_region_id, 
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
