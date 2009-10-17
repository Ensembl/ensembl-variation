use strict;
use warnings;
#object that contains the specific methods to dump data when there are chromosome coordinates from dbSNP (not contigs, as usual). 
#So far, this is the case for rat and chicken
package dbSNP::GenericChromosome;

use dbSNP::GenericContig;
use vars qw(@ISA);
use ImportUtils qw(debug load dumpSQL create_and_load);

@ISA = ('dbSNP::GenericContig');

sub variation_feature{
    my $self = shift;

     debug("Dumping seq_region data");

     #only take toplevel coordinates
     dumpSQL($self->{'dbCore'}->dbc()->db_handle, qq{SELECT sr.seq_region_id, 
  				      if (sr.name like "E%", CONCAT("LG",sr.name),sr.name) ##add LG for chicken
  				      FROM   seq_region_attrib sra, attrib_type at, seq_region sr
  				      WHERE sra.attrib_type_id=at.attrib_type_id 
 	                              AND at.code="toplevel" 
                                       AND sr.seq_region_id = sra.seq_region_id 
 				    });


     debug("Loading seq_region data");
     load($self->{'dbVar'}, "seq_region", "seq_region_id", "name");

     debug("Dumping SNPLoc data");
    
     my ($tablename1,$tablename2,$row);

     print "assembly_version is ",$self->{'assembly_version'},"\n";
     my ($assembly_version) =  $self->{'assembly_version'} =~ /^[a-zA-Z]+\_?(\d+)\.*.*$/;
     $assembly_version=$1 if $self->{'assembly_version'} =~ /RGSC\d\.(\d+)/;

     my $sth = $self->{'dbSNP'}->prepare(qq{SHOW TABLES LIKE 
 					   '$self->{'dbSNP_version'}\_SNPContigLoc\_$assembly_version\_%'});
	 
	 print qq{SHOW TABLES LIKE 
 					   '$self->{'dbSNP_version'}\_SNPContigLoc\_$assembly_version\_%'};
	 
     $sth->execute();

     while($row = $sth->fetchrow_arrayref()) {
       $tablename1 = $row->[0];
     }

     my $sth1 = $self->{'dbSNP'}->prepare(qq{SHOW TABLES LIKE 
 					   '$self->{'dbSNP_version'}\_ContigInfo\_$assembly_version\_%'});
     $sth1->execute();

     while($row = $sth1->fetchrow_arrayref()) {
       $tablename2 = $row->[0];
     }
     print "table_name1 is $tablename1 table_name2 is $tablename2\n";

    if (!$tablename1) {
      debug("core db has assembly version : $assembly_version, which is different from dbsnp");
      my $tablename1_ref = $self->{'dbSNP'}->selectall_arrayref(qq{SHOW TABLES LIKE 
 					   '$self->{'dbSNP_version'}\_SNPContigLoc\_%\_1'});
      $tablename1 = $tablename1_ref->[0][0] if $tablename1_ref;
      my $tablename2_ref = $self->{'dbSNP'}->selectall_arrayref(qq{SHOW TABLES LIKE 
 					   '$self->{'dbSNP_version'}\_ContigInfo\_%\_1'});
      $tablename2 = $tablename1_ref->[0][0] if $tablename2_ref;
    }

     #my $tablename = $self->{'species_prefix'} . 'SNPContigLoc';
     dumpSQL($self->{'dbSNP'}, qq{SELECT t1.snp_id, t2.contig_acc,t1.lc_ngbr+2,t1.rc_ngbr,
 				 IF(t2.group_term like "ref_%",t2.contig_chr,t2.contig_label), 
 				 IF(t1.loc_type = 3, t1.phys_pos_from+2, t1.phys_pos_from+1),
 				 IF(t1.loc_type = 3,  t1.phys_pos_from+1, t1.phys_pos_from+length(t1.allele)),
 				 IF(t1.orientation, -1, 1)
 				 FROM $tablename1 t1, $tablename2 t2 
 				 WHERE t1.ctg_id = t2.ctg_id
 				 #AND t2.group_term like "ref_%"
                                  $self->{'limit'}});
    
    
    debug("Loading SNPLoc data");
    
     create_and_load($self->{'dbVar'}, "tmp_contig_loc_chrom", "snp_id i*", "ctg *", "ctg_start i", "ctg_end i", "chr *", "start i", "end i", "strand i");

    debug("Creating genotyped variations");
    #creating the temporary table with the genotyped variations

    my $gtype_ref = $self->{'dbVar'}->db_handle->selectall_arrayref(qq{SELECT COUNT(*) FROM  tmp_individual_genotype_single_bp});
    my $gtype_row = $gtype_ref->[0][0] if $gtype_ref;
    if ($gtype_row) {
      $self->{'dbVar'}->do(qq{CREATE TABLE tmp_genotyped_var SELECT DISTINCT variation_id FROM tmp_individual_genotype_single_bp});
      $self->{'dbVar'}->do(qq{CREATE UNIQUE INDEX variation_idx ON tmp_genotyped_var (variation_id)});
      $self->{'dbVar'}->do(qq{INSERT IGNORE INTO tmp_genotyped_var SELECT DISTINCT variation_id FROM individual_genotype_multiple_bp});
    }
    debug("Creating tmp_variation_feature_chrom data");
    #if tcl.end>1, this means we have coordinates for chromosome, we will use it
    dumpSQL($self->{'dbVar'},qq{SELECT v.variation_id, ts.seq_region_id, 
                                      tcl.start,tcl.end,
                                      tcl.strand, v.name, v.source_id, v.validation_status
				      FROM variation v, tmp_contig_loc_chrom tcl, seq_region ts
				      WHERE v.snp_id = tcl.snp_id
				      AND tcl.start>2 #to get rid of lots of start=1
                                      AND tcl.chr = ts.name
    });

    create_and_load($self->{'dbVar'},'tmp_variation_feature_chrom',"variation_id *","seq_region_id", "seq_region_start", "seq_region_end", "seq_region_strand", "variation_name", "source_id", "validation_status");
    
    debug("Creating tmp_variation_feature_ctg data");
    #if tcl.start = 1 or tcl.end=1, this means we don't have mappings on chromosome, we take ctg coordinates if it is in toplevel
    dumpSQL($self->{'dbVar'},qq{SELECT v.variation_id, ts.seq_region_id, 
                                      tcl.ctg_start,tcl.ctg_end,
                                      tcl.strand, v.name, v.source_id, v.validation_status
				      FROM variation v, tmp_contig_loc_chrom tcl, seq_region ts
				      WHERE v.snp_id = tcl.snp_id
				      AND (tcl.start = 1 or tcl.end=1)
                                      AND tcl.ctg = ts.name
   });

    create_and_load($self->{'dbVar'},'tmp_variation_feature_ctg',"variation_id *","seq_region_id", "seq_region_start", "seq_region_end", "seq_region_strand", "variation_name", "source_id", "validation_status");

    debug("Dumping data into variation_feature table");
    if ($gtype_row) {
      foreach my $table ("tmp_variation_feature_chrom","tmp_variation_feature_ctg") {
	$self->{'dbVar'}->do(qq{INSERT INTO variation_feature (variation_id, seq_region_id,seq_region_start, seq_region_end, seq_region_strand,variation_name, flags, source_id, validation_status)
				  SELECT tvf.variation_id, tvf.seq_region_id, tvf.seq_region_start, tvf.seq_region_end, tvf.seq_region_strand,tvf.variation_name,IF(tgv.variation_id,'genotyped',NULL), tvf.source_id, tvf.validation_status
				  FROM $table tvf LEFT JOIN tmp_genotyped_var tgv ON tvf.variation_id = tgv.variation_id
				  });
      }
      #last fill in flags with genotyped
      $self->{'dbVar'}->do(qq{UPDATE variation_feature vf, tmp_genotyped_var tgv
                              SET vf.flags = "genotyped"
                              WHERE vf.variation_id = tgv.variation_id
                             });

    }
    else {

      debug("Dumping data into variation_feature table only used if table tmp_genotyped_var is not ready");
      foreach my $table ("tmp_variation_feature_chrom","tmp_variation_feature_ctg") {
	$self->{'dbVar'}->do(qq{INSERT INTO variation_feature (variation_id, seq_region_id,seq_region_start, seq_region_end, seq_region_strand,variation_name, flags, source_id, validation_status)
 				  SELECT tvf.variation_id, tvf.seq_region_id, tvf.seq_region_start, tvf.seq_region_end, tvf.seq_region_strand,tvf.variation_name,NULL, tvf.source_id, tvf.validation_status
 				  FROM $table tvf
 				  });
      }
    }


    #$self->{'dbVar'}->do("DROP TABLE tmp_contig_loc_chrom");
    #$self->{'dbVar'}->do("DROP TABLE tmp_genotyped_var");
    #$self->{'dbVar'}->do("DROP TABLE tmp_variation_feature_chrom");
    #$self->{'dbVar'}->do("DROP TABLE tmp_variation_feature_ctg");
    #for the chicken, delete 13,000 SNPs that cannot be mapped to EnsEMBL coordinate
    if ($self->{'dbCore'}->species =~ /gga/i){
	$self->{'dbVar'}->do("DELETE FROM variation_feature WHERE seq_region_end = -1");
    }
}

1;
