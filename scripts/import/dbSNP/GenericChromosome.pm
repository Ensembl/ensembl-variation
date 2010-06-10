use strict;
use warnings;
#object that contains the specific methods to dump data when there are chromosome coordinates from dbSNP (not contigs, as usual). 
#So far, this is the case for rat and chicken
package dbSNP::GenericChromosome;

use dbSNP::GenericContig;
use vars qw(@ISA);
use ImportUtils qw(debug load dumpSQL create_and_load);
use Progress;

@ISA = ('dbSNP::GenericContig');

sub variation_feature{
    my $self = shift;

     debug(localtime() . "\tDumping seq_region data");

     #only take toplevel coordinates
     dumpSQL($self->{'dbCore'}->dbc()->db_handle, qq{SELECT sr.seq_region_id, 
  				      if (sr.name like "E%", CONCAT("LG",sr.name),sr.name) ##add LG for chicken
  				      FROM   seq_region_attrib sra, attrib_type at, seq_region sr
  				      WHERE sra.attrib_type_id=at.attrib_type_id 
 	                              AND at.code="toplevel" 
                                       AND sr.seq_region_id = sra.seq_region_id 
 				    });


     debug(localtime() . "\tLoading seq_region data");
     load($self->{'dbVar'}, "seq_region", "seq_region_id", "name");
  print Progress::location();
    
     debug(localtime() . "\tDumping SNPLoc data");
    
     my ($tablename1,$tablename2,$row);

     print "assembly_version is ",$self->{'assembly_version'},"\n";
     my ($assembly_version) =  $self->{'assembly_version'} =~ /^[a-zA-Z]+\_?(\d+)\.*.*$/;
     $assembly_version=$1 if $self->{'assembly_version'} =~ /RGSC\d\.(\d+)/;

     my $stmt = qq{
                   SELECT 
                     name 
                   FROM 
                     $self->{'snp_dbname'}..sysobjects 
                   WHERE 
                     name LIKE '$self->{'dbSNP_version'}\_SNPContigLoc\_$assembly_version\_%'
                  };
     my $sth = $self->{'dbSNP'}->prepare($stmt);
     $sth->execute();

     while($row = $sth->fetchrow_arrayref()) {
       $tablename1 = $row->[0];
     }

     $stmt = qq{
                SELECT 
                  name 
                FROM 
                  $self->{'snp_dbname'}..sysobjects 
                WHERE 
                  name LIKE '$self->{'dbSNP_version'}\_ContigInfo\_$assembly_version\_%'
               };
     my $sth1 = $self->{'dbSNP'}->prepare($stmt);
     $sth1->execute();

     while($row = $sth1->fetchrow_arrayref()) {
       $tablename2 = $row->[0];
     }
     print "table_name1 is $tablename1 table_name2 is $tablename2\n";
    if (!$tablename1) {
      debug("core db has assembly version : $assembly_version, which is different from dbsnp");
      $stmt = qq{
	SELECT 
	    name 
        FROM 
	    $self->{'snp_dbname'}..sysobjects 
	WHERE 
	    name LIKE '$self->{'dbSNP_version'}\_SNPContigLoc\_%\_1'
      };
      my $tablename1_ref = $self->{'dbSNP'}->selectall_arrayref($stmt);
      $tablename1 = $tablename1_ref->[0][0] if $tablename1_ref;
      $stmt = qq{
	SELECT 
	    name 
        FROM 
	    $self->{'snp_dbname'}..sysobjects 
	WHERE 
	    name LIKE '$self->{'dbSNP_version'}\_ContigInfo\_%\_1'
      };
      my $tablename2_ref = $self->{'dbSNP'}->selectall_arrayref($stmt);
      $tablename2 = $tablename1_ref->[0][0] if $tablename2_ref;
    }

     #my $tablename = $self->{'species_prefix'} . 'SNPContigLoc';

    # In the query below, the pre-131 syntax was ref-assembly. In 131 it is GRCh37 for human. What is it for other species??
    my $group_term = 'ref_';
    my ($release) = $self->{'dbSNP_version'} =~ m/^b?(\d+)$/;
    $group_term = 'GRCh' if ($self->{'dbCore'}->species =~ m/homo|human/i && $release > 130);
    
	#	     t2.group_term LIKE 'ref_%'
     $stmt = "SELECT ";
     if ($self->{'limit'}) {
       $stmt .= "TOP $self->{'limit'} ";
     }
     $stmt .= qq{
                   t1.snp_id AS sorting_id, 
                   t2.contig_acc,
                   t1.lc_ngbr+2,t1.rc_ngbr,
		   CASE WHEN
		     t2.group_term LIKE '$group_term%'
		   THEN
		     t2.contig_chr
		   ELSE
		     t2.contig_label
		   END, 
		   CASE WHEN
		     t1.loc_type = 3
		   THEN
		     t1.phys_pos_from+2
		   ELSE
		     t1.phys_pos_from+1
		   END,
		   CASE WHEN
		     t1.loc_type = 3
		   THEN
		     t1.phys_pos_from+1
		   ELSE
		     t1.phys_pos_from+LEN(t1.allele)
		   END,
		   CASE WHEN
		     t1.orientation = 1
		   THEN
		     -1
		   ELSE
		     1
		   END
		 FROM 
		   $tablename1 t1, 
		   $tablename2 t2 
		 WHERE 
		   t1.ctg_id = t2.ctg_id
	        };
     if ($self->{'limit'}) {
       $stmt .= qq{    
		   ORDER BY
		     sorting_id ASC  
	          };
     }
 				 #AND t2.group_term like "ref_%"
     dumpSQL($self->{'dbSNP'},$stmt);
    
    
    debug(localtime() . "\tLoading SNPLoc data");
    
     create_and_load($self->{'dbVar'}, "tmp_contig_loc_chrom", "snp_id i*", "ctg *", "ctg_start i", "ctg_end i", "chr *", "start i", "end i", "strand i");
  print Progress::location();

    debug(localtime() . "\tCreating genotyped variations");
    #creating the temporary table with the genotyped variations

    my $gtype_ref = $self->{'dbVar'}->db_handle->selectall_arrayref(qq{SELECT COUNT(*) FROM  tmp_individual_genotype_single_bp});
    my $gtype_row = $gtype_ref->[0][0] if $gtype_ref;
    if ($gtype_row) {
      $self->{'dbVar'}->do(qq{CREATE TABLE tmp_genotyped_var SELECT DISTINCT variation_id FROM tmp_individual_genotype_single_bp});
  print Progress::location();
      $self->{'dbVar'}->do(qq{CREATE UNIQUE INDEX variation_idx ON tmp_genotyped_var (variation_id)});
  print Progress::location();
      $self->{'dbVar'}->do(qq{INSERT IGNORE INTO tmp_genotyped_var SELECT DISTINCT variation_id FROM individual_genotype_multiple_bp});
  print Progress::location();
    }
    debug(localtime() . "\tCreating tmp_variation_feature_chrom data");
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
  print Progress::location();
    
    debug(localtime() . "\tCreating tmp_variation_feature_ctg data");
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
  print Progress::location();

    debug(localtime() . "\tDumping data into variation_feature table");
    if ($gtype_row) {
      foreach my $table ("tmp_variation_feature_chrom","tmp_variation_feature_ctg") {
	$self->{'dbVar'}->do(qq{INSERT INTO variation_feature (variation_id, seq_region_id,seq_region_start, seq_region_end, seq_region_strand,variation_name, flags, source_id, validation_status)
				  SELECT tvf.variation_id, tvf.seq_region_id, tvf.seq_region_start, tvf.seq_region_end, tvf.seq_region_strand,tvf.variation_name,IF(tgv.variation_id,'genotyped',NULL), tvf.source_id, tvf.validation_status
				  FROM $table tvf LEFT JOIN tmp_genotyped_var tgv ON tvf.variation_id = tgv.variation_id
				  });
  print Progress::location();
      }
      #last fill in flags with genotyped
      $self->{'dbVar'}->do(qq{UPDATE variation_feature vf, tmp_genotyped_var tgv
                              SET vf.flags = "genotyped"
                              WHERE vf.variation_id = tgv.variation_id
                             });
  print Progress::location();

    }
    else {

      debug(localtime() . "\tDumping data into variation_feature table only used if table tmp_genotyped_var is not ready");
      foreach my $table ("tmp_variation_feature_chrom","tmp_variation_feature_ctg") {
	$self->{'dbVar'}->do(qq{INSERT INTO variation_feature (variation_id, seq_region_id,seq_region_start, seq_region_end, seq_region_strand,variation_name, flags, source_id, validation_status)
 				  SELECT tvf.variation_id, tvf.seq_region_id, tvf.seq_region_start, tvf.seq_region_end, tvf.seq_region_strand,tvf.variation_name,NULL, tvf.source_id, tvf.validation_status
 				  FROM $table tvf
 				  });
  print Progress::location();
      }
    }


    #$self->{'dbVar'}->do("DROP TABLE tmp_contig_loc_chrom");
    #$self->{'dbVar'}->do("DROP TABLE tmp_genotyped_var");
    #$self->{'dbVar'}->do("DROP TABLE tmp_variation_feature_chrom");
    #$self->{'dbVar'}->do("DROP TABLE tmp_variation_feature_ctg");
    #for the chicken, delete 13,000 SNPs that cannot be mapped to EnsEMBL coordinate
    if ($self->{'dbCore'}->species =~ /gga/i){
	$self->{'dbVar'}->do("DELETE FROM variation_feature WHERE seq_region_end = -1");
  print Progress::location();
    }
}

1;
