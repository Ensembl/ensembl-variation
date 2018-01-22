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

    ## new mysql version errors with empty not null columns
    ## switch to null allowable here then back to non null in QC process 
    $self->{'dbVar'}->do("alter table variation_feature modify column map_weight int default null");

     debug(localtime() . "\tDumping seq_region data using GenericChromosome");

     #only take toplevel coordinates
     dumpSQL($self->{'dbCore'}->db_handle, qq{SELECT sr.seq_region_id, sr.name
  				       	      FROM seq_region_attrib sra, attrib_type at, seq_region sr
  				      WHERE sra.attrib_type_id=at.attrib_type_id 
 	                              AND at.code="toplevel" 
                                       AND sr.seq_region_id = sra.seq_region_id 
 				    });


     debug(localtime() . "\tLoading seq_region data");
     load($self->{'dbVar'}, "seq_region", "seq_region_id", "name");
     print Progress::location();

     my $version_number = $self->{'dbSNP_version'};
    $version_number =~ s/^b//;
     debug(localtime() . "\tDumping SNPLoc data for dbSNP version $version_number " );
    
     my ($tablename1,$tablename2,$row);

    if(  $version_number  < 137){ ## table rename for dbSNP - keeping this temporarily for backwards comparibility


     print "assembly_version is ",$self->{'assembly_version'},"\n";
     my ($assembly_version) =  $self->{'assembly_version'} =~ /^[a-zA-Z]+\_?(\d+)\.*.*$/;
     $assembly_version=$1 if $self->{'assembly_version'} =~ /RGSC\d\.(\d+)/;

     my $stmt = qq{ SELECT  name 
                    FROM    $self->{'snp_dbname'}..sysobjects 
                    WHERE   name LIKE '$self->{'dbSNP_version'}\_SNPContigLoc\__%'
                    ORDER BY name DESC };
                  
     my $sth = $self->{'dbSNP'}->prepare($stmt);
     $sth->execute();

    my $genome_build;
    my @genome_builds;
     while($row = $sth->fetchrow_arrayref()) {
	($genome_build) = $row->[0] =~ m/SNPContigLoc\_(.+)$/;
	push(@genome_builds,$genome_build) if ($genome_build);
     }

    if (scalar(@genome_builds) != 1) {
	if (!scalar(@genome_builds)) {
	    die("Could not find the " . $self->{'dbSNP_version'} . "_SNPContigLoc_NNN table!");
	}
	warn("SNPContigLoc tables for multiple builds found, guessing that the one to use is " . $genome_builds[0]);
    }
    $genome_build = shift(@genome_builds);
    $tablename1 = $self->{'dbSNP_version'} . "_SNPContigLoc_" . $genome_build;
    
     $stmt = qq{ SELECT  name 
                 FROM    $self->{'snp_dbname'}..sysobjects 
                 WHERE name LIKE '$self->{'dbSNP_version'}\_ContigInfo\_$assembly_version\_%'};
               
     my $sth1 = $self->{'dbSNP'}->prepare($stmt);
     $sth1->execute();

     while($row = $sth1->fetchrow_arrayref()) {
       $tablename2 = $row->[0];
     }
     
    $tablename2 = $self->{'dbSNP_version'} . "_ContigInfo_" . $genome_build;
    
     print "SNPContigLoc table is $tablename1 and ContigInfo table is $tablename2\n";
   

    }

     else{
         ## tables names now _SNPContigLoc_105/ _SNPContigLoc_106 for human
         $tablename1 = $self->{'dbSNP_version'} . "_SNPContigLoc";
         $tablename2 = $self->{'dbSNP_version'} . "_ContigInfo";

	 $tablename1 = 'b' .  $tablename1 unless  $tablename1 =~/^b/;
	 $tablename2 = 'b' .  $tablename2 unless  $tablename2 =~/^b/;
     }
    ## hack for human multi-build support
    $tablename1 = $tablename1 . "_108" if $self->{'group_label'} =~ /GRCh38/;
    $tablename2 = $tablename2 . "_108" if $self->{'group_label'} =~ /GRCh38/;

    $tablename1 = $tablename1 . "_105" if $self->{'group_label'} =~ /GRCh37/;
    $tablename2 = $tablename2 . "_105" if $self->{'group_label'} =~ /GRCh37/;

    my $stmt;
    #The group term (the name of the reference assembly in the dbSNP b[version]_SNPContigInfo_[assembly]_[assembly version] table) is either specified via the config file or, if not, attempted to automatically determine from the data
    my $group_term = $self->{'group_term'};
    my $group_label = $self->{'group_label'};

    if (defined($group_term) && defined($group_label)) {
	warn "Using group_term:$group_term and group_label:$group_label to extract mappings \n";
    }
    else{        
        #If no group term was specified, use the one with the most entries in the dbSNP table. This may be wrong though so warn about it.
        $stmt = qq{
            SELECT
                ctg.group_term,
                ctg.group_label,
                COUNT(*) AS N
            FROM
                $tablename1 loc JOIN 
		        $tablename2 ctg ON (
		            ctg.ctg_id = loc.ctg_id
		        )
		    GROUP BY
 		        ctg.group_term,
		        ctg.group_label
            ORDER BY
                N DESC
        };
        my $result = $self->{'dbSNP'}->db_handle->selectall_arrayref($stmt);
        $group_term = $result->[0][0];
        $group_label = $result->[0][1];
        print Progress::location();
        
        #Warn about the group_term we settled for
        debug(
            qq{
                There was no 'group_term' specified in the config file. 
                What this means is that dbSNP maps variations to multiple assemblies. These
                are indicated by the 'group_term' column in the $tablename2 table. We only import
                chromosome labels for the mappings to the reference assembly we have in Ensembl and
                for the other mappings we import the contig label. For example, in human dbSNP132, the
                group_term for GRCh37 is 'GRCh37'. You can specify this with the 'group_term' option in
                the configuration file.
                However, in order to save you the trouble, for now we'll just take the group_term with the
                most entries in the current dbSNP table and hope that it is the correct one. If it's not you'll
                have to truncate the variation_feature table and re-run the variation_feature subroutine.
                
                Based on this approach, we'll use the group term '$group_term' and group label '$group_label'.
            }
        );
    }


    # In the query below, the pre-131 syntax was ref-assembly. In 131 it is GRCh37 for human. What is it for other species??
    #my $group_term = 'ref_';
    #my ($release) = $self->{'dbSNP_version'} =~ m/^b?(\d+)$/;
    #$group_term = 'GRCh' if ($self->{'dbm'}->dbCore()->species =~ m/homo|human/i && $release > 130);
    
	#	     t2.group_term LIKE 'ref_%'
 
    ###Extract mappings for Mitochondria as well as named reference 
    my $extract_mappings_for = qq['$group_term', 'non-nuclear'];
  
## Type 3: DelOnCtg  Deletion on the contig: part of the snp flanking sequence including 
##                   the snp was absent on the contig sequence in the alignment

     $stmt = qq{ SELECT
                   loc.snp_id AS sorting_id, 
                   ctg.contig_acc,
		   ctg.contig_gi,
                   loc.lc_ngbr+2,
		   loc.rc_ngbr,
		   ctg.contig_chr, 
		   CASE WHEN
		     loc.loc_type = 3
		   THEN
		     loc.phys_pos_from+2
		   ELSE
		     loc.phys_pos_from+1
		   END,
		   CASE WHEN
		     loc.loc_type = 3
		   THEN
		     loc.phys_pos_from+1
		   ELSE
		     (loc.phys_pos_from + 1 + loc.asn_to - loc.asn_from )
		   END,
		   CASE WHEN
		     loc.orientation = 1
		   THEN
		     -1
		   ELSE
		     1
		   END, 
                   loc.aln_quality,
                    CASE WHEN
                     ctg.orient = 1
                   THEN
                     -1
                   ELSE
                     1
                   END 
		 FROM 
		   $tablename1 loc JOIN 
		   $tablename2 ctg ON (
		     ctg.ctg_id = loc.ctg_id
		   )
      WHERE
        ctg.group_term in($extract_mappings_for) AND
        ctg.group_label LIKE '$group_label'
       ORDER BY  sorting_id ASC  
	          };


     dumpSQL($self->{'dbSNP'},$stmt, $self->{source_engine} );
    
    
    debug(localtime() . "\tLoading SNPLoc data");

     create_and_load($self->{'dbVar'}, "tmp_contig_loc_chrom", "snp_id i* not_null", "ctg * not_null", "ctg_gi i", "ctg_start i not_null", "ctg_end i", "chr *", "start i", "end i", "strand i", "aln_quality d", "ctg_strand i");
  print Progress::location();


    debug(localtime() . "\tCreating genotyped variations");
    #creating the temporary table with the genotyped variations

    my $gtype_ref = $self->{'dbVar'}->db_handle->selectall_arrayref(qq{SELECT COUNT(*) FROM  tmp_sample_genotype_single_bp});
    my $gtype_row = $gtype_ref->[0][0] if $gtype_ref;
    if ($gtype_row) {
	$self->{'dbVar'}->do(qq{CREATE TABLE tmp_genotyped_var SELECT DISTINCT variation_id FROM tmp_sample_genotype_single_bp});
	print Progress::location();
	$self->{'dbVar'}->do(qq{CREATE UNIQUE INDEX variation_idx ON tmp_genotyped_var (variation_id)});
	print Progress::location();
	$self->{'dbVar'}->do(qq{INSERT IGNORE INTO tmp_genotyped_var SELECT DISTINCT variation_id FROM sample_genotype_multiple_bp});
        print Progress::location();
    }
    debug(localtime() . "\tCreating tmp_variation_feature_chrom data in GenericChromosome");
    #if tcl.end>1, this means we have coordinates for chromosome, we will use it

    ## take account of variant to contig and contig to reference orientation
    dumpSQL($self->{'dbVar'},qq{SELECT v.variation_id, ts.seq_region_id, 
                                      tcl.start,tcl.end,
                                      (tcl.strand * tcl.ctg_strand), v.name, v.source_id, tcl.aln_quality
				      FROM variation v, tmp_contig_loc_chrom tcl, seq_region ts
				      WHERE v.snp_id = tcl.snp_id
				      AND tcl.start>2 
                                      AND tcl.chr = ts.name
    });#to get rid of lots of start=1

    create_and_load($self->{'dbVar'},'tmp_variation_feature_chrom',"variation_id i* not_null","seq_region_id i", "seq_region_start i", "seq_region_end i", "seq_region_strand", "variation_name", "source_id", "aln_quality d");
  print Progress::location();
    
    debug(localtime() . "\tCreating tmp_variation_feature_ctg data in GenericChromosome");
    #if tcl.start = 1 or tcl.end=1, this means we don't have mappings on chromosome, we take ctg coordinates if it is in toplevel
    dumpSQL($self->{'dbVar'},qq{SELECT v.variation_id, ts.seq_region_id, 
                                      tcl.ctg_start,tcl.ctg_end,
                                      tcl.strand, v.name, v.source_id, tcl.aln_quality
				      FROM variation v, tmp_contig_loc_chrom tcl, seq_region ts
				      WHERE v.snp_id = tcl.snp_id
				      AND (
					tcl.start = 1 OR
					tcl.end = 1 OR
					tcl.start IS NULL OR
					tcl.end IS NULL
				      )
                                      AND tcl.ctg = ts.name
   });

    create_and_load($self->{'dbVar'},'tmp_variation_feature_ctg',"variation_id i* not_null","seq_region_id i ", "seq_region_start i", "seq_region_end i", "seq_region_strand", "variation_name", "source_id",  "aln_quality d");
  print Progress::location();

    debug(localtime() . "\tDumping data into variation_feature table in GenericChromosome");
    ## the copy process is slow for large dbs, so run 1M at a time
    my $get_max_sth = $self->{'dbVar'}->prepare(qq[select max(variation_id) from variation]);
    $get_max_sth->execute()||die;
    my $max = $get_max_sth->fetchall_arrayref();

    if ($gtype_row) {
      foreach my $table ("tmp_variation_feature_chrom","tmp_variation_feature_ctg") {

     	my $vf_ins_sth = $self->{'dbVar'}->prepare(qq{INSERT INTO variation_feature (variation_id, seq_region_id,seq_region_start, seq_region_end, seq_region_strand,variation_name, flags, source_id, alignment_quality, somatic)
				  SELECT tvf.variation_id, tvf.seq_region_id, tvf.seq_region_start, tvf.seq_region_end, tvf.seq_region_strand,tvf.variation_name,IF(tgv.variation_id,'genotyped',NULL), tvf.source_id, tvf.aln_quality,  v.somatic
				  FROM $table tvf LEFT JOIN tmp_genotyped_var tgv ON tvf.variation_id = tgv.variation_id
                                  LEFT JOIN variation v on tvf.variation_id = v.variation_id
				  WHERE tvf.variation_id between ? and ?
				  });

       my $batch = 1000000;
       my $start = 1;
       debug(localtime() . "\tStarting binned variation_feature copy");
       while( $start  < $max->[0]->[0] ){

           my $end = $start + $batch ;
           $vf_ins_sth->execute($start, $end)||die;
           $start = $end + 1;
       }
       

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
	 my $vf_ins_sth =$self->{'dbVar'}->prepare(qq{INSERT INTO variation_feature (variation_id, seq_region_id,seq_region_start, seq_region_end, seq_region_strand,variation_name, flags, source_id, alignment_quality, somatic)
 				  SELECT tvf.variation_id, tvf.seq_region_id, tvf.seq_region_start, tvf.seq_region_end, tvf.seq_region_strand,tvf.variation_name,NULL, tvf.source_id, tvf.aln_quality, v.somatic
 				  FROM $table tvf, variation v
                                  WHERE v.variation_id = tvf.variation_id
				  AND tvf.variation_id between ? and ?
 				  });
  
	my $batch = 1000000;
	my $start = 1;
	debug(localtime() . "\tStarting binned variation_feature copy");
	while( $start  < $max->[0]->[0] ){	
	    
	    my $end = $start + $batch ;   
	    $vf_ins_sth->execute($start, $end)||die;
	    $start = $end + 1;
	}
	print Progress::location();
      }
    }

    if ($self->{'dbm'}->dbCore()->species =~ /homo/i){
	debug(localtime() . "\tDumping non-primary assembly data for human");
	$self->extract_haplotype_mappings($tablename1, $tablename2, $group_label);
    }
   
    $self->{'dbVar'}->do("DROP TABLE IF EXISTS tmp_contig_loc_chrom");
    $self->{'dbVar'}->do("DROP TABLE IF EXISTS tmp_genotyped_var");
    $self->{'dbVar'}->do("DROP TABLE IF EXISTS tmp_variation_feature_chrom");
    $self->{'dbVar'}->do("DROP TABLE IF EXISTS tmp_variation_feature_ctg");
    

}

## haplotype & patch mappings are held against the contigs rather than chromosomes

sub extract_haplotype_mappings{
    
    my $self         = shift;
    my $SNPContigLoc = shift;
    my $ContigInfo   = shift;
    my $group_label  = shift;
        

    ## copy synonyms & fake chrom seq_region_ids from core db to temp table
    my $syn_ext_stmt = qq[ select sr2.seq_region_id, sr2.name, srs.synonym, assembly.asm_start, assembly.asm_end
                           from  seq_region sr1, seq_region_synonym srs, seq_region sr2, assembly
                           where sr1.seq_region_id = srs.seq_region_id 
                           and sr2.seq_region_id = assembly.asm_seq_region_id
                           and sr1.seq_region_id = assembly.cmp_seq_region_id
                         ];


    dumpSQL($self->{'dbCore'},$syn_ext_stmt);       
    
    create_and_load($self->{'dbVar'}, "tmp_hap_synonym", "seq_region_id i* not_null", "name * not_null", "synonym * not_null", "asm_start i", "asm_end i"   );

    ### check for reverse strand patches - not expecting any, so merely warn
    my $strand_ch_sth = $self->{'dbVar'}->prepare(qq[ select seq_region_id, name from tmp_hap_synonym where asm_start > asm_end]);
    $strand_ch_sth->execute();
    my $problem = $strand_ch_sth->fetchall_arrayref();
    warn "DANGER: Patch id $problem->[0]->[0] name $problem->[0]->[1] on reverse strand\n" if defined $problem->[0]->[0];


 my $concat_syntax ;
    if($self->{source_engine} =~/mssql/){
	$concat_syntax = qq[ ctg.genbank_acc + '.' + cast(ctg.genbank_ver as nvarchar(10)) ];
    }
    elsif($self->{source_engine} =~/postgreSQL/i){
	warn "setting search path for postgreSQL: " .$self->{'schema_name'} . "\n";
	my $sth = "SET search_path TO $self->{'schema_name'},dbsnp_main,public";
	$concat_syntax = qq[ ctg.genbank_acc ||  '.' || ctg.genbank_ver   ];
    }
    else{
	$concat_syntax = qq[  CONCAT(ctg.genbank_acc, '.', ctg.genbank_ver)  ];
    }
    

    ## extract variant mappings on patches/haplotypes
    my $dat_ext_stmt =  qq[ SELECT loc.snp_id,                            
                             $concat_syntax,
                            loc.lc_ngbr+2, 
                            loc.rc_ngbr,                            
                            CASE WHEN loc.orientation = 1
                            THEN   -1
                            ELSE  1
                            END,
                            loc.aln_quality FROM  $SNPContigLoc loc JOIN  $ContigInfo ctg ON ( ctg.ctg_id = loc.ctg_id )
                            WHERE  ctg.group_term !=   'Primary_Assembly' and ctg.group_term !=   'non-nuclear'
                            and ctg.group_label LIKE '$group_label'                
                           ];

    
    dumpSQL($self->{'dbSNP'},$dat_ext_stmt);
        
    debug(localtime() . "\tLoading SNPLoc data for haplotypes");

    create_and_load($self->{'dbVar'}, "tmp_contig_loc_hap", "snp_id i* not_null", "seq_av * not_null",  "seq_start i not_null", "seq_end i", "strand i", "aln_quality d");
    print Progress::location();

    debug(localtime() . "\tAdding VariationFeature data for haplotypes");
    ### copy to variation_feature table
    $self->{'dbVar'}->do(qq[ INSERT INTO variation_feature (variation_id, seq_region_id,seq_region_start, seq_region_end, seq_region_strand,
                                                            variation_name, source_id, alignment_quality, somatic)
                                                                   SELECT v.variation_id, ths.seq_region_id, tvf.seq_start+ths.asm_start-1, tvf.seq_end+ths.asm_start-1, tvf.strand, 
                                         v.name, v.source_id, tvf.aln_quality,  v.somatic
                                  FROM tmp_contig_loc_hap tvf  
                                  LEFT JOIN variation v on tvf.snp_id = v.snp_id
                                   LEFT JOIN tmp_hap_synonym ths on ths.synonym = tvf.seq_av 
                                  WHERE ths.seq_region_id is not null            
                                  ]);

    debug(localtime() . "\tDone VariationFeature data for haplotypes");

}

1;
