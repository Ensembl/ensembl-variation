-- Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

To limit dbsnp sample_ids not in thses wild derived strains, use this sql:
mysql> create table tmp_strain_gtype_final like tmp_strain_39_dbsnp;
mysql> insert into tmp_strain_gtype_final select s.* from tmp_strain_39_dbsnp s where s.sample_id1 not in (174,215,198,230,178) and s.sample_id2 not in (174,215,198,230,178);

Five wild strain names are: CAST/EiJ,PWD/PhJ,MOLF/EiJ,WSB/EiJ,CZECHll/Ei
next time, we could limit this 5 strains when we generate table:tmp_strain_39_dbsnp

to validate repeats_SNPs are they occur in dbSNP, by using sql:
mysql> select count(*) from tmp_ids1 t, var_name_repeat vn, variation v where vn.variation_name=v.name and v.variation_id=t.variation_id1;
to validate non-repeats_SNPs are they occur in dbSNP, by using sql:
select count(*) from variation v left join var_name_repeat vn on v.name=vn.variation_name, tmp_ids1 t where v.source_id=2 and v.variation_id=t.variation_id2 and vn.variation_name is null;
+----------+
| count(*) |
+----------+
|  1868743 |
+----------+
1 row in set (1 min 55.75 sec)

then to merge dbSNP and sanger strain_gtype table by using:
insert ignore into tmp_strain_gtype_final select * from tmp_strain_39_sanger;
because variation_id in sanger and dbsnp is different, so the effect is adding rather than merge.

count how many SNPs between a pair of sample_ids, merge same variation between sanger and dbsnp:
CREATE TABLE tmp_strain_gtype_count_final 
 				  SELECT sg1.sample_id1 as sample_id1,sg2.sample_id2 as sample_id2,count(*) as count 
 				  FROM tmp_strain_gtype_final sg1 
                                        LEFT JOIN tmp_ids1 t1 on sg1.variation_id = t1.variation_id1,  
                                      tmp_strain_gtype_final sg2
 				  WHERE sg1.variation_id=sg2.variation_id 
 				  AND sg1.sample_id1 != sg2.sample_id2 
 				  AND sg1.sample_id1 = sg2.sample_id1 
 				  AND sg1.sample_id2 = sg2.sample_id2 
                                  AND t1.variation_id1 is NULL #to avoid duplicate variation from both dbsnp and sanger
 				  GROUP BY sg1.sample_id1,sg2.sample_id2});
     $self->{'dbVariation'}->do("ALTER TABLE tmp_strain_gtype_count_final ADD INDEX sample_idx(sample_id1,sample_id2)");
Query OK, 6578 rows affected (29 min 33.50 sec)
Records: 6578  Duplicates: 0  Warnings: 0




mysql> create table tmp_ids1 select vf1.variation_id as variation_id1, vf2.variation_id as variation_id2 from variation_feature vf1, variation_feature vf2 where vf1.seq_region_id=vf2.seq_region_id and vf1.seq_region_start=vf2.seq_region_start and vf1.seq_region_end=vf2.seq_region_end and vf1.source_id=1 and vf2.source_id=2;
Query OK, 2601215 rows affected (1 min 59.73 sec)
Records: 2601215  Duplicates: 0  Warnings: 0

CREATE TABLE tmp_strain_gtype_final (
                                        variation_id int not null,
                                        allele_string char(3) not null,
                                        sample_id1 int(11) not null,
                                        sample_id2 int(11) not null,
                                        sample_name1 varchar(50) not null,
                                        sample_name2 varchar(50) not null,
                                        unique key variation_idx(variation_id,sample_id1,sample_id2),
                                        key sample_idx(sample_id1,sample_id2)) MAX_ROWS = 100000000
                                   });
do this separately for dbsnp and sanger snps, 
for dbsnps:
mysql> insert ignore into tmp_strain_39_dbsnp SELECT ig1.variation_id,concat(ig1.allele_1,"/",ig2.allele_1) as allele_string,ip1.population_sample_id as sample_id1,ip2.population_sample_id as sample_id2, s1.name as sample_name1, s2.name as sample_name2 FROM population p1, population p2,tmp_individual_genotype_single_bp ig1, tmp_individual_genotype_single_bp ig2,individual_population ip1,individual_population ip2, sample s1, sample s2 WHERE ig1.variation_id=ig2.variation_id and ig1.allele_1=ig1.allele_2 AND ig2.allele_1 = ig2.allele_2 AND ig1.sample_id=ip1.individual_sample_id AND ig2.sample_id=ip2.individual_sample_id AND p1.sample_id = ip1.population_sample_id AND p2.sample_id = ip2.population_sample_id AND ig1.allele_1 != ig2.allele_1 and p1.sample_id != p2.sample_id and p1.is_strain=1 and p2.is_strain=1 and s1.sample_id = p1.sample_id AND s2.sample_id = p2.sample_id and ig1.variation_id <=6491546 and ig2.variation_id <=6491546;
Query OK, 141999022 rows affected (5-9 hours 58 min 35.82 sec depends on how busy is the matchine)
Records: 142921258  Duplicates: 922236  Warnings: 0

for sanger snps:
for sanger against sanger, it should be as the one above:
mysql> insert ignore into tmp_strain_39_sanger SELECT ig1.variation_id,concat(ig1.allele_1,"/",ig2.allele_1) as allele_string, ip1.population_sample_id as sample_id1,ip2.population_sample_id as sample_id2, s1.name as sample_name1, s2.name as sample_name2 FROM population p1, population p2,tmp_individual_genotype_single_bp ig1, tmp_individual_genotype_single_bp ig2,individual_population ip1,individual_population ip2, sample s1, sample s2 WHERE ig1.variation_id=ig2.variation_id and ig1.allele_1=ig1.allele_2 AND ig2.allele_1 = ig2.allele_2 AND ig1.sample_id=ip1.individual_sample_id AND ig2.sample_id=ip2.individual_sample_id  AND p1.sample_id = ip1.population_sample_id AND p2.sample_id = ip2.population_sample_id AND ig1.allele_1 != ig2.allele_1 and  p1.sample_id != p2.sample_id and p1.is_strain=1 and p2.is_strain=1 and s1.sample_id = p1.sample_id AND s2.sample_id = p2.sample_id and ig1.variation_id >6491546 and ig2.variation_id >6491546;
Query OK, 13530554 rows affected (1 hour 30 min 59.93 sec)
Records: 20549764  Duplicates: 7019210  Warnings: 0

for sanger against dbsnp one way:
mysql> insert ignore into tmp_strain_39_sanger SELECT ig1.variation_id,concat(ig1.allele_1,"/",ig2.allele_1) as allele_string,ip1.population_sample_id as sample_id1,ip2.population_sample_id as sample_id2, s1.name as sample_name1, s2.name as sample_name2 FROM population p1, population p2,tmp_individual_genotype_single_bp ig1, tmp_individual_genotype_single_bp ig2,individual_population ip1,individual_population ip2, sample s1, sample s2, tmp_ids1 ids WHERE ig1.variation_id=ids.variation_id2 and ids.variation_id1=ig2.variation_id and ig1.allele_1=ig1.allele_2 AND ig2.allele_1 = ig2.allele_2 AND ig1.sample_id=ip1.individual_sample_id AND ig2.sample_id=ip2.individual_sample_id AND p1.sample_id = ip1.population_sample_id AND p2.sample_id = ip2.population_sample_id AND ig1.allele_1 != ig2.allele_1 and p1.sample_id != p2.sample_id and p1.is_strain=1 and p2.is_strain=1 and s1.sample_id = p1.sample_id AND s2.sample_id = p2.sample_id ;

for sanger against dbsnp another way back:
mysql> insert ignore into tmp_strain_39_sanger SELECT ig2.variation_id , concat(ig1.allele_1,"/",ig2.allele_1) as allele_string,ip1.population_sample_id as sample_id1,ip2.population_sample_id as sample_id2, s1.name as sample_name1, s2.name as sample_name2 FROM population p1, population p2,tmp_individual_genotype_single_bp ig1, tmp_individual_genotype_single_bp ig2,individual_population ip1,individual_population ip2, sample s1, sample s2, tmp_ids1 ids WHERE ig2.variation_id=ids.variation_id2 and ids.variation_id1=ig1.variation_id and ig1.allele_1=ig1.allele_2 AND ig2.allele_1 = ig2.allele_2 AND ig1.sample_id=ip1.individual_sample_id AND ig2.sample_id=ip2.individual_sample_id AND p1.sample_id = ip1.population_sample_id AND p2.sample_id = ip2.population_sample_id AND ig1.allele_1 != ig2.allele_1 and p1.sample_id != p2.sample_id and p1.is_strain=1 and p2.is_strain=1 and s1.sample_id = p1.sample_id AND s2.sample_id = p2.sample_id ;
Query OK, 11753142 rows affected (4 hours 31 min 54.28 sec)
Records: 17276574  Duplicates: 5523432  Warnings: 0

sanger_strain.sample_id in (139,143,147,165,168,182,202,232);
