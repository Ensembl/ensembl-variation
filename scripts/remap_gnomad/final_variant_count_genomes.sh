#!/bin/bash
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2023] EMBL-European Bioinformatics Institute
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
dir=/variation/data/gnomAD/v2.1/
for chr in {1..22} 'X'
do
  i_37=`bcftools index --nrecords $dir/grch37/genomes/gnomad.genomes.r2.1.sites.chr${chr}_noVEP.vcf.gz`
  i_unmap=`bcftools index --nrecords /hps/nobackup2/production/ensembl/anja/gnomad/Genomes/mapping_results/unmapped/gnomad.genomes.r2.1.sites.grch38.chr${chr}_noVEP.vcf.unmap.gz`
  i_crossmap_mapped=`bcftools index --nrecords /hps/nobackup2/production/ensembl/anja/gnomad/Genomes/mapping_results/unmapped/gnomad.genomes.r2.1.sites.grch38.chr${chr}_noVEP.vcf.gz`
  i_ensembl_unique_map=`wc -l < /hps/nobackup2/production/ensembl/anja/release_96/human/remap_gnomad/genomes/unique_mappings/chrom${chr}.txt`
  i_successfully_mapped=$((i_crossmap+i_ensembl_unique_map))
  i_end_file=`bcftools index --nrecords $dir/grch38/genomes/gnomad.genomes.r2.1.sites.grch38.chr${chr}_noVEP.vcf.gz`
  i_qc_add=$((i_crossmap+i_unmap))
  i_gnomad_unmapped=$((i_unmap-i_ensembl_unique_map))
  i_final_count=$((i_end_file+i_gnomad_unmapped))
  echo $chr $i_37 $i_crossmap $i_crossmap_mapped
done
