#!/bin/bash
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2025] EMBL-European Bioinformatics Institute
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
i=$1
working_dir=$2
CrossMap.py vcf \
$working_dir/GRCh37_to_GRCh38.chain.gz \
$working_dir/Genomes/gnomad.genomes.r2.1.sites.chr${i}_noVEP.vcf.gz \
$working_dir/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
$working_dir/Genomes/mapping_results/gnomad.genomes.r2.1.sites.grch38.chr${i}_noVEP.vcf
