#!/bin/bash
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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
dir=$2
vcf_file=${dir}/gnomad.genomes.r2.1.sites.grch38.chr${i}_noVEP.vcf.unmap
vcf_file_gz=${dir}/gnomad.genomes.r2.1.sites.grch38.chr${i}_noVEP.vcf.unmap.gz
vcf-sort -t ${dir} < $vcf_file | bgzip > $vcf_file_gz
