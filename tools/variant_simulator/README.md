Variant Simulator
=====================================
The script generates all single base substitutions (SNPs) for protein_coding gene(s) given a specific species or chromosome or gene. One can restrict the SNPs to be generated only for the introns, exons or only coding exons and a specific number of bases around each of them.

##### Table of contents
* [Usage](#usage)
* [Examples](#examples)
---
<a name="usage"></a>
### Usage
The script generates all single base substitutions for protein coding genes given a specific species or chromosome or gene. The script uses the Ensembl API and Ensembl databases.

#### Running the script
```bash
./simulate_variation -gene BRCA2
```

#### Script options
* `-species [species]` : species to use (default: human)
* `-assembly [assembly]` : assembly to use if species is human (default: grch38)
* `-refseq`: use RefSeq genes/transcripts if species is human
* `-registry [registry_file]` : [registry_file](https://www.ensembl.org/info/docs/api/registry.html)
* `-chrom [chromosome]` : generate SNPs only for specified chromosome
* `-gene [gene_ensembl_stable_id | gene_symbol ]` : generate SNPs only for specified gene (if chromosome and gene are both specified but don't agree, the gene options is used)
* `-spliceai`: generate files for SpliceAI tool
* `-exonsOnly` : generate SNPs only for all exons of the protein_coding genes
* `-intronsOnly` : generate SNPs only for all introns of the protein_coding genes
* `-codingOnly` : generate SNPs only for translatable exons of protein_coding transcripts
* `-onlyMane` : generate SNPs only for MANE transcripts
* `-edge` : upstream/downstream bp of each feature (default: 0)
* `-output|-o [output_file]` : output vcf file (default: simulated.vcf)
* `-help` : print usage message

<a name="examples"></a>
### Examples 
#### Generate SNPs for a chromosome
```bash
./simulate_variation -chrom 2
```
First 7 rows of the output:
```txt
##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol or Gene stable id">
##INFO=<ID=FEATURE,Number=1,Type=String,Description="Feature id">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
2	38814	2-38814-T-A	T	A	.	.	GENE=FAM110C;FEATURE=ENSG00000184731
2	38814	2-38814-T-C	T	C	.	.	GENE=FAM110C;FEATURE=ENSG00000184731
2	38814	2-38814-T-G	T	G	.	.	GENE=FAM110C;FEATURE=ENSG00000184731
```

#### Generate SNPs for a gene
```bash
./simulate_variation -gene ENSG00000139618
```
or
```bash
./simulate_variation -gene BRCA2
```
First 7 rows of the output:
```txt
##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol or Gene stable id">
##INFO=<ID=FEATURE,Number=1,Type=String,Description="Feature id">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
13	32315474	13-32315474-G-A	G	A	.	.	GENE=BRCA2;FEATURE=ENSG00000139618
13	32315474	13-32315474-G-C	G	C	.	.	GENE=BRCA2;FEATURE=ENSG00000139618
13	32315474	13-32315474-G-T	G	T	.	.	GENE=BRCA2;FEATURE=ENSG00000139618
```

#### Generate SNPs for a gene using exonsOnly
```bash
./simulate_variation -gene BRCA2 -exonsOnly
```
First 7 rows of the output:
```txt
##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol or Gene stable id">
##INFO=<ID=FEATURE,Number=1,Type=String,Description="Feature id">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
13	32357742	13-32357742-C-A	C	A	.	.	GENE=BRCA2;FEATURE=ENSE00003719469
13	32357742	13-32357742-C-T	C	T	.	.	GENE=BRCA2;FEATURE=ENSE00003719469
13	32357742	13-32357742-C-G	C	G	.	.	GENE=BRCA2;FEATURE=ENSE00003719469
```

#### Generate SNPs for a gene using codingOnly exons
```bash
./simulate_variation -gene BRCA2 -codingOnly
```
First 7 rows of the output:
```txt
##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol or Gene stable id">
##INFO=<ID=FEATURE,Number=1,Type=String,Description="Feature id">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
13	32325076	13-32325076-G-A	G	A	.	.	GENE=BRCA2;FEATURE=ENSE00003659301
13	32325076	13-32325076-G-C	G	C	.	.	GENE=BRCA2;FEATURE=ENSE00003659301
13	32325076	13-32325076-G-T	G	T	.	.	GENE=BRCA2;FEATURE=ENSE00003659301
```

#### Generate SNPs for a gene using codingOnly exons with 5bp upstream/downstream of each exon
```bash
./simulate_variation -gene BRCA2 -codingOnly -edge 5
```
First 7 rows of the output:
```txt
##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol or Gene stable id">
##INFO=<ID=FEATURE,Number=1,Type=String,Description="Feature id">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
13	32325071	13-32325071-T-A	T	A	.	.	GENE=BRCA2;FEATURE=ENSE00003659301
13	32325071	13-32325071-T-C	T	C	.	.	GENE=BRCA2;FEATURE=ENSE00003659301
13	32325071	13-32325071-T-G	T	G	.	.	GENE=BRCA2;FEATURE=ENSE00003659301
```

#### Generate SNPs for a specific species
```bash
./simulate_variation -species pig -chrom 12
```
First 7 rows of the output:
```txt
##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol or Gene stable id">
##INFO=<ID=FEATURE,Number=1,Type=String,Description="Feature id">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
12	128254	12-128254-T-A	T	A	.	.	GENE=ENSSSCG00000034696;FEATURE=ENSSSCG00000034696
12	128254	12-128254-T-C	T	C	.	.	GENE=ENSSSCG00000034696;FEATURE=ENSSSCG00000034696
12	128254	12-128254-T-G	T	G	.	.	GENE=ENSSSCG00000034696;FEATURE=ENSSSCG00000034696
```
