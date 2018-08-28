Variant Simulator
=====================================
The script generates all single base substitutions (SNPs) for protein_coding gene(s) given a specific species or chromosome or gene. One can restrict the SNPs to be generated only for the introns, exons or only coding exons and a specific number of bases around each of them.

##### Table of contents
* [Usage](#usage)
* [Examples](#examples)

---
<a name="usage"></a>
### Usage
The script generates all single base substitutions for protein coding genes given a specific species or chromosome or gene. The script uses the Esembl API and Ensembl databases.

#### Running the script
```bash
./simulate_variation.pl -registry registry_file
```
* `-registry [registry_file]` : [registry_file](https://www.ensembl.org/info/docs/api/registry.html)

#### Script options
* `-species [species]` : species to use (default: human)
* `-chrom [chromosome]` : generate SNPs only for specified chromosome
* `-gene [gene_ensembl_stable_id | gene_symbol ]` : generate SNPs only for specified gene (if chromosome and gene are both specified but don't agree, the gene options is used)
* `-exonsOnly` : generate SNPs only for all exons of the protein_coding genes
* `-intronsOnly` : generate SNPs only for all introns of the protein_coding genes
* `-codingOnly` : generate SNPs only for translatable exons of protein_coding transcripts
* `-edge` : upstream/downstream bp of each feature (default: 0)
* `-output|-o [output_file]` : output vcf file (default: simulated.vcf)
* `-help` : print usage message


<a name="examples"></a>
### Examples 
#### Generate SNPs for a chromosome
```bash
./simulate_variation.pl -registry registry_file -chrom 2
```
First 6 rows of the output:

```txt
##fileformat=VCFv4.1
##INFO=gene_name:feature_id
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
2	38814	2-38814-T-A	T	A	.	.	FAM110C:ENSG00000184731
2	38814	2-38814-T-C	T	C	.	.	FAM110C:ENSG00000184731
2	38814	2-38814-T-G	T	G	.	.	FAM110C:ENSG00000184731
```

#### Generate SNPs for a gene
```bash
./simulate_variation.pl -registry registry_file -gene ENSG00000139618
```
or

```bash
./simulate_variation.pl -registry registry_file -gene BRCA2
```
First 6 rows of the output:

```txt
##fileformat=VCFv4.1
##INFO=gene_name:feature_id
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
13	32315474	13-32315474-G-A	G	A	.	.	BRCA2:ENSG00000139618
13	32315474	13-32315474-G-C	G	C	.	.	BRCA2:ENSG00000139618
13	32315474	13-32315474-G-T	G	T	.	.	BRCA2:ENSG00000139618
```

#### Generate SNPs for a gene using exonsOnly
```bash
./simulate_variation.pl -registry registry_file -gene BRCA2 -exonsOnly
```
First 6 rows of the output:

```txt
##fileformat=VCFv4.1
##INFO=gene_name:feature_id
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
13	32357742	13-32357742-C-A	C	A	.	.	BRCA2:ENSE00003719469
13	32357742	13-32357742-C-T	C	T	.	.	BRCA2:ENSE00003719469
13	32357742	13-32357742-C-G	C	G	.	.	BRCA2:ENSE00003719469
```

#### Generate SNPs for a gene using codingOnly exons
```bash
./simulate_variation.pl -registry registry_file -gene BRCA2 -codingOnly
```
First 6 rows of the output:

```txt
##fileformat=VCFv4.1
##INFO=gene_name:feature_id
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
13	32325076	13-32325076-G-A	G	A	.	.	BRCA2:ENSE00003659301
13	32325076	13-32325076-G-C	G	C	.	.	BRCA2:ENSE00003659301
13	32325076	13-32325076-G-T	G	T	.	.	BRCA2:ENSE00003659301
```

#### Generate SNPs for a gene using codingOnly exons with 5bp upstream/downstream of each exon
```bash
./simulate_variation.pl -registry registry_file -gene BRCA2 -codingOnly -edge 5
```
First 6 rows of the output:

```txt
##fileformat=VCFv4.1
##INFO=gene_name:feature_id
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
13	32325071	13-32325071-T-A	T	A	.	.	BRCA2:ENSE00003659301
13	32325071	13-32325071-T-C	T	C	.	.	BRCA2:ENSE00003659301
13	32325071	13-32325071-T-G	T	G	.	.	BRCA2:ENSE00003659301
```
