Phenotype Annotation
=====================================
The script retrieves phenotype annotations from Ensembl for a given list of genes OR phenotypes OR rsIDs OR VCF file OR VEP output file.

##### Table of contents
* [Usage](#usage)
* [Examples](#examples)
---
<a name="usage"></a>
### Usage
The script is a wrapper of the Ensembl Variation API to retrieve the phenotype annotations.

#### Running the script
```bash
./phenotype_annotation -input list_genes.txt -inputType gene
```

#### Script options
* `-input_file | i [input_file]`: input file (list or VCF or VEP output)
* `-inputType [type]`: type of input data: [variant, gene, phenotype, vep] (Default: gene)
* `-species [species]`: species to use (default: human)
* `-assembly [assembly]`: assembly to use if species is human (default: grch38)
* `-registry [registry_file]`: [registry_file](https://www.ensembl.org/info/docs/api/registry.html)
* `-associated`: if inputType is `gene`, include GWAS variant-phenotype associations associated where the gene is reported as a likely candidate
* `-overlap`: if inputType is `gene`, include phenotypes overlapping the input data
* `-all`: if inputType is `gene`, annotate all phenotypes (directly assigned, associated and overlapping)
* `-output_file|-o [output_file]`: output file (default: phenotype_annotation_output.txt)
* `-output_format [format]`: format of the output data: [txt, bed, vcf] (default: txt)
* `-help`: print usage message

<a name="examples"></a>
### Examples 
#### Annotate a list of genes 
##### Annotate a list of genes with directly assigned phenotypes
```bash
./phenotype_annotation -input_file list_genes.txt -inputType gene
```

##### Annotate a list of genes including overlaps
```bash
./phenotype_annotation -input_file list_genes.txt -inputType gene -overlap
```

##### Annotate a list of genes including associated phenotypes
```bash
./phenotype_annotation -input_file list_genes.txt -inputType gene -overlap
```

##### Annotate a list of genes with all linked phenotypes (direct, associated and overlapping)
```bash
./phenotype_annotation -input_file list_genes.txt -inputType gene -all
```

##### Annotate a list of genes and write output BED
```bash
./phenotype_annotation -input_file list_genes.txt -inputType gene -o output.bed -output_format bed
```


#### Annotate a list of phenotypes
```bash
./phenotype_annotation -input_file list_phenotypes.txt -inputType phenotype
```


#### Annotate a list of variants
##### Annotate a list of rsIDs
```bash
./phenotype_annotation -input_file list_rsIDs.txt -inputType variant
```
##### Annotate a VCF
```bash
./phenotype_annotation -input_file input.vcf -inputType vcf
```
##### Annotate a VEP output
```bash
./phenotype_annotation -input_file vep_output.txt -inputType vep
```
