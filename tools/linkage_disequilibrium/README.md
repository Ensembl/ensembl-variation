Linkage Disequilibrium (LD) Calculator
=====================================
The script calculates LD between variants using genotypes from a selected population. We only support LD calculation for variants for which we have genotypes from at least 40 samples in the selected population. At the moment we only have sufficient amounts of genotype data from the 1000 Genomes project for human.

There is a [web tool](http://www.ensembl.org/Multi/Tools/LD) version with an interface for easy parameter selection. 

##### Table of contents
* [Usage](#usage)
* [Show all populations](#show_all_populations)
* [LD calculations](#ld_calculations)
* [Input file formats](#input_file_formats)
* [Output file formats](#output_file_formats)
* [Examples](#examples)
---
<a name="usage"></a>
### Usage
The script wraps around the LD calculation code in the Ensembl Variation API. Before running the ld_tool script follow the [installation instructions](https://github.com/Ensembl/ensembl-variation/blob/master/C_code/README.txt).

#### Running the script
```bash
./ld_tool -i input_file --calculation [center|pairwise|region] 
```
* `-i|--input_file [input_file]` : [input_file](#input_file_formats)
* `-c|--calculation [center|pairwise|region]` : [calculation](#ld_calculation)

#### Script options
* `--species [species]` : species to use (default: human)
* `--assembly [assembly]` : assembly to use (If human is used default assembly is grch38)
* `--output_file [output_file]` : If not provided output_file name is set to [input_file].out or if input_file is not provided to [calculation]_year_month_day_hour_minutes.out
* `--warnings_file [warnings_file]` : If not provided warnings_file name is set to [input_file].warnings or if input_file is not provided to [calculation]_year_month_day_hour_minutes.warnings
* `--population [population_name]` : List of population(s) for LD calculation
* `--region [region]` : List of region(s) as input for region calculation. A region is defined as chromosome:start-end
* `--variant [variant]` : List of variant(s) as input for either center or pairwise calculations
* `--r2 [r2]` : Only include variants to the result whose r<sup>2</sup> value is greater than or equal to the given value. r<sup>2</sup> needs to be in the range of 0.0 and 1.0.
* `--d_prime [d_prime]` : Only include variants to the result whose D<sup>'</sup> value is greater than or equal to the given value. D<sup>'</sup> needs to be in the range of 0.0 and 1.0.
* `--add_variant_attribs` : Add variant attributes (evidence values and consequence type) to the output

<a name="show_all_populations"></a>
### Show all populations
Show all populations that can be used for the given species. 
```bash
./ld_tool --show_all_populations --species species
```
<a name="ld_calculations"></a>
### LD calculations
The script supports three different types of calculations.
#### Center
Compute all pairwise LD values for a given variant with all variants within a window of the specified size, centered on the input variant.
#### Pairwise
Compute all pairwise LD values for a list of variants.
#### Region
Compute all pairwise LD values for all variants in a given region. 

<a name="input_file_formats"></a>
### Input file formats
#### Variant identifiers
Provide a list of variant identifiers for center and pairwise calculations.
```txt
rs6778
rs66778
```
#### Regions
Provide a list of regions for the region calculation.
```txt
12 10824960 11171544
```
<a name="output_file_formats"></a>
### Output file format
Default output format is a list of tab-separated values:
```txt
first variant name
first variant location
second variant name
second variant location
r2
d_prime
population name
```
If `--add_variant_attributes` variant evidence and variant consequence are added to the output:
```txt
first variant name
first variant location
second variant name
second variant location
r2
d_prime
population name
first variant consequence
first variant evidence
second variant consequence
second variant evidence
```
<a name="examples"></a>
### Examples 
#### Calculate LD in a region
```bash
./ld_tool -calculation region --region 12:10824960-11171544 --population 1000GENOMES:phase_3:CEU
```
First 5 rows of the output
```txt
rs36161147  12:11102008-11102007  rs10700124  12:11104454 0.181909  0.999963  1000GENOMES:phase_3:CEU
rs2600335   12:11099321           rs67181394  12:11121206 1.000000  1.000000  1000GENOMES:phase_3:CEU
rs570660067 12:10827629           rs669503    12:10830475 0.066715  0.999998  1000GENOMES:phase_3:CEU
rs67181394  12:11121206-11121209  rs66547886  12:11124355 1.000000  1.000000  1000GENOMES:phase_3:CEU
rs10845219  12:10825055           rs762496999 12:10829216 0.073683  0.999991  1000GENOMES:phase_3:CEU
```
#### Calculate pairwise LD for a list of variants
```bash
./ld_tool -calculation pairwise --variant rs5961768 rs112287390 rs7058648 --population 1000GENOMES:phase_3:ASW
```
```txt
rs112287390 X:5636653 rs5961768 X:5637210 1.000000  1.000000  1000GENOMES:phase_3:ASW
rs112287390 X:5636653 rs7058648 X:5637986 1.000000  1.000000  1000GENOMES:phase_3:ASW
rs5961768   X:5637210 rs7058648 X:5637986 1.000000  1.000000  1000GENOMES:phase_3:ASW
```
#### Calculate pairwise LD for a variant and all variants within a 100 kb window of that variant
```bash
./ld_tool -calculation center -variant rs2708377 --population 1000GENOMES:phase_3:PEL --d_prime 1.0
```
First 5 rows of the output
```txt
rs2708377 12:11063716 rs2708336   12:11120807 1.000000  1.000000  1000GENOMES:phase_3:PEL
rs2708377 12:11063716 rs2255418   12:11064373 1.000000  1.000000  1000GENOMES:phase_3:PEL
rs2708377 12:11063716 rs3911150   12:11049649 0.562500  1.000000  1000GENOMES:phase_3:PEL
rs2708377 12:11063716 rs67679362  12:11080572 0.479167  1.000000  1000GENOMES:phase_3:PEL
rs2708377 12:11063716 rs35318883  12:11157151 1.000000  1.000000  1000GENOMES:phase_3:PEL
```

