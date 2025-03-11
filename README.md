# DupGen_finder

The DupGen_finder was developed to identify different modes of duplicated gene pairs. [MCScanX](http://chibba.pgml.uga.edu/mcscan2/) algorithm was incorporated in this pipeline.

| | |
| --- | --- |
| Authors | Xin Qiao ([Xin Qiao](https://github.com/qiao-xin)) |
| | Qionghou Li ([Qionghou Li](https://github.com/LQHHHHH)) |
| | Yupeng Wang ([Yupeng Wang](https://github.com/wyp1125)) |
| | Andrew H. Paterson ([PGML](http://www.plantgenome.uga.edu)) |
| Email   | <qiaoxin@njau.edu.cn> |

## The schematic diagram of DupGen_finder pipeline
<p align="center">
<img src="https://github.com/LQHHHHH/DupGen_finder/blob/master/data/The%20flowchart%20of%20DupGen_finder%20pipeline.png"  height="500" width="300">
<p align="center">

<t align="center">
  Figure 1: The flowchart of DupGen_finder pipeline
</t>

## Contents
* [Dependencies](#dependencies)
* [Installation](#installation)
* [Preparing input files](#preparing-input-files)
* [Running](#running)
* [Result Files](#result-files)
* [Citation](#citation)

## Dependencies

- [Perl](https://www.perl.org)

## Installation

```bash
cd ~/software  # or any directory of your choice
git clone https://github.com/qiao-xin/DupGen_finder.git
cd DupGen_finder
make
chmod 775 DupGen_finder.pl
chmod 775 DupGen_finder-unique.pl
chmod 775 set_PATH.sh
source set_PATH.sh
```

Test you can run DupGen_finder:
```bash
DupGen_finder.pl
```
DupGen_finder should print its 'help' text.

**\*\*Note\*\***

DupGen_finder software package includes a custom MCScanX algorithm which can output sorted `gff` files such as `Ath.gff.sorted`, and is slightly different from the original MCScanX algorithm implemented in [MCScanX](http://chibba.pgml.uga.edu/mcscan2/) software package.

## Preparing input files

Pre-computed BLAST results (-outfmt 6) and gene location information (GFF format) are required for running DupGen_finder successfully.

1. For the target genome in which gene duplicaiton modes will be classified, please prepare two input files:
	- ```target_species.gff```, a gene position file for the target species, following a tab-delimited format. For example, "Ath.gff".
	- ```target_species.blast```, a blastp output file (-outfmt 6) for the target species (self-genome comparison). For example, "Ath.blast".

2. For the outgroup genome, please prepare two input files:
	- ```[target_species]_[outgroup_species].gff```, a gene position file for the target_species and outgroup_species, following a tab-delimited format.
	- ```[target_species]_[outgroup_species].blast```, a blastp output file (-outfmt 6) between the target and outgroup species (cross-genome comparison).

For example, assuming that you are going to classify gene duplication modes in *Arabidopsis thaliana* (Ath), using *Nelumbo nucifera* (Nnu) as outgroup, you need to prepare 4 input files: ```Ath.gff```,```Ath.blast```, ```Ath_Nnu.gff```, ```Ath_Nnu.blast```

**gff file**

```Ath.gff``` is in the following format (tab separated):
```
Species_abbrev-Chr_ID	gene_ID	start_position	end_position
```

The data in ```Ath.gff``` looks like this (tab separated):
```
Ath-Chr1	AT1G01010.1	3631	5899
Ath-Chr1	AT1G01020.1	5928	8737
Ath-Chr1	AT1G01030.1	11649	13714
Ath-Chr1	AT1G01040.2	23416	31120
Ath-Chr1	AT1G01050.1	31170	33153
```

The below command can be used to creat ```Ath_Nnu.gff```:
```bash
cat Ath.gff Nnu.gff >Ath_Nnu.gff
```

The data in ```Ath_Nnu.gff``` looks like this (tab separated):
```
Ath-Chr1	AT1G01010.1	3631	5899
Ath-Chr1	AT1G01020.1	5928	8737
...
Nnu-megascaffold_32	NNU_00001-RA	57481	63788
Nnu-megascaffold_32	NNU_00002-RA	32491	41125
...
```

**blast file**

```Ath.blast``` is in the following format:
```
query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
```

The data in ```Ath.blast``` looks like this (tab separated):
```
ATCG00500.1	ATCG00500.1	100.00	488	0	0	1	488	1	488	0.0	 932
ATCG00510.1	ATCG00510.1	100.00	37	0	0	1	37	1	37	2e-19	73.9
ATCG00280.1	ATCG00280.1	100.00	473	0	0	1	473	1	473	0.0	 876
ATCG00890.1	ATCG01250.1	100.00	389	0	0	1	389	1	389	0.0	 660
ATCG00890.1	ATCG00890.1	100.00	389	0	0	1	389	1	389	0.0	 660
```

Here is the typical parameter setting for generating the xyz.blast file:

The example file ```Ath.pep``` contains the whole genome protein sequences (FASTA format) of Arabidopsis.

- For BLAST software
```bash
# Create a reference database
makeblastdb -in Ath.pep -dbtype prot -title Ath -parse_seqids -out Ath
# Align protein query sequences against the reference database
blastp -query query_file -db database -evalue 1e-10 -max_target_seqs 5 -outfmt 6 -out xyz.blast
# For example
blastp -query Ath.pep -db Ath -evalue 1e-10 -max_target_seqs 5 -outfmt 6 -out Ath.blast
```
- For DIAMOND software
```bash
# Create a reference database
diamond makedb --in Ath.pep -d Ath
# Align protein query sequences against the reference database
diamond blastp -d Ath -q Ath.pep -o Ath.blast -p 20 --sensitive --max-target-seqs 5 --evalue 1e-10 --quiet
```
---
**NOTE: All above input files should be stored under the same folder (the "-i" option). For more parameters please see below.**

## Running

Run the following command to get help information about **DupGen_finder**:

```bash
DupGen_finder.pl
```

This command will print a full list of options:
```
Usage: DupGen_finder.pl -i data_directory -t target_species -c outgroup_species -o output_directory
#####################
Optional:
-a 1 or 0(are segmental duplicates ancestral loci or not? default: 1, yes)
-d number_of_genes(maximum distance to call proximal, default: 10)
#####################
The following are optional MCScanX parameters:
-k match_score(cutoff score of collinear blocks for MCScanX, default: 50)
-g gap_penalty(gap penalty for MCScanX, default: -1)
-s match_size(number of genes required to call a collinear block for MCScanX, default: 5)
-e e_value(alignment significance for MCScanX, default: 1e-05)
-m max_gaps(maximum gaps allowed for MCScanX, default: 25)
-w overlap_window(maximum distance in terms of gene number, to collapse BLAST matches for MCScanX, default: 5)
```

A typical command to identify different modes of duplicated gene pairs in a given species could look like this:
```bash
DupGen_finder.pl -i data -t Ath -c Nnu -o results
```
Here, **DupGen_finder** attempts to identify the different modes of duplicated gene pairs in *A.thaliana* by using *N.nucifera* as outgroup. All required data files should be stored under this directory ```data```. The output files will be stored under this directory ```results```. For more details please see below. Ath: *A.thaliana*, Nnu: *N.nucifera*.

**Note**: We recommend that the "data_directory" or "output_directory" should be given a full path. For example, ```/home/the_path_to_your_data_directory/```

### *DupGen_finder-unique*
Moreover, to eliminate redundant duplicate genes among different modes, we provide a stricter version of **DupGen_finder** named **DupGen_finder-unique** by which each duplicate gene was assigned to a unique mode after all of the duplicated gene pairs were classified into different gene duplication types. The priority of the duplicate genes is as follows: WGD > tandem > proximal > transposed > dispersed.

```bash
DupGen_finder-unique.pl -i data -t Ath -c Nnu -o results
```

## Result Files
### 1 - Duplicated gene pairs: 
- Ath.wgd.pairs
- Ath.tandem.pairs
- Ath.proximal.pairs
- Ath.transposed.pairs
- Ath.dispersed.pairs
 
These files includes duplicated gene pairs derived from five modes of gene duplication, including WGD (**Ath.wgd.pairs**), tandem duplication (**Ath.tandem.pairs**), proximal duplication (**Ath.proximal.pairs**), transposed duplication (**Ath.transposed.pairs**), dispersed duplication (**Ath.dispersed.pairs**). The gene pairs contained in these files looks like this (tab separated):
```
Duplicate 1	Location	Duplicate 2	Location	E-value
AT1G01010.1	Ath-Chr1:3631	AT4G01550.1	Ath-Chr4:673862	5e-52
AT1G01020.1	Ath-Chr1:5928	AT4G01510.1	Ath-Chr4:642733	5e-74
AT1G01030.1	Ath-Chr1:11649	AT1G13260.1	Ath-Chr1:4542168	3e-45
AT1G01050.1	Ath-Chr1:31170	AT3G53620.1	Ath-Chr3:19880504	6e-126
AT1G01060.1	Ath-Chr1:33666	AT5G17300.1	Ath-Chr5:5690227	3e-30
```

### 2 - Ath.singletons

It includes genes that have no homologous genes within target species.
```
GeneID	Location
AT2G32600.1	Ath-Chr2:13833545
AT5G11810.1	Ath-Chr5:3808720
AT4G16610.1	Ath-Chr4:9354321
AT1G66520.1	Ath-Chr1:24816128
AT1G51110.1	Ath-Chr1:18935329
```

### 3 - Ath.stats

The number of duplicated gene pairs derived from different modes.
```
Types	NO. of gene pairs
WGD-pairs	4352
TD-pairs	2063
PD-pairs	788
TRD-pairs	4447
DSD-pairs	16130
```

### 4 - Collinearity files
- Ath.collinearity
- Ath_Nnu.collinearity

## Citation
*[Qiao X, Li Q, Yin H, Qi K, Li L, Wang R, Zhang S, Paterson AH: Gene duplication and evolution in recurring polyploidization–diploidization cycles in plants. Genome Biology 2019, 20:38.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1650-2)*
