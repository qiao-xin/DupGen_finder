# DupGen_finder

The DupGen_finder was developed to identify different modes of duplicated gene pairs. [MCScanX](http://chibba.pgml.uga.edu/mcscan2/) algorithm was incorporated in this pipeline.

| | |
| --- | --- |
| Authors | Xin Qiao ([qiaoxin](https://github.com/qiao-xin)) |
| | Andrew Paterson ([PGML](http://www.plantgenome.uga.edu)) |
| Email   | <qiaoxinqx2011@126.com> |

## Dependencies

- [Perl](https://www.perl.org)

## Installation

```bash
git clone https://github.com/qiao-xin/DupGen_finder.git
cd DupGen_finder
make
```

## Prepeations

Pre-computed BLAST results and gene location information (GFF format) are required for running DupGen_finder successfully.

1. For the target genome in which gene duplicaiton modes will be classified, please prepare two input files:
a. "[target_species].gff", a gene position file for the target species, following a tab-delimited format: "sp&chr_NO      gene    starting_position       ending_position". For example, "Ath.gff".

```
Ath-Chr1	AT1G01010.1	3631	5899
Ath-Chr1	AT1G01020.1	5928	8737
Ath-Chr1	AT1G01030.1	11649	13714
Ath-Chr1	AT1G01040.2	23416	31120
Ath-Chr1	AT1G01050.1	31170	33153
```

b. "[target_species].blast", a blastp output file (m8 format) for the target species (self-genome comparison). For example, "Ath.blast".

```
ATCG00500.1	ATCG00500.1	100.00	488	0	0	1	488	1	488	0.0	 932
ATCG00510.1	ATCG00510.1	100.00	37	0	0	1	37	1	37	2e-19	73.9
ATCG00280.1	ATCG00280.1	100.00	473	0	0	1	473	1	473	0.0	 876
ATCG00890.1	ATCG01250.1	100.00	389	0	0	1	389	1	389	0.0	 660
ATCG00890.1	ATCG00890.1	100.00	389	0	0	1	389	1	389	0.0	 660
```

2. For the outgroup genome, please prepare two input files:
   - a. "[target_species]_[outgroup_species].gff", a gene position file for the target_species and outgroup_species, following a tab-delimited format:"sp&chr_NO      gene    starting_position       ending_position"
   - b. "[target_species]_[outgroup_species].blast", a blastp output file (m8 format) between the target and outgroup species (cross-genome comparison).
3. For example, assuming that you are going to classify gene duplication modes in Arabidopsis thaliana (ID: Ath), using Nelumbo nucifera (ID: Nnu)as outgroups, you need to prepare 6 input files: "Ath.gff","Ath.blast", "Ath_Nnu.gff", "Ath_Nnu.blast"

**NOTE**: All input files should be stored under ONE folder (the "data_directory" parameter)

## Running

You can simply run the following command to get help information about **DupGen_finder**:

```bash
$ perl DupGen_finder.pl
```

Help information:

```
  Usage: perl DupGen_finder.pl -i data_directory -t target_species -c outgroup_species -o output_directory
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
