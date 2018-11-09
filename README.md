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
   - a. "[target_species].gff", a gene position file for the target species, following a tab-delimited format: "sp&chr_NO      gene    starting_position       ending_position"
   - b. "[target_species].blast", a blastp output file (m8 format) for the target species (self-genome comparison).
2. For each outgroup genome, please prepare two input files:
   - a. "[target_species]_[outgroup_species].gff", a gene position file for the target_species and outgroup_species, following a tab-delimited format:"sp&chr_NO      gene    starting_position       ending_position"
   - b. "[target_species]_[outgroup_species].blast", a blastp output file (m8 format) between the target and outgroup species (cross-genome comparison).
3. For example, assuming that you are going to classify gene duplication modes in Arabidopsis thaliana (ID: Ath), using Nelumbo nucifera (ID: Nnu)as outgroups, you need to prepare 6 input files: "Ath.gff","Ath.blast", "Ath_Nnu.gff", "Ath_Nnu.blast"

**NOTE**: All input files should be stored under ONE folder(the "data_directory" parameter)

## Running

You can simply run the following command to get help information about **DupGen_finder**:

```bash
$ perl DupGen_finder.pl
```

Help information

```
  Usage: perl DupGen_finder.pl -i data_directory -t target_species -c outgroup_species(comma_delimited) -o output_directory
  #####################
  Optional:
  -x number_of_different_epoches (if specified, outgroup species must be provided in the order of divergence from the target species(most recent first), default: 1, only consider the transposed duplications that occurred after the divergence between target species and all outgroups )
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
