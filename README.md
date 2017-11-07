# IT Miner

A Pipeline to predict the **I**ntrinsic **T**erminators  sites of bacterial genomes with low false positive rate. 

_IT_Miner_ relies on the **[RNIE](https://github.com/ppgardne/RNIE/tree/master/paper)** algorithm to produce a list of predicted intrinsic terminators and filter out this primary output based on several genomic criteria to retain only the highly probable terminators.

## Output
IT Miner produces a gff file containing the coordinates of all predicted intrinsic terminators in the genome of interest.

## Inputs
IT Miner requires two input files:
 - a genome sequence in fasta.
 - an annotation of the genome in either gff, gff3 or gtf format.

## Installation and Usage
If you want to install and use this program, please take a look at the [wiki](https://github.com/charles-bernard/IT_Miner/wiki) of the repository.
