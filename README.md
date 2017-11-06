# IT Miner

A Pipeline to predict the Intrinsic Terminators of your bacterial genome with a low false positive rate.

## Table of contents
- [**Output**](#output)
- [**Inputs**](#input)
- [**Installation**](#installation)
	- [Clone](#clone)
	- [Dependencies](#dependencies)
	- [Requirements](#requirements)
- [**Usage**](#usage)
- [**Command Line Example**](#command-line-example)


## Output
IT Miner produces a gff file containing the coordinates of all predicted intrinsic terminators in the genome of interest.

## Inputs
IT Miner requires two input files:
 - a genome sequence in fasta.
 - an annotation of the genome in either gff, gff3 or gtf format.

## Installation

### Clone
    IT_MINER_DIR=<your_path>
    cd "$IT_MINE_DIR"
    git clone https://github.com/charles-bernard/IT_Miner.git

### Dependencies
IT Miner depends on RNIE, which itself depends on Infernal 1.02

Installing Infernal 1.02 (Bash)

    # Define your Installation Directory
    INFERNAL_DIR=<your_path>
    cd "$INFERNAL_DIR"

    # Download Infernal 1.02
    wget http://eddylab.org/software/infernal/infernal-1.0.2.tar.gz
    uncompress infernal-1.0.2.tar.gz
    tar xf infernal-1.0.2.tar

    # Configure and Install Infernal
    cd infernal-1.0.2
    ./configure
    make
    make check
    make install

Installing RNIE

NB: This step is optional since the RNIE Directory is by default
in the repository of IT Miner.

By default RNIE_DIR="$IT_MINER_DIR"/Subscripts/RNIE

    git clone https://github.com/ppgardne/RNIE
    cp -r RNIE "$IT_MINER_DIR"/Subscripts 
    RNIE_DIR="$IT_MINER_DIR"/Subscripts/RNIE

### Requirements
RNIE needs some binaries from Infernal to run. One has
to create symlinks of these binaries in the path of RNIE.

    ln -s "$INFERNAL_DIR"/src/cmsearch "$RNIE_DIR"
    ln -s "$INFERNAL_DIR"/src/cmalign "$RNIE_DIR"
    ln -s "$INFERNAL_DIR"/easel/miniapps/esl-sfetch "$RNIE_DIR"

## Usage

To run IT_Miner, you can open terminal and type:

    bash IT_Miner.sh

with the following arguments:
* _**-o/--output-dir**_ Output Directory
* _**-g/--genome**_ Path to genome sequence in fasta
* _**-a/--annotation**_ Path to genome annotation in gff/gff3/gtf
* _**-c/--cutoff**_ specifies the type of cutoff distance from upstream gene
to apply for filering. Can be either 'conservative', 'inclusive', or 'average' btw the two.
* _**-l/--log**_ weither you want to specify a path for the log file

## Command Line Example
    # Escherichia coli 
    bash "<my_IT_Miner_path>"/IT_Miner.sh \
     --output-dir "/home/charles/clones/IT_Miner/Out/Escherichia_coli/U00096.2" \
     --genome "/home/charles/Documents/Genomes/Escherichia_coli_K12/release2_U00096.2.fasta" \
     --annotation "/home/charles/Documents/Genomes/Escherichia_coli_K12/release2_U00096.2.gff3" \
     --cutoff conservative
