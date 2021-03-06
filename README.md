# metemirge
Bioinformatics pipeline for reconstructing marker genes from metagenome sequence data, using EMIRGE.

## Overview

This pipeline is written in snakemake and designed to automate and control the submission of processes to the Synergy server at the University of Calgary in the lab of Dr. Laura Sycuro. Developed by Shalabh Thakur, Laura Sycuro, and Alana Schick.

EMIRGE is a tool for reconstructing full-length marker gene sequences from microbial community short-read sequencing data (Miller et al, 2011, Genome Biology). 

Input: filtered and cleaned fastq files. 

Output: ???

## Installation

To use this pipeline, clone this repository into your project directory using the following command:

```
git clone https://github.com/alanaschick/metemirge.git metemirge
```

Note: you need to have **conda** and **snakemake** installed in order to run this. To install conda, see the instructions [here](https://github.com/ucvm/synergy/wiki). 

To install snakemake using conda, run the following line:

```
conda install -c bioconda -c conda-forge snakemake
```

See the snakemake installation [webpage](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for further details.

## Config file

This pipeline requires a config file, written in yaml, to run. Enter any custom parameters in `config.yaml`. This is the only file that should be modified before running the pipeline. 

## Running the pipeline on Synergy

Test the pipeline by running `snakemake -np`. This command prints out the jobs that will be submitted to the Synergy compute cluster without actually submitting them.

To run the pipeline on the cluster, enter the following command from the project directory:

```
snakemake --cluster-config cluster.json --cluster 'bsub -n {cluster.n} -R {cluster.resources} -W {cluster.walllim} -We {cluster.time} -M {cluster.maxmem} -oo {cluster.output} -e {cluster.error}' --jobs 100 --use-conda
```

Note: the file `cluster.json` contains the parameters for the LSF job submission system. They are by default the same for each process, but can be modified in this file.

## Results and log files

All output files will be placed in the `results` directory and logs of each rule/process will be written to the `logs` directory.

## Pipeline summary

**1) Reconstruct marker sequences using EMIRGE.**

This step runs EMIRGE with the `emirge.py` script for reconstructing marking sequences using a reference sequence database. 

**2) Blast.**

**3) Select best blast hit.**

**4) Optional metaphlan? Compare to EMIRGE output.**
