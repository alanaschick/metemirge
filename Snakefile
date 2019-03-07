# **********************************
# * Snakefile for metemirge pipeline *
# **********************************

# **** Variables ****

configfile: "config.yaml"

# Not sure if this is needed for metemirge, it specifies a list of files to run the pipeline on
import pandas as pd
SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# **** Rules ****

rule all:
    input:
        "results/result1file",
        "results/result2file"

rule reconstruct:
    input:
        r1 = "data/filtdata/{sample}_1.fastq.gz",
        r2 = "data/filtdata/{sample}_2.fastq.gz"
    output: "data/emirge/{sample}"
    conda: "envs/reconstruct_env.yaml"
    shell: "python scripts/runEmirgeCluster.py <options>"

rule blast:
    input: "data/emirge/{sample}"
    output: "results/result1file"
    conda: "envs/blast_env.yaml"
    shell: "python scripts/runParseBlast.py <options>"

rule selectbest:
    input: "results/result1file"
    output: "results/result2file"
    conda: "envs/selectbest_env.yaml"
    shell: "python scripts/runEmirgeSelectBlastHit.py <options>"
