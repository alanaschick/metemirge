# ************************************
# * Snakefile for metemirge pipeline *
# ************************************

# **** Variables ****

configfile: "config.yaml"

# Specify the list of files to run the pipeline on
import pandas as pd
SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# **** Rules ****

rule all:
    input:
        "results/{sample}/summary.csv",
        "results/{sample}/{sample}_best.tsv"

rule runemirge:
    input:
        r1 = "{params.dir}{sample}_read1.fastq",
        r2 = "{params.dir}{sample}_read2.fastq"
    output: "data/emirge/{sample}/iter.{config[max_num_iter]}/iter.{config[max_num_iter]}.cons.fasta}"
    params:
        dir = {config[emirge_dir]}
    conda: "envs/emirge_env.yaml"
    shell:
        "emirge.py {config[emirge_dir]} -1 {input.r1} -2 {input.r2} "
        "-f {config[fasta_db]} -b {config[bowtie_db]} -l {config[max_read_length]} "
        "-i {config[insert_mean]} -s {config[insert_stddev]} -n {config[num_iter]} "
        "-a {config[num_threads]}"

#rule rename:
#    input: "data/emirge/{sample}/iter.{config[max_num_iter]}/iter.{config[max_num_iter]}.cons.fasta}"
#    output: "data/emirge/{sample}/{sample}_renamed.fasta"
#    conda: "envs/emirge_env.yaml"
#    shell: "python scripts/emirge_rename_fasta.py <options>"

#rule blast:
#    input: "data/emirge/{sample}/{sample}_renamed.fasta"
#    output: "data/emirge/{sample}/{sample}_raw_blast_table.tsv"
#    params:
#        db = "{config[blast_db]}"
#    conda: "envs/blast_env.yaml"
#    shell: "blast <options>"

#rule parseblast:
#    input: "data/emirge/{sample}/{sample}_raw_blast_table.tsv"
#    output: "data/emirge/{sample}/{sample}_parsed.tsv"
#    conda:
#    shell: "python scripts/parsescript.py"

#rule selectbest:
#    input: expand("data/emirge/{sample}/{database}/{sample}_parsed.tsv", database = DATABASE)
#    output: "results/{sample}/{sample}_best.tsv"
#    conda:
#    shell: "python scripts/selecbest.py"
