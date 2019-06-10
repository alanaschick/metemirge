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
    input: expand("data/emirge/{sample}/{sample}_emirge.fasta", sample=SAMPLES)
        #"results/{sample}/summary.csv",
        #"results/{sample}/{sample}_best.tsv"

rule runemirge:
    input:
        r1 = config["input_dir"]+"{sample}_read1.fastq",
        r2 = config["input_dir"]+"{sample}_read2.fastq"
    output: "data/emirge/{sample}/iter."+config["num_iter_str"]
    conda: "envs/emirge_env.yaml"
    params:
        outdir = "data/emirge/{sample}/"
    shell:
        "emirge.py {params.outdir} -1 {input.r1} -2 {input.r2} "
        "-f {config[fasta_db]} -b {config[bowtie_db]} -l {config[max_read_length]} "
        "-i {config[insert_mean]} -s {config[insert_stddev]} -n {config[num_iter]} "
        "-a {config[num_threads]} --phred33"

rule rename:
    input: "data/emirge/{sample}/iter."+config["num_iter_str"]
    output: "data/emirge/{sample}/{sample}_emirge.fasta"
    conda: "envs/emirge_env.yaml"
    shell: "emirge_rename_fasta.py {input} > {output}"

rule blast:
    input: "data/emirge/{sample}/{sample}_emirge.fasta"
    output: expand("data/emirge/{{sample}}/blast/{database}/{{sample}}_raw_blast_table.tsv", database = config["database"])
    conda: "envs/blast_env.yaml"
    params:
        outdir = "data/emirge/{sample}/blast"
    shell:
        "python scripts/runBlast.py --query_file {input} --db_dir {config[blast_db_dir]} "
        "--db_names {config[database]} --output_dir {params.out_dir} --evalue {config[evalue]} "
        "--max_target_seq {config[max_target_seqs]} --output_format {config[outfmt]} "
        "--output_columns {config[out_columns]} --prefix {wildcards.sample}"

#rule parseblast:
#    input: "data/emirge/{sample}_raw_blast_table.tsv"
#    output: "data/emirge/{sample}_parsed.tsv"
#    conda: "envs/py_env.yaml"
#    shell: "python scripts/parsescript.py"

#rule selectbest:
#    input: expand("data/emirge/{sample}/{database}/{sample}_parsed.tsv", database = DATABASE)
#    output: "results/{sample}/{sample}_best.tsv"
#    conda: " envs/py_env.yaml"
#    shell: "python scripts/selecbest.py"
