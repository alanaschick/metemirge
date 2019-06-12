# ************************************
# * Snakefile for metemirge pipeline *
# ************************************

# **** Variables ****
configfile: "config.yaml"
import os
os.environ["BLASTDB"]=config["path_to_taxdb"]
os.environ["PATH"]+=os.pathsep+"/export/common/bin"

# Specify the list of files to run the pipeline on
import pandas as pd
SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# **** Rules ****

rule all:
    input: expand("data/blast/{database}/{sample}_raw_blast_table.tsv", sample=SAMPLES, database=config["list_blastdb_names"])
        #"results/{sample}/summary.csv",
        #"results/{sample}/{sample}_best.tsv"

rule runemirge:
    input:
        r1 = config["input_dir"]+"{sample}_read1.fastq",
        r2 = config["input_dir"]+"{sample}_read2.fastq"
    output: directory("data/emirge/{sample}/iter."+config["num_iter_str"])
    conda: "metemirge_files/envs/emirge_env.yaml"
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
    conda: "metemirge_files/envs/emirge_env.yaml"
    shell: "emirge_rename_fasta.py {input} > {output}"

rule blast:
    input: "data/emirge/{sample}/{sample}_emirge.fasta"
    output: "data/blast/{database}/{sample}_raw_blast_table.tsv"
    conda: "metemirge_files/envs/blast_env.yaml"
    params:
        db = config["path_to_blastdbs"]+"{database}"+"/"+"{database}"
    shell:
        "blastn -query {input} -db {params.db} -out {output} "
        "-evalue {config[evalue]} -max_target_seqs {config[max_target_seqs]} "
        "-outfmt '6 qseqid qacc qlen sseqid sacc slen stitle qstart qend sstart send length evalue bitscore pident qcovs nident mismatch positive gaps qframe sframe staxids sskingdoms'"

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
