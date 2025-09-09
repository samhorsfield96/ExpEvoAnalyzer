import glob
import os

configfile: "config.yaml"

def get_tokenized_genomes(wildcards):
    checkpoint_output = checkpoints.list_pyrodigal_full.get(**wildcards).output[0]
    return expand("{out_dir}/tokenised_genomes/tokenized_genomes_batch_{batch_ID}.txt", out_dir=config['output_dir'], batch_ID=range(config['mmseqs2_num_batches']))
    
# Define the final output
rule all:
   input:
        f"{config['output_dir']}/merged_clusters/all_clusters.pkl",
        f"{config['output_dir']}/merged_clusters/all_clusters.tsv",
        f"{config['output_dir']}/merged_clusters/reps.pkl",
        get_tokenized_genomes

rule shovill:
  input:
    r1=f"{config['ref_reads']}_R1.fq.gz",
    r2=f"{config['ref_reads']}_R2.fq.gz",
  output:
    raw_assembly=f"{config['output_dir']}/ref_assembly.fa",
    contigs=f"{config['output_dir']}/ref_contigs.fa"
  params:
    extra=""
  log:
    "logs/shovill/{sample}.{assembler}.log"
  threads: 1
  wrapper:
    "v6.0.1/bio/shovill"

rule ska_build:
    input:
        r1=f"{config['read_dir']}/{{sample}}_R1_{id}.fastq.gz",
        r2=f"{config['read_dir']}/{{sample}}_R2_{id}.fastq.gz"
    output:
    params:
        ksize=f"{config['ksize']}"
    threads: 1
    shell:
    """
    ska build
    """

rule ska_map:
    input:
        contigs=f"{config['output_dir']}/ref_contigs.fa",
        r1=f"{config['read_dir']}/{{sample}}_R1_{id}.fastq.gz",
        r2=f"{config['read_dir']}/{{sample}}_R2_{id}.fastq.gz"
    output:
    params:
        ksize=f"{config['ksize']}"
    threads: 1
    shell:
    """
    ska build
    ska map -f vcf --threads {threads} {input.contigs} 
    """

