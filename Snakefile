import glob
import os

configfile: "config.yaml"

def get_ska_map(wildcards):
   checkpoint_output = checkpoints.ska_map.get(**wildcards).output[0]
   # read reads_list and get sample ids
   sample_ids = []
   with open(f"{config['output_dir']}/reads_list.tsv", "r") as f:
    for line in f:
      split_line = f.rstrip().split("\t")
      sample_ids.append(split_line[0])
   return expand("{out_dir}/{sample}.vcf", out_dir=f"{config['output_dir']}, batch_ID=sample_ids)
    
# Define the final output
rule all:
   input:
        get_ska_map

rule create_file_list:
  input:
    indir=f"{config['input_dir']}",
  output:
    output_file=f"{config['output_dir']}/reads_list.tsv",
  params:
    reference_reads_R1=f"{config['reference_reads_R1']}",
    reference_reads_R2=f"{config['reference_reads_R2']}",
  script:
    "scripts/create_file_list.py"

rule shovill:
  input:
    r1=f"{config['reference_reads_R1']}",
    r2=f"{config['reference_reads_R1']}",
  output:
    raw_assembly=f"{config['output_dir']}/ref_assembly.fa",
    contigs=f"{config['output_dir']}/ref_contigs.fa"
  params:
    extra=""
  log:
    "logs/shovill/reference.log"
  threads: 1
  wrapper:
    "v7.2.0/bio/shovill"

checkpoint ska_build:
    input:
      infile=f"{config['output_dir']}/reads_list.tsv"
    output:
      outpref=f"{config['output_dir']}/ska_"
    params:
        ksize=f"{config['ksize']}"
    threads: 40
    shell:
    """
    echo "Building from {input.infile}"
    ska build -o {output.outpref} -k {params.ksize} --threads {threads} -f {input.infile} 
    """

def get_ska_build(wildcards):
   checkpoint_output = checkpoints.ska_build.get(**wildcards).output[0]
   # read reads_list and get sample ids
   sample_ids = []
   with open(f"{config['output_dir']}/reads_list.tsv", "r") as f:
    for line in f:
      split_line = f.rstrip().split("\t")
      sample_ids.append(split_line[0])
   return expand("{out_dir}/ska_{sample}.skf", out_dir=f"{config['output_dir']}, batch_ID=sample_ids)

checkpoint ska_map:
    input:
      checkpoint=get_ska_build
      contigs=f"{config['output_dir']}/ref_contigs.fa",
      infile=f"{config['output_dir']}/ska_{sample}.skf"
    output:
      outfile=f"{config['output_dir']}/{sample}.vcf"
    params:
        ksize=f"{config['ksize']}"
    threads: 1
    shell:
    """
    echo "Mapping {input.checkpoint}"
    ska map -f vcf --threads {threads} {input.contigs} {input.infile} 
    """

