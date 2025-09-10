import glob
import os

configfile: "config.yaml"

def get_ska_map(wildcards):
    checkpoint_output = checkpoints.create_file_list.get(**wildcards).output[0]
    SAMPLES = [l.rstrip().split("\t")[0] for l in open(f"{config['output_dir']}/reads_list.tsv")]      
    return expand("{out_dir}/ska_map/{sample}.vcf", out_dir=config['output_dir'], sample=SAMPLES)
    
# Define the final output
rule all:
   input:
      f"{config['output_dir']}/reads_list.tsv",
      f"{config['output_dir']}/shovill/contigs.fa",
      get_ska_map

checkpoint create_file_list:
  input:
    indir=f"{config['input_dir']}",
  output:
    output_file=f"{config['output_dir']}/reads_list.tsv",
  params:
    reference_reads_R1=f"{config['reference_reads_R1']}",
    reference_reads_R2=f"{config['reference_reads_R2']}"
  script:
    "scripts/create_file_list.py"

rule shovill:
  input:
    r1=f"{config['reference_reads_R1']}",
    r2=f"{config['reference_reads_R2']}"
  output:
    outdir=directory(f"{config['output_dir']}/shovill"),
    contigs=f"{config['output_dir']}/shovill/contigs.fa"
  threads: 40
  log:
    f"{config['output_dir']}/logs/shovill.txt"
  params:
    tmp_dir=f"{config['output_dir']}/tmp"
  shell:
    """
    shovill --assembler spades --outdir {output.outdir} --force --R1 {input.r1} --R2 {input.r2} --cpus {threads} --tmpdir {params.tmp_dir} &> {log}
    echo "Assembly: {output.contigs}"
    """

checkpoint ska_build:
  input:
    infile=f"{config['output_dir']}/reads_list.tsv"
  output:
    outpref=f"{config['output_dir']}/ska_build/ska_"
  params:
    ksize=f"{config['ksize']}"
  threads: 40
  log:
    f"{config['output_dir']}/logs/ska_build.txt"
  shell:
    """
    echo "Building from {input.infile}"
    ska build -v -o {output.outpref} -k {params.ksize} --threads {threads} -f {input.infile} &> {log}
    """

def get_ska_build(wildcards):
  checkpoint_output = checkpoints.ska_build.get(**wildcards).output[0]
  # For a given sample, return just that one file
  return f"{config['output_dir']}/ska_build/ska_{wildcards.sample}.skf"

checkpoint ska_map:
  input:
    infile=get_ska_build,
    contigs=f"{config['output_dir']}/ref_contigs.fa"
  output:
    outfile=f"{config['output_dir']}/ska_map/{{sample}}.vcf"
  threads: 1
  log:
    f"{config['output_dir']}/logs/ska_map.txt"
  shell:
    """
    ska map -f vcf -v -o {output.outfile} --threads {threads} {input.contigs} {input.infile} &> {log}
    """


