import glob
import os

configfile: "config.yaml"

def get_ska_map(wildcards):
    ckpt1 = checkpoints.create_file_list.get(**wildcards).output[0] 
    ckpt2 = checkpoints.ska_build.get(**wildcards)
    reads_list = f"{ckpt1.output.outdir}/../reads_list.tsv"
    with open(reads_list) as f:
        samples = [line.rstrip().split("\t")[0] for line in f]
    return expand(
        "{outdir}/ska_map/{sample}.vcf",
        outdir=config["output_dir"],
        sample=samples
    )

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
    rmdir {params.tmp_dir}
    echo "Assembly: {output.contigs}"
    """

checkpoint ska_build:
  input:
    infile=f"{config['output_dir']}/reads_list.tsv"
  output:
    outdir=directory(f"{config['output_dir']}/ska_build")
  params:
    ksize=f"{config['ksize']}"
  threads: 40
  log:
    f"{config['output_dir']}/logs/ska_build.txt"
  shell:
    """
    echo "Building from {input.infile}"
    ska build -v -o {output.outdir}/ska_ -k {params.ksize} --threads {threads} -f {input.infile} &> {log}
    """

def get_ska_build(wildcards):
    ckpt = checkpoints.ska_build.get(**wildcards)
    outdir = ckpt.output.outdir

    # discover all .skf files in the directory
    return [os.path.join(outdir, f) for f in os.listdir(outdir) if f.endswith(".skf")]

rule ska_map:
  input:
    infile=get_ska_build,
    contigs=f"{config['output_dir']}/shovill/contigs.fa"
  output:
    outfile=f"{config['output_dir']}/ska_map/{{sample}}.vcf"
  threads: 1
  log:
    f"{config['output_dir']}/logs/ska_map_{{sample}}.txt"
  shell:
    """
    ska map -f vcf -v -o {output.outfile} --threads {threads} {input.contigs} {input.infile} &> {log}
    """
