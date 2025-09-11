import glob
import os

configfile: "config.yaml"

def get_snpeff_run(wildcards):
    # Wait for create_file_list to finish
    ckpt1 = checkpoints.create_file_list.get(**wildcards)
    ckpt2 = checkpoints.snpeff_run.get(**wildcards)

    # Grab the reads list it produced
    reads_list = ckpt1.output.outfile

    # Collect sample names from the file
    with open(reads_list) as f:
        samples = [line.rstrip().split("\t")[0] for line in f]

    # Return all snpEff-annotated VCFs expected for these samples
    return expand(
        f"{config['output_dir']}/snpeff/{{sample}}.ann.vcf",
        sample=samples
    )

# Define the final output
rule all:
   input:
      f"{config['output_dir']}/reads_list/reads_list_all.tsv",
      f"{config['output_dir']}/shovill/contigs.fa",
      f"{config['output_dir']}/bakta/genes.gbk",
      f"{config['output_dir']}/bakta/snpEffectPredictor.bin",
      get_snpeff_run,
      f"{config['output_dir']}/annotated_variants.csv"

checkpoint create_file_list:
  input:
    indir=f"{config['input_dir']}",
  output:
    output_dir=directory(f"{config['output_dir']}/reads_list"),
    outfile=f"{config['output_dir']}/reads_list/reads_list_all.tsv"
  params:
    reference_reads_R1=f"{config['reference_reads_R1']}",
    reference_reads_R2=f"{config['reference_reads_R2']}"
  script:
    "scripts/create_file_list.py"

def get_create_file_list(wildcards):
    ckpt = checkpoints.create_file_list.get(**wildcards)
    outdir = ckpt.output.output_dir

    # discover all .txt files in the directory
    return [os.path.join(outdir, f) for f in os.listdir(outdir) if f.endswith(".txt")]

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
    tmp_dir=f"{config['output_dir']}/tmp",
    shovill_params=f"{config['shovill_params']}"
  shell:
    """
    shovill {params.shovill_params} --assembler spades --outdir {output.outdir} --force --R1 {input.r1} --R2 {input.r2} --cpus {threads} --tmpdir {params.tmp_dir} &> {log}
    rmdir {params.tmp_dir}
    echo "Assembly: {output.contigs}"
    """

rule bakta:
    input:
      contigs=f"{config['output_dir']}/shovill/contigs.fa"
    output:
      ann_dir = directory(f"{config['output_dir']}/bakta"),
      gbk_file = f"{config['output_dir']}/bakta/genes.gbk",
      gbff_file = f"{config['output_dir']}/bakta/genes.gbff",
    threads: 40
    params:
        DB = directory(f"{config['bakta_db']}"),
        tmp_dir=f"{config['output_dir']}/tmp"
    log:
        f"{config['output_dir']}/logs/bakta.log"
    shell:
        """
        mkdir -p {params.tmp_dir}
        bakta {input.contigs} --tmp-dir {params.tmp_dir} --force --keep-contig-headers --db {params.DB} --prefix genes --translation-table 11 --skip-plot --threads {threads} --output {output.ann_dir} &> {log}
        rmdir {params.tmp_dir}
        cp {output.gbff_file} {output.gbk_file}
        """

rule snpeff_build:
  input:
    gbk_file = f"{config['output_dir']}/bakta/genes.gbk"
  output:
    sequence=f"{config['output_dir']}/bakta/sequence.bin",
    snpEffectPredictor=f"{config['output_dir']}/bakta/snpEffectPredictor.bin"
  params:
    config=f"{config['snpEFF_config']}",
    indir=f"{config['output_dir']}",
  threads: 1
  log:
    f"{config['output_dir']}/logs/snpeff_build.txt"
  shell:
    """
    echo "Reading from {input.gbk_file}"
    snpEff build -genbank -nodownload -dataDir {params.indir} -c {params.config} bakta &> {log}
    echo "Built {output.snpEffectPredictor} and {output.sequence}"
    """

rule ska_build:
  input:
    infile=f"{config['output_dir']}/reads_list/{{sample}}.txt",
  output:
    outfile=f"{config['output_dir']}/ska_build/{{sample}}.skf"
  params:
    outdir=f"{config['output_dir']}/ska_build",
    outpref=f"{config['output_dir']}/ska_build/{{sample}}",
    ksize=f"{config['ksize']}"
  threads: 8
  log:
    f"{config['output_dir']}/logs/ska_build_{{sample}}.txt"
  shell:
    """
    echo "Building from {input.infile}"
    mkdir -p {params.outdir}
    ska build -v -o {params.outpref} -k {params.ksize} --threads {threads} -f {input.infile} &> {log}
    """

rule ska_map:
  input:
    infile=f"{config['output_dir']}/ska_build/{{sample}}.skf",
    contigs=f"{config['output_dir']}/shovill/contigs.fa"
  output:
    outfile=f"{config['output_dir']}/ska_map/{{sample}}.vcf"
  params:
    outdir=f"{config['output_dir']}/ska_map",
  threads: 1
  log:
    f"{config['output_dir']}/logs/ska_map_{{sample}}.txt"
  shell:
    """
    mkdir -p {params.outdir}
    ska map -f vcf -v -o {output.outfile} --threads {threads} {input.contigs} {input.infile} &> {log}
    """

checkpoint snpeff_run:
  input:
    sequence=f"{config['output_dir']}/bakta/sequence.bin",
    snpEffectPredictor=f"{config['output_dir']}/bakta/snpEffectPredictor.bin",
    vcf=f"{config['output_dir']}/ska_map/{{sample}}.vcf"
  output:
    ann_vcf=f"{config['output_dir']}/snpeff/{{sample}}.ann.vcf"
  params:
    config=f"{config['snpEFF_config']}",
    indir=f"{config['output_dir']}",
    outdir=f"{config['output_dir']}/snpeff"
  threads: 1
  log:
    f"{config['output_dir']}/logs/snpeff_eff_{{sample}}.txt"
  shell:
    """
    mkdir -p {params.outdir}
    snpEff eff -i vcf -o vcf -noStats -lof -noLog -v -nodownload -dataDir {params.indir} -c {params.config} bakta {input.vcf} > {output.ann_vcf} 2> {log}
    echo "Built {input.snpEffectPredictor} and {input.sequence}"
    """

rule read_ann_vcf:
  input:
    ann_vcf=get_snpeff_run,
    indir=f"{config['output_dir']}/snpeff/"
  output:
    outfile=f"{config['output_dir']}/annotated_variants.csv"
  script:
    "scripts/read_ann_vcf.py"
