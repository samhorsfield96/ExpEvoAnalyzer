import glob
import os
import pandas as pd

configfile: "config.yaml"

def get_snpeff_run(wildcards):
    # Wait for create_file_list to finish
    ckpt1 = checkpoints.create_file_list.get(**wildcards)

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
    reference_reads_R2=f"{config['reference_reads_R2']}",
    reference_assembly=f"{config['reference_assembly']}"
  script:
    "scripts/create_file_list.py"

def get_create_file_list(wildcards):
    ckpt = checkpoints.create_file_list.get(**wildcards)
    outdir = ckpt.output.output_dir

    # discover all .txt files in the directory
    return [os.path.join(outdir, f) for f in os.listdir(outdir) if f.endswith(".txt")]

if config['reference_reads_R1'] != "" and config['reference_reads_R2'] != "" and config['reference_assembly'] == "":
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
      """

# rename reference assembly
elif config['reference_assembly'] != "":
  def get_assembly(wildcards):
    base = f"{config['reference_assembly']}"
    if os.path.exists(base):
      return base
    else:
        raise ValueError(f"No assembly file found for sample {config['reference_assembly']}")

  rule copy_ref_assembly:
    input:
      assembly=get_assembly
    output:
      contigs=f"{config['output_dir']}/shovill/contigs.fa"
    threads: 1
    params:
      outdir=f"{config['output_dir']}/shovill"
    log:
      f"{config['output_dir']}/logs/copy_reference_assembly.txt"
    shell:
      r"""
      mkdir -p {params.outdir}
      if [[ "{input.assembly}" == *.gz ]]; then
          gunzip -c "{input.assembly}" > "{output.contigs}"
      else
          cp "{input.assembly}" "{output.contigs}"
      fi
      """

# only run bakta if not already annotated
if config['reference_assembly_gbk'] == "":
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

elif config['reference_assembly_gbk'] != "":
  def get_assembly_gbk(wildcards):
    base = f"{config['reference_assembly_gbk']}"
    if os.path.exists(base):
      return base
    else:
        raise ValueError(f"No assembly file found for sample {config['reference_assembly']}")

  rule copy_ref_assembly_gbk:
    input:
      assembly=get_assembly_gbk
    output:
      gbk=f"{config['output_dir']}/bakta/genes.gbk",
    threads: 1
    params:
      outdir=f"{config['output_dir']}/bakta"
    log:
      f"{config['output_dir']}/logs/copy_reference_assembly_gbk.txt"
    shell:
      r"""
      mkdir -p {params.outdir}
      cp "{input.assembly}" "{output.gbk}"
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
    snpEff build -genbank -nodownload -dataDir {params.indir} -c {params.config} bakta &> {log}
    """

if config['alignment_method'] in ["ska"]:
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
      mkdir -p {params.outdir}
      ska build -v -o {params.outpref} -k {params.ksize} --threads {threads} -f {input.infile} &> {log}
      """

  rule ska_map:
    input:
      infile=f"{config['output_dir']}/ska_build/{{sample}}.skf",
      contigs=f"{config['output_dir']}/shovill/contigs.fa"
    output:
      outfile=f"{config['output_dir']}/vcfs/{{sample}}.vcf"
    params:
      outdir=f"{config['output_dir']}/vcfs",
    threads: 1
    log:
      f"{config['output_dir']}/logs/ska_map_{{sample}}.txt"
    shell:
      """
      mkdir -p {params.outdir}
      ska map -f vcf -v -o {output.outfile} --threads {threads} {input.contigs} {input.infile} &> {log}
      """

elif config['alignment_method'] in ["bwa"]:
  # indexing
  rule index_reference:
    input:
        f"{config['output_dir']}/shovill/contigs.fa"
    output:
        idx=f"{config['output_dir']}/shovill/contigs.fa.bwt"
    log:
      f"{config['output_dir']}/logs/bwa_index.txt"
    shell:
        """
        bwa index {input} &> {log}
        """

  # Run BWA MEM per sample
  rule bwa_map:
    input:
      contigs=f"{config['output_dir']}/shovill/contigs.fa",
      sheet=f"{config['output_dir']}/reads_list/{{sample}}.txt",
      idx=f"{config['output_dir']}/shovill/contigs.fa.bwt"
    output:
      outfile=f"{config['output_dir']}/bwa_mem/{{sample}}.bam"
    threads: 8
    params:
      outdir=f"{config['output_dir']}/bwa_mem"
    log:
      f"{config['output_dir']}/logs/bwa_mem_{{sample}}.txt"
    run:
      os.makedirs(params.outdir, exist_ok=True)
      df = pd.read_csv(input.sheet, sep="\t", header=None, names=["sample", "r1", "r2"])
      r1, r2 = df.loc[0, ["r1", "r2"]]
      shell(f"bwa mem -t {threads} {input.contigs} {r1} {r2} 2> {log} | samtools view -Sb - > {output.outfile}")

  # SORT BAM
  rule sort_bam:
    input:
        f"{config['output_dir']}/bwa_mem/{{sample}}.bam"
    output:
        f"{config['output_dir']}/bwa_mem/{{sample}}.sorted.bam"
    log:
      f"{config['output_dir']}/logs/samtools_sort_{{sample}}.txt"
    shell:
        "samtools sort -o {output} {input} &> {log}"

  # Mark duplicates
  rule mark_duplicates:
      input:
          f"{config['output_dir']}/bwa_mem/{{sample}}.sorted.bam"
      output:
          f"{config['output_dir']}/bwa_mem/{{sample}}.sorted.dedup.bam"
      threads: 8
      params:
        tmp_pref=f"{config['output_dir']}/bwa_mem/{{sample}}.sorted.dedup.bam.tmp"
      shell:
        """
        samtools markdup -r -@ {threads} -T {params.tmp_pref} {input} {output}
        """

  # INDEX BAM
  rule index_bam:
    input:
        f"{config['output_dir']}/bwa_mem/{{sample}}.sorted.dedup.bam"
    output:
        f"{config['output_dir']}/bwa_mem/{{sample}}.sorted.dedup.bam.bai"
    log:
      f"{config['output_dir']}/logs/samtools_index_{{sample}}.txt"
    shell:
        "samtools index -o {output} {input} &> {log}"

  # CALL VARIANTS
  rule call_variants:
    input:
        contigs=f"{config['output_dir']}/shovill/contigs.fa",
        bam=f"{config['output_dir']}/bwa_mem/{{sample}}.sorted.dedup.bam",
        bai=f"{config['output_dir']}/bwa_mem/{{sample}}.sorted.dedup.bam.bai"
    output:
        outfile=f"{config['output_dir']}/vcfs/{{sample}}.vcf"
    params:
      outdir=f"{config['output_dir']}/vcfs",
    log:
      f"{config['output_dir']}/logs/bcftools_{{sample}}.txt"
    shell:
        """
        mkdir -p {params.outdir}
        bcftools mpileup -Ou -f {input.contigs} {input.bam} | \
        bcftools call -mv -Ov -o {output.outfile} &> {log}
        """

checkpoint snpeff_run:
  input:
    sequence=f"{config['output_dir']}/bakta/sequence.bin",
    snpEffectPredictor=f"{config['output_dir']}/bakta/snpEffectPredictor.bin",
    vcf=f"{config['output_dir']}/vcfs/{{sample}}.vcf"
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
    """

rule read_ann_vcf:
  input:
    ann_vcf=get_snpeff_run
  output:
    outfile=f"{config['output_dir']}/annotated_variants.csv"
  params:
    indir=f"{config['output_dir']}/snpeff"
  script:
    "scripts/read_ann_vcf.py"
