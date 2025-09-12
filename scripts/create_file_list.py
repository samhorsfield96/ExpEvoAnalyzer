import os
import re

def detect_paired_fastq(input_dir, output_dir, reference_reads_R1=None, reference_reads_R2=None, reference_assembly=None):
    """
    Detect paired FASTQ files in a directory and write them to a tab-delimited file.
    
    Args:
        input_dir (str): Path to directory containing FASTQ files.
        output_dir (str): Output directory.
        reference_reads_R1 (str): Reference reads R2 to ignore.
        reference_reads_R2 (str): Reference reads R2 to ignore.
    """
    
    # Regex pattern for R1/R2 fastq files
    fastq_pattern = re.compile(r"(.+?)(?:_R|_)([12])(?:[^/]*)\.(?:fastq|fq)(?:\.gz)?$")
    fasta_pattern = re.compile(r"(.+?)\.(?:fa|fna|fasta)(?:\.gz)?$")

    pairs = {}
    assemblies = {}

    for fname in os.listdir(input_dir):
        reference = False
        full_path = os.path.join(input_dir, fname)
        if reference_reads_R1 and reference_reads_R1 in full_path:
            reference = True
        if reference_reads_R2 and reference_reads_R2 in full_path:
            reference = True
        if reference_assembly and reference_assembly in full_path:
            reference = True

        fastq_match = fastq_pattern.match(fname)
        fasta_match = fasta_pattern.match(fname)

        # Case 1: FASTQ paired reads
        if fastq_match:
            sample_id, read_direction = fastq_match.groups()
            if sample_id not in pairs:
                pairs[sample_id] = {}
            if read_direction == "1":
                pairs[sample_id]["R1"] = (full_path, reference)
            elif read_direction == "2":
                pairs[sample_id]["R2"] = (full_path, reference)

        # Case 2: FASTA assemblies
        elif fasta_match:
            sample_id = fasta_match.group(1)
            assemblies[sample_id] = (full_path, reference)

    # Write paired reads file
    output_file = os.path.join(output_dir, "reads_list_all.tsv")
    with open(output_file, "w") as out:
        for sample_id, reads in sorted(pairs.items()):
            indiv_output_file = os.path.join(output_dir, sample_id + ".txt")
            if "R1" in reads and "R2" in reads:  # Only keep properly paired
                full_path_R1, reference = reads['R1']
                full_path_R2, reference = reads['R2']
                ref_str = ".REF" if reference else ""
                out.write(f"{sample_id + ref_str}\t{full_path_R1}\t{full_path_R2}\n")
                with open(indiv_output_file, "w") as indiv_out:
                    indiv_out.write(f"{sample_id + ref_str}\t{full_path_R1}\t{full_path_R2}\n")

        for sample_id, fasta in sorted(assemblies.items()):
            full_path, reference = fasta
            ref_str = ".REF" if reference else ""
            out.write(f"{sample_id + ref_str}\t{full_path}\n")
            indiv_output_file = os.path.join(output_dir, sample_id + ".txt")
            with open(indiv_output_file, "w") as indiv_out:
                indiv_out.write(f"{sample_id + ref_str}\t{full_path}\n")


detect_paired_fastq(snakemake.input.indir, snakemake.output.output_dir, snakemake.params.reference_reads_R1, snakemake.params.reference_reads_R2, snakemake.params.reference_assembly)

# indir = "/nfs/research/jlees/shorsfield/McClean_group_analysis/Illumina_Seq"    
# reference_reads_R1 = "/nfs/research/jlees/shorsfield/McClean_group_analysis/Illumina_Seq/Ancestral_Strain_020_S28_R1_001.fastq.gz"
# reference_reads_R2 = "/nfs/research/jlees/shorsfield/McClean_group_analysis/Illumina_Seq/Ancestral_Strain_020_S28_R2_001.fastq.gz"
# output_file = "/nfs/research/jlees/shorsfield/McClean_group_analysis/paired_reads.txt"

# detect_paired_fastq(indir, output_file, reference_reads_R1, reference_reads_R2)


    


