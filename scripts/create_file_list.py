import os
import re

def detect_paired_fastq(input_dir, output_file, reference_reads_R1=None, reference_reads_R2=None):
    """
    Detect paired FASTQ files in a directory and write them to a tab-delimited file.
    
    Args:
        input_dir (str): Path to directory containing FASTQ files.
        output_file (str): Path to output TSV file.
        reference_reads_R1 (str): Reference reads R2 to ignore.
        reference_reads_R2 (str): Reference reads R2 to ignore.
    """
    
    # Regex pattern for R1/R2 fastq files
    fastq_pattern = re.compile(r"(.+?)(?:_R|_)([12])(?:[^/]*)\.(?:fastq|fq)(?:\.gz)?$")
    
    pairs = {}
    
    for fname in os.listdir(input_dir):
        full_path = os.path.join(input_dir, fname)
        if reference_reads_R1 and reference_reads_R1 in full_path:
            continue
        if reference_reads_R2 and reference_reads_R2 in full_path:
            continue
        
        match = fastq_pattern.match(fname)
        if not match:
            continue
        
        sample_id, read_direction = match.groups()
        
        if sample_id not in pairs:
            pairs[sample_id] = {}
        
        if "1" in read_direction:
            pairs[sample_id]["R1"] = full_path
        elif "2" in read_direction:
            pairs[sample_id]["R2"] = full_path
    
    with open(output_file, "w") as out:
        for sample_id, reads in sorted(pairs.items()):
            if "R1" in reads and "R2" in reads:  # Only keep properly paired
                out.write(f"{sample_id}\t{reads['R1']}\t{reads['R2']}\n")

detect_paired_fastq(snakemake.input.indir, snakemake.output.output_file, snakemake.params.reference_reads_R1, snakemake.params.reference_reads_R2)

# indir = "/nfs/research/jlees/shorsfield/McClean_group_analysis/Illumina_Seq"    
# reference_reads_R1 = "/nfs/research/jlees/shorsfield/McClean_group_analysis/Illumina_Seq/Ancestral_Strain_020_S28_R1_001.fastq.gz"
# reference_reads_R2 = "/nfs/research/jlees/shorsfield/McClean_group_analysis/Illumina_Seq/Ancestral_Strain_020_S28_R2_001.fastq.gz"
# output_file = "/nfs/research/jlees/shorsfield/McClean_group_analysis/paired_reads.txt"

# detect_paired_fastq(indir, output_file, reference_reads_R1, reference_reads_R2)


    


