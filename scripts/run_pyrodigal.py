from Bio import SeqIO
import pyrodigal
import gzip
import os
import sys

def read_fasta(input_file):
    """Reads a FASTA file, gzipped or not, and returns the sequences."""
    with open(input_file, "rb") as f:
        magic = f.read(2)
    try:
        if magic == b"\x1f\x8b":
            with gzip.open(input_file, "rt") as handle:
                sequences = list(SeqIO.parse(handle, "fasta"))
        else:
            with open(input_file, "r") as handle:
                sequences = list(SeqIO.parse(handle, "fasta"))
    # if error encountered with opening, return empty list
    except:
        sequences = []
    return sequences

def get_basename(file_path):
    """Gets base name of file"""
    base_name = os.path.basename(file_path)
    file_name_without_extension, _ = os.path.splitext(base_name)
    return file_name_without_extension

def run_pyrodigal(file_name, output_dir, basename=None):
    """Trains and runs Pyrodigal on input sequences."""
    
    if not os.path.exists(output_dir):
            os.mkdir(output_dir)

    sequences = read_fasta(file_name)

    # train the orf_finder
    orf_finder = pyrodigal.GeneFinder()
    try:
        orf_finder.train("TTAATTAATTAA".join([str(record.seq) for record in sequences]))
    # pass error with pyrodigal for given file
    except ValueError as e:
        print(f"Error with pyrodigal for file {file_name}: {str(e)}")
        sys.exit(0)

    # get file basename
    if basename == None:
        basename = get_basename(file_name)

    # generate output filename
    output_path = os.path.join(output_dir, basename)

    # iterate over records and generate gene sequences, writing sequences to file
    protein_records = []
    DNA_records = []
    with open(output_path + ".gff", "w") as o_gff, open(output_path + ".faa", "w") as o_faa,  open(output_path + ".ffn", "w") as o_ffn, open(output_path + ".gbk", "w") as o_gbk:               
        for record in sequences:
            # find genes in contig
            orfs = orf_finder.find_genes(str(record.seq))

            # write to open files
            orfs.write_genes(o_ffn, basename + "_" + record.id, full_id=True)
            orfs.write_translations(o_faa, basename + "_" + record.id, full_id=True)
            orfs.write_gff(o_gff, basename + "_" + record.id, include_translation_table=False, full_id=True)
            orfs.write_genbank(o_gbk, basename + "_" + record.id)
        
run_pyrodigal(snakemake.input.contigs, snakemake.output.outdir, snakemake.params.basename)
