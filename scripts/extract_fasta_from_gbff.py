from Bio import SeqIO
import argparse

def get_options():
    description = 'Gets a nearesty neighbour for each tip on a tree.'
    
    parser = argparse.ArgumentParser(description=description,
                                     prog='python get_nearest_tips.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--infile',
                    required=True,
                    help='gbff file.')
    IO.add_argument('--outfile',
                    required=True,
                    help='Output file.')

    return parser.parse_args()

def main():
    options = get_options()
    input_file = options.infile
    output_file = options.outfile

    with open(output_file, "w") as out:
        for record in SeqIO.parse(input_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    # Skip if translation or locus tag is missing
                    if "locus_tag" not in feature.qualifiers:
                        continue
                    if "translation" not in feature.qualifiers:
                        continue

                    locus_tag = feature.qualifiers["locus_tag"][0]
                    protein_seq = feature.qualifiers["translation"][0]

                    out.write(f">{locus_tag} {locus_tag}~~~{locus_tag}~~~\n{protein_seq}\n")

if __name__ == "__main__":
    main()