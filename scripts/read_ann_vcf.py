import os
import argparse
import pandas as pd
import vcf  # pip install PyVCF

#Define impact priority
IMPACT_ORDER = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}

def parse_vcfs(input_dir, ignore_n=False):
    all_sites = {}
    sample_variants = {}

    for file in os.listdir(input_dir):
        if file.endswith(".ann.vcf"):
            sample_name = file.replace(".ann.vcf", "")
            sample_variants[sample_name] = {}

            reader = vcf.Reader(filename=os.path.join(input_dir, file))

            for record in reader:
                # Optionally ignore N bases
                if ignore_n and ('N' in record.REF or any('N' in str(alt) for alt in record.ALT)):
                    continue

                ann_field = record.INFO.get("ANN")
                if not ann_field:
                    continue

                # Select the annotation with the highest impact
                best_ann = None
                best_score = -1

                for ann_entry in ann_field:
                    parts = ann_entry.split('|')
                    if len(parts) < 4:
                        continue
                    effect = parts[1]
                    impact = parts[2]
                    gene = parts[3]

                    score = IMPACT_ORDER.get(impact, -1)
                    if score > best_score:
                        best_ann = (gene, effect, impact)
                        best_score = score

                if best_ann is None:
                    continue  # skip if no valid annotation

                gene, effect, impact = best_ann
                site_id = f"{record.CHROM}:{record.POS}|{gene}|{impact}|{effect}"

                # Store reference allele (same across samples)
                if site_id not in all_sites:
                    all_sites[site_id] = record.REF

                # Store ALT allele(s) for this sample
                alt_alleles = ",".join(str(a) for a in record.ALT)
                sample_variants[sample_name][site_id] = alt_alleles

    return all_sites, sample_variants


def build_matrix(all_sites, sample_variants):
    sites_sorted = sorted(all_sites.keys())
    matrix = []

    # Reference row first
    ref_row = ["REF"] + [all_sites[site] for site in sites_sorted]
    matrix.append(ref_row)

    # Then samples
    for sample, vars in sample_variants.items():
        row = [sample]
        for site in sites_sorted:
            if site in vars:
                row.append(vars[site])  # ALT allele(s)
            else:
                row.append(all_sites[site])  # matches reference
        matrix.append(row)

    df = pd.DataFrame(matrix, columns=["Sample"] + sites_sorted)
    return df


def read_ann_vcf(indir, outfile, ignore_n=True):
    all_variants, sample_variants = parse_vcfs(indir, ignore_n=ignore_n)
    df = build_matrix(all_variants, sample_variants)
    df.to_csv(outfile, index=False)

read_ann_vcf(snakemake.params.indir, snakemake.output.outfile)

def main():
    parser = argparse.ArgumentParser(description="Summarize SNPeff-annotated VCFs into Excel matrix")
    parser.add_argument("--indir", help="Directory containing *.ann.vcf files")
    parser.add_argument("--outfile", help="Output Excel file")
    parser.add_argument("--ignore-n", action="store_true", help="Ignore variants with N bases")
    args = parser.parse_args()

    all_sites, sample_variants = parse_vcfs(args.indir, ignore_n=args.ignore_n)
    df = build_matrix(all_sites, sample_variants)

    df.to_csv(args.outfile, index=False)

if __name__ == "__main__":
    main()

