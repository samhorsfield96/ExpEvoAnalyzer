#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import vcf  # pip install PyVCF

# Define impact priority (unused in build logic but kept for annotation selection)
IMPACT_ORDER = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}


def parse_vcf_file(vcf_path, ignore_n=False):
    """
    Parse one annotated VCF and return:
      variants: dict { "chr:pos" : "ALT1,ALT2|GENE|IMPACT|EFFECT" }
    Only sites that have ANN entries are returned (consistent with previous behaviour).
    """
    variants = {}
    reader = vcf.Reader(filename=vcf_path)

    for record in reader:
        if ignore_n and ('N' in record.REF or any('N' in str(alt) for alt in record.ALT)):
            continue

        pos_id = f"{record.CHROM}:{record.POS}"
        ann_field = record.INFO.get("ANN")
        if not ann_field:
            continue

        # Pick best annotation (highest impact) among ANN entries
        best_ann = None
        best_score = -1
        for ann_entry in ann_field:
            parts = ann_entry.split('|')
            if len(parts) < 4:
                continue
            effect = parts[1]
            impact = parts[2]
            gene = parts[3]
            locus = parts[4]
            score = IMPACT_ORDER.get(impact, -1)
            if score > best_score:
                best_ann = (gene, locus, effect, impact)
                best_score = score

        if best_ann is None:
            continue

        gene, locus, effect, impact = best_ann
        alt_alleles = ",".join(str(a) for a in record.ALT)
        variants[pos_id] = f"{alt_alleles}|{gene}|{locus}|{impact}|{effect}"

    return variants


def parse_vcfs(input_dir, ignore_n=False):
    """
    Parse all *.ann.vcf in a directory.
    - auto-detects a single *.REF.ann.vcf (first one found) as assembly reference
    - returns:
        all_sites: dict { pos_id : canonical_ref_base }
        sample_variants: dict { sample_name : { pos_id : "ALT|GENE|IMPACT|EFFECT" } }
        ref_variants: dict for the REF.vcf in the same format as sample_variants
    """
    all_sites = {}
    sample_variants = {}
    ref_variants = {}
    ref_vcf = None

    # detect ref vcf
    for fname in os.listdir(input_dir):
        if fname.endswith(".REF.ann.vcf"):
            ref_vcf = os.path.join(input_dir, fname)
            break

    # parse sample files (exclude REF file)
    for fname in os.listdir(input_dir):
        if not fname.endswith(".ann.vcf"):
            continue
        if fname.endswith(".REF.ann.vcf"):
            continue
        path = os.path.join(input_dir, fname)
        sample_name = fname.replace(".ann.vcf", "")
        variants = parse_vcf_file(path, ignore_n=ignore_n)
        sample_variants[sample_name] = variants

        # collect canonical REF base for positions seen in this file
        reader = vcf.Reader(filename=path)
        for record in reader:
            pos_id = f"{record.CHROM}:{record.POS}"
            if pos_id not in all_sites:
                all_sites[pos_id] = record.REF

    # if no sample provided a site but ref_vcf has it, we should still capture the reference base from ref_vcf
    if ref_vcf:
        # collect reference bases from ref_vcf if missing in all_sites
        reader = vcf.Reader(filename=ref_vcf)
        for record in reader:
            pos_id = f"{record.CHROM}:{record.POS}"
            if pos_id not in all_sites:
                all_sites[pos_id] = record.REF

        # parse annotations from ref_vcf (assembly alleles)
        ref_variants = parse_vcf_file(ref_vcf, ignore_n=ignore_n)

    return all_sites, sample_variants, ref_variants, ref_vcf


def build_matrix(all_sites, sample_variants, ref_variants, ignore_ref=True, debug=False):
    """
    Build matrix where columns are strictly positions (CHROM:POS).
    For each site:
      - if ignore_ref and ref_variants has assembly alleles for that site:
          * any sample allele equal to an assembly allele is replaced with canonical reference base
      - after that replacement, if ALL samples equal the canonical reference -> drop the site (assembly-only)
      - otherwise keep site; cell value is either canonical_ref or "ALT|GENE|IMPACT|EFFECT"
    """
    sites_sorted = sorted(all_sites.keys())
    samples = sorted(sample_variants.keys())  # deterministic order

    keep_sites = []
    per_site_final_values = {}  # site -> { sample: value }
    dropped_sites = []

    for site in sites_sorted:
        canonical_ref = all_sites[site]
        assembly_alleles = set()
        if ignore_ref and site in ref_variants:
            assembly_str = ref_variants[site].split('|')[0]
            assembly_alleles = set(a for a in assembly_str.split(',') if a != '')
        # Build final values for each sample
        final_values = {}
        for sample in samples:
            vars_for_sample = sample_variants[sample]
            if site in vars_for_sample:
                allele_full = vars_for_sample[site]        # e.g. "T|GENE|MODERATE|missense_variant"
                allele_base = allele_full.split('|')[0]   # e.g. "T" or "T,G"
                # if multi-allelic allele_base contains commas, treat as set
                allele_bases_set = set(allele_base.split(','))
                # if assembly alleles present and all alleles of this sample are in assembly_alleles -> replace
                if assembly_alleles and allele_bases_set.issubset(assembly_alleles):
                    final_values[sample] = canonical_ref
                else:
                    final_values[sample] = allele_full
            else:
                # no variant reported in sample -> canonical reference
                final_values[sample] = canonical_ref

        # decide whether to drop: drop only if every sample equals canonical_ref after replacement
        if all(val == canonical_ref for val in final_values.values()):
            dropped_sites.append(site)
            continue
        # else keep
        per_site_final_values[site] = final_values
        keep_sites.append(site)

    if debug:
        print(f"Total sites found: {len(sites_sorted)}")
        print(f"Reference VCF provided: {'yes' if ref_variants else 'no'}")
        print(f"Sites dropped due to REF-matching assembly alleles: {len(dropped_sites)}")
        if dropped_sites:
            print("Example dropped sites (up to 10):", dropped_sites[:10])
        print(f"Sites kept: {len(keep_sites)}")

    # build matrix rows
    matrix = []
    # Reference row (canonical refs for kept sites)
    ref_row = ["REF"] + [all_sites[site] for site in keep_sites]
    matrix.append(ref_row)

    # Sample rows
    for sample in samples:
        row = [sample] + [per_site_final_values[site][sample] for site in keep_sites]
        matrix.append(row)

    # Build dataframe with columns ["Sample"] + keep_sites
    columns = ["Sample"] + keep_sites
    df = pd.DataFrame(matrix, columns=columns)
    return df, dropped_sites

def read_ann_vcf(indir, outfile, ignore_n=True, ignore_ref=True, debug=False): 
    all_sites, sample_variants, ref_variants, ref_vcf = parse_vcfs(indir, ignore_n=ignore_n)
    df, dropped_sites = build_matrix(all_sites, sample_variants, ref_variants, ignore_ref=ignore_ref, debug=debug)

    df.to_csv(outfile, index=False)

read_ann_vcf(snakemake.params.indir, snakemake.output.outfile)

# def main():
#     parser = argparse.ArgumentParser(description="Summarize SNPeff-annotated VCFs into a matrix (CHROM:POS columns).")
#     parser.add_argument("--indir", required=True, help="Directory containing *.ann.vcf files")
#     parser.add_argument("--outfile", required=True, help="Output CSV file (or .xlsx if you prefer df.to_excel)")
#     parser.add_argument("--ignore-n", action="store_true", help="Ignore variants with N bases")
#     parser.add_argument("--ignore-ref", action="store_true", help="Use *.REF.ann.vcf to filter assembly-matching alleles")
#     parser.add_argument("--debug", action="store_true", help="Print debug info (counts, example dropped sites)")
#     args = parser.parse_args()

#     all_sites, sample_variants, ref_variants, ref_vcf = parse_vcfs(args.indir, ignore_n=args.ignore_n)
#     df, dropped_sites = build_matrix(all_sites, sample_variants, ref_variants, ignore_ref=args.ignore_ref, debug=args.debug)

#     df.to_csv(args.outfile, index=False)

#     if args.debug:
#         print(f"Wrote output to {args.outfile}")
#         if ref_vcf:
#             print(f"Used REF VCF: {ref_vcf}")
#         if len(dropped_sites) == 0:
#             print("No sites were dropped.")
#         else:
#             print(f"Dropped {len(dropped_sites)} sites (example first 10): {dropped_sites[:10]}")


# if __name__ == "__main__":
#     main()


