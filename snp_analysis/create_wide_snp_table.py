#!/usr/bin/env python3
"""
Create a wide-format SNP table showing which nucleotide each strain has at each position.
"""

import re
import pandas as pd
from collections import defaultdict

def parse_snps_by_sample(snp_file):
    """Parse SNP file and return dict of {sample: {position: alt_allele}}"""
    sample_snps = defaultdict(dict)
    current_sample = None

    with open(snp_file, 'r') as f:
        for line in f:
            line = line.strip()

            # Skip header lines
            if line.startswith('SNP Details') or line.startswith('Core genome') or line.startswith('===') or not line:
                continue

            # Check for sample name
            sample_match = re.match(r'^(.+):\s+(\d+)\s+SNPs', line)
            if sample_match:
                current_sample = sample_match.group(1)
                continue

            # Check for separator
            if line.startswith('---'):
                continue

            # Parse SNP line: Position  550880:   T>G (Ref: T, Alt: G)
            snp_match = re.match(r'Position\s+(\d+):\s+\S+>\S+\s+\(Ref:\s+(\S+),\s+Alt:\s+(\S+)\)', line)
            if snp_match:
                position = int(snp_match.group(1))
                ref = snp_match.group(2)
                alt = snp_match.group(3)

                # Skip N (ambiguous)
                if alt == 'N':
                    continue

                # Handle multiple alternatives (e.g., "C,G") - take first one
                if ',' in alt:
                    alt = alt.split(',')[0]

                sample_snps[current_sample][position] = alt

    return sample_snps

def load_annotations(annot_file):
    """Load the SNP annotations we created earlier."""
    df = pd.read_csv(annot_file, sep='\t')

    # Create a lookup dict: position -> annotation info
    annot_dict = {}
    for _, row in df.iterrows():
        pos = row['Position']
        if pos not in annot_dict:  # Take first annotation for each position
            annot_dict[pos] = row

    return annot_dict

def format_effect(row):
    """Format the effect string like 'Gly202Val' or blank for intergenic."""
    if row['Effect'] == 'intergenic':
        return ''
    elif row['Effect'] == 'synonymous':
        aa_pos = int(row['AA_Pos']) if pd.notna(row['AA_Pos']) else ''
        return f"{row['Ref_AA']}{aa_pos}{row['Alt_AA']}" if aa_pos else ''
    elif row['Effect'] == 'missense':
        aa_pos = int(row['AA_Pos']) if pd.notna(row['AA_Pos']) else ''
        return f"{row['Ref_AA']}{aa_pos}{row['Alt_AA']}" if aa_pos else ''
    elif row['Effect'] == 'nonsense':
        aa_pos = int(row['AA_Pos']) if pd.notna(row['AA_Pos']) else ''
        return f"{row['Ref_AA']}{aa_pos}*" if aa_pos else ''
    elif row['Effect'] == 'readthrough':
        aa_pos = int(row['AA_Pos']) if pd.notna(row['AA_Pos']) else ''
        return f"*{aa_pos}{row['Alt_AA']}" if aa_pos else ''
    else:
        return row['Effect']

def create_wide_table(sample_snps, annot_dict, output_file):
    """Create wide-format table."""

    # Get all samples, excluding EcN RefSeq and any *.ref duplicates
    all_samples = sorted([s for s in sample_snps.keys()
                         if s != 'EcN RefSeq' and not s.endswith('.ref')])

    # Get all unique positions from these samples only
    all_positions = set()
    for sample in all_samples:
        all_positions.update(sample_snps[sample].keys())

    all_positions = sorted(all_positions)

    # Reference is Stock-1
    ref_name = "Stock-1"

    # Build table
    rows = []

    for pos in all_positions:
        # Get annotation for this position
        if pos in annot_dict:
            annot = annot_dict[pos]
            ref_nt = annot['Ref']

            # Determine protein name
            if annot['Gene'] and annot['Gene'] != 'intergenic':
                protein = annot['Product']
                if annot['Gene']:
                    protein = f"{annot['Product']};{annot['Gene']}"
            else:
                protein = 'intergenic'

            # Determine change
            change = f"{ref_nt}>{annot['Alt']}"

            # Format effect
            effect_str = format_effect(annot)
        else:
            # Shouldn't happen, but handle gracefully
            ref_nt = '?'
            protein = 'unknown'
            change = '?'
            effect_str = ''

        # Build row
        row = {
            'Ref': ref_name,
            'POS': pos,
            'Protein': protein,
            'REF': ref_nt,
            'Change (nt)': change,
            'Effect': effect_str
        }

        # Add each sample's nucleotide at this position
        for sample in all_samples:
            if pos in sample_snps[sample]:
                # This sample has the alternate allele
                row[sample] = sample_snps[sample][pos]
            else:
                # This sample has the reference allele
                row[sample] = ref_nt

        rows.append(row)

    # Create DataFrame
    columns = ['Ref', 'POS', 'Protein', 'REF', 'Change (nt)', 'Effect'] + all_samples
    df = pd.DataFrame(rows, columns=columns)

    # Write to file
    df.to_csv(output_file, sep='\t', index=False)

    return df

def main():
    snp_file = '/fastpool/active_data/ecoli_genomics/parsnp_with_outgroups/parsnp_output/snp_details_by_sample_filtered.txt'
    annot_file = '/fastpool/active_data/ecoli_genomics/snp_annotations_with_effects.tsv'
    output_file = '/fastpool/active_data/ecoli_genomics/snp_table_wide_format.tsv'

    print("Parsing SNPs by sample...")
    sample_snps = parse_snps_by_sample(snp_file)

    print(f"Found {len(sample_snps)} samples:")
    for sample, snps in sorted(sample_snps.items()):
        print(f"  {sample}: {len(snps)} SNPs")

    print("\nLoading annotations...")
    annot_dict = load_annotations(annot_file)
    print(f"Loaded annotations for {len(annot_dict)} positions")

    print("\nCreating wide-format table...")
    df = create_wide_table(sample_snps, annot_dict, output_file)

    print(f"\nTable created with {len(df)} rows and {len(df.columns)} columns")
    print(f"Saved to: {output_file}")

    print("\nFirst 10 rows:")
    print(df.head(10).to_string())

if __name__ == '__main__':
    main()
