#!/usr/bin/env python3
"""
Extract Murein lipoprotein (Lpp/MLP) sequences from E. coli strains.

Strains to analyze:
- EcN (Nissle 1917)
- UTI89
- RS218
- SCB12, SCB29, SCB61, SCB34, SCB58, SCB37, SCB60
"""

import os
import subprocess
from Bio import SeqIO

# Known reference sequences from UniProt/NCBI
# Lpp (Murein lipoprotein) is highly conserved at 78 aa
REFERENCE_LPP = {
    'K12_MG1655': 'MKATKLVLGAVILGSTLLAGCSSNAKIDQLSSDVQTLNAKVDQLSNDVNAMRSDVQAAKDDAARANQRLDNMATKYRK',
}

# N2SKTQ sample annotations directory
N2SKTQ_BASE = "N2SKTQ_results"

def extract_lpp_from_faa(faa_file, strain_name):
    """Extract Lpp sequence from a protein FASTA file."""
    if not os.path.exists(faa_file):
        print(f"  File not found: {faa_file}")
        return None

    for record in SeqIO.parse(faa_file, 'fasta'):
        desc = record.description.lower()
        # Look for murein lipoprotein or lpp
        if 'murein lipoprotein' in desc or ' lpp ' in desc or desc.endswith(' lpp'):
            return str(record.seq)

    return None


def extract_lpp_from_tsv(tsv_file, fna_file, strain_name):
    """Extract Lpp gene coordinates from TSV and get sequence from FNA."""
    if not os.path.exists(tsv_file) or not os.path.exists(fna_file):
        return None

    # Parse TSV for lpp gene location
    lpp_info = None
    with open(tsv_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                gene = parts[6] if len(parts) > 6 else ""
                product = parts[7] if len(parts) > 7 else ""
                if gene.lower() == 'lpp' or 'murein lipoprotein' in product.lower():
                    lpp_info = {
                        'contig': parts[0],
                        'start': int(parts[2]),
                        'end': int(parts[3]),
                        'strand': parts[4],
                        'gene': gene,
                        'product': product
                    }
                    break

    return lpp_info


def main():
    print("="*80)
    print("MUREIN LIPOPROTEIN (Lpp/MLP) SEQUENCE EXTRACTION")
    print("="*80)
    print()

    # Collect all Lpp sequences
    lpp_sequences = {}

    # Reference sequence
    lpp_sequences['K12_MG1655'] = REFERENCE_LPP['K12_MG1655']
    print(f"Reference K12_MG1655: {len(REFERENCE_LPP['K12_MG1655'])} aa")

    # Extract from N2SKTQ samples
    print("\nExtracting from N2SKTQ samples:")
    for i in range(1, 15):
        sample_name = f"N2SKTQ_{i}_{i}"
        faa_file = f"{N2SKTQ_BASE}/{sample_name}/ONT-only/annotation/{sample_name}.faa"

        if os.path.exists(faa_file):
            seq = extract_lpp_from_faa(faa_file, sample_name)
            if seq:
                lpp_sequences[sample_name] = seq
                print(f"  {sample_name}: {len(seq)} aa - Found")
            else:
                print(f"  {sample_name}: Lpp not found")
        else:
            print(f"  {sample_name}: No annotation file")

    # Try to get reference strain sequences from NCBI dataset
    ncbi_protein = "ncbi_dataset/data/GCF_000800845.1/protein.faa"
    if os.path.exists(ncbi_protein):
        print(f"\nExtracting from NCBI reference:")
        for record in SeqIO.parse(ncbi_protein, 'fasta'):
            desc = record.description.lower()
            if 'murein lipoprotein lpp' in desc:
                # This is EcN (Nissle 1917)
                lpp_sequences['EcN'] = str(record.seq)
                print(f"  EcN: {len(record.seq)} aa - Found")
                break

    # Write output FASTA
    output_file = "lpp_sequences.fasta"
    with open(output_file, 'w') as f:
        for name, seq in sorted(lpp_sequences.items()):
            f.write(f">{name}\n")
            f.write(f"{seq}\n")

    print(f"\nWrote {len(lpp_sequences)} sequences to {output_file}")

    # Show sequences for comparison
    print("\n" + "="*80)
    print("LPP SEQUENCES COMPARISON")
    print("="*80)

    ref_seq = REFERENCE_LPP['K12_MG1655']
    print(f"\nReference (K12): {ref_seq}")
    print(f"Length: {len(ref_seq)} aa")

    print("\nComparing all sequences to K12 reference:")
    for name, seq in sorted(lpp_sequences.items()):
        if name == 'K12_MG1655':
            continue

        if seq == ref_seq:
            status = "IDENTICAL"
        else:
            # Count differences
            diffs = sum(1 for a, b in zip(seq, ref_seq) if a != b)
            len_diff = abs(len(seq) - len(ref_seq))
            status = f"{diffs} aa differences, length diff: {len_diff}"

        print(f"  {name}: {status}")

    print("\n" + "="*80)


if __name__ == "__main__":
    main()
