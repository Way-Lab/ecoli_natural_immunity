#!/usr/bin/env python3
"""
Generate full multi-FASTA alignments showing the peptide region from ALL sequences.
Each peptide will have one FASTA file with all 130 sequences (or however many contain it).
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import os

# Define peptide sequences (only the ones found)
PEPTIDES = {
    'Peptide_1': 'WSQYHDTGFINNNGPTHEN',
    'Peptide_2': 'GRMPYKGSVENGAYKA',
    'Peptide_3a': 'KSNVYGKN',
    'Peptide_3b': 'KANVPGGASFKD',
    'Peptide_4': 'TNNIGDAHTIGTRPDNGM',
    'Peptide_2b': 'GRMPYKGDNINGAYKA',
}


def find_peptide_in_sequence(peptide, sequence, allow_mismatches=3):
    """
    Find peptide in sequence allowing for some mismatches.
    Returns list of (position, num_mismatches, matched_subsequence) tuples.
    Default allows up to 3 mismatches.
    """
    matches = []
    peptide = peptide.upper()
    sequence = str(sequence).upper().replace('-', '')

    pep_len = len(peptide)
    seq_len = len(sequence)

    for i in range(seq_len - pep_len + 1):
        subseq = sequence[i:i + pep_len]
        mismatches = sum(1 for a, b in zip(peptide, subseq) if a != b)

        if mismatches <= allow_mismatches:
            matches.append((i, mismatches, subseq))

    return matches


def extract_species_from_header(header):
    """Extract clean species/strain name from FASTA header."""
    header_lower = header.lower()

    # Check for specific E. coli strain identifiers
    if 'scb' in header_lower:
        # Extract SCB number
        parts = header.split()
        for part in parts:
            if 'scb' in part.lower():
                return f"E_coli_{part.strip()}"

    if 'nissle' in header_lower or 'ecn' in header_lower:
        return 'E_coli_Nissle'
    elif 'o157' in header_lower or 'sakai' in header_lower:
        return 'E_coli_O157H7'
    elif 'ec958' in header_lower:
        return 'E_coli_EC958'
    elif 'uti189' in header_lower:
        return 'E_coli_UTI189'
    elif 'cft073' in header_lower:
        return 'E_coli_CFT073'
    elif 'dh5' in header_lower:
        return 'E_coli_DH5alpha'
    elif 'k-12' in header_lower or 'mg1655' in header_lower:
        return 'E_coli_K12'
    elif 'coli' in header_lower:
        # Generic E. coli - use accession
        acc = header.split()[0].split('|')[0]
        return f"E_coli_{acc}"
    elif 'klebsiella' in header_lower:
        acc = header.split()[0].split('|')[0]
        if 'pneumoniae' in header_lower:
            return f'Klebsiella_pneumoniae_{acc}'
        elif 'oxytoca' in header_lower:
            return f'Klebsiella_oxytoca_{acc}'
        else:
            return f'Klebsiella_{acc}'
    elif 'enterobacter' in header_lower:
        acc = header.split()[0].split('|')[0]
        if 'cloacae' in header_lower:
            return f'Enterobacter_cloacae_{acc}'
        elif 'aerogenes' in header_lower:
            return f'Enterobacter_aerogenes_{acc}'
        else:
            return f'Enterobacter_{acc}'
    elif 'citrobacter' in header_lower:
        acc = header.split()[0].split('|')[0]
        if 'freundii' in header_lower:
            return f'Citrobacter_freundii_{acc}'
        elif 'koseri' in header_lower:
            return f'Citrobacter_koseri_{acc}'
        else:
            return f'Citrobacter_{acc}'
    elif 'salmonella' in header_lower:
        acc = header.split()[0].split('|')[0]
        return f'Salmonella_{acc}'
    elif 'shigella' in header_lower:
        acc = header.split()[0].split('|')[0]
        if 'flexneri' in header_lower:
            return f'Shigella_flexneri_{acc}'
        elif 'sonnei' in header_lower:
            return f'Shigella_sonnei_{acc}'
        else:
            return f'Shigella_{acc}'
    elif 'proteus' in header_lower:
        acc = header.split()[0].split('|')[0]
        return f'Proteus_mirabilis_{acc}'
    elif 'serratia' in header_lower:
        acc = header.split()[0].split('|')[0]
        return f'Serratia_marcescens_{acc}'
    else:
        # Use first accession/ID
        return header.split()[0].split('|')[0][:30]


def create_full_alignment_fasta(peptide_name, peptide_seq, fasta_file, output_file):
    """
    Create a multi-FASTA file with ALL sequences that contain this peptide region.
    """
    print(f"\nProcessing {peptide_name}...")

    sequences = list(SeqIO.parse(fasta_file, 'fasta'))

    # Collect all sequences with this peptide
    found_sequences = []
    not_found = []

    for seq_record in sequences:
        seq_id = seq_record.id
        sequence = str(seq_record.seq).upper()

        matches = find_peptide_in_sequence(peptide_seq, sequence, allow_mismatches=3)

        if matches:
            # Use best match
            best_match = min(matches, key=lambda x: x[1])
            position, mismatches, matched_seq = best_match

            # Create clean species name
            species_name = extract_species_from_header(seq_record.description)

            found_sequences.append({
                'species': species_name,
                'original_id': seq_id,
                'sequence': matched_seq,
                'position': position,
                'mismatches': mismatches
            })
        else:
            not_found.append(seq_record.description)

    # Sort by species name
    found_sequences.sort(key=lambda x: x['species'])

    # Write FASTA file
    with open(output_file, 'w') as f:
        # Write header
        f.write(f"# Multi-sequence alignment for {peptide_name}\n")
        f.write(f"# Reference: {peptide_seq}\n")
        f.write(f"# Found in: {len(found_sequences)}/{len(sequences)} sequences\n")
        f.write(f"# Not found in: {len(not_found)} sequences\n")
        f.write("#\n")
        f.write(f"# Format: >Species_Name | original_id | position:X | mismatches:Y\n")
        f.write("#\n\n")

        # Write reference sequence first
        f.write(f">REFERENCE_{peptide_name}\n")
        f.write(f"{peptide_seq}\n\n")

        # Write all found sequences
        for idx, seq_info in enumerate(found_sequences, 1):
            header = f">{seq_info['species']} | {seq_info['original_id']} | pos:{seq_info['position']} | mm:{seq_info['mismatches']}"
            f.write(f"{header}\n")
            f.write(f"{seq_info['sequence']}\n")

    print(f"  Found in: {len(found_sequences)}/{len(sequences)} sequences")
    print(f"  Written to: {output_file}")

    # Also create a summary text file
    summary_file = output_file.replace('.fasta', '_summary.txt')
    with open(summary_file, 'w') as f:
        f.write(f"{'='*80}\n")
        f.write(f"ALIGNMENT SUMMARY: {peptide_name}\n")
        f.write(f"{'='*80}\n\n")
        f.write(f"Reference sequence: {peptide_seq}\n")
        f.write(f"Length: {len(peptide_seq)} amino acids\n\n")
        f.write(f"Found in: {len(found_sequences)}/{len(sequences)} sequences ({100*len(found_sequences)/len(sequences):.1f}%)\n\n")

        if not_found:
            f.write(f"NOT FOUND IN ({len(not_found)} sequences):\n")
            for desc in not_found[:10]:
                f.write(f"  - {desc}\n")
            if len(not_found) > 10:
                f.write(f"  ... and {len(not_found)-10} more\n")
            f.write("\n")

        # Group by species
        species_counts = defaultdict(list)
        for seq_info in found_sequences:
            species = seq_info['species'].split('_')[0]
            species_counts[species].append(seq_info)

        f.write(f"SPECIES DISTRIBUTION:\n")
        for species in sorted(species_counts.keys()):
            f.write(f"  {species}: {len(species_counts[species])} sequences\n")
        f.write("\n")

        # Count unique sequences
        unique_seqs = defaultdict(list)
        for seq_info in found_sequences:
            unique_seqs[seq_info['sequence']].append(seq_info['species'])

        f.write(f"UNIQUE VARIANTS: {len(unique_seqs)}\n")
        for seq, species_list in sorted(unique_seqs.items(), key=lambda x: len(x[1]), reverse=True):
            mismatches = sum(1 for a, b in zip(peptide_seq, seq) if a != b)
            f.write(f"\n  {seq} ({len(species_list)} sequences, {mismatches} mm)\n")
            for species in species_list[:5]:
                f.write(f"    - {species}\n")
            if len(species_list) > 5:
                f.write(f"    ... and {len(species_list)-5} more\n")

    print(f"  Summary written to: {summary_file}")

    return len(found_sequences), len(sequences)


def main():
    import sys

    if len(sys.argv) < 2:
        print("Usage: python generate_full_alignments.py <fasta_file> [output_dir]")
        sys.exit(1)

    fasta_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "full_peptide_alignments"

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    print(f"Generating full multi-FASTA alignments for all peptides...")
    print(f"Input: {fasta_file}")
    print(f"Output directory: {output_dir}/\n")

    results = []
    for peptide_name, peptide_seq in PEPTIDES.items():
        output_file = f"{output_dir}/{peptide_name}_all_sequences.fasta"
        found, total = create_full_alignment_fasta(peptide_name, peptide_seq, fasta_file, output_file)
        results.append((peptide_name, found, total))

    # Create master index
    with open(f"{output_dir}/README.txt", 'w') as f:
        f.write("FULL PEPTIDE ALIGNMENTS\n")
        f.write("="*80 + "\n\n")
        f.write("This directory contains multi-FASTA alignments showing the peptide region\n")
        f.write("from ALL sequences analyzed (not just unique variants).\n\n")
        f.write("For each peptide:\n")
        f.write("  - {peptide}_all_sequences.fasta = Multi-FASTA with all sequences\n")
        f.write("  - {peptide}_all_sequences_summary.txt = Summary statistics\n\n")
        f.write("SUMMARY:\n")
        f.write("-"*80 + "\n")
        for peptide_name, found, total in results:
            f.write(f"{peptide_name:15} Found in {found:3}/{total} sequences ({100*found/total:5.1f}%)\n")
        f.write("\n")
        f.write("These files can be opened in:\n")
        f.write("  - Jalview, MEGA, Geneious, SnapGene (FASTA visualization)\n")
        f.write("  - Text editors (human-readable format)\n")
        f.write("  - ClustalW/MUSCLE for re-alignment if needed\n")

    print(f"\n{'='*80}")
    print("COMPLETE!")
    print(f"{'='*80}")
    print(f"\nAll alignment files written to: {output_dir}/")
    print(f"See {output_dir}/README.txt for details\n")


if __name__ == "__main__":
    main()
