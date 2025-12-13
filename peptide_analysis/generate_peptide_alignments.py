#!/usr/bin/env python3
"""
Generate multi-FASTA alignments for each peptide showing all variants found.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import os

# Define peptide sequences
PEPTIDES = {
    'Peptide_1': 'WSQYHDTGFINNNGPTHEN',
    'Peptide_2': 'GRMPYKGSVENGAYKA',
    'Peptide_3a': 'KSNVYGKN',
    'Peptide_3b': 'KANVPGGASFKD',
    'Peptide_4': 'TNNIGDAHTIGTRPDNGM',
    'Peptide_2b': 'GRMPYKGDNINGAYKA',

    # These won't be found but we'll try
    'Peptide_1S': 'TQPNSHNDINTNHGYGWEF',
    'Peptide_2S': 'RSMAGKGPKNAGYVYE',
    'Peptide_3aS': 'KYNGVNSK',
    'Peptide_3bS': 'FKAKDNGVGSPA',
    'Peptide_4S': 'MTNTHNTAIGDRNIGPGD',
    'Peptide_3BS2': 'GVSFDKAKNAPG',
    'Peptide_2bS': 'ANGGGYIADYNKRKPM',
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
    """Extract species/strain information from FASTA header."""
    header_lower = header.lower()

    # Check for specific identifiers
    if 'ecn' in header_lower or 'nissle' in header_lower:
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
    elif 'scb' in header_lower:
        scb_match = header.split()[0]
        return f'E_coli_{scb_match}'
    elif 'coli' in header_lower:
        return 'E_coli'
    elif 'klebsiella' in header_lower:
        if 'pneumoniae' in header_lower:
            return 'K_pneumoniae'
        elif 'oxytoca' in header_lower:
            return 'K_oxytoca'
        else:
            return 'Klebsiella_sp'
    elif 'enterobacter' in header_lower:
        if 'cloacae' in header_lower:
            return 'E_cloacae'
        elif 'aerogenes' in header_lower:
            return 'E_aerogenes'
        else:
            return 'Enterobacter_sp'
    elif 'citrobacter' in header_lower:
        if 'freundii' in header_lower:
            return 'C_freundii'
        elif 'koseri' in header_lower:
            return 'C_koseri'
        else:
            return 'Citrobacter_sp'
    elif 'salmonella' in header_lower:
        return 'Salmonella'
    elif 'shigella' in header_lower:
        if 'flexneri' in header_lower:
            return 'S_flexneri'
        elif 'sonnei' in header_lower:
            return 'S_sonnei'
        else:
            return 'Shigella_sp'
    elif 'proteus' in header_lower:
        return 'P_mirabilis'
    elif 'serratia' in header_lower:
        return 'S_marcescens'
    else:
        # Try to get first part of ID
        return header.split()[0].split('|')[0][:20]


def create_aligned_fasta_with_reference(peptide_name, peptide_seq, variants_dict, output_file):
    """
    Create a FASTA file with reference sequence and all variants aligned.
    Show differences with lowercase and special markers.
    """
    with open(output_file, 'w') as f:
        # Write header
        f.write(f"# Alignment for {peptide_name}\n")
        f.write(f"# Reference sequence: {peptide_seq}\n")
        f.write(f"# Total unique variants: {len(variants_dict)}\n")
        f.write(f"# Variable positions marked with lowercase in variants\n")
        f.write("#\n")

        # Write reference sequence first
        f.write(f">REFERENCE_{peptide_name}\n")
        f.write(f"{peptide_seq}\n")
        f.write("\n")

        # Sort variants by number of sequences (most common first)
        sorted_variants = sorted(variants_dict.items(),
                                key=lambda x: len(x[1]),
                                reverse=True)

        # Write each variant
        for variant_seq, seq_records in sorted_variants:
            # Calculate mismatches
            mismatches = sum(1 for a, b in zip(peptide_seq, variant_seq) if a != b)

            # Create a version showing differences
            highlighted_seq = ""
            for ref_aa, var_aa in zip(peptide_seq, variant_seq):
                if ref_aa == var_aa:
                    highlighted_seq += var_aa
                else:
                    highlighted_seq += var_aa.lower()  # lowercase = different

            # Get list of species/strains
            species_list = [extract_species_from_header(rec['id']) for rec in seq_records]
            species_counts = defaultdict(int)
            for sp in species_list:
                species_counts[sp] += 1

            species_summary = ", ".join([f"{sp}({count})" for sp, count in
                                        sorted(species_counts.items(),
                                              key=lambda x: x[1],
                                              reverse=True)])

            # Write variant header
            f.write(f">VARIANT_{mismatches}mm_n{len(seq_records)} | {species_summary}\n")
            f.write(f"{highlighted_seq}\n")

            # Write actual sequence (all uppercase) on next line for reference
            f.write(f"# {variant_seq} (mismatches: {mismatches})\n")

            # List individual sequences
            for rec in seq_records[:5]:  # Show first 5
                f.write(f"#   - {rec['id']}\n")
            if len(seq_records) > 5:
                f.write(f"#   - ... and {len(seq_records) - 5} more\n")
            f.write("\n")


def create_simple_alignment(peptide_name, peptide_seq, variants_dict, output_file):
    """
    Create a simple visual alignment showing all variants.
    """
    with open(output_file, 'w') as f:
        f.write(f"{'='*80}\n")
        f.write(f"PEPTIDE ALIGNMENT: {peptide_name}\n")
        f.write(f"{'='*80}\n\n")

        if not variants_dict:
            f.write("*** PEPTIDE NOT FOUND IN ANY SEQUENCES ***\n\n")
            f.write(f"Reference: {peptide_seq}\n")
            return

        # Position numbers
        f.write("Position:  ")
        for i in range(len(peptide_seq)):
            f.write(f"{(i+1):>3}")
        f.write("\n")
        f.write("           " + "---" * len(peptide_seq) + "\n")

        # Reference sequence
        f.write(f"REFERENCE: ")
        for aa in peptide_seq:
            f.write(f"  {aa}")
        f.write(f"   (reference)\n")
        f.write("           " + "---" * len(peptide_seq) + "\n")

        # Sort variants by frequency
        sorted_variants = sorted(variants_dict.items(),
                                key=lambda x: len(x[1]),
                                reverse=True)

        # Write each variant
        for idx, (variant_seq, seq_records) in enumerate(sorted_variants, 1):
            # Calculate mismatches
            mismatches = sum(1 for a, b in zip(peptide_seq, variant_seq) if a != b)

            # Count sequences
            count = len(seq_records)

            # Show sequence with markers for differences
            f.write(f"Variant {idx}:  ")
            for ref_aa, var_aa in zip(peptide_seq, variant_seq):
                if ref_aa == var_aa:
                    f.write(f"  {var_aa}")
                else:
                    f.write(f" *{var_aa}")  # * marks differences

            f.write(f"   ({count} seqs, {mismatches} mm)\n")

            # Show which positions differ
            if mismatches > 0:
                diff_positions = [(i+1, peptide_seq[i], variant_seq[i])
                                 for i in range(len(peptide_seq))
                                 if peptide_seq[i] != variant_seq[i]]

                f.write("           ")
                for pos, ref, var in diff_positions:
                    f.write(f"Pos{pos}:{ref}→{var}  ")
                f.write("\n")

            # Show species distribution
            species_list = [extract_species_from_header(rec['id']) for rec in seq_records]
            species_counts = defaultdict(int)
            for sp in species_list:
                species_counts[sp] += 1

            f.write("           Species: ")
            species_str = ", ".join([f"{sp}({count})" for sp, count in
                                    sorted(species_counts.items(),
                                          key=lambda x: x[1],
                                          reverse=True)[:5]])
            f.write(f"{species_str}\n")
            f.write("\n")

        # Summary
        f.write("           " + "---" * len(peptide_seq) + "\n")
        f.write(f"\nSummary:\n")
        f.write(f"  Total unique variants: {len(variants_dict)}\n")
        f.write(f"  Total sequences with peptide: {sum(len(v) for v in variants_dict.values())}\n")
        f.write(f"  Conservation: ")

        # Calculate per-position conservation
        pos_conservation = []
        for pos in range(len(peptide_seq)):
            ref_aa = peptide_seq[pos]
            total = sum(len(v) for v in variants_dict.values())
            ref_count = sum(len(v) for var_seq, v in variants_dict.items()
                          if var_seq[pos] == ref_aa)
            cons = (ref_count / total * 100) if total > 0 else 0
            pos_conservation.append(cons)

        avg_cons = sum(pos_conservation) / len(pos_conservation) if pos_conservation else 0
        f.write(f"{avg_cons:.1f}% average\n")

        # Most variable positions
        variable_positions = [(i+1, pos_conservation[i])
                             for i in range(len(pos_conservation))
                             if pos_conservation[i] < 90]

        if variable_positions:
            f.write(f"  Variable positions (<90% conserved):\n")
            for pos, cons in variable_positions:
                f.write(f"    Position {pos}: {cons:.1f}%\n")

        f.write("\n")


def generate_all_alignments(fasta_file, output_dir="peptide_alignments"):
    """
    Generate alignments for all peptides.
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    print(f"Loading sequences from {fasta_file}...")
    sequences = list(SeqIO.parse(fasta_file, 'fasta'))
    print(f"Loaded {len(sequences)} sequences\n")

    # Analyze each peptide
    for peptide_name, peptide_seq in PEPTIDES.items():
        print(f"Processing {peptide_name}...")

        # Find all variants
        variants = defaultdict(list)

        for seq_record in sequences:
            seq_id = seq_record.id
            sequence = str(seq_record.seq).upper()

            matches = find_peptide_in_sequence(peptide_seq, sequence, allow_mismatches=3)

            if matches:
                # Use best match
                best_match = min(matches, key=lambda x: x[1])
                position, mismatches, matched_seq = best_match

                variants[matched_seq].append({
                    'id': seq_id,
                    'description': seq_record.description,
                    'position': position,
                    'mismatches': mismatches
                })

        # Generate output files
        base_name = f"{output_dir}/{peptide_name}"

        # FASTA format
        create_aligned_fasta_with_reference(
            peptide_name,
            peptide_seq,
            variants,
            f"{base_name}.fasta"
        )

        # Text alignment
        create_simple_alignment(
            peptide_name,
            peptide_seq,
            variants,
            f"{base_name}_alignment.txt"
        )

        print(f"  Found {len(variants)} unique variants in {sum(len(v) for v in variants.values())} sequences")
        print(f"  Written to: {base_name}.fasta and {base_name}_alignment.txt\n")

    print(f"\nAll alignments written to directory: {output_dir}/")

    # Create index file
    with open(f"{output_dir}/INDEX.txt", 'w') as f:
        f.write("PEPTIDE ALIGNMENT FILES INDEX\n")
        f.write("="*80 + "\n\n")
        f.write("For each peptide, two files are generated:\n\n")

        for peptide_name in PEPTIDES.keys():
            f.write(f"{peptide_name}:\n")
            f.write(f"  - {peptide_name}.fasta (FASTA format with variants)\n")
            f.write(f"  - {peptide_name}_alignment.txt (visual alignment)\n\n")

    print(f"Index written to: {output_dir}/INDEX.txt")


def main():
    import sys

    if len(sys.argv) < 2:
        print("Usage: python generate_peptide_alignments.py <fasta_file> [output_dir]")
        sys.exit(1)

    fasta_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "peptide_alignments"

    generate_all_alignments(fasta_file, output_dir)


if __name__ == "__main__":
    main()
