#!/usr/bin/env python3
"""
Analyze variability of specific peptide sequences within OmpA proteins
across E. coli and Enterobacteriaceae strains.

This script:
1. Loads peptide sequences of interest
2. Maps them to OmpA protein sequences
3. Analyzes conservation/variability at those positions
4. Reports substitutions and their frequencies
"""

from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from collections import defaultdict, Counter
import argparse
import sys

# Define peptide sequences
PEPTIDES = {
    'Peptide_1': 'WSQYHDTGFINNNGPTHEN',
    'Peptide_2': 'GRMPYKGSVENGAYKA',
    'Peptide_3a': 'KSNVYGKN',
    'Peptide_3b': 'KANVPGGASFKD',
    'Peptide_4': 'TNNIGDAHTIGTRPDNGM',

    'Peptide_1S': 'TQPNSHNDINTNHGYGWEF',
    'Peptide_2S': 'RSMAGKGPKNAGYVYE',
    'Peptide_3aS': 'KYNGVNSK',
    'Peptide_3bS': 'FKAKDNGVGSPA',
    'Peptide_4S': 'MTNTHNTAIGDRNIGPGD',

    'Peptide_3BS2': 'GVSFDKAKNAPG',
    'Peptide_2b': 'GRMPYKGDNINGAYKA',
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
    sequence = str(sequence).upper().replace('-', '')  # Remove gaps for searching

    pep_len = len(peptide)
    seq_len = len(sequence)

    for i in range(seq_len - pep_len + 1):
        subseq = sequence[i:i + pep_len]
        mismatches = sum(1 for a, b in zip(peptide, subseq) if a != b)

        if mismatches <= allow_mismatches:
            matches.append((i, mismatches, subseq))

    return matches


def analyze_aligned_sequences(alignment_file, peptides, output_prefix):
    """
    Analyze aligned sequences for peptide variability.
    """
    print(f"Loading alignment from {alignment_file}...")

    # Try different formats
    alignment = None
    for fmt in ['fasta', 'clustal']:
        try:
            alignment = AlignIO.read(alignment_file, fmt)
            print(f"Successfully loaded alignment in {fmt} format")
            break
        except:
            continue

    if alignment is None:
        # Try loading as unaligned FASTA
        print("Loading as unaligned FASTA...")
        sequences = list(SeqIO.parse(alignment_file, 'fasta'))
        if not sequences:
            print(f"Error: Could not load sequences from {alignment_file}")
            return
        print(f"Loaded {len(sequences)} unaligned sequences")
    else:
        sequences = list(alignment)
        print(f"Loaded {len(sequences)} aligned sequences")

    # Analyze each peptide
    results = {}

    for peptide_name, peptide_seq in peptides.items():
        print(f"\n{'='*80}")
        print(f"Analyzing {peptide_name}: {peptide_seq}")
        print(f"{'='*80}")

        peptide_results = {
            'peptide': peptide_seq,
            'length': len(peptide_seq),
            'found_in': [],
            'not_found_in': [],
            'variants': defaultdict(list),
            'position_variability': defaultdict(Counter),
        }

        # Search for peptide in each sequence
        # Allow up to 3 mismatches
        for seq_record in sequences:
            seq_id = seq_record.id
            sequence = str(seq_record.seq).upper()

            matches = find_peptide_in_sequence(peptide_seq, sequence, allow_mismatches=3)

            if matches:
                # Use best match (fewest mismatches)
                best_match = min(matches, key=lambda x: x[1])
                position, mismatches, matched_seq = best_match

                peptide_results['found_in'].append({
                    'id': seq_id,
                    'position': position,
                    'mismatches': mismatches,
                    'sequence': matched_seq
                })

                peptide_results['variants'][matched_seq].append(seq_id)

                # Track position-specific variability
                for i, (ref_aa, obs_aa) in enumerate(zip(peptide_seq, matched_seq)):
                    peptide_results['position_variability'][i][obs_aa] += 1
            else:
                peptide_results['not_found_in'].append(seq_id)

        results[peptide_name] = peptide_results

    # Generate report
    report_file = f"{output_prefix}_variability_report.txt"
    with open(report_file, 'w') as f:
        f.write("OmpA PEPTIDE VARIABILITY ANALYSIS\n")
        f.write("="*80 + "\n\n")
        f.write(f"Total sequences analyzed: {len(sequences)}\n\n")

        for peptide_name, data in results.items():
            f.write("\n" + "="*80 + "\n")
            f.write(f"PEPTIDE: {peptide_name}\n")
            f.write(f"Sequence: {data['peptide']}\n")
            f.write(f"Length: {data['length']} amino acids\n")
            f.write("="*80 + "\n\n")

            found_count = len(data['found_in'])
            not_found_count = len(data['not_found_in'])

            f.write(f"Found in: {found_count}/{len(sequences)} sequences ({100*found_count/len(sequences):.1f}%)\n")
            f.write(f"Not found in: {not_found_count}/{len(sequences)} sequences\n\n")

            if data['not_found_in']:
                f.write("Sequences where peptide was NOT found:\n")
                for seq_id in data['not_found_in']:
                    f.write(f"  - {seq_id}\n")
                f.write("\n")

            # Report variants
            f.write(f"\nUnique variants found: {len(data['variants'])}\n")
            f.write("-"*80 + "\n")

            for variant_seq, seq_ids in sorted(data['variants'].items(),
                                               key=lambda x: len(x[1]),
                                               reverse=True):
                # Calculate mismatches
                mismatches = sum(1 for a, b in zip(data['peptide'], variant_seq) if a != b)
                percent = 100 * len(seq_ids) / found_count

                f.write(f"\nVariant: {variant_seq}\n")
                f.write(f"  Mismatches: {mismatches}\n")
                f.write(f"  Frequency: {len(seq_ids)}/{found_count} ({percent:.1f}%)\n")
                f.write(f"  Found in:\n")
                for seq_id in seq_ids:
                    f.write(f"    - {seq_id}\n")

            # Position-specific variability
            f.write("\n" + "-"*80 + "\n")
            f.write("Position-specific amino acid variability:\n")
            f.write("-"*80 + "\n\n")

            for pos in range(data['length']):
                ref_aa = data['peptide'][pos]
                aa_counts = data['position_variability'][pos]

                # Calculate conservation
                total = sum(aa_counts.values())
                if total == 0:
                    continue

                conservation = max(aa_counts.values()) / total * 100

                f.write(f"Position {pos+1} (reference: {ref_aa}):\n")
                f.write(f"  Conservation: {conservation:.1f}%\n")

                # Show all amino acids found at this position
                for aa, count in sorted(aa_counts.items(), key=lambda x: x[1], reverse=True):
                    percent = count / total * 100
                    marker = " (REFERENCE)" if aa == ref_aa else ""
                    f.write(f"    {aa}: {count}/{total} ({percent:.1f}%){marker}\n")

                # Highlight if variable (< 90% conservation)
                if conservation < 90:
                    f.write(f"  *** VARIABLE POSITION ***\n")

                f.write("\n")

    print(f"\nReport written to: {report_file}")

    # Generate summary CSV
    csv_file = f"{output_prefix}_variability_summary.csv"
    with open(csv_file, 'w') as f:
        f.write("Peptide,Position,Reference_AA,Observed_AA,Count,Frequency_Percent,Conservation_Percent\n")

        for peptide_name, data in results.items():
            found_count = len(data['found_in'])
            if found_count == 0:
                continue

            for pos in range(data['length']):
                ref_aa = data['peptide'][pos]
                aa_counts = data['position_variability'][pos]
                total = sum(aa_counts.values())

                if total == 0:
                    continue

                conservation = max(aa_counts.values()) / total * 100

                for aa, count in sorted(aa_counts.items(), key=lambda x: x[1], reverse=True):
                    percent = count / total * 100
                    f.write(f"{peptide_name},{pos+1},{ref_aa},{aa},{count},{percent:.2f},{conservation:.2f}\n")

    print(f"Summary CSV written to: {csv_file}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Analyze peptide sequence variability in OmpA proteins'
    )
    parser.add_argument('alignment_file',
                        help='FASTA file with OmpA protein sequences (aligned or unaligned)')
    parser.add_argument('--output-prefix', default='peptide_analysis',
                        help='Prefix for output files')
    parser.add_argument('--custom-peptides',
                        help='File with custom peptides (one per line: name,sequence)')

    args = parser.parse_args()

    # Load custom peptides if provided
    peptides = PEPTIDES.copy()
    if args.custom_peptides:
        with open(args.custom_peptides, 'r') as f:
            for line in f:
                line = line.strip()
                if line and ',' in line:
                    name, seq = line.split(',', 1)
                    peptides[name.strip()] = seq.strip()

    # Run analysis
    results = analyze_aligned_sequences(args.alignment_file, peptides, args.output_prefix)


if __name__ == "__main__":
    main()
