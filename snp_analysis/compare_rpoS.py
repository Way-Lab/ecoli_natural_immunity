#!/usr/bin/env python3
"""
Extract and compare rpoS gene sequences between two E. coli strains.
Includes upstream regulatory region for promoter analysis.
"""

from Bio import SeqIO
from Bio.Seq import Seq
import argparse


def parse_annotation_for_rpoS(tsv_file):
    """Parse annotation TSV file to find rpoS gene location(s)."""
    rpoS_entries = []

    with open(tsv_file, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue
            if 'rpos' in line.lower():
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    contig = parts[0]
                    start = int(parts[2])
                    end = int(parts[3])
                    strand = parts[4]
                    locus_tag = parts[5] if len(parts) > 5 else "unknown"
                    product = parts[7] if len(parts) > 7 else "unknown"
                    rpoS_entries.append({
                        'contig': contig,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'locus_tag': locus_tag,
                        'product': product
                    })

    return rpoS_entries


def extract_sequence_with_upstream(fna_file, contig_name, start, end, strand, upstream_bp=500):
    """Extract gene sequence with upstream regulatory region."""
    # Parse the genome file
    genome = {}
    for record in SeqIO.parse(fna_file, "fasta"):
        genome[record.id] = record.seq

    if contig_name not in genome:
        raise ValueError(f"Contig {contig_name} not found in genome file")

    contig_seq = genome[contig_name]

    # Calculate extraction coordinates with upstream region
    if strand == '+':
        extract_start = max(0, start - 1 - upstream_bp)  # -1 for 0-based indexing
        extract_end = end
        extracted_seq = contig_seq[extract_start:extract_end]
    else:  # negative strand
        extract_start = start - 1
        extract_end = min(len(contig_seq), end + upstream_bp)
        extracted_seq = contig_seq[extract_start:extract_end].reverse_complement()

    return extracted_seq, extract_start + 1, extract_end  # +1 to return 1-based coords


def compare_sequences(seq1, seq2, name1, name2):
    """Compare two sequences and identify differences."""
    print(f"\n{'='*80}")
    print(f"SEQUENCE COMPARISON: {name1} vs {name2}")
    print(f"{'='*80}\n")

    print(f"{name1} length: {len(seq1)} bp")
    print(f"{name2} length: {len(seq2)} bp")
    print(f"Length difference: {abs(len(seq1) - len(seq2))} bp\n")

    # Check if sequences are identical
    if seq1 == seq2:
        print("✓ Sequences are IDENTICAL\n")
        return

    print("✗ Sequences are DIFFERENT\n")

    # Find differences
    min_len = min(len(seq1), len(seq2))
    differences = []

    for i in range(min_len):
        if seq1[i] != seq2[i]:
            differences.append(i)

    if differences:
        print(f"Found {len(differences)} nucleotide differences in aligned region:")
        print(f"First 10 differences:")
        for i, pos in enumerate(differences[:10]):
            print(f"  Position {pos+1}: {name1}={seq1[pos]} vs {name2}={seq2[pos]}")
        if len(differences) > 10:
            print(f"  ... and {len(differences)-10} more differences")

    # Check for insertions/deletions
    if len(seq1) != len(seq2):
        print(f"\nLength difference suggests insertion/deletion event")
        if len(seq1) > len(seq2):
            print(f"  {name1} has {len(seq1) - len(seq2)} extra bp")
        else:
            print(f"  {name2} has {len(seq2) - len(seq1)} extra bp")

    print()


def main():
    parser = argparse.ArgumentParser(
        description='Compare rpoS gene sequences between two strains'
    )
    parser.add_argument('--strain1-name', default='N2SKTQ_8_8',
                        help='Name of first strain')
    parser.add_argument('--strain2-name', default='N2SKTQ_11_11',
                        help='Name of second strain')
    parser.add_argument('--upstream', type=int, default=500,
                        help='Number of base pairs upstream to include (default: 500)')
    parser.add_argument('--output-fasta', default='rpoS_comparison.fasta',
                        help='Output FASTA file with extracted sequences')

    args = parser.parse_args()

    # File paths
    base_path = "N2SKTQ_results"

    strain1_tsv = f"{base_path}/{args.strain1_name}/ONT-only/annotation/{args.strain1_name}.tsv"
    strain1_fna = f"{base_path}/{args.strain1_name}/ONT-only/annotation/{args.strain1_name}.fna"

    strain2_tsv = f"{base_path}/{args.strain2_name}/ONT-only/annotation/{args.strain2_name}.tsv"
    strain2_fna = f"{base_path}/{args.strain2_name}/ONT-only/annotation/{args.strain2_name}.fna"

    print(f"Analyzing rpoS gene in {args.strain1_name} and {args.strain2_name}")
    print(f"Including {args.upstream} bp upstream for regulatory region analysis\n")

    # Parse annotations
    print(f"Parsing {args.strain1_name} annotations...")
    strain1_rpoS = parse_annotation_for_rpoS(strain1_tsv)
    print(f"Found {len(strain1_rpoS)} rpoS entry(ies):")
    for entry in strain1_rpoS:
        print(f"  {entry['contig']}: {entry['start']}-{entry['end']} ({entry['strand']}) - {entry['product']}")

    print(f"\nParsing {args.strain2_name} annotations...")
    strain2_rpoS = parse_annotation_for_rpoS(strain2_tsv)
    print(f"Found {len(strain2_rpoS)} rpoS entry(ies):")
    for entry in strain2_rpoS:
        print(f"  {entry['contig']}: {entry['start']}-{entry['end']} ({entry['strand']}) - {entry['product']}")

    # Extract sequences
    print(f"\n{'='*80}")
    print("EXTRACTING SEQUENCES")
    print(f"{'='*80}\n")

    # For strain 1, if multiple entries, combine them
    if len(strain1_rpoS) > 1:
        print(f"WARNING: {args.strain1_name} has {len(strain1_rpoS)} rpoS fragments!")
        print("This suggests a frameshift or split gene. Extracting region covering all fragments...")
        strain1_start = min(e['start'] for e in strain1_rpoS)
        strain1_end = max(e['end'] for e in strain1_rpoS)
        strain1_contig = strain1_rpoS[0]['contig']
        strain1_strand = strain1_rpoS[0]['strand']
    else:
        strain1_start = strain1_rpoS[0]['start']
        strain1_end = strain1_rpoS[0]['end']
        strain1_contig = strain1_rpoS[0]['contig']
        strain1_strand = strain1_rpoS[0]['strand']

    seq1, s1_start, s1_end = extract_sequence_with_upstream(
        strain1_fna, strain1_contig, strain1_start, strain1_end,
        strain1_strand, args.upstream
    )
    print(f"{args.strain1_name}: Extracted {len(seq1)} bp from {strain1_contig}:{s1_start}-{s1_end}")

    # For strain 2
    strain2_start = strain2_rpoS[0]['start']
    strain2_end = strain2_rpoS[0]['end']
    strain2_contig = strain2_rpoS[0]['contig']
    strain2_strand = strain2_rpoS[0]['strand']

    seq2, s2_start, s2_end = extract_sequence_with_upstream(
        strain2_fna, strain2_contig, strain2_start, strain2_end,
        strain2_strand, args.upstream
    )
    print(f"{args.strain2_name}: Extracted {len(seq2)} bp from {strain2_contig}:{s2_start}-{s2_end}")

    # Compare sequences
    compare_sequences(seq1, seq2, args.strain1_name, args.strain2_name)

    # Write to FASTA
    with open(args.output_fasta, 'w') as f:
        f.write(f">{args.strain1_name}_rpoS_{strain1_contig}:{s1_start}-{s1_end}_strand{strain1_strand}\n")
        # Write sequence in 80 bp lines
        for i in range(0, len(seq1), 80):
            f.write(str(seq1[i:i+80]) + "\n")

        f.write(f">{args.strain2_name}_rpoS_{strain2_contig}:{s2_start}-{s2_end}_strand{strain2_strand}\n")
        for i in range(0, len(seq2), 80):
            f.write(str(seq2[i:i+80]) + "\n")

    print(f"Sequences written to {args.output_fasta}")
    print("\nTo view alignment, you can use tools like:")
    print(f"  mafft {args.output_fasta} > rpoS_aligned.fasta")
    print(f"  Or use NCBI BLAST to compare: https://blast.ncbi.nlm.nih.gov/")


if __name__ == "__main__":
    main()
