#!/usr/bin/env python3
"""
Analyze the gap region between the two rpoS fragments in N2SKTQ_8_8
"""

from Bio import SeqIO

# Read the genome
fna_file = "N2SKTQ_results/N2SKTQ_8_8/ONT-only/annotation/N2SKTQ_8_8.fna"

for record in SeqIO.parse(fna_file, "fasta"):
    if record.id == "contig_3":
        # Extract the gap region and surrounding sequence
        gap_start = 279169  # End of first fragment
        gap_end = 279250    # Start of second fragment

        # Get sequence with some context
        context = 50
        region_start = gap_start - context
        region_end = gap_end + context

        sequence = record.seq[region_start-1:region_end]  # -1 for 0-based indexing

        print(f"Region around the gap in N2SKTQ_8_8 (contig_3:{region_start}-{region_end})")
        print(f"Fragment 1 ends at: {gap_start}")
        print(f"Fragment 2 starts at: {gap_end}")
        print(f"Gap size: {gap_end - gap_start - 1} bp")
        print(f"\nSequence (with 50 bp context on each side):")
        print(f"Position {region_start}:")

        # Print with position markers
        gap_relative_start = context
        gap_relative_end = gap_relative_start + (gap_end - gap_start - 1)

        seq_str = str(sequence)
        for i in range(0, len(seq_str), 80):
            chunk = seq_str[i:i+80]
            print(f"{region_start+i:6d}  {chunk}")

        print(f"\nThe gap region ({gap_start+1}-{gap_end-1}):")
        gap_seq = record.seq[gap_start:gap_end-1]  # The actual gap
        print(str(gap_seq))
        print(f"Length: {len(gap_seq)} bp")

        # Check for stop codons in all three frames
        print("\nTranslation of gap region in all frames:")
        for frame in range(3):
            if frame < len(gap_seq):
                protein = gap_seq[frame:].translate()
                print(f"Frame {frame}: {protein}")
                if '*' in str(protein):
                    print(f"  -> Contains STOP codon(s)!")

        break
