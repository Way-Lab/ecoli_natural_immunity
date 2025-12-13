#!/usr/bin/env python3
"""
Analyze which species/strains have which peptide variants.
This helps identify species-specific or pathotype-specific patterns.
"""

from Bio import SeqIO
from collections import defaultdict
import re

# Define peptide sequences
PEPTIDES = {
    'Peptide_1': 'WSQYHDTGFINNNGPTHEN',
    'Peptide_2': 'GRMPYKGSVENGAYKA',
    'Peptide_3a': 'KSNVYGKN',
    'Peptide_3b': 'KANVPGGASFKD',
    'Peptide_4': 'TNNIGDAHTIGTRPDNGM',
    'Peptide_2b': 'GRMPYKGDNINGAYKA',
}


def extract_species_from_header(header):
    """
    Extract species/strain information from FASTA header.
    """
    # Try to identify organism
    if 'coli' in header.lower():
        return 'Escherichia_coli'
    elif 'klebsiella' in header.lower():
        if 'pneumoniae' in header.lower():
            return 'Klebsiella_pneumoniae'
        elif 'oxytoca' in header.lower():
            return 'Klebsiella_oxytoca'
        else:
            return 'Klebsiella_sp'
    elif 'enterobacter' in header.lower():
        if 'cloacae' in header.lower():
            return 'Enterobacter_cloacae'
        elif 'aerogenes' in header.lower():
            return 'Enterobacter_aerogenes'
        else:
            return 'Enterobacter_sp'
    elif 'citrobacter' in header.lower():
        if 'freundii' in header.lower():
            return 'Citrobacter_freundii'
        elif 'koseri' in header.lower():
            return 'Citrobacter_koseri'
        else:
            return 'Citrobacter_sp'
    elif 'salmonella' in header.lower():
        return 'Salmonella_enterica'
    elif 'shigella' in header.lower():
        if 'flexneri' in header.lower():
            return 'Shigella_flexneri'
        elif 'sonnei' in header.lower():
            return 'Shigella_sonnei'
        else:
            return 'Shigella_sp'
    elif 'proteus' in header.lower():
        return 'Proteus_mirabilis'
    elif 'serratia' in header.lower():
        return 'Serratia_marcescens'
    else:
        return 'Unknown'


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


def analyze_species_variants(fasta_file, output_file):
    """
    Analyze which species have which variants.
    """
    print(f"Loading sequences from {fasta_file}...")
    sequences = list(SeqIO.parse(fasta_file, 'fasta'))
    print(f"Loaded {len(sequences)} sequences\n")

    # Group sequences by species
    species_seqs = defaultdict(list)
    for seq in sequences:
        species = extract_species_from_header(seq.description)
        species_seqs[species].append(seq)

    print(f"Found {len(species_seqs)} species/groups:")
    for species, seqs in sorted(species_seqs.items(), key=lambda x: -len(x[1])):
        print(f"  {species}: {len(seqs)} sequences")

    # Analyze each peptide
    with open(output_file, 'w') as f:
        f.write("# Species-Specific Peptide Variant Analysis\n\n")
        f.write(f"Total sequences: {len(sequences)}\n")
        f.write(f"Total species/groups: {len(species_seqs)}\n\n")

        for peptide_name, peptide_seq in PEPTIDES.items():
            f.write(f"\n{'='*80}\n")
            f.write(f"## {peptide_name}: {peptide_seq}\n")
            f.write(f"{'='*80}\n\n")

            # Track variants by species
            species_variants = defaultdict(lambda: defaultdict(int))
            species_total = defaultdict(int)

            for species, seqs in species_seqs.items():
                for seq in seqs:
                    species_total[species] += 1
                    sequence = str(seq.seq).upper()
                    matches = find_peptide_in_sequence(peptide_seq, sequence, allow_mismatches=3)

                    if matches:
                        best_match = min(matches, key=lambda x: x[1])
                        variant = best_match[2]
                        species_variants[species][variant] += 1

            # Report results
            f.write("### Species Distribution\n\n")
            f.write("| Species | Total Seqs | Found | Percentage | Variants |\n")
            f.write("|---------|-----------|-------|------------|----------|\n")

            for species in sorted(species_total.keys()):
                total = species_total[species]
                found = sum(species_variants[species].values())
                percent = (found / total * 100) if total > 0 else 0
                num_variants = len(species_variants[species])

                f.write(f"| {species} | {total} | {found} | {percent:.1f}% | {num_variants} |\n")

            # Show variants for each species
            f.write("\n### Variants by Species\n\n")

            for species in sorted(species_variants.keys()):
                if not species_variants[species]:
                    continue

                f.write(f"\n#### {species}\n\n")

                for variant, count in sorted(species_variants[species].items(),
                                            key=lambda x: x[1], reverse=True):
                    mismatches = sum(1 for a, b in zip(peptide_seq, variant) if a != b)
                    total = species_total[species]
                    percent = count / total * 100

                    f.write(f"- **{variant}** ({count}/{total}, {percent:.1f}%)\n")
                    if mismatches > 0:
                        f.write(f"  - {mismatches} mismatch(es) from reference\n")
                        # Show specific differences
                        diffs = [(i+1, peptide_seq[i], variant[i])
                                for i in range(len(peptide_seq))
                                if peptide_seq[i] != variant[i]]
                        for pos, ref, obs in diffs:
                            f.write(f"    - Position {pos}: {ref} → {obs}\n")

            # Identify species-specific variants
            f.write("\n### Species-Specific Patterns\n\n")

            # Find variants unique to one species
            variant_species = defaultdict(set)
            for species, variants in species_variants.items():
                for variant in variants:
                    variant_species[variant].add(species)

            species_specific = {v: list(s) for v, s in variant_species.items() if len(s) == 1}

            if species_specific:
                f.write("Variants found in only one species:\n\n")
                for variant, species_list in sorted(species_specific.items()):
                    species = species_list[0]
                    count = species_variants[species][variant]
                    f.write(f"- **{variant}** - unique to {species} ({count} sequences)\n")
            else:
                f.write("No species-specific variants found (all variants present in multiple species)\n")

            f.write("\n")

    print(f"\nSpecies analysis written to: {output_file}")


def main():
    import sys

    if len(sys.argv) < 2:
        print("Usage: python analyze_species_variants.py <fasta_file> [output_file]")
        sys.exit(1)

    fasta_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'species_variant_analysis.md'

    analyze_species_variants(fasta_file, output_file)


if __name__ == "__main__":
    main()
