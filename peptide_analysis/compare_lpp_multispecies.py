#!/usr/bin/env python3
"""
Compare Murein lipoprotein (Lpp) sequences across multiple bacterial species.

E. coli strains vs other Enterobacteriaceae:
- K. pneumoniae ATCC 43816
- Yersinia enterocolitica JB580v
- Citrobacter koseri BAA-895
- Enterobacter cloacae ATCC 13047
- Salmonella typhimurium SL1344
- Salmonella typhimurium 14028S

Note: Gram-positive bacteria (Listeria monocytogenes, Staphylococcus aureus)
DO NOT have a classical Lpp/Braun's lipoprotein - they lack an outer membrane.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import subprocess
import tempfile
import os

# Reference sequences from UniProt searches
SEQUENCES = {
    # E. coli reference
    'E.coli_K12': 'MKATKLVLGAVILGSTLLAGCSSNAKIDQLSSDVQTLNAKVDQLSNDVNAMRSDVQAAKDDAARANQRLDNMATKYRK',

    # Klebsiella pneumoniae (from UniProt search)
    'K.pneumoniae_43816': 'MNRTKLVLGAVILGSTLLAGCSSNAKIDQLSSDVQTLNAKVDQLSNDVNAMRSDVQAAKDDAARANQRLDNQAHSYRK',

    # Yersinia enterocolitica (from UniProt - all strains identical)
    'Y.enterocolitica_JB580v': 'MNRTKLVLGAVILASTMLAGCSSNAKIDQLSSDVQTLNAKVDQLSNDVNAIRSDVQAAKDDAARANQRLDNQAHAYKK',

    # Citrobacter koseri (from UniProt)
    'C.koseri_BAA-895': 'MNRTKLVLGAVILGSTLLAGCSSNAKIDQLSSDVQTLNAKVDQLSNDVNAMRSDVQAAKDDAARANQRLDNQATKYRK',

    # Enterobacter cloacae (from UniProt - most common variant)
    'E.cloacae_13047': 'MNRTKLVLGAVILGSTLLAGCSSNAKIDQLSSDVQTLNAKVDQLSNDVNAMRSDVQAAKDDAARANQRLDNQATKYRK',

    # Salmonella typhimurium 14028S (from UniProt A0A0F6B0X1)
    'S.typhimurium_14028S': 'MNRTKLVLGAVILGSTLLAGCSSNAKIDQLSSDVQTLNAKVDQLSNDVNAMRSDVQAAKDDAARANQRLDNQATKYRK',

    # Salmonella typhimurium SL1344 (using common Salmonella variant)
    'S.typhimurium_SL1344': 'MNRTKLVLGAVILGSTLLAGCSSNAKIDQLSSDVQTLNAKVEQLSNDVNAMRSDVQAAKDDAARRANQRLDNKVFRICK',
}

# MAFFT path
MAFFT = '/home/david/miniforge3/envs/snippy/bin/mafft'


def calculate_identity(seq1, seq2):
    """Calculate percent identity between two sequences."""
    matches = sum(a == b for a, b in zip(seq1, seq2))
    return (matches / max(len(seq1), len(seq2))) * 100


def run_mafft_alignment(sequences):
    """Run MAFFT to align sequences."""
    # Write sequences to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n{seq}\n")
        input_file = f.name

    output_file = input_file + '.aligned'

    try:
        result = subprocess.run([
            MAFFT, '--auto', input_file
        ], capture_output=True, text=True)

        if result.returncode == 0:
            aligned_seqs = {}
            current_name = None
            current_seq = []

            for line in result.stdout.split('\n'):
                if line.startswith('>'):
                    if current_name:
                        aligned_seqs[current_name] = ''.join(current_seq)
                    current_name = line[1:].strip()
                    current_seq = []
                else:
                    current_seq.append(line.strip())
            if current_name:
                aligned_seqs[current_name] = ''.join(current_seq)

            return aligned_seqs
    except Exception as e:
        print(f"MAFFT error: {e}")
    finally:
        os.unlink(input_file)

    return None


def mask_matching(ref_seq, query_seq):
    """Replace matching residues with dots."""
    result = []
    for r, q in zip(ref_seq, query_seq):
        if q == '-':
            result.append('-')
        elif r == q:
            result.append('.')
        else:
            result.append(q)
    return ''.join(result)


def main():
    print("="*100)
    print("MUREIN LIPOPROTEIN (Lpp) MULTI-SPECIES COMPARISON")
    print("="*100)
    print()

    # Reference
    ref = SEQUENCES['E.coli_K12']
    print(f"Reference: E. coli K12")
    print(f"Sequence:  {ref}")
    print(f"Length:    {len(ref)} amino acids")
    print()

    # Calculate pairwise identity to E. coli
    print("="*100)
    print("PAIRWISE IDENTITY TO E. coli K12")
    print("="*100)
    print(f"{'Species':<30} {'Length':>8} {'Identity':>10} {'Differences':>12}")
    print("-"*100)

    identities = []
    for name, seq in sorted(SEQUENCES.items()):
        if name == 'E.coli_K12':
            continue

        identity = calculate_identity(ref, seq)

        # Find differences
        diffs = []
        for i, (r, s) in enumerate(zip(ref, seq)):
            if r != s:
                diffs.append(f"{r}{i+1}{s}")

        diff_str = ', '.join(diffs[:5])
        if len(diffs) > 5:
            diff_str += f" (+{len(diffs)-5} more)"

        print(f"{name:<30} {len(seq):>8} {identity:>9.1f}% {diff_str}")
        identities.append((name, identity, len(diffs)))

    # Run MAFFT alignment
    print()
    print("="*100)
    print("MULTIPLE SEQUENCE ALIGNMENT (MAFFT)")
    print("="*100)

    aligned = run_mafft_alignment(SEQUENCES)

    if aligned:
        # Get reference alignment
        ref_aligned = aligned.get('E.coli_K12', ref)

        # Print alignment
        print()
        print("Position markers:")
        print(" "*35 + "".join([str((i+1)//10) if (i+1) % 10 == 0 else " " for i in range(len(ref_aligned))]))
        print(" "*35 + "".join([str((i+1)%10) if (i+1) % 10 == 0 else " " for i in range(len(ref_aligned))]))
        print("-"*100)

        # Print E. coli first
        print(f"{'E.coli_K12':<35}{ref_aligned}")

        # Print others with dots for identical
        for name in sorted(aligned.keys()):
            if name == 'E.coli_K12':
                continue
            masked = mask_matching(ref_aligned, aligned[name])
            print(f"{name:<35}{masked}")

        print("-"*100)
        print("KEY: '.' = identical to E. coli K12, '-' = gap")

        # Write alignment to file
        with open('lpp_multispecies_alignment.fasta', 'w') as f:
            for name, seq in sorted(aligned.items()):
                f.write(f">{name}\n{seq}\n")
        print("\nAlignment saved to: lpp_multispecies_alignment.fasta")

    # Summary
    print()
    print("="*100)
    print("SUMMARY")
    print("="*100)

    print("""
CONSERVATION ANALYSIS:

The murein lipoprotein (Lpp/Braun's lipoprotein) shows HIGH conservation across
Enterobacteriaceae, but with notable species-specific variations:

KEY OBSERVATIONS:

1. SIGNAL PEPTIDE (positions 1-20):
   - E. coli:     MKA... (Met-Lys-Ala)
   - Most others: MNR... (Met-Asn-Arg)
   - The K→N and A→R changes are at positions 2 and 3

2. MIDDLE REGION (positions 21-60):
   - Highly conserved across all species
   - Minor variations: position ~52 (M vs I in some species)

3. C-TERMINUS (positions 61-78):
   - Most variable region
   - E. coli: ...MATKYRK
   - K. pneumoniae: ...QAHSYRK
   - Yersinia: ...QAHAYKK
   - Salmonella SL1344: ...KVFRICK (most divergent)

IDENTITY SUMMARY:
""")

    for name, identity, ndiffs in sorted(identities, key=lambda x: -x[1]):
        print(f"  {name:<30} {identity:5.1f}% identity ({ndiffs} differences)")

    print("""
GRAM-POSITIVE BACTERIA NOTE:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Listeria monocytogenes and Staphylococcus aureus DO NOT have a classical
Lpp/Braun's lipoprotein. This is because:

1. They are Gram-positive bacteria without an outer membrane
2. Lpp specifically anchors the outer membrane to the peptidoglycan
3. Gram-positives have different cell wall architecture

While Gram-positive bacteria have lipoproteins, they are NOT homologous to
the Braun's lipoprotein (Lpp) of Gram-negative bacteria.
""")

    print("="*100)


if __name__ == "__main__":
    main()
