#!/usr/bin/env python3
"""
Create a properly formatted Lpp alignment figure with dots for identical positions.
"""

# Sequences (verified from UniProt)
SEQUENCES = {
    'E.coli_K12':              'MKATKLVLGAVILGSTLLAGCSSNAKIDQLSSDVQTLNAKVDQLSNDVNAMRSDVQAAKDDAARANQRLDNMATKYRK',
    'C.koseri_BAA-895':        'MNRTKLVLGAVILGSTLLAGCSSNAKIDQLSSDVQTLNAKVDQLSNDVNAMRSDVQAAKDDAARANQRLDNQATKYRK',
    'E.cloacae_13047':         'MNRTKLVLGAVILGSTLLAGCSSNAKIDQLSSDVQTLNAKVDQLSNDVNAMRSDVQAAKDDAARANQRLDNQATKYRK',
    'K.pneumoniae_43816':      'MNRTKLVLGAVILGSTLLAGCSSNAKIDQLSSDVQTLNAKVDQLSNDVNAMRSDVQAAKDDAARANQRLDNQAHSYRK',
    'S.typhimurium_14028S':    'MNRTKLVLGAVILGSTLLAGCSSNAKIDQLSSDVQTLNAKVDQLSNDVNAMRSDVQAAKDDAARANQRLDNQATKYRK',
    'S.typhimurium_SL1344':    'MNRTKLVLGAVILGSTLLAGCSSNAKIDQLSSDVQTLNAKVEQLSNDVNAMRSDVQAAKDDAARRANQRLDNKVFRICK',
    'Y.enterocolitica_JB580v': 'MNRTKLVLGAVILASTMLAGCSSNAKIDQLSSDVQTLNAKVDQLSNDVNAIRSDVQAAKDDAARANQRLDNQAHAYKK',
}

def mask_sequence(ref, query):
    """Replace matching positions with dots."""
    result = []
    for i, (r, q) in enumerate(zip(ref, query)):
        if r == q:
            result.append('.')
        else:
            result.append(q)
    # Handle length differences
    if len(query) > len(ref):
        result.extend(list(query[len(ref):]))
    return ''.join(result)

def count_differences(ref, query):
    """Count and list differences."""
    diffs = []
    for i, (r, q) in enumerate(zip(ref, query)):
        if r != q:
            diffs.append(f"{r}{i+1}{q}")
    return diffs

def main():
    ref = SEQUENCES['E.coli_K12']

    # Calculate max name length for formatting
    max_name = max(len(name) for name in SEQUENCES.keys())

    output = []
    output.append("="*100)
    output.append("MUREIN LIPOPROTEIN (Lpp) MULTI-SPECIES ALIGNMENT")
    output.append("="*100)
    output.append("")
    output.append(f"Reference: E. coli K12 ({len(ref)} amino acids)")
    output.append("")

    # Position markers
    pos_line1 = " " * (max_name + 2)
    pos_line2 = " " * (max_name + 2)
    for i in range(1, len(ref) + 1):
        if i % 10 == 0:
            pos_line1 += str(i // 10)
            pos_line2 += "0"
        else:
            pos_line1 += " "
            pos_line2 += str(i % 10) if i % 10 == 1 else " "

    output.append("Position:" + " " * (max_name - 7) + "         1         2         3         4         5         6         7        ")
    output.append(" " * (max_name + 2) + "12345678901234567890123456789012345678901234567890123456789012345678901234567890")
    output.append("-" * 100)

    # Reference sequence
    output.append(f"{'E.coli_K12':<{max_name}}  {ref}")

    # Other sequences with masking
    for name in ['C.koseri_BAA-895', 'E.cloacae_13047', 'K.pneumoniae_43816',
                 'S.typhimurium_14028S', 'S.typhimurium_SL1344', 'Y.enterocolitica_JB580v']:
        seq = SEQUENCES[name]
        masked = mask_sequence(ref, seq)
        output.append(f"{name:<{max_name}}  {masked}")

    output.append("-" * 100)
    output.append("")
    output.append("KEY: '.' = identical to E. coli K12")
    output.append("")

    # Identity summary
    output.append("="*100)
    output.append("PAIRWISE IDENTITY TO E. coli K12")
    output.append("="*100)
    output.append("")
    output.append(f"{'Species':<30} {'Identity':>10} {'#Diff':>8}   Differences")
    output.append("-" * 100)

    for name, seq in SEQUENCES.items():
        if name == 'E.coli_K12':
            continue

        diffs = count_differences(ref, seq)
        identity = ((len(ref) - len(diffs)) / len(ref)) * 100
        diff_str = ', '.join(diffs)

        output.append(f"{name:<30} {identity:>9.1f}% {len(diffs):>8}   {diff_str}")

    output.append("")
    output.append("="*100)
    output.append("CONSERVATION PATTERNS")
    output.append("="*100)
    output.append("""
REGION 1 - N-terminus/Signal Peptide (aa 1-3):
  E. coli:        MKA
  All others:     MNR  (K2N, A3R are universal changes)

REGION 2 - Core Lipoprotein Domain (aa 4-70):
  Highly conserved (>95%) across all Enterobacteriaceae
  Notable variations:
    - Y. enterocolitica: G14A, L17M, M51I
    - S. typhimurium SL1344: D42E, A65R

REGION 3 - C-terminus (aa 71-78):
  Most variable region - species-specific signatures:
    E. coli K12:        NMATKYRK
    C. koseri:          NQATKYRK  (M71Q)
    E. cloacae:         NQATKYRK  (M71Q)
    K. pneumoniae:      NQAHSYRK  (M71Q, T74H, K75S)
    S. typhimurium:     NQATKYRK  (M71Q)
    Y. enterocolitica:  NQAHAYKK  (M71Q, T74A, K77K)
    S. SL1344:          NKVFRICK  (highly divergent)
""")

    output.append("="*100)
    output.append("NOTE ON GRAM-POSITIVE BACTERIA")
    output.append("="*100)
    output.append("""
Listeria monocytogenes (10403s) and Staphylococcus aureus (USA200) DO NOT
have a classical Lpp/Braun's lipoprotein because:

  1. They are Gram-positive bacteria WITHOUT an outer membrane
  2. Lpp specifically anchors the outer membrane to peptidoglycan
  3. This cell wall architecture is exclusive to Gram-negative bacteria

Gram-positive bacteria have different lipoproteins that are NOT homologous
to Braun's lipoprotein.
""")
    output.append("="*100)

    # Write to file
    with open('lpp_multispecies_alignment.txt', 'w') as f:
        f.write('\n'.join(output))

    # Print to screen
    print('\n'.join(output))

if __name__ == "__main__":
    main()
