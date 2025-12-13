#!/usr/bin/env python3
"""
Analyze selected ompA sequences and compare peptide loop regions

This script:
1. Filters sequences to keep only specified strains
2. Performs multiple sequence alignment with MAFFT
3. Extracts peptide loop regions
4. Creates comprehensive comparison output with specified ordering

Strains included:
- E. coli: EcN, RS218, UTI189
- Other species: C.koseri_BAA-895, E.cloacae_ATCC13047, K.pneumoniae_ATCC43816,
                 S.Typhimurium_14028S, S.Typhimurium_SL1344, Y.enterocolitica_JB580v

Usage:
    python analyze_selected_ompA_strains.py
"""

import subprocess
import sys
from Bio import SeqIO
from Bio.Seq import Seq


# Define peptide sequences to look for (from previous analysis)
PEPTIDES = {
    'Peptide_1': 'WSQYHDTGFINNNGPTHEN',
    'Peptide_2': 'GRMPYKGSVENGAYKA',
    'Peptide_2b': 'GRMPYKGDNINGAYKA',
    'Peptide_3a': 'KSNVYGKN',
    'Peptide_3b': 'KANVPGGASFKD',
    'Peptide_4': 'TNNIGDAHTIGTRPDNGM',
}

# Define which strains to keep (in order)
STRAIN_ORDER = [
    'EcN',
    'RS218',
    'UTI189',
    'C.koseri_BAA-895',
    'E.cloacae_ATCC13047',
    'K.pneumoniae_ATCC43816',
    'S.Typhimurium_14028S',
    'S.Typhimurium_SL1344',
    'Y.enterocolitica_JB580v'
]


def get_strain_name(record):
    """Extract a clean strain name from a SeqRecord."""
    desc = record.description.lower()
    seq_id = record.id

    # Check E. coli strains FIRST (before other genera)
    # This prevents false matches like "Enterobacteriaceae" matching "enterobacter"
    if 'rs218' in desc or seq_id == 'RS218':
        return 'RS218'
    elif seq_id == 'CFT073' or 'cft073' in desc:
        return 'CFT073'
    elif seq_id == 'EcN' or 'nissle' in desc or 'ecn' in desc:
        return 'EcN'
    elif 'eco57' in desc or ('o157' in desc and 'p0a911' in seq_id.lower()):
        return 'O157:H7'
    elif 'sakai' in desc or ('o157' in desc and 'np_309068' in seq_id.lower()):
        return 'O157:H7_Sakai'
    elif 'ec958' in desc:
        return 'EC958'
    elif 'uti189' in desc:
        return 'UTI189'
    elif 'dh5' in desc:
        return 'DH5alpha'
    elif 'k-12' in desc or 'mg1655' in desc:
        return 'K-12'
    elif 'scb12' in desc:
        return 'SCB12'
    elif 'scb29' in desc:
        return 'SCB29'
    elif 'scb34' in desc:
        return 'SCB34'
    elif 'scb37' in desc:
        return 'SCB37'
    elif 'scb58' in desc:
        return 'SCB58'
    elif 'scb60' in desc or seq_id == 'SCB60':
        return 'SCB60'
    elif 'scb61' in desc or seq_id == 'SCB61':
        return 'SCB61'

    # Other species (check AFTER E. coli strains)
    # Yersinia: rename Y11 to JB580v as requested
    elif 'yersinia' in desc:
        return 'Y.enterocolitica_JB580v'
    elif 'klebsiella' in desc and 'atcc 43816' in desc:
        return 'K.pneumoniae_ATCC43816'
    elif 'klebsiella' in desc:
        return 'K.pneumoniae'
    elif 'salmonella' in desc and 'sl1344' in desc:
        return 'S.Typhimurium_SL1344'
    elif 'salmonella' in desc and '14028s' in desc:
        return 'S.Typhimurium_14028S'
    elif 'enterobacter' in desc and 'atcc 13047' in desc:
        return 'E.cloacae_ATCC13047'
    elif 'enterobacter' in desc:
        return 'E.cloacae'
    elif 'citrobacter' in desc and ('baa-895' in desc or 'baa 895' in desc):
        return 'C.koseri_BAA-895'
    elif 'citrobacter' in desc:
        return 'C.koseri'
    else:
        return seq_id[:25]


def mask_matching_residues(ref_seq, query_seq):
    """Replace matching residues with dots to highlight differences.

    Dashes (gaps) are preserved to show insertions/deletions clearly.
    """
    result = []
    for r, q in zip(ref_seq, query_seq):
        if q == '-':
            # Keep dashes as dashes (gaps should be visible)
            result.append('-')
        elif r == q:
            # Matching residues become dots
            result.append('.')
        else:
            # Different residues stay as is
            result.append(q)
    return ''.join(result)


def map_ungapped_to_gapped(aligned_seq, ungapped_start, ungapped_end):
    """Map ungapped sequence positions to gapped alignment positions."""
    ungapped_pos = 0
    gapped_start = None
    gapped_end = None

    for i, aa in enumerate(aligned_seq):
        if aa != '-':
            if ungapped_pos == ungapped_start:
                gapped_start = i
            if ungapped_pos == ungapped_end:
                gapped_end = i
                break
            ungapped_pos += 1

    return gapped_start, gapped_end


def find_peptide_in_sequence(ungapped_seq, peptide_seq):
    """Find a peptide sequence in an ungapped sequence."""
    pos = ungapped_seq.find(peptide_seq)
    if pos != -1:
        return pos, pos + len(peptide_seq)
    return None, None


def main():
    print("="*100)
    print("OmpA PEPTIDE LOOP REGION ANALYSIS - SELECTED STRAINS")
    print("Comparing selected E. coli and other Enterobacteriaceae strains")
    print("="*100)
    print()

    # Step 1: Filter and combine sequences
    print("Step 1: Filtering sequences to keep only selected strains...")

    # Read existing sequences (use the most recent alignment file with SCB60/61)
    existing_file = 'mafft_aligned_with_SCB60_61.fasta'
    new_file = 'new_ompA_sequences.fasta'
    filtered_file = 'selected_ompA_sequences.fasta'

    # Extract sequences from existing file (skip MAFFT log at beginning)
    all_seqs = []

    for record in SeqIO.parse(existing_file, 'fasta'):
        # Skip if the record ID looks like MAFFT output
        if not record.id.startswith('outputhat') and len(str(record.seq)) > 100:
            # Remove gaps from existing alignment for re-alignment
            ungapped_seq = str(record.seq).replace('-', '')
            record.seq = Seq(ungapped_seq)
            all_seqs.append(record)

    # Read new sequences
    new_seqs = list(SeqIO.parse(new_file, 'fasta'))
    all_seqs.extend(new_seqs)

    print(f"  Loaded {len(all_seqs)} total sequences")

    # Filter sequences to keep only those in STRAIN_ORDER
    filtered_seqs = []
    strain_counts = {strain: 0 for strain in STRAIN_ORDER}

    for record in all_seqs:
        strain = get_strain_name(record)
        if strain in STRAIN_ORDER:
            filtered_seqs.append((strain, record))
            strain_counts[strain] += 1

    print(f"  Filtered to {len(filtered_seqs)} sequences:")
    for strain in STRAIN_ORDER:
        count = strain_counts[strain]
        if count > 0:
            print(f"    {strain}: {count}")
        else:
            print(f"    {strain}: MISSING!")

    # Write filtered sequences
    SeqIO.write([rec for _, rec in filtered_seqs], filtered_file, 'fasta')
    print(f"  Saved filtered sequences -> {filtered_file}")
    print()

    # Step 2: Align with MAFFT
    print("Step 2: Aligning selected sequences with MAFFT...")
    aligned_file = 'selected_ompA_aligned.fasta'
    mafft_cmd = [
        '/home/david/miniforge3/envs/snippy/bin/mafft',
        '--auto',
        filtered_file
    ]

    with open(aligned_file, 'w') as out_f:
        with open('selected_mafft.log', 'w') as log_f:
            result = subprocess.run(mafft_cmd, stdout=out_f, stderr=log_f)

    if result.returncode == 0:
        print(f"  Alignment complete -> {aligned_file}")
    else:
        print(f"  ERROR: MAFFT alignment failed!")
        sys.exit(1)
    print()

    # Step 3: Extract peptide regions
    print("Step 3: Extracting peptide loop regions...")

    # Read aligned sequences
    aligned_seqs = list(SeqIO.parse(aligned_file, 'fasta'))
    print(f"  Loaded {len(aligned_seqs)} aligned sequences")

    # Find RS218 and EcN as references
    rs218_record = None
    ecn_record = None

    for record in aligned_seqs:
        strain = get_strain_name(record)
        if strain == 'RS218':
            rs218_record = record
        elif strain == 'EcN':
            ecn_record = record

    if not rs218_record or not ecn_record:
        print("  ERROR: Could not find RS218 or EcN reference!")
        sys.exit(1)

    rs218_aligned = str(rs218_record.seq)
    rs218_ungapped = rs218_aligned.replace('-', '')
    ecn_aligned = str(ecn_record.seq)
    ecn_ungapped = ecn_aligned.replace('-', '')

    print(f"  Found reference sequences: RS218 and EcN")
    print()

    # Find peptide positions
    peptide_positions = {}

    print("  Locating peptide regions:")
    for pep_name, pep_seq in PEPTIDES.items():
        # Try RS218 first
        start, end = find_peptide_in_sequence(rs218_ungapped, pep_seq)
        if start is not None:
            ref_seq = rs218_aligned
            ref_name = 'RS218'
            peptide_positions[pep_name] = {
                'reference': ref_name,
                'ref_aligned': ref_seq,
                'ungapped_start': start,
                'ungapped_end': end,
                'sequence': pep_seq
            }
            print(f"    {pep_name}: found in RS218 at position {start}")
        else:
            # Try EcN
            start, end = find_peptide_in_sequence(ecn_ungapped, pep_seq)
            if start is not None:
                ref_seq = ecn_aligned
                ref_name = 'EcN'
                peptide_positions[pep_name] = {
                    'reference': ref_name,
                    'ref_aligned': ref_seq,
                    'ungapped_start': start,
                    'ungapped_end': end,
                    'sequence': pep_seq
                }
                print(f"    {pep_name}: found in EcN at position {start}")
            else:
                print(f"    {pep_name}: NOT FOUND in either RS218 or EcN!")
    print()

    # Step 4: Generate output
    print("Step 4: Generating peptide alignment output...")
    output_file = 'peptide_alignment_selected_strains.txt'

    with open(output_file, 'w') as out:
        out.write("="*100 + "\n")
        out.write("OmpA PEPTIDE LOOP REGIONS - SELECTED STRAINS COMPARISON\n")
        out.write("="*100 + "\n\n")
        out.write(f"Total sequences: {len(aligned_seqs)}\n")
        out.write(f"E. coli strains: EcN, RS218, UTI189\n")
        out.write(f"Other species: C. koseri, E. cloacae, K. pneumoniae, S. Typhimurium (2 strains), Y. enterocolitica\n")
        out.write(f"Alignment: MAFFT\n")
        out.write(f"Format: '.' = match to reference, letter = difference, '-' = gap\n\n")

        # Process each peptide
        for pep_name in ['Peptide_1', 'Peptide_2', 'Peptide_2b', 'Peptide_3a', 'Peptide_3b', 'Peptide_4']:
            if pep_name not in peptide_positions:
                continue

            pep_info = peptide_positions[pep_name]
            ref_aligned = pep_info['ref_aligned']
            ref_name = pep_info['reference']

            # Map to gapped positions
            start, end = map_ungapped_to_gapped(
                ref_aligned,
                pep_info['ungapped_start'],
                pep_info['ungapped_end']
            )

            if start is None or end is None:
                print(f"    WARNING: Could not map {pep_name} to alignment")
                continue

            # Extract peptide region from reference
            ref_peptide = ref_aligned[start:end]

            # Collect all sequences in a dictionary
            strain_to_peptide = {}
            for record in aligned_seqs:
                strain = get_strain_name(record)
                aligned_seq = str(record.seq)
                peptide_region = aligned_seq[start:end]
                strain_to_peptide[strain] = peptide_region

            # Write output in specified order
            out.write("="*100 + "\n")
            out.write(f"{pep_name}: {pep_info['sequence']} (reference: {ref_name})\n")
            out.write("="*100 + "\n")
            out.write(f"{'Strain':<30} {ref_peptide}\n")
            out.write(f"{'-'*30} {'-'*len(ref_peptide)}\n")

            for strain in STRAIN_ORDER:
                if strain in strain_to_peptide:
                    peptide_seq = strain_to_peptide[strain]
                    if strain == ref_name:
                        out.write(f"{strain:<30} {peptide_seq}\n")
                    else:
                        masked = mask_matching_residues(ref_peptide, peptide_seq)
                        out.write(f"{strain:<30} {masked}\n")

            out.write("\n")

        out.write("="*100 + "\n")
        out.write("\nKEY:\n")
        out.write("  '.' = Matches reference sequence\n")
        out.write("  Letter = Differs from reference\n")
        out.write("  '-' = Gap in alignment\n")
        out.write("="*100 + "\n")

    print(f"  Complete! Output saved to: {output_file}")
    print()
    print("="*100)
    print("ANALYSIS COMPLETE")
    print("="*100)
    print()
    print(f"Files generated:")
    print(f"  1. {filtered_file} - Filtered sequences (selected strains only)")
    print(f"  2. {aligned_file} - MAFFT alignment of selected sequences")
    print(f"  3. {output_file} - Peptide loop region comparison")
    print()


if __name__ == "__main__":
    main()
