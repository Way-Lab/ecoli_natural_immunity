#!/usr/bin/env python3
"""
Analyze new ompA sequences and compare peptide loop regions

This script:
1. Combines new sequences with existing E. coli sequences
2. Performs multiple sequence alignment with MAFFT
3. Extracts peptide loop regions
4. Creates comprehensive comparison output

Usage:
    python analyze_new_ompA_with_loops.py
"""

import subprocess
import sys
from Bio import SeqIO


# Define peptide sequences to look for (from previous analysis)
PEPTIDES = {
    'Peptide_1': 'WSQYHDTGFINNNGPTHEN',
    'Peptide_2': 'GRMPYKGSVENGAYKA',
    'Peptide_2b': 'GRMPYKGDNINGAYKA',
    'Peptide_3a': 'KSNVYGKN',
    'Peptide_3b': 'KANVPGGASFKD',
    'Peptide_4': 'TNNIGDAHTIGTRPDNGM',
}


def get_strain_name(record):
    """Extract a clean strain name from a SeqRecord."""
    desc = record.description.lower()
    seq_id = record.id

    # New strains from new_ompA_sequences.fasta
    if 'yersinia' in desc:
        return 'Y.enterocolitica'
    elif 'klebsiella' in desc and 'atcc 43816' in desc:
        return 'K.pneumoniae_ATCC43816'
    elif 'klebsiella' in desc and 'xvp22030' in seq_id.lower():
        return 'K.pneumoniae'
    elif 'salmonella' in desc and 'sl1344' in desc:
        return 'S.Typhimurium_SL1344'
    elif 'salmonella' in desc and '14028s' in desc:
        return 'S.Typhimurium_14028S'
    elif 'enterobacter' in desc and 'atcc 13047' in desc:
        return 'E.cloacae_ATCC13047'
    elif 'enterobacter' in desc and 'kgb06032' in seq_id.lower():
        return 'E.cloacae'
    elif 'citrobacter' in desc:
        return 'C.koseri'

    # Existing E. coli strains
    elif 'rs218' in desc or seq_id == 'RS218':
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
    else:
        return seq_id[:25]


def mask_matching_residues(ref_seq, query_seq):
    """Replace matching residues with dots to highlight differences."""
    result = []
    for r, q in zip(ref_seq, query_seq):
        if r == q:
            result.append('.')
        else:
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
    print("OmpA PEPTIDE LOOP REGION ANALYSIS")
    print("Comparing new sequences with existing E. coli strains")
    print("="*100)
    print()

    # Step 1: Combine sequences
    print("Step 1: Combining new sequences with existing E. coli sequences...")

    # Read existing sequences (use the most recent alignment file with SCB60/61)
    existing_file = 'mafft_aligned_with_SCB60_61.fasta'
    new_file = 'new_ompA_sequences.fasta'
    combined_file = 'combined_ompA_all_species.fasta'

    # Extract sequences from existing file (skip MAFFT log at beginning)
    existing_seqs = []
    in_sequences = False
    for record in SeqIO.parse(existing_file, 'fasta'):
        # Skip if the record ID looks like MAFFT output
        if not record.id.startswith('outputhat') and len(str(record.seq)) > 100:
            # Remove gaps from existing alignment for re-alignment
            from Bio.Seq import Seq
            ungapped_seq = str(record.seq).replace('-', '')
            record.seq = Seq(ungapped_seq)
            existing_seqs.append(record)

    print(f"  Loaded {len(existing_seqs)} existing E. coli sequences")

    # Read new sequences
    new_seqs = list(SeqIO.parse(new_file, 'fasta'))
    print(f"  Loaded {len(new_seqs)} new sequences")

    # Combine and write
    all_seqs = existing_seqs + new_seqs
    SeqIO.write(all_seqs, combined_file, 'fasta')
    print(f"  Combined {len(all_seqs)} sequences -> {combined_file}")
    print()

    # Step 2: Align with MAFFT
    print("Step 2: Aligning all sequences with MAFFT...")
    aligned_file = 'combined_ompA_all_species_aligned.fasta'
    mafft_cmd = [
        '/home/david/miniforge3/envs/snippy/bin/mafft',
        '--auto',
        combined_file
    ]

    with open(aligned_file, 'w') as out_f:
        with open('combined_mafft.log', 'w') as log_f:
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
    output_file = 'peptide_alignment_with_new_species.txt'

    with open(output_file, 'w') as out:
        out.write("="*100 + "\n")
        out.write("OmpA PEPTIDE LOOP REGIONS - CROSS-SPECIES COMPARISON\n")
        out.write("="*100 + "\n\n")
        out.write(f"Total sequences: {len(aligned_seqs)}\n")
        out.write(f"Species included: E. coli, Klebsiella, Salmonella, Enterobacter, Citrobacter, Yersinia\n")
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

            # Collect all sequences
            all_peptides = []
            for record in aligned_seqs:
                strain = get_strain_name(record)
                aligned_seq = str(record.seq)
                peptide_region = aligned_seq[start:end]
                all_peptides.append((strain, peptide_region))

            # Sort: E. coli first, then others
            ecoli_strains = []
            other_strains = []

            for strain, peptide in all_peptides:
                if strain.startswith(('RS218', 'EcN', 'O157', 'EC958', 'UTI189',
                                      'DH5', 'K-12', 'CFT073', 'SCB')):
                    ecoli_strains.append((strain, peptide))
                else:
                    other_strains.append((strain, peptide))

            # Sort alphabetically within each group
            ecoli_strains.sort()
            other_strains.sort()
            sorted_peptides = ecoli_strains + other_strains

            # Write output
            out.write("="*100 + "\n")
            out.write(f"{pep_name}: {pep_info['sequence']} (reference: {ref_name})\n")
            out.write("="*100 + "\n")
            out.write(f"{'Strain':<30} {ref_peptide}\n")
            out.write(f"{'-'*30} {'-'*len(ref_peptide)}\n")

            for strain, peptide_seq in sorted_peptides:
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
    print(f"  1. {combined_file} - Combined sequences")
    print(f"  2. {aligned_file} - MAFFT alignment of all sequences")
    print(f"  3. {output_file} - Peptide loop region comparison")
    print()


if __name__ == "__main__":
    main()
