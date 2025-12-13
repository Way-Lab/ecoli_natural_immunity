#!/usr/bin/env python3
"""
Set up parsnp analysis including divergent outgroup samples.
These will be included to define the core genome but excluded from the heatmap.
"""

import os
import shutil

# Target Nissle strains for heatmap (same 13 as before)
TARGET_STRAINS = [
    'Nissle Stock-1',
    'Nissle Stock-2',
    'Nissle Control-1',
    'Nissle Control-2',
    'Nissle Control-3',
    'Nissle Control-4',
    'Nissle Mark stock',
    'Nissle Preg-1',
    'Nissle Preg-2',
    'Nissle Preg-5',
    'Nissle Preg-6',
    'Nissle Preg-13',
    'Nissle Preg-14',
]

# Divergent outgroups to include in analysis but exclude from heatmap
OUTGROUPS = [
    'CFT073/rpoS',
    'MS8707',
    'MS25509',
    'Nissle Preg-16',
    'Nissle Preg-18',
]

# All samples for parsnp
ALL_SAMPLES = TARGET_STRAINS + OUTGROUPS

# Reverse mapping: display name to sequence file name
NAME_TO_FILE = {
    'Nissle Stock-1': 'N2SKTQ_1_1',
    'Nissle Stock-2': 'N2SKTQ_2_2',
    'Nissle Control-1': 'N2SKTQ_3_3',
    'Nissle Control-2': 'N2SKTQ_4_4',
    'Nissle Control-3': 'N2SKTQ_5_5',
    'Nissle Control-4': 'N2SKTQ_6_6',
    'Nissle Mark stock': 'N2SKTQ_7_7',
    'CFT073/rpoS': 'N2SKTQ_8_8',
    'Nissle Preg-1': 'N2SKTQ_9_9',
    'Nissle Preg-2': 'N2SKTQ_10_10',
    'MS8707': 'N2SKTQ_11_11',
    'MS25509': 'N2SKTQ_12_12',
    'Nissle Preg-5': 'Z6M7Y5_1',
    'Nissle Preg-6': 'Z6M7Y5_2',
    'Nissle Preg-13': 'Z6M7Y5_3',
    'Nissle Preg-14': 'Z6M7Y5_4',
    'Nissle Preg-16': 'Z6M7Y5_5',
    'Nissle Preg-18': 'Z6M7Y5_6',
}

# Possible locations for assembly files
SEARCH_PATHS = [
    'N2SKTQ_results/parsnp_results_nissle_refseq/assemblies',
    'Z6M7Y5_comparative_analysis/genomes',
    'Z6M7Y5_comparative_analysis_with_refs/genomes',
    'Z6M7Y5_assemblies',
]

# Reference genome
REFERENCE_GENOME = 'N2SKTQ_results/parsnp_results_nissle_refseq/assemblies/reference_Nissle_RefSeq.fasta'

def find_assembly_file(file_base):
    """Find assembly file in search paths."""
    for search_path in SEARCH_PATHS:
        # Try different extensions and subdirectories
        for ext in ['.fasta', '_hybrid_assembly.fasta']:
            # Try direct path
            candidate = os.path.join(search_path, f"{file_base}{ext}")
            if os.path.exists(candidate):
                return candidate
            # Try subdirectory
            candidate = os.path.join(search_path, file_base, f"{file_base}{ext}")
            if os.path.exists(candidate):
                return candidate
    return None

def main():
    print("="*80)
    print("Setting up Parsnp Analysis with Outgroups")
    print("="*80)
    print()

    # Create output directory
    output_dir = 'parsnp_with_outgroups'
    assemblies_dir = os.path.join(output_dir, 'assemblies')
    os.makedirs(assemblies_dir, exist_ok=True)
    print(f"Created directory: {assemblies_dir}\n")

    # Find and copy assembly files
    print("Target strains for heatmap (13):")
    print("-" * 80)
    for strain in TARGET_STRAINS:
        print(f"  • {strain}")

    print("\nOutgroups (for core genome definition only, 5):")
    print("-" * 80)
    for strain in OUTGROUPS:
        print(f"  • {strain}")

    print("\n" + "="*80)
    print("Locating and copying assembly files:")
    print("-" * 80)

    found_files = []
    missing_files = []

    for strain in ALL_SAMPLES:
        file_base = NAME_TO_FILE[strain]
        source_path = find_assembly_file(file_base)

        if source_path:
            # Copy with standardized naming
            dest_filename = f"{file_base}.fasta"
            dest_path = os.path.join(assemblies_dir, dest_filename)
            shutil.copy2(source_path, dest_path)
            found_files.append((strain, file_base, source_path))
            marker = "✓" if strain in TARGET_STRAINS else "○"
            print(f"{marker} {strain:25s} ({file_base:15s}) <- {source_path}")
        else:
            missing_files.append((strain, file_base))
            print(f"✗ {strain:25s} ({file_base:15s}) NOT FOUND")

    print()

    # Copy reference genome
    if os.path.exists(REFERENCE_GENOME):
        ref_dest = os.path.join(assemblies_dir, 'reference_Nissle_RefSeq.fasta')
        shutil.copy2(REFERENCE_GENOME, ref_dest)
        print(f"✓ EcN RefSeq reference genome copied")
    else:
        print(f"✗ ERROR: Reference genome not found at {REFERENCE_GENOME}")
        return

    print()
    print("="*80)
    print("Summary")
    print("="*80)
    print(f"Found and copied: {len(found_files)} assembly files")
    print(f"  - Target strains: {len([s for s in found_files if s[0] in TARGET_STRAINS])}")
    print(f"  - Outgroups: {len([s for s in found_files if s[0] in OUTGROUPS])}")
    print(f"Missing: {len(missing_files)} files")

    if missing_files:
        print("\nMissing files:")
        for strain, file_base in missing_files:
            print(f"  - {strain} ({file_base})")

    print()
    print("="*80)
    print("Ready to run parsnp!")
    print("="*80)
    print()
    print("Note: Outgroups will be included in analysis to define core genome,")
    print("      but excluded from final heatmap visualization.")

if __name__ == "__main__":
    main()
