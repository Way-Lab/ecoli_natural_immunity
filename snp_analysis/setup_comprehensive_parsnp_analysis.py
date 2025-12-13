#!/usr/bin/env python3
"""
Set up a comprehensive parsnp analysis including:
- All samples from yesterday's heatmap (12 strains)
- Mark Nissle Stock (N2SKTQ_7_7)
- EcN reference genome
"""

import os
import shutil
import csv

# Target strains from yesterday's heatmap
TARGET_STRAINS = [
    'Nissle Stock-1',
    'Nissle Stock-2',
    'Nissle Control-1',
    'Nissle Control-2',
    'Nissle Control-3',
    'Nissle Control-4',
    'Nissle Preg-1',
    'Nissle Preg-2',
    'Nissle Preg-5',
    'Nissle Preg-6',
    'Nissle Preg-13',
    'Nissle Preg-14',
    'Nissle Mark stock'  # Additional sample
]

# Reverse mapping: display name to sequence file name
NAME_TO_FILE = {
    'Nissle Stock-1': 'N2SKTQ_1_1',
    'Nissle Stock-2': 'N2SKTQ_2_2',
    'Nissle Control-1': 'N2SKTQ_3_3',
    'Nissle Control-2': 'N2SKTQ_4_4',
    'Nissle Control-3': 'N2SKTQ_5_5',
    'Nissle Control-4': 'N2SKTQ_6_6',
    'Nissle Mark stock': 'N2SKTQ_7_7',
    'Nissle Preg-1': 'N2SKTQ_9_9',
    'Nissle Preg-2': 'N2SKTQ_10_10',
    'Nissle Preg-5': 'Z6M7Y5_1',
    'Nissle Preg-6': 'Z6M7Y5_2',
    'Nissle Preg-13': 'Z6M7Y5_3',
    'Nissle Preg-14': 'Z6M7Y5_4',
}

# Possible locations for assembly files
SEARCH_PATHS = [
    'N2SKTQ_results/parsnp_results_nissle_refseq/assemblies',
    'Z6M7Y5_comparative_analysis/genomes',
    'Z6M7Y5_comparative_analysis_with_refs/genomes',
    'Z6M7Y5_assemblies',
    'Z6M7Y5_comparative_analysis/parsnp_results/assemblies',
    'Z6M7Y5_comparative_analysis_with_refs/parsnp_results/assemblies',
    'N2SKTQ_results/parsnp_results/assemblies',
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
            # Try subdirectory (for Z6M7Y5_assemblies structure)
            candidate = os.path.join(search_path, file_base, f"{file_base}{ext}")
            if os.path.exists(candidate):
                return candidate
    return None

def main():
    print("="*80)
    print("Setting up Comprehensive Parsnp Analysis")
    print("="*80)
    print()

    # Create output directory
    output_dir = 'comprehensive_parsnp_analysis'
    assemblies_dir = os.path.join(output_dir, 'assemblies')
    os.makedirs(assemblies_dir, exist_ok=True)
    print(f"Created directory: {assemblies_dir}\n")

    # Find and copy assembly files
    print("Locating and copying assembly files:")
    print("-" * 80)

    found_files = []
    missing_files = []

    for strain in TARGET_STRAINS:
        file_base = NAME_TO_FILE[strain]
        source_path = find_assembly_file(file_base)

        if source_path:
            # Copy with standardized naming
            dest_filename = f"{file_base}.fasta"
            dest_path = os.path.join(assemblies_dir, dest_filename)
            shutil.copy2(source_path, dest_path)
            found_files.append((strain, file_base, source_path))
            print(f"✓ {strain:25s} ({file_base:15s}) <- {source_path}")
        else:
            missing_files.append((strain, file_base))
            print(f"✗ {strain:25s} ({file_base:15s}) NOT FOUND")

    print()

    # Copy reference genome
    if os.path.exists(REFERENCE_GENOME):
        ref_dest = os.path.join(assemblies_dir, 'reference_Nissle_RefSeq.fasta')
        shutil.copy2(REFERENCE_GENOME, ref_dest)
        print(f"✓ EcN Reference genome copied from: {REFERENCE_GENOME}")
    else:
        print(f"✗ ERROR: Reference genome not found at {REFERENCE_GENOME}")
        return

    print()
    print("="*80)
    print("Summary")
    print("="*80)
    print(f"Found and copied: {len(found_files)} assembly files")
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
    print("Command to run:")
    print(f"parsnp -r {output_dir}/assemblies/reference_Nissle_RefSeq.fasta \\")
    print(f"       -d {output_dir}/assemblies \\")
    print(f"       -o {output_dir}/parsnp_output \\")
    print(f"       -c -v")
    print()

if __name__ == "__main__":
    main()
