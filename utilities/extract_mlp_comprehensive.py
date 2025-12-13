#!/usr/bin/env python3
"""
Comprehensive extraction and comparison of Murein lipoprotein (Lpp/MLP) sequences
from E. coli strains using tBLASTn.

Strains to analyze:
- EcN (Nissle 1917)
- UTI89
- RS218
- SCB strains: SCB12, SCB29, SCB61, SCB34, SCB58, SCB37, SCB60
"""

import os
import subprocess
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from io import StringIO

# Reference Lpp protein sequence (E. coli K12)
LPP_REFERENCE = "MKATKLVLGAVILGSTLLAGCSSNAKIDQLSSDVQTLNAKVDQLSNDVNAMRSDVQAAKDDAARANQRLDNMATKYRK"

# Genome files for reference strains
GENOME_FILES = {
    'UTI89': 'HaslamNanoporeSeq/phylogenetic_analysis/reference_genomes/ecoli_UTI89.fna',
    'RS218': 'HaslamNanoporeSeq/phylogenetic_analysis/reference_genomes/ecoli_RS218.fna.gz',
    'EcN': 'ncbi_dataset/data/GCF_000800845.1/GCF_000800845.1_ASM80084v1_genomic.fna',  # May need to find
}

# SCB scaffold/genome files
SCB_GENOME_FILES = {
    'SCB29': 'SCB29_scaffolds.fasta',
    'SCB37': 'SCB37_scaffolds.fasta',
    'SCB58': 'SCB58_scaffolds.fasta',
    'SCB60': 'SCB60_genome.contigs.fasta',
    'SCB61': 'SCB61_genome.contigs.fasta',
}

# N2SKTQ annotation files (already have protein sequences)
N2SKTQ_FAA = {
    f'N2SKTQ_{i}_{i}': f'N2SKTQ_results/N2SKTQ_{i}_{i}/ONT-only/annotation/N2SKTQ_{i}_{i}.faa'
    for i in range(1, 15)
}


def extract_lpp_from_faa(faa_file, strain_name):
    """Extract Lpp sequence from a protein FASTA file."""
    if not os.path.exists(faa_file):
        return None

    for record in SeqIO.parse(faa_file, 'fasta'):
        desc = record.description.lower()
        if 'murein lipoprotein' in desc and 'lpp' in desc:
            return str(record.seq)
    return None


def run_tblastn(protein_seq, genome_file, strain_name):
    """Run tBLASTn to find protein in genome and return translated hit."""
    # Handle gzipped files
    if genome_file.endswith('.gz'):
        import gzip
        with gzip.open(genome_file, 'rt') as f:
            genome_content = f.read()
        # Write to temp file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fna', delete=False) as tmp:
            tmp.write(genome_content)
            genome_file = tmp.name

    if not os.path.exists(genome_file):
        print(f"  {strain_name}: Genome file not found: {genome_file}")
        return None

    # Create temp query file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as query_file:
        query_file.write(f">lpp_query\n{protein_seq}\n")
        query_path = query_file.name

    # Run tBLASTn
    try:
        result = subprocess.run([
            '/home/david/miniforge3/envs/bakta/bin/tblastn',
            '-query', query_path,
            '-subject', genome_file,
            '-outfmt', '6 sseqid sstart send sframe evalue pident sseq',
            '-max_target_seqs', '1',
            '-evalue', '1e-10'
        ], capture_output=True, text=True, timeout=60)

        os.unlink(query_path)

        if result.returncode != 0 or not result.stdout.strip():
            return None

        # Parse result
        line = result.stdout.strip().split('\n')[0]
        parts = line.split('\t')
        if len(parts) >= 7:
            pident = float(parts[5])
            sseq = parts[6].replace('-', '')  # Remove gaps
            return {'identity': pident, 'sequence': sseq}

    except Exception as e:
        print(f"  {strain_name}: BLAST error - {e}")
        os.unlink(query_path)
        return None

    return None


def main():
    print("="*80)
    print("MUREIN LIPOPROTEIN (Lpp/MLP) COMPREHENSIVE EXTRACTION")
    print("="*80)
    print()

    lpp_sequences = {}

    # Reference
    lpp_sequences['K12_reference'] = LPP_REFERENCE
    print(f"Reference K12: {len(LPP_REFERENCE)} aa")

    # Extract from N2SKTQ annotation files
    print("\n--- N2SKTQ Samples (from protein annotations) ---")
    for strain, faa in sorted(N2SKTQ_FAA.items()):
        seq = extract_lpp_from_faa(faa, strain)
        if seq:
            lpp_sequences[strain] = seq
            print(f"  {strain}: {len(seq)} aa")
        elif os.path.exists(faa):
            print(f"  {strain}: Lpp not found in annotation")
        else:
            print(f"  {strain}: No annotation file")

    # Extract from SCB scaffold files using tBLASTn
    print("\n--- SCB Strains (from genome scaffolds via tBLASTn) ---")
    for strain, genome in sorted(SCB_GENOME_FILES.items()):
        result = run_tblastn(LPP_REFERENCE, genome, strain)
        if result:
            lpp_sequences[strain] = result['sequence']
            print(f"  {strain}: {len(result['sequence'])} aa ({result['identity']:.1f}% identity)")
        else:
            print(f"  {strain}: Not found or BLAST failed")

    # Extract from reference genomes via tBLASTn
    print("\n--- Reference Strains (from genomes via tBLASTn) ---")
    for strain, genome in sorted(GENOME_FILES.items()):
        if os.path.exists(genome) or os.path.exists(genome.replace('.gz', '')):
            genome_path = genome if os.path.exists(genome) else genome.replace('.gz', '')
            result = run_tblastn(LPP_REFERENCE, genome_path, strain)
            if result:
                lpp_sequences[strain] = result['sequence']
                print(f"  {strain}: {len(result['sequence'])} aa ({result['identity']:.1f}% identity)")
            else:
                print(f"  {strain}: Not found via BLAST")
        else:
            print(f"  {strain}: Genome file not found")

    # Write output FASTA
    output_file = "mlp_sequences_all.fasta"
    with open(output_file, 'w') as f:
        for name, seq in sorted(lpp_sequences.items()):
            f.write(f">{name}\n{seq}\n")

    print(f"\nWrote {len(lpp_sequences)} sequences to {output_file}")

    # Comparison
    print("\n" + "="*80)
    print("SEQUENCE COMPARISON TO K12 REFERENCE")
    print("="*80)
    print(f"\nK12 Reference: {LPP_REFERENCE}")
    print(f"Length: {len(LPP_REFERENCE)} aa\n")

    # Group by identity
    identical = []
    different = []

    for name, seq in sorted(lpp_sequences.items()):
        if name == 'K12_reference':
            continue

        if seq == LPP_REFERENCE:
            identical.append(name)
        else:
            # Find differences
            diffs = []
            for i, (a, b) in enumerate(zip(seq, LPP_REFERENCE)):
                if a != b:
                    diffs.append(f"{b}{i+1}{a}")

            len_diff = len(seq) - len(LPP_REFERENCE)
            different.append((name, seq, diffs, len_diff))

    print(f"IDENTICAL to K12 ({len(identical)} strains):")
    for name in identical:
        print(f"  - {name}")

    if different:
        print(f"\nDIFFERENT from K12 ({len(different)} strains):")
        for name, seq, diffs, len_diff in different:
            diff_str = ', '.join(diffs[:5])
            if len(diffs) > 5:
                diff_str += f"... (+{len(diffs)-5} more)"
            print(f"  - {name}: {len(diffs)} substitutions, length diff: {len_diff}")
            if diffs:
                print(f"    Changes: {diff_str}")
            print(f"    Sequence: {seq}")
    else:
        print("\nAll strains have IDENTICAL Lpp sequences!")

    print("\n" + "="*80)


if __name__ == "__main__":
    main()
