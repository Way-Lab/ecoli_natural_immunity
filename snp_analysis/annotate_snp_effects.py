#!/usr/bin/env python3
"""
Annotate SNPs from parsnp with gene information from bakta annotation
and determine if SNPs cause amino acid changes.
"""

import re
from Bio import SeqIO
from Bio.Seq import Seq
import sys

# Genetic code
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def parse_gff3(gff_file):
    """Parse GFF3 file and return list of gene features."""
    genes = []
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attributes = parts

            # Only process CDS features
            if feature_type != 'CDS':
                continue

            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value

            gene_name = attr_dict.get('gene', '')
            product = attr_dict.get('Name', attr_dict.get('product', 'Unknown'))
            locus_tag = attr_dict.get('locus_tag', '')

            genes.append({
                'seqid': seqid,
                'start': int(start),
                'end': int(end),
                'strand': strand,
                'gene': gene_name,
                'product': product,
                'locus_tag': locus_tag,
                'phase': int(phase) if phase != '.' else 0
            })

    return genes

def parse_snps(snp_file):
    """Parse SNP file and return list of SNP positions."""
    snps = []
    current_sample = None

    with open(snp_file, 'r') as f:
        for line in f:
            line = line.strip()

            # Skip header lines
            if line.startswith('SNP Details') or line.startswith('Core genome') or line.startswith('===') or not line:
                continue

            # Check for sample name
            sample_match = re.match(r'^(\S+):\s+(\d+)\s+SNPs', line)
            if sample_match:
                current_sample = sample_match.group(1)
                continue

            # Check for separator
            if line.startswith('---'):
                continue

            # Parse SNP line: Position  550880:   T>G (Ref: T, Alt: G)
            snp_match = re.match(r'Position\s+(\d+):\s+(\S+)>(\S+)\s+\(Ref:\s+(\S+),\s+Alt:\s+(\S+)\)', line)
            if snp_match:
                position = int(snp_match.group(1))
                ref = snp_match.group(4)
                alt = snp_match.group(5)

                # Handle multiple alternatives (e.g., "C,G")
                alts = alt.split(',')

                for a in alts:
                    snps.append({
                        'sample': current_sample,
                        'position': position,
                        'ref': ref,
                        'alt': a
                    })

    return snps

def find_gene_for_position(position, genes):
    """Find gene that contains the given position."""
    for gene in genes:
        if gene['start'] <= position <= gene['end']:
            return gene
    return None

def annotate_snp(snp, gene, genome_seq):
    """
    Determine the effect of SNP on protein coding.

    Returns a dictionary with annotation information.
    """
    annotation = {
        'position': snp['position'],
        'ref': snp['ref'],
        'alt': snp['alt'],
        'gene': gene['gene'] if gene else 'intergenic',
        'locus_tag': gene['locus_tag'] if gene else '',
        'product': gene['product'] if gene else '',
        'effect': 'intergenic',
        'ref_codon': '',
        'alt_codon': '',
        'ref_aa': '',
        'alt_aa': '',
        'codon_position': '',
        'aa_position': ''
    }

    if not gene:
        return annotation

    # Calculate position within gene (1-based)
    if gene['strand'] == '+':
        pos_in_gene = snp['position'] - gene['start'] + 1
    else:
        pos_in_gene = gene['end'] - snp['position'] + 1

    # Adjust for phase (reading frame offset)
    pos_in_gene_adjusted = pos_in_gene - gene['phase']

    if pos_in_gene_adjusted < 1:
        annotation['effect'] = 'upstream'
        return annotation

    # Determine codon position (0, 1, or 2 within codon)
    codon_pos = (pos_in_gene_adjusted - 1) % 3
    codon_number = (pos_in_gene_adjusted - 1) // 3 + 1

    # Get the codon sequence
    codon_start_in_gene = (codon_number - 1) * 3 + gene['phase']

    if gene['strand'] == '+':
        genome_codon_start = gene['start'] + codon_start_in_gene - 1
        genome_codon_end = genome_codon_start + 3

        if genome_codon_end > len(genome_seq):
            annotation['effect'] = 'out_of_bounds'
            return annotation

        ref_codon = str(genome_seq[genome_codon_start-1:genome_codon_end-1]).upper()
        alt_codon = list(ref_codon)
        alt_codon[codon_pos] = snp['alt']
        alt_codon = ''.join(alt_codon)
    else:
        # For negative strand, work in reverse
        genome_codon_end = gene['end'] - codon_start_in_gene + 1
        genome_codon_start = genome_codon_end - 2

        if genome_codon_start < 1:
            annotation['effect'] = 'out_of_bounds'
            return annotation

        # Get sequence and reverse complement
        ref_codon = str(genome_seq[genome_codon_start-1:genome_codon_end]).upper()
        ref_codon = reverse_complement(ref_codon)

        # Create alt codon
        alt_codon = list(ref_codon)
        alt_codon[codon_pos] = snp['alt']
        alt_codon = ''.join(alt_codon)

    # Translate codons
    ref_aa = GENETIC_CODE.get(ref_codon, '?')
    alt_aa = GENETIC_CODE.get(alt_codon, '?')

    # Determine effect
    if ref_aa == alt_aa:
        effect = 'synonymous'
    elif alt_aa == '*':
        effect = 'nonsense'
    elif ref_aa == '*':
        effect = 'readthrough'
    else:
        effect = 'missense'

    annotation.update({
        'effect': effect,
        'ref_codon': ref_codon,
        'alt_codon': alt_codon,
        'ref_aa': ref_aa,
        'alt_aa': alt_aa,
        'codon_position': codon_pos + 1,
        'aa_position': codon_number
    })

    return annotation

def main():
    # File paths
    gff_file = '/fastpool/active_data/ecoli_genomics/N2SKTQ_results/N2SKTQ_1_1_bakta_annotation/N2SKTQ_1_1.gff3'
    snp_file = '/fastpool/active_data/ecoli_genomics/parsnp_with_outgroups/parsnp_output/snp_details_by_sample_filtered.txt'
    genome_file = '/fastpool/active_data/ecoli_genomics/N2SKTQ_results/N2SKTQ_1_1_bakta_annotation/N2SKTQ_1_1.fna'
    output_file = '/fastpool/active_data/ecoli_genomics/snp_annotations_with_effects.tsv'

    print("Loading reference genome...")
    genome_seq = None
    for record in SeqIO.parse(genome_file, 'fasta'):
        if genome_seq is None:
            genome_seq = record.seq
        else:
            # Handle multi-contig genome
            print(f"Warning: Multiple contigs found. Using first contig only.")
            break

    print(f"Genome length: {len(genome_seq)} bp")

    print("Parsing GFF3 annotation...")
    genes = parse_gff3(gff_file)
    print(f"Found {len(genes)} CDS features")

    print("Parsing SNPs...")
    snps = parse_snps(snp_file)
    print(f"Found {len(snps)} SNPs")

    print("Annotating SNPs...")
    annotations = []

    for snp in snps:
        # Skip SNPs with 'N' as alternative (ambiguous bases)
        if snp['alt'] == 'N':
            continue

        gene = find_gene_for_position(snp['position'], genes)
        annotation = annotate_snp(snp, gene, genome_seq)
        annotation['sample'] = snp['sample']
        annotations.append(annotation)

    print(f"Writing results to {output_file}...")
    with open(output_file, 'w') as f:
        # Write header
        f.write('\t'.join([
            'Sample', 'Position', 'Ref', 'Alt', 'Gene', 'Locus_Tag', 'Product',
            'Effect', 'Codon_Pos', 'AA_Pos', 'Ref_Codon', 'Alt_Codon', 'Ref_AA', 'Alt_AA'
        ]) + '\n')

        # Write data
        for ann in annotations:
            f.write('\t'.join([
                str(ann.get('sample', '')),
                str(ann['position']),
                ann['ref'],
                ann['alt'],
                ann['gene'],
                ann['locus_tag'],
                ann['product'],
                ann['effect'],
                str(ann.get('codon_position', '')),
                str(ann.get('aa_position', '')),
                ann.get('ref_codon', ''),
                ann.get('alt_codon', ''),
                ann.get('ref_aa', ''),
                ann.get('alt_aa', '')
            ]) + '\n')

    print("\nSummary of effects:")
    effect_counts = {}
    for ann in annotations:
        effect = ann['effect']
        effect_counts[effect] = effect_counts.get(effect, 0) + 1

    for effect, count in sorted(effect_counts.items()):
        print(f"  {effect}: {count}")

    print(f"\nDone! Results written to {output_file}")

if __name__ == '__main__':
    main()
