#!/usr/bin/env python3
"""
Annotate SNPs with gene information from GFF3 file.
"""

import pandas as pd
from collections import defaultdict

def parse_gff3(gff_file):
    """Parse GFF3 file to get gene annotations."""

    genes = []
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != 'CDS':  # Focus on protein-coding genes
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]

            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value

            # Get gene name or locus tag
            gene_name = attr_dict.get('Name', attr_dict.get('gene', ''))
            locus_tag = attr_dict.get('locus_tag', '')
            product = attr_dict.get('product', 'hypothetical protein')

            genes.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'gene_name': gene_name,
                'locus_tag': locus_tag,
                'product': product
            })

    return genes

def find_gene_for_position(position, genes):
    """Find which gene (if any) contains the given position."""

    for gene in genes:
        if gene['start'] <= position <= gene['end']:
            return gene
    return None

def annotate_vcf_snps(vcf_file, genes, sample_names, output_file):
    """Annotate all SNPs in VCF with gene information."""

    # Read VCF
    snps = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                header = line.strip().split('\t')
                vcf_samples = header[9:]
                continue
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            genotypes = fields[9:]

            # Find gene
            gene = find_gene_for_position(pos, genes)

            # Determine which samples have this SNP
            samples_with_snp = []
            for i, gt in enumerate(genotypes):
                if gt == '1':  # Has alternate allele
                    sample_name = sample_names.get(vcf_samples[i], vcf_samples[i])
                    samples_with_snp.append(sample_name)

            if gene:
                gene_info = f"{gene['gene_name'] or gene['locus_tag']}"
                product = gene['product']
                location = f"{gene['start']}-{gene['end']}"
            else:
                gene_info = "intergenic"
                product = "N/A"
                location = "N/A"

            snps.append({
                'position': pos,
                'ref': ref,
                'alt': alt,
                'change': f"{ref}>{alt}",
                'gene': gene_info,
                'product': product,
                'gene_location': location,
                'samples': ', '.join(samples_with_snp) if samples_with_snp else 'reference'
            })

    # Create DataFrame and save
    df = pd.DataFrame(snps)
    df.to_csv(output_file, sep='\t', index=False)

    return df

# Mapping from file names to display names
FILE_TO_NAME = {
    'N2SKTQ_1_1.fasta.ref': 'Nissle Stock-1 (ref)',
    'N2SKTQ_1_1.fasta': 'Nissle Stock-1',
    'N2SKTQ_2_2.fasta': 'Nissle Stock-2',
    'N2SKTQ_3_3.fasta': 'Nissle Control-1',
    'N2SKTQ_4_4.fasta': 'Nissle Control-2',
    'N2SKTQ_5_5.fasta': 'Nissle Control-3',
    'N2SKTQ_6_6.fasta': 'Nissle Control-4',
    'N2SKTQ_7_7.fasta': 'Nissle Mark stock',
    'N2SKTQ_9_9.fasta': 'Nissle Preg-1',
    'N2SKTQ_10_10.fasta': 'Nissle Preg-2',
    'Z6M7Y5_1.fasta': 'Nissle Preg-5',
    'Z6M7Y5_2.fasta': 'Nissle Preg-6',
    'Z6M7Y5_3.fasta': 'Nissle Preg-13',
    'Z6M7Y5_4.fasta': 'Nissle Preg-14',
    'reference_Nissle_RefSeq.fasta': 'EcN RefSeq',
}

def main():
    print("="*80)
    print("Annotating SNPs with Gene Information")
    print("="*80)
    print()

    # Parse GFF3
    gff_file = 'N2SKTQ_results/N2SKTQ_1_1_bakta_annotation/N2SKTQ_1_1.gff3'
    print(f"Loading gene annotations from: {gff_file}")
    genes = parse_gff3(gff_file)
    print(f"Loaded {len(genes)} gene annotations\n")

    # Annotate SNPs
    vcf_file = 'comprehensive_parsnp_analysis/parsnp_output/parsnp.vcf'
    output_file = 'comprehensive_parsnp_analysis/parsnp_output/snps_annotated.tsv'

    print("Annotating SNPs...")
    df = annotate_vcf_snps(vcf_file, genes, FILE_TO_NAME, output_file)

    print(f"✓ Annotated SNPs saved to: {output_file}\n")

    # Print summary
    print("="*80)
    print("SNP Annotation Summary")
    print("="*80)
    print(f"Total SNPs: {len(df)}")
    print(f"SNPs in genes: {len(df[df['gene'] != 'intergenic'])}")
    print(f"Intergenic SNPs: {len(df[df['gene'] == 'intergenic'])}")
    print()

    # Show first 20 annotated SNPs
    print("First 20 annotated SNPs:")
    print("-" * 80)
    for idx, row in df.head(20).iterrows():
        print(f"Position {row['position']:>7}: {row['change']:>5} | "
              f"Gene: {row['gene']:20s} | {row['product'][:40]:40s} | "
              f"Samples: {row['samples'][:30]}")

    print()
    print("="*80)
    print("Complete! Full annotations available in the TSV file.")
    print("="*80)

if __name__ == "__main__":
    main()
