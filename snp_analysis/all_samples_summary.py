#!/usr/bin/env python3
"""
Generate a comprehensive summary of all sequencing results across all 14 samples.
"""

import os
import glob

def parse_summary_file(summary_file):
    """Parse summary.txt file for key metrics."""
    stats = {}
    if not os.path.exists(summary_file):
        return None

    with open(summary_file, 'r') as f:
        for line in f:
            line = line.strip()
            if '|' in line and not line.startswith('|:'):
                parts = [p.strip() for p in line.split('|') if p.strip()]
                if len(parts) == 2:
                    key = parts[0]
                    value = parts[1]
                    stats[key] = value
    return stats

def parse_stats_file(stats_file):
    """Parse stats.tsv file for numeric metrics."""
    stats = {}
    if not os.path.exists(stats_file):
        return None

    with open(stats_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                stats[parts[0]] = parts[1]
    return stats

def main():
    base_path = "N2SKTQ_results"
    samples = []

    # Collect all sample directories
    for i in range(1, 15):
        sample_name = f"N2SKTQ_{i}_{i}"
        summary_file = f"{base_path}/{sample_name}/ONT-only/{sample_name}_summary.txt"
        stats_file = f"{base_path}/{sample_name}/ONT-only/{sample_name}_stats.tsv"

        summary = parse_summary_file(summary_file)
        stats = parse_stats_file(stats_file)

        samples.append({
            'name': sample_name,
            'summary': summary,
            'stats': stats,
            'number': i
        })

    # Print comprehensive table
    print("="*150)
    print("COMPREHENSIVE SEQUENCING RESULTS SUMMARY")
    print("="*150)
    print()

    # Header
    print(f"{'Sample':<15} {'Total Reads':<12} {'Total bp':<12} {'Longest':<10} {'Genome':<10} {'Contigs':<10} {'Genes':<10} {'Status':<15}")
    print("-"*150)

    # Data rows
    failed_samples = []
    for sample in samples:
        name = sample['name']
        stats = sample['stats']
        summary = sample['summary']

        if stats is None:
            print(f"{name:<15} {'N/A':<12} {'N/A':<12} {'N/A':<10} {'N/A':<10} {'N/A':<10} {'N/A':<10} {'FAILED':<15}")
            failed_samples.append(name)
            continue

        total_reads = stats.get('total_num_reads', 'N/A')
        total_bp = stats.get('total_seq_mbp', 'N/A')
        longest = stats.get('longest_read', 'N/A')

        if summary:
            genome_size = summary.get('Genome size (Mb)', 'N/A')
            contigs = summary.get('Number of contigs', 'N/A')
            genes = summary.get('Number of genes annotated', 'N/A')

            # Determine status
            if genome_size and 'Mb' in genome_size:
                size_mb = float(genome_size.replace(' Mb', ''))
                if size_mb >= 4.5:
                    status = "✓ SUCCESS"
                elif size_mb >= 1.0:
                    status = "⚠ PARTIAL"
                else:
                    status = "✗ FAILED"
                    failed_samples.append(name)
            else:
                status = "✗ FAILED"
                failed_samples.append(name)
        else:
            genome_size = stats.get('genome_size_mbp', 'N/A')
            contigs = 'N/A'
            genes = 'N/A'
            status = "✗ FAILED"
            failed_samples.append(name)

        print(f"{name:<15} {total_reads:<12} {total_bp:<12} {longest:<10} {genome_size:<10} {contigs:<10} {genes:<10} {status:<15}")

    # Summary statistics
    print()
    print("="*150)
    print("SUMMARY")
    print("="*150)
    print(f"\nTotal samples: 14")
    print(f"Successful assemblies: {14 - len(failed_samples)}")
    print(f"Failed assemblies: {len(failed_samples)}")
    print(f"Success rate: {((14 - len(failed_samples))/14)*100:.1f}%")

    if failed_samples:
        print(f"\nFailed samples: {', '.join(failed_samples)}")

    print("\n" + "="*150)
    print("NOTES")
    print("="*150)
    print("""
Expected E. coli genome characteristics:
- Genome size: ~5.0-5.5 Mb
- Number of contigs: 50-200 (depends on assembly quality)
- Number of genes: ~4,500-5,500
- Coverage: 30-100x

Failed samples analysis:
- N2SKTQ_13_13: Assembly collapse - only 4 tiny contigs (11.7 kb total)
- N2SKTQ_14_14: Complete assembly failure - no genome produced

Recommended action:
- Re-extract DNA from N2SKTQ_13_13 and N2SKTQ_14_14
- Verify sample identity before re-sequencing
- Perform quality control on DNA before library preparation
    """)

if __name__ == "__main__":
    main()
