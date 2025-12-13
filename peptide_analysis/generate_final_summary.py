#!/usr/bin/env python3
"""
Generate a comprehensive summary of peptide variability analysis.
"""

import pandas as pd
import re

def parse_variability_report(report_file):
    """
    Parse the detailed variability report and extract key information.
    """
    with open(report_file, 'r') as f:
        content = f.read()

    # Extract total sequences
    total_match = re.search(r'Total sequences analyzed: (\d+)', content)
    total_seqs = int(total_match.group(1)) if total_match else 0

    # Split by peptides
    peptide_sections = re.split(r'={80,}\nPEPTIDE: ', content)

    results = []

    for section in peptide_sections[1:]:
        lines = section.split('\n')
        peptide_name = lines[0].strip()

        # Extract sequence
        seq_match = re.search(r'Sequence: ([A-Z]+)', section)
        sequence = seq_match.group(1) if seq_match else ""

        # Extract found/not found counts
        found_match = re.search(r'Found in: (\d+)/\d+ sequences \(([\d.]+)%\)', section)
        found_count = int(found_match.group(1)) if found_match else 0
        found_percent = float(found_match.group(2)) if found_match else 0.0

        # Extract unique variants
        variants_match = re.search(r'Unique variants found: (\d+)', section)
        num_variants = int(variants_match.group(1)) if variants_match else 0

        # Extract variable positions
        variable_positions = []
        for match in re.finditer(r'Position (\d+) \(reference: ([A-Z])\):[^*]*\*\*\* VARIABLE POSITION \*\*\*', section):
            pos = int(match.group(1))
            ref_aa = match.group(2)
            variable_positions.append((pos, ref_aa))

        results.append({
            'Peptide': peptide_name,
            'Sequence': sequence,
            'Length': len(sequence),
            'Found_Count': found_count,
            'Found_Percent': found_percent,
            'Num_Variants': num_variants,
            'Num_Variable_Positions': len(variable_positions),
            'Variable_Positions': variable_positions
        })

    return results, total_seqs


def create_markdown_report(results, total_seqs, output_file):
    """
    Create a markdown-formatted summary report.
    """
    with open(output_file, 'w') as f:
        f.write("# OmpA Peptide Variability Analysis Report\n\n")
        f.write(f"**Total sequences analyzed:** {total_seqs}\n\n")
        f.write("## Summary\n\n")

        f.write("This analysis examined the conservation of specific peptide sequences within OmpA proteins across ")
        f.write("130 sequences from E. coli and other Enterobacteriaceae species (Klebsiella, Enterobacter, ")
        f.write("Citrobacter, Salmonella, Shigella, Proteus, and Serratia).\n\n")

        # Categorize peptides
        highly_conserved = [r for r in results if r['Found_Percent'] > 0 and r['Num_Variable_Positions'] <= 1]
        moderately_variable = [r for r in results if r['Found_Percent'] > 0 and 1 < r['Num_Variable_Positions'] <= 3]
        highly_variable = [r for r in results if r['Found_Percent'] > 0 and r['Num_Variable_Positions'] > 3]
        not_found = [r for r in results if r['Found_Percent'] == 0]

        f.write("### Peptide Categories\n\n")
        f.write(f"- **Highly conserved** (≤1 variable position): {len(highly_conserved)} peptides\n")
        f.write(f"- **Moderately variable** (2-3 variable positions): {len(moderately_variable)} peptides\n")
        f.write(f"- **Highly variable** (>3 variable positions): {len(highly_variable)} peptides\n")
        f.write(f"- **Not found in any sequences**: {len(not_found)} peptides\n\n")

        f.write("---\n\n")

        f.write("## Detailed Results\n\n")

        for result in sorted(results, key=lambda x: (x['Found_Count'] == 0, -x['Found_Percent'])):
            f.write(f"### {result['Peptide']}\n\n")
            f.write(f"**Sequence:** `{result['Sequence']}`  \n")
            f.write(f"**Length:** {result['Length']} amino acids  \n")
            f.write(f"**Found in:** {result['Found_Count']}/{total_seqs} sequences ({result['Found_Percent']:.1f}%)  \n")

            if result['Found_Count'] > 0:
                f.write(f"**Unique variants:** {result['Num_Variants']}  \n")
                f.write(f"**Variable positions:** {result['Num_Variable_Positions']}  \n")

                if result['Variable_Positions']:
                    f.write(f"\n**Variable positions:**\n")
                    for pos, ref_aa in result['Variable_Positions']:
                        f.write(f"- Position {pos} (reference: {ref_aa})\n")

                # Conservation assessment
                if result['Num_Variable_Positions'] == 0:
                    conservation = "Perfectly conserved"
                elif result['Num_Variable_Positions'] <= 1:
                    conservation = "Highly conserved"
                elif result['Num_Variable_Positions'] <= 3:
                    conservation = "Moderately conserved"
                else:
                    conservation = "Highly variable"

                f.write(f"\n**Conservation assessment:** {conservation}\n")
            else:
                f.write(f"\n**Status:** Not found in any sequences (likely not present in OmpA)\n")

            f.write("\n---\n\n")

        # Key findings
        f.write("## Key Findings\n\n")

        f.write("### Most Conserved Peptides\n\n")
        conserved = sorted([r for r in results if r['Found_Count'] > 0],
                          key=lambda x: x['Num_Variable_Positions'])[:5]
        for i, r in enumerate(conserved, 1):
            f.write(f"{i}. **{r['Peptide']}** - {r['Num_Variable_Positions']} variable positions, ")
            f.write(f"found in {r['Found_Percent']:.1f}% of sequences\n")

        f.write("\n### Most Variable Peptides\n\n")
        variable = sorted([r for r in results if r['Found_Count'] > 0],
                         key=lambda x: x['Num_Variable_Positions'], reverse=True)[:5]
        for i, r in enumerate(variable, 1):
            f.write(f"{i}. **{r['Peptide']}** - {r['Num_Variable_Positions']} variable positions, ")
            f.write(f"found in {r['Found_Percent']:.1f}% of sequences\n")

        f.write("\n### Peptides Not Found\n\n")
        if not_found:
            f.write("The following peptides were not found in any OmpA sequences:\n\n")
            for r in not_found:
                f.write(f"- **{r['Peptide']}** (`{r['Sequence']}`)\n")
            f.write("\nThese peptides may represent:\n")
            f.write("- Synthetic variants designed for testing\n")
            f.write("- Sequences from non-OmpA proteins\n")
            f.write("- Significantly modified peptides not present in natural OmpA variants\n")
        else:
            f.write("All peptides were found in at least some sequences.\n")

        f.write("\n## Recommendations\n\n")
        f.write("Based on this analysis:\n\n")

        if highly_conserved:
            f.write("1. **Highly conserved peptides** ")
            f.write("(e.g., " + ", ".join([r['Peptide'] for r in highly_conserved[:3]]) + ") ")
            f.write("are excellent candidates for:\n")
            f.write("   - Broad-spectrum vaccines or therapeutics targeting Enterobacteriaceae\n")
            f.write("   - Diagnostic assays with high specificity\n")
            f.write("   - Antibody development with cross-reactivity across species\n\n")

        if highly_variable:
            f.write("2. **Highly variable peptides** ")
            f.write("(e.g., " + ", ".join([r['Peptide'] for r in highly_variable[:3]]) + ") ")
            f.write("may be useful for:\n")
            f.write("   - Strain-specific identification\n")
            f.write("   - Understanding evolutionary divergence\n")
            f.write("   - Targeted therapies for specific pathogenic variants\n\n")

        f.write("3. **Species distribution analysis** would benefit from examining which sequences ")
        f.write("contain which variants to identify species-specific or pathotype-specific patterns.\n")


def main():
    import sys

    if len(sys.argv) < 2:
        print("Usage: python generate_final_summary.py <variability_report.txt>")
        sys.exit(1)

    report_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'peptide_analysis_summary.md'

    results, total_seqs = parse_variability_report(report_file)
    create_markdown_report(results, total_seqs, output_file)

    print(f"Summary report written to: {output_file}")


if __name__ == "__main__":
    main()
