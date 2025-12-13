#!/usr/bin/env python3

import sys
import os
import subprocess
import tempfile
import time
import argparse
from pathlib import Path
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import numpy as np


class SequenceAligner:
    def __init__(self, method='mafft', conda_env='snippy'):
        """
        Initialize the aligner with specified method.

        Args:
            method: Alignment tool to use ('mafft', 'clustalo', 'muscle', 'web')
            conda_env: Conda environment containing the alignment tools
        """
        self.method = method
        self.conda_env = conda_env
        self.available_methods = self._check_available_methods()

    def _check_available_methods(self):
        """Check which alignment methods are available."""
        methods = {}

        # Check for local tools in conda environment
        for tool in ['mafft', 'clustalo', 'muscle']:
            try:
                result = subprocess.run(
                    f"conda run -n {self.conda_env} which {tool}",
                    shell=True, capture_output=True, text=True
                )
                if result.returncode == 0:
                    methods[tool] = ('conda', self.conda_env)
            except:
                pass

        # Check for tools in system PATH
        for tool in ['mafft', 'clustalo', 'muscle']:
            if tool not in methods:
                try:
                    result = subprocess.run(
                        f"which {tool}", shell=True, capture_output=True, text=True
                    )
                    if result.returncode == 0:
                        methods[tool] = ('system', None)
                except:
                    pass

        # Web service is always available
        methods['web'] = ('web', None)

        return methods

    def align_sequences(self, fasta_file, output_file=None, mark_divergent=False):
        """
        Perform multiple sequence alignment using the specified method.

        Args:
            fasta_file: Path to input FASTA file
            output_file: Path to save aligned sequences (optional)
            mark_divergent: If True, mark divergent residues as lowercase
        """
        if self.method not in self.available_methods:
            print(f"Warning: {self.method} not available. Available methods: {list(self.available_methods.keys())}")
            # Fall back to first available method
            if self.available_methods:
                self.method = list(self.available_methods.keys())[0]
                print(f"Using {self.method} instead.")
            else:
                raise RuntimeError("No alignment methods available!")

        print(f"Using {self.method} for alignment...")
        start_time = time.time()

        method_type, env = self.available_methods[self.method]

        if method_type == 'web':
            alignment = self._align_with_web(fasta_file)
        else:
            alignment = self._align_with_local_tool(fasta_file, self.method, method_type, env)

        elapsed_time = time.time() - start_time
        print(f"Alignment completed in {elapsed_time:.2f} seconds")

        if output_file:
            if mark_divergent:
                # Save with divergent residues marked as lowercase
                self._save_with_divergent_marked(alignment, output_file)
            else:
                AlignIO.write(alignment, output_file, "fasta")
            print(f"Alignment saved to {output_file}")

        return alignment

    def _save_with_divergent_marked(self, alignment, output_file):
        """
        Save alignment with divergent residues marked as lowercase.

        Args:
            alignment: The alignment object
            output_file: Path to save the modified alignment
        """
        if len(alignment) == 0:
            AlignIO.write(alignment, output_file, "fasta")
            return

        # Get reference sequence (first sequence)
        ref_seq = str(alignment[0].seq).upper()

        # Create modified records
        modified_records = []

        # First sequence stays uppercase
        modified_records.append(SeqRecord(
            Seq(ref_seq),
            id=alignment[0].id,
            description=alignment[0].description
        ))

        # Process remaining sequences
        for record in alignment[1:]:
            seq_str = str(record.seq)
            modified_seq = []

            for i, (ref_aa, seq_aa) in enumerate(zip(ref_seq, seq_str)):
                if seq_aa == '-':
                    # Keep gaps as they are
                    modified_seq.append('-')
                elif ref_aa == '-':
                    # If reference has gap, mark as divergent (lowercase)
                    modified_seq.append(seq_aa.lower())
                elif seq_aa.upper() == ref_aa.upper():
                    # Matching residue - uppercase
                    modified_seq.append(seq_aa.upper())
                else:
                    # Divergent residue - lowercase
                    modified_seq.append(seq_aa.lower())

            modified_records.append(SeqRecord(
                Seq(''.join(modified_seq)),
                id=record.id,
                description=record.description
            ))

        # Write modified alignment
        with open(output_file, 'w') as f:
            for record in modified_records:
                f.write(f">{record.id}")
                if record.description:
                    f.write(f" {record.description}")
                f.write("\n")
                # Write sequence in chunks of 60 characters
                seq_str = str(record.seq)
                for i in range(0, len(seq_str), 60):
                    f.write(seq_str[i:i+60] + "\n")

        print(f"Alignment with divergent residues marked (lowercase) saved to {output_file}")

    def _align_with_local_tool(self, fasta_file, tool, method_type, env):
        """Use local alignment tool."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_out:
            temp_out_name = temp_out.name

        try:
            if tool == 'mafft':
                # MAFFT with auto strategy (fastest for most cases)
                cmd = f"mafft --auto --thread -1 {fasta_file} > {temp_out_name}"
            elif tool == 'clustalo':
                # Clustal Omega with auto threads
                cmd = f"clustalo -i {fasta_file} -o {temp_out_name} --threads=0 --auto"
            elif tool == 'muscle':
                # MUSCLE v5 if available, otherwise v3
                cmd = f"muscle -in {fasta_file} -out {temp_out_name}"

            if method_type == 'conda':
                cmd = f"conda run -n {env} {cmd}"

            print(f"Running: {cmd}")
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if result.returncode != 0:
                print(f"Error running {tool}: {result.stderr}")
                raise RuntimeError(f"Alignment with {tool} failed")

            # Read the alignment
            alignment = AlignIO.read(temp_out_name, "fasta")
            return alignment

        finally:
            # Clean up temp file
            if os.path.exists(temp_out_name):
                os.remove(temp_out_name)

    def _align_with_web(self, fasta_file):
        """Use Clustal Omega web service (fallback option)."""
        try:
            import requests
        except ImportError:
            raise ImportError("requests library required for web service. Install with: pip install requests")

        print("Using Clustal Omega web service (this may take a while)...")

        # Read sequences
        records = list(SeqIO.parse(fasta_file, "fasta"))
        sequences = "\n".join([f">{rec.id}\n{rec.seq}" for rec in records])

        # Submit job
        url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run"
        data = {
            "email": "test@example.com",
            "sequence": sequences,
            "stype": "protein"
        }

        response = requests.post(url, data=data)
        if response.status_code != 200:
            raise RuntimeError(f"Failed to submit job: {response.status_code}")

        job_id = response.text
        print(f"Job submitted with ID: {job_id}")

        # Poll for results
        status_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}"
        while True:
            response = requests.get(status_url)
            status = response.text
            if status == "FINISHED":
                break
            elif status in ["ERROR", "FAILURE"]:
                raise RuntimeError(f"Job failed with status: {status}")
            print(f"Status: {status}, waiting...")
            time.sleep(5)

        # Get results
        result_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/aln-fasta"
        response = requests.get(result_url)

        if response.status_code != 200:
            raise RuntimeError(f"Failed to get results: {response.status_code}")

        # Parse alignment
        from io import StringIO
        alignment = AlignIO.read(StringIO(response.text), "fasta")
        return alignment


def calculate_conservation(alignment):
    """Calculate conservation score for each position in the alignment."""
    conservation_scores = []
    alignment_array = np.array([list(str(record.seq)) for record in alignment])

    for col in range(alignment_array.shape[1]):
        column = alignment_array[:, col]
        # Count unique amino acids (excluding gaps)
        unique_aa = set(column) - {'-'}
        if len(unique_aa) == 0:
            conservation_scores.append(0)
        elif len(unique_aa) == 1:
            conservation_scores.append(1)
        else:
            # Calculate Shannon entropy
            from collections import Counter
            counts = Counter(column)
            total = sum(counts.values())
            entropy = -sum((count/total) * np.log2(count/total)
                          for count in counts.values() if count > 0)
            # Convert to conservation score (1 - normalized entropy)
            max_entropy = np.log2(len(counts))
            conservation_scores.append(1 - (entropy / max_entropy if max_entropy > 0 else 0))

    return conservation_scores


def calculate_statistics(alignment):
    """Calculate alignment statistics."""
    stats = {
        'num_sequences': len(alignment),
        'alignment_length': alignment.get_alignment_length(),
        'identities': [],
        'gaps': []
    }

    # Calculate pairwise identities with first sequence
    first_seq = str(alignment[0].seq)
    for i in range(1, len(alignment)):
        seq = str(alignment[i].seq)
        matches = sum(1 for a, b in zip(first_seq, seq) if a == b and a != '-')
        non_gaps = sum(1 for a, b in zip(first_seq, seq) if a != '-' and b != '-')
        identity = (matches / non_gaps * 100) if non_gaps > 0 else 0
        gaps = seq.count('-') / len(seq) * 100

        stats['identities'].append({
            'id': alignment[i].id,
            'identity': identity,
            'gap_percentage': gaps
        })

    return stats


def plot_alignment_heatmap(alignment, output_file="alignment_heatmap.png"):
    """Create an improved heatmap visualization of the alignment."""
    # Create numerical representation
    amino_acids = "ACDEFGHIKLMNPQRSTVWY-"
    aa_to_int = {aa: i for i, aa in enumerate(amino_acids)}

    # Limit to first 200 positions for visibility
    max_len = min(200, alignment.get_alignment_length())

    num_alignment = []
    for record in alignment:
        num_seq = [aa_to_int.get(aa, -1) for aa in str(record.seq[:max_len])]
        num_alignment.append(num_seq)

    num_alignment = np.array(num_alignment)

    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, len(alignment) * 0.3 + 2),
                                    gridspec_kw={'height_ratios': [1, len(alignment)]})

    # Plot conservation scores
    conservation = calculate_conservation(alignment)[:max_len]
    ax1.plot(conservation, color='blue', linewidth=1)
    ax1.fill_between(range(len(conservation)), conservation, alpha=0.3)
    ax1.set_xlim(0, max_len)
    ax1.set_ylabel('Conservation')
    ax1.set_title('Sequence Conservation')
    ax1.grid(True, alpha=0.3)

    # Plot alignment heatmap
    im = ax2.imshow(num_alignment, cmap="viridis", interpolation="none", aspect="auto")

    # Set labels
    seq_labels = [f"{rec.id[:30]}..." if len(rec.id) > 30 else rec.id
                  for rec in alignment]
    ax2.set_yticks(range(len(alignment)))
    ax2.set_yticklabels(seq_labels, fontsize=8)
    ax2.set_xlabel('Position')
    ax2.set_title(f'Multiple Sequence Alignment (first {max_len} positions)')

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax2, fraction=0.02, pad=0.01)
    cbar.set_label('Amino Acid')

    # Set colorbar ticks for amino acids
    cbar.set_ticks(range(0, len(amino_acids), 2))
    cbar.set_ticklabels([amino_acids[i] for i in range(0, len(amino_acids), 2)])

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Alignment heatmap saved to {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Perform multiple sequence alignment with visualization')
    parser.add_argument('fasta_file', help='Input FASTA file with sequences')
    parser.add_argument('--method', choices=['mafft', 'clustalo', 'muscle', 'web', 'auto'],
                        default='auto', help='Alignment method (default: auto)')
    parser.add_argument('--output', help='Output file for aligned sequences')
    parser.add_argument('--conda-env', default='snippy',
                        help='Conda environment with alignment tools (default: snippy)')
    parser.add_argument('--no-plot', action='store_true', help='Skip visualization')
    parser.add_argument('--plot-output', default='alignment_heatmap.png',
                        help='Output file for alignment plot')

    args = parser.parse_args()

    # Check if input file exists
    if not os.path.exists(args.fasta_file):
        print(f"Error: Input file {args.fasta_file} not found")
        sys.exit(1)

    # Count sequences
    num_sequences = len(list(SeqIO.parse(args.fasta_file, "fasta")))
    print(f"Input file contains {num_sequences} sequences")

    # Initialize aligner
    if args.method == 'auto':
        # Try methods in order of preference
        for method in ['mafft', 'clustalo', 'muscle', 'web']:
            aligner = SequenceAligner(method=method, conda_env=args.conda_env)
            if method in aligner.available_methods and method != 'web':
                break
    else:
        aligner = SequenceAligner(method=args.method, conda_env=args.conda_env)

    print(f"Available alignment methods: {list(aligner.available_methods.keys())}")
    print(f"Selected method: {aligner.method}")

    # Perform alignment
    try:
        # Default output file name if not specified
        output_file = args.output if args.output else "aligned_sequences.fasta"
        alignment = aligner.align_sequences(args.fasta_file, output_file, mark_divergent=True)

        # Calculate and display statistics
        stats = calculate_statistics(alignment)

        print(f"\n{'='*60}")
        print("ALIGNMENT STATISTICS")
        print(f"{'='*60}")
        print(f"Number of sequences: {stats['num_sequences']}")
        print(f"Alignment length: {stats['alignment_length']} positions")

        if stats['identities']:
            print(f"\nPairwise identities with {alignment[0].id}:")
            print(f"{'Sequence':<30} {'Identity':>10} {'Gaps':>10}")
            print("-" * 52)
            for item in stats['identities']:
                seq_id = item['id'][:27] + "..." if len(item['id']) > 30 else item['id']
                print(f"{seq_id:<30} {item['identity']:>9.1f}% {item['gap_percentage']:>9.1f}%")

        # Create visualization
        if not args.no_plot:
            plot_alignment_heatmap(alignment, args.plot_output)

        print(f"\n{'='*60}")
        print("Alignment completed successfully!")

    except Exception as e:
        print(f"Error during alignment: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()