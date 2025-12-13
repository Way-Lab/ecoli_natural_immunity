#!/usr/bin/env python3
"""
Download OmpA sequences from Enterobacteriaceae species using NCBI Entrez.
Requires Biopython and internet connection.
"""

from Bio import Entrez, SeqIO
import time
import sys

# Always tell NCBI who you are
Entrez.email = "your_email@example.com"  # Replace with your email

def search_and_download_ompA(organism, max_seqs=20):
    """
    Search for OmpA protein sequences from a specific organism.
    """
    print(f"\nSearching for OmpA from {organism}...")

    # Search query
    query = f'"{organism}"[Organism] AND (ompA[Gene] OR "outer membrane protein A"[Protein Name])'

    try:
        # Search
        handle = Entrez.esearch(db="protein", term=query, retmax=max_seqs)
        record = Entrez.read(handle)
        handle.close()

        id_list = record["IdList"]
        print(f"  Found {len(id_list)} sequences")

        if not id_list:
            return []

        # Fetch sequences
        time.sleep(0.5)  # Be nice to NCBI
        handle = Entrez.efetch(db="protein", id=id_list, rettype="fasta", retmode="text")
        sequences = list(SeqIO.parse(handle, "fasta"))
        handle.close()

        print(f"  Downloaded {len(sequences)} sequences")
        return sequences

    except Exception as e:
        print(f"  Error: {e}")
        return []


def main():
    # Organisms to search
    organisms = [
        "Escherichia coli",
        "Klebsiella pneumoniae",
        "Klebsiella oxytoca",
        "Enterobacter cloacae",
        "Enterobacter aerogenes",
        "Citrobacter freundii",
        "Citrobacter koseri",
        "Salmonella enterica",
        "Shigella flexneri",
        "Shigella sonnei",
        "Proteus mirabilis",
        "Serratia marcescens",
    ]

    all_sequences = []

    for organism in organisms:
        seqs = search_and_download_ompA(organism, max_seqs=10)
        all_sequences.extend(seqs)
        time.sleep(1)  # Be nice to NCBI servers

    # Write to file
    output_file = "enterobacteriaceae_ompA.fasta"
    with open(output_file, 'w') as f:
        SeqIO.write(all_sequences, f, "fasta")

    print(f"\n{'='*80}")
    print(f"Total sequences downloaded: {len(all_sequences)}")
    print(f"Saved to: {output_file}")
    print(f"{'='*80}")


if __name__ == "__main__":
    main()
