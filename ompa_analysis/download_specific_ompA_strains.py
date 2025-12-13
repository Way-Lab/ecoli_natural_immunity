#!/usr/bin/env python3
"""
Download OmpA sequences from specific bacterial strains using NCBI Entrez.
Requires Biopython and internet connection.
"""

from Bio import Entrez, SeqIO
import time
import sys

# Always tell NCBI who you are
Entrez.email = "your_email@example.com"  # Replace with your email

def search_strain_ompA(organism, strain, gene="ompA"):
    """
    Search for OmpA protein sequences from a specific strain.
    """
    print(f"\nSearching for OmpA from {organism} strain {strain}...")

    # Try multiple query approaches
    queries = [
        f'"{organism}"[Organism] AND {strain}[Strain] AND (ompA[Gene] OR "outer membrane protein A"[Protein Name])',
        f'"{organism} {strain}"[Organism] AND (ompA[Gene] OR "outer membrane protein A"[Protein Name])',
        f'{organism}[Organism] AND {strain} AND ompA[Gene]',
        f'"{organism}"[Organism] AND "{strain}" AND ompA',
    ]

    for i, query in enumerate(queries):
        try:
            print(f"  Trying query {i+1}: {query}")

            # Search protein database
            handle = Entrez.esearch(db="protein", term=query, retmax=50)
            record = Entrez.read(handle)
            handle.close()

            id_list = record["IdList"]
            print(f"  Found {len(id_list)} sequences")

            if id_list:
                # Fetch sequences
                time.sleep(0.5)  # Be nice to NCBI
                handle = Entrez.efetch(db="protein", id=id_list, rettype="fasta", retmode="text")
                sequences = list(SeqIO.parse(handle, "fasta"))
                handle.close()

                print(f"  Downloaded {len(sequences)} sequences")

                # Also try to get more details
                if sequences:
                    return sequences

            time.sleep(0.5)

        except Exception as e:
            print(f"  Error with query {i+1}: {e}")
            continue

    print(f"  No sequences found for {organism} {strain}")
    return []


def search_genome_for_ompA(organism, strain):
    """
    Alternative approach: search for complete genomes and extract OmpA
    """
    print(f"\nSearching genome database for {organism} strain {strain}...")

    queries = [
        f'"{organism}"[Organism] AND {strain}[Strain] AND "complete genome"[Title]',
        f'"{organism} {strain}"[Organism] AND genome',
    ]

    for query in queries:
        try:
            print(f"  Trying: {query}")
            handle = Entrez.esearch(db="nuccore", term=query, retmax=10)
            record = Entrez.read(handle)
            handle.close()

            if record["IdList"]:
                print(f"  Found {len(record['IdList'])} genome records")
                return record["IdList"]

            time.sleep(0.5)

        except Exception as e:
            print(f"  Error: {e}")

    return []


def main():
    # Specific strains to search
    strains_to_search = [
        ("Klebsiella pneumoniae", "43816"),
        ("Salmonella typhimurium", "SL1344"),
        ("Salmonella enterica serovar Typhimurium", "SL1344"),  # Alternative name
        ("Salmonella typhimurium", "ST14028"),
        ("Salmonella enterica serovar Typhimurium", "14028S"),  # Alternative name
        ("Enterobacter cloacae", "13047"),
        ("Citrobacter koseri", "BAA-895"),
    ]

    all_sequences = []
    found_strains = []

    for organism, strain in strains_to_search:
        seqs = search_strain_ompA(organism, strain)

        if seqs:
            all_sequences.extend(seqs)
            found_strains.append(f"{organism} {strain}")
        else:
            # Try genome search as backup
            genome_ids = search_genome_for_ompA(organism, strain)
            if genome_ids:
                print(f"  Found genome(s) but need to extract OmpA manually")

        time.sleep(1)  # Be nice to NCBI servers

    # Write to file
    if all_sequences:
        output_file = "new_ompA_sequences.fasta"
        with open(output_file, 'w') as f:
            SeqIO.write(all_sequences, f, "fasta")

        print(f"\n{'='*80}")
        print(f"Total sequences downloaded: {len(all_sequences)}")
        print(f"Found from strains: {', '.join(found_strains)}")
        print(f"Saved to: {output_file}")
        print(f"{'='*80}")
    else:
        print(f"\n{'='*80}")
        print("No sequences were found. You may need to:")
        print("1. Check if the strain names are correct in NCBI")
        print("2. Search manually on NCBI website")
        print("3. Try alternative strain designations")
        print(f"{'='*80}")


if __name__ == "__main__":
    main()
