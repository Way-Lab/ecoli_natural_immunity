from Bio import SeqIO

def fix_duplicate_headers(input_file, output_file):
    records = []
    for record in SeqIO.parse(input_file, "fasta"):
        # The ID appears to be duplicated in the description
        # Keep only the ID, remove the description
        record.description = ""
        records.append(record)

    SeqIO.write(records, output_file, "fasta")
    print(f"Fixed headers for {len(records)} sequences")
    print(f"Output saved to {output_file}")

# Fix the aligned file headers
fix_duplicate_headers("OmpA_aligned_clean_with_lowercase.fasta",
                     "OmpA_aligned_final.fasta")