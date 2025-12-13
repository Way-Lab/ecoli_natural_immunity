from Bio import SeqIO

def clean_fasta_headers(input_file, output_file):
    records = []
    seen_ids = set()
    count = 1
    for record in SeqIO.parse(input_file, "fasta"):
        new_id = record.id
        original_id = new_id
        while new_id in seen_ids:
            new_id = f"{original_id}_{count}"
            count += 1
        record.id = new_id
        record.description = ""
        seen_ids.add(new_id)
        records.append(record)
    SeqIO.write(records, output_file, "fasta")
    print(f"Cleaned {len(records)} sequences")
    print(f"Output saved to {output_file}")

clean_fasta_headers("OmpA_peptide_sequences_processed.fasta", "OmpA_peptide_sequences_processed_clean.fasta")