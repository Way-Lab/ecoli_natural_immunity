from Bio import SeqIO

# List of sequence IDs to remove (using the accession numbers)
sequences_to_remove = {
    'A0A1X7PQ46',
    'A0A2X1JB46',
    'A0A2X3M0D2',
    'Q1RDQ7',
    'A0AAN3V584',
    'A0A8E0FTW8',
    'A0A1X3JK18',
    'W1F721',
    'A0ABC7ZPI0',
    'A0A0H2YXT7',
    'A0A0H3EH07'
}

# Read the aligned file and filter out unwanted sequences
input_file = 'OmpA_peptide_sequences_aligned.fasta'
output_file = 'OmpA_peptide_sequences_aligned_filtered.fasta'

sequences_kept = []
sequences_removed = []

with open(output_file, 'w') as out_handle:
    for record in SeqIO.parse(input_file, 'fasta'):
        # Check if any of the IDs to remove are in the record ID or description
        should_remove = False
        for seq_id in sequences_to_remove:
            if seq_id in record.id or seq_id in record.description:
                should_remove = True
                sequences_removed.append(record.description)
                break

        if not should_remove:
            SeqIO.write(record, out_handle, 'fasta')
            sequences_kept.append(record.description)

print(f"Original file had sequences from: {input_file}")
print(f"Filtered file written to: {output_file}")
print(f"Sequences kept: {len(sequences_kept)}")
print(f"Sequences removed: {len(sequences_removed)}")
print("\nRemoved sequences:")
for seq in sequences_removed:
    print(f"  - {seq[:80]}...")  # Print first 80 chars of each removed sequence header