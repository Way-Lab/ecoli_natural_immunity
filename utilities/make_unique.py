
from Bio import SeqIO

def make_fasta_headers_unique(input_file, output_file):
    records = []
    ids = set()
    count = 1
    for record in SeqIO.parse(input_file, "fasta"):
        original_id = record.id
        while record.id in ids:
            record.id = f"{original_id}_{count}"
            count += 1
        ids.add(record.id)
        records.append(record)
    SeqIO.write(records, output_file, "fasta")

make_fasta_headers_unique('aligned_sequences_with_SCB.fasta', 'aligned_sequences_unique.fasta')
