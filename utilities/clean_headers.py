
from Bio import SeqIO
import re

def clean_headers(input_file, output_file):
    records = []
    seen_ids = set()
    counter = 1
    for record in SeqIO.parse(input_file, "fasta"):
        description = record.description
        description = description.replace(record.id, '').strip()
        parts = description.split(' ')
        new_id = parts[0]
        new_id = re.sub(r'[^a-zA-Z0-9_.-]', '', new_id)

        original_new_id = new_id
        while new_id in seen_ids:
            new_id = f"{original_new_id}_{counter}"
            counter += 1
        
        record.id = new_id
        record.description = ''
        seen_ids.add(new_id)
        records.append(record)

    with open(output_file, "w") as f_out:
        SeqIO.write(records, f_out, "fasta")

clean_headers('aligned_sequences_with_SCB.fasta', 'aligned_sequences_strains.fasta')
