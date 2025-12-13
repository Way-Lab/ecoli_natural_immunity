from Bio import SeqIO

def format_fasta(input_file, output_file):
    with open(output_file, 'w') as f_out:
        for record in SeqIO.parse(input_file, "fasta"):
            f_out.write(f">{record.id}\n")
            seq = str(record.seq)
            for i in range(0, len(seq), 60):
                line = seq[i:i+60]
                spaced_line = ' '.join([line[j:j+10] for j in range(0, len(line), 10)])
                f_out.write(spaced_line + '\n')

format_fasta('aligned_sequences_cleaned.fasta', 'aligned_sequences_formatted.fasta')
