from Bio import SeqIO

def linearize_fasta(input_file, output_file):
    with open(output_file, 'w') as f_out:
        for record in SeqIO.parse(input_file, "fasta"):
            f_out.write(f">{record.id}\n")
            f_out.write(f"{str(record.seq)}\n")

linearize_fasta('aligned_sequences_corrected.fasta', 'aligned_sequences_linearized.fasta')