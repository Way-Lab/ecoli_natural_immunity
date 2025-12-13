from Bio import SeqIO

def process_fasta(input_file, output_file):
    with open(output_file, 'w') as f_out:
        for record in SeqIO.parse(input_file, "fasta"):
            f_out.write(f">{record.description}\n")
            f_out.write(f"{str(record.seq)}\n")

process_fasta('OmpA_peptide_seqeunces.fasta', 'OmpA_peptide_sequences_processed.fasta')