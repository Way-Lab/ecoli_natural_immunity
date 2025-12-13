from Bio import SeqIO

def linearize_fasta(input_file, output_file):
    with open(output_file, 'w') as f_out:
        for record in SeqIO.parse(input_file, "fasta"):
            f_out.write(f">{record.id}\n")
            f_out.write(f"{str(record.seq)}\n")

# Linearize SCB60 peptide
linearize_fasta('SCB60_ompA_peptide.fasta', 'SCB60_ompA_peptide_linear.fasta')
print("Linearized SCB60_ompA_peptide.fasta to SCB60_ompA_peptide_linear.fasta")

# Linearize SCB61 peptide
linearize_fasta('SCB61_ompA_peptide.fasta', 'SCB61_ompA_peptide_linear.fasta')
print("Linearized SCB61_ompA_peptide.fasta to SCB61_ompA_peptide_linear.fasta")