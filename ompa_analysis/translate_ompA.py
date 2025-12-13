from Bio.Seq import Seq
import glob

def translate_sequences():
    fasta_files = glob.glob('*_ompA.fasta')
    for fasta_file in fasta_files:
        with open(fasta_file, 'r') as f:
            content = f.read()
        
        # Assuming single sequence per file
        parts = content.split('\n')
        header = parts[0]
        nucleotide_sequence = "".join(parts[1:])
        
        # Translate
        seq_obj = Seq(nucleotide_sequence)
        peptide_sequence = seq_obj.translate()
        
        # Write output
        output_filename = fasta_file.replace('.fasta', '_peptide.fasta')
        with open(output_filename, 'w') as out_f:
            out_f.write(f'>{output_filename.replace(".fasta", "")}\n')
            out_f.write(str(peptide_sequence) + '\n')
        print(f"Translated {fasta_file} to {output_filename}")

translate_sequences()
