from Bio import SeqIO

def format_fasta(input_file, output_file):
    with open(output_file, 'w') as f_out:
        first_record = True
        for record in SeqIO.parse(input_file, "fasta"):
            # Add blank line between sequences (but not before the first one)
            if not first_record:
                f_out.write('\n')
            first_record = False

            # Write header
            f_out.write(f">{record.id}\n")

            # Format sequence: 60 residues per line, space every 10
            seq = str(record.seq)
            for i in range(0, len(seq), 60):
                line = seq[i:i+60]
                spaced_line = ' '.join([line[j:j+10] for j in range(0, len(line), 10)])
                f_out.write(spaced_line + '\n')

    print(f"Formatted {input_file} -> {output_file}")
    print("Added spaces every 10 residues, 60 residues per line")
    print("Added blank lines between sequences")

# Format the aligned OmpA file
format_fasta('OmpA_aligned_final.fasta', 'OmpA_aligned_formatted.fasta')