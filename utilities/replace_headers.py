
def replace_accession_with_strain(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                parts = line[1:].split(' ', 1)
                if len(parts) > 1:
                    f_out.write('>' + parts[1])
                else:
                    f_out.write(line) # Keep original if no space
            else:
                f_out.write(line)

replace_accession_with_strain('aligned_sequences_unique.fasta', 'aligned_sequences_strains.fasta')
