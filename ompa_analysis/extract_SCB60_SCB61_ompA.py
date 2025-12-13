import re
from Bio.Seq import Seq

def extract_sequences(scaffold_file, scaffold_id, start, end, strand, output_file):
    with open(scaffold_file, 'r') as f:
        content = f.read()

    # Find the specific scaffold
    match = re.search(f'>{scaffold_id}[^\n]*\n([ATGC\n]+)', content)
    if match:
        sequence = match.group(1).replace('\n', '')

        # Adjust for 1-based indexing from BLAST
        start -= 1

        # Extract the subsequence
        subsequence = sequence[start:end]

        if strand == '-':
            seq_obj = Seq(subsequence)
            subsequence = str(seq_obj.reverse_complement())

        # Write nucleotide sequence
        with open(output_file, 'w') as out_f:
            out_f.write(f'>{output_file.replace(".fasta", "")}\n')
            out_f.write(subsequence + '\n')
        print(f"Extracted nucleotide sequence to {output_file}")

        # Translate and write peptide sequence
        seq_obj = Seq(subsequence)
        peptide = str(seq_obj.translate())

        peptide_file = output_file.replace('_ompA.fasta', '_ompA_peptide.fasta')
        with open(peptide_file, 'w') as out_f:
            out_f.write(f'>{output_file.replace("_ompA.fasta", "")} ompA peptide\n')
            out_f.write(peptide + '\n')
        print(f"Translated peptide sequence to {peptide_file}")
    else:
        print(f"Scaffold {scaffold_id} not found in {scaffold_file}")

# SCB60
extract_sequences(
    'SCB60_genome.contigs.fasta',
    'scaffold_1',
    151505,
    152542,
    '+',
    'SCB60_ompA.fasta'
)

# SCB61
extract_sequences(
    'SCB61_genome.contigs.fasta',
    'scaffold_11',
    19043,
    20080,
    '-',
    'SCB61_ompA.fasta'
)