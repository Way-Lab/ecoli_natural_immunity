
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

        with open(output_file, 'w') as out_f:
            out_f.write(f'>{output_file.replace(".fasta", "")}\n')
            out_f.write(subsequence + '\n')
        print(f"Extracted sequence to {output_file}")
    else:
        print(f"Scaffold {scaffold_id} not found in {scaffold_file}")

# SCB58
extract_sequences(
    'SCB58_scaffolds.fasta',
    'JAHUUC010000003.1',
    173295,
    174344,
    '+',
    'SCB58_ompA.fasta'
)

# SCB29
extract_sequences(
    'SCB29_scaffolds.fasta',
    'JAHUUB010000010.1',
    19043,
    20080,
    '-',
    'SCB29_ompA.fasta'
)

# SCB37
extract_sequences(
    'SCB37_scaffolds.fasta',
    'JBPXEM010000002.1',
    151505,
    152542,
    '+',
    'SCB37_ompA.fasta'
)

