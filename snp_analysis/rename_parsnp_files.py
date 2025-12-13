#!/usr/bin/env python3
"""
Script to rename parsnp output files using FileNameConversion.csv mapping.
Replaces sequence result file names with Our ID names.
"""

import csv
import re
import os
import shutil

def load_conversion_mapping(csv_file):
    """Load the name conversion mapping from CSV file."""
    mapping = {}
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Skip empty rows
            if not row['Our ID'] or not row['Sequence Result Files']:
                continue
            # Map from sequence result files to our ID
            seq_name = row['Sequence Result Files'].strip()
            our_id = row['Our ID'].strip()
            mapping[seq_name] = our_id
            # Also handle with .fasta extension and .fasta.ref
            mapping[f"{seq_name}.fasta"] = f"{our_id}.fasta"
            mapping[f"{seq_name}.fasta.ref"] = f"{our_id}.fasta.ref"
    return mapping

def replace_names_in_content(content, mapping):
    """Replace all occurrences of sequence names in content."""
    for old_name, new_name in mapping.items():
        # Replace exact matches (with word boundaries when possible)
        content = re.sub(r'\b' + re.escape(old_name) + r'\b', new_name, content)
    return content

def process_file(input_file, output_file, mapping):
    """Process a single file and replace sequence names."""
    try:
        # Try to read as text file
        with open(input_file, 'r', encoding='utf-8') as f:
            content = f.read()

        # Replace names
        new_content = replace_names_in_content(content, mapping)

        # Write to output file
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(new_content)
        return True
    except UnicodeDecodeError:
        # If it's a binary file, just copy it
        print(f"  Note: Binary file detected, copying without modification")
        shutil.copy2(input_file, output_file)
        return False

def main():
    # Load the conversion mapping
    csv_file = 'FileNameConversion.csv'
    mapping = load_conversion_mapping(csv_file)

    print(f"Loaded {len(mapping)} name mappings")
    print("\nMapping examples:")
    for old, new in list(mapping.items())[:5]:
        print(f"  {old} -> {new}")

    # Create output directory
    input_dir = 'parsnp_results'
    output_dir = 'renamed_parsnp_output'
    os.makedirs(output_dir, exist_ok=True)

    # Process all files in parsnp_results
    files_processed = 0
    files_renamed = 0
    for filename in os.listdir(input_dir):
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename)

        if os.path.isfile(input_path):
            print(f"\nProcessing: {filename}")
            try:
                is_text = process_file(input_path, output_path, mapping)
                files_processed += 1
                if is_text:
                    files_renamed += 1
                print(f"  ✓ Written to: {output_path}")
            except Exception as e:
                print(f"  ✗ Error processing {filename}: {e}")

    print(f"\n{'='*60}")
    print(f"Completed! Processed {files_processed} files.")
    print(f"  - Text files with names replaced: {files_renamed}")
    print(f"  - Binary files copied: {files_processed - files_renamed}")
    print(f"Output directory: {output_dir}")

if __name__ == '__main__':
    main()
