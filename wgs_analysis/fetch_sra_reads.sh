#!/usr/bin/bash
# --- fetch_sra_reads.sh ---
# Download the deposited read files from SRA (BioProject PRJNA1510228) and put
# them where the analysis scripts in this repository expect to find them.
#
# The files are deposited under public isolate names (Nissle_Stock-1,
# Nissle_Preg-5, ...) while the analysis scripts refer to the internal
# sequencing-run IDs (N2SKTQ_1_1, Z6M7Y5_1, ...). SRA_file_map.tsv in the
# repository root carries that correspondence; this script reads it and writes
# the internal names directly.
#
# Layout produced (relative to the directory this is run from):
#   N2SKTQ_results/N2SKTQ_<i>_<i>_illumina_R1.fastq.gz, _R2, _nanopore
#   Z6M7Y5_illumina_prepared/Z6M7Y5_<i>_R1.fastq.gz, _R2
#   Z6M7Y5_nanopore_prepared/Z6M7Y5_<i>_nanopore.fastq.gz
#
# This replaces prepare_Z6M7Y5_data.sh for SRA-derived data. That script globs
# on the vendor accession IDs embedded in the original filenames, which are not
# part of the public deposit, so it cannot match downloaded files.
#
# Requires the SRA Toolkit (prefetch, fasterq-dump) on PATH:
#   conda install -c bioconda sra-tools
#
# Usage: ./fetch_sra_reads.sh [threads]

set -euo pipefail

THREADS=${1:-8}
MAP="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/SRA_file_map.tsv"

if [ ! -f "$MAP" ]; then
    echo "ERROR: mapping file not found: $MAP" >&2
    exit 1
fi

for tool in prefetch fasterq-dump; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: $tool not on PATH. Install the SRA Toolkit:" >&2
        echo "  conda install -c bioconda sra-tools" >&2
        exit 1; }
done

# The run accessions are filled into SRA_file_map.tsv once NCBI processes the
# submission. Refuse to run against placeholders rather than emit confusing
# prefetch errors for every row.
if cut -f6 "$MAP" | tail -n +2 | grep -qx TBD; then
    echo "ERROR: SRA_file_map.tsv still lists run_accession=TBD." >&2
    echo "The SRR accessions are assigned when the submission finishes" >&2
    echo "processing; fill in column 6 before running this script." >&2
    exit 1
fi

echo "=========================================="
echo "Fetching PRJNA1510228 reads from SRA"
echo "=========================================="
echo "Mapping:  $MAP"
echo "Threads:  $THREADS"
echo ""

mkdir -p N2SKTQ_results Z6M7Y5_illumina_prepared Z6M7Y5_nanopore_prepared sra_tmp

n_done=0
n_skip=0

# Columns: isolate internal_run_id platform biosample library_id run_accession
#          sra_filename_1 sra_filename_2 analysis_path_1 analysis_path_2
#
# Single-end rows carry "-" in the second-file columns rather than an empty
# field. Tab is an IFS whitespace character, so `read` collapses runs of tabs
# and an empty column would shift every later value one position left.
while IFS=$'\t' read -r isolate run_id platform biosample lib srr f1 f2 a1 a2; do
    [ -n "$a1" ] && [ -n "$srr" ] || {
        echo "ERROR: malformed row in $MAP for '${run_id:-?}' (${platform:-?})" >&2
        echo "Expected 10 tab-separated columns with '-' for absent values." >&2
        exit 1; }

    # Already present from an earlier run: leave it alone.
    if [ -f "$a1" ] && { [ "$a2" = "-" ] || [ -f "$a2" ]; }; then
        echo "  = $run_id ($platform) - already present, skipping"
        n_skip=$((n_skip + 1))
        continue
    fi

    echo ""
    echo "--- $isolate / $run_id ($platform) - $srr ---"

    prefetch --max-size u --output-directory sra_tmp "$srr"
    fasterq-dump --threads "$THREADS" --split-3 --outdir sra_tmp "sra_tmp/$srr"

    if [ "$platform" = "illumina" ]; then
        # --split-3 writes ${SRR}_1.fastq / ${SRR}_2.fastq for paired runs.
        [ -f "sra_tmp/${srr}_1.fastq" ] && [ -f "sra_tmp/${srr}_2.fastq" ] || {
            echo "ERROR: expected paired output for $srr" >&2; exit 1; }
        gzip -c "sra_tmp/${srr}_1.fastq" > "$a1"
        gzip -c "sra_tmp/${srr}_2.fastq" > "$a2"
        rm -f "sra_tmp/${srr}_1.fastq" "sra_tmp/${srr}_2.fastq"
        echo "  -> $a1"
        echo "  -> $a2"
    else
        # Single-end nanopore run: --split-3 writes ${SRR}.fastq.
        [ -f "sra_tmp/${srr}.fastq" ] || {
            echo "ERROR: expected single-end output for $srr" >&2; exit 1; }
        gzip -c "sra_tmp/${srr}.fastq" > "$a1"
        rm -f "sra_tmp/${srr}.fastq"
        echo "  -> $a1"
    fi

    rm -rf "sra_tmp/$srr"
    n_done=$((n_done + 1))
done < <(tail -n +2 "$MAP")

rmdir sra_tmp 2>/dev/null || true

echo ""
echo "=========================================="
echo "Download complete"
echo "=========================================="
echo "Fetched: $n_done runs"
echo "Skipped: $n_skip runs (already present)"
echo ""
echo "NOTE: this deposit covers the 7 isolates analysed in the paper."
echo "Scripts that iterate over the full sequencing batches will report the"
echo "other samples as missing and skip them; see the 'Reproducing from SRA'"
echo "section of README.md for which analyses that affects."
echo ""
echo "Next steps:"
echo "  Illumina/Nanopore hybrid assemblies (N2SKTQ): ./run_parallel_assemblies.sh"
echo "  Hybrid assemblies (Z6M7Y5):  ./run_Z6M7Y5_hybrid_assemblies.sh $THREADS"
