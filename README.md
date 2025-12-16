# *E. coli* Natural Immunity

Scripts supporting the publication: **"Natural maternal immunity protects neonates from *Escherichia coli* sepsis"**

Diep, R.E., et al. *Nature* (2025).

This repository contains bioinformatics scripts used for whole genome sequencing (WGS) analysis, OmpA protein analysis, peptide homology studies, and statistical analysis of maternal antibody protection against neonatal *E. coli* sepsis.

**Repository:** https://github.com/Way-Lab/ecoli_natural_immunity

---

## System Requirements

### Operating System
- **Tested on:** Ubuntu 24.04 LTS (Linux)
- **Compatible with:** Linux distributions (Debian, Ubuntu, CentOS), macOS (with Homebrew/conda)
- **Note:** Windows users should use WSL2 (Windows Subsystem for Linux)

### Hardware Requirements
- **Minimum:** 8 GB RAM, 4 CPU cores
- **Recommended for WGS analysis:** 32 GB RAM, 16 CPU cores
- **Storage:** ~10 GB for software dependencies; additional space for sequence data
- No specialized hardware (GPU) required

### Software Dependencies

#### Python (version 3.10+)
Tested with Python 3.12.10

| Package | Version Tested | Purpose |
|---------|----------------|---------|
| pandas | 2.3.3 | Data manipulation |
| numpy | 2.3.5 | Numerical computing |
| matplotlib | 3.10.3 | Visualization |
| seaborn | 0.13.2 | Statistical graphics |
| scipy | 1.16.3 | Scientific computing |
| statsmodels | 0.14.5 | Statistical models |
| biopython | 1.85 | Sequence analysis |

#### Bioinformatics Tools (for WGS analysis)
| Tool | Purpose | Installation |
|------|---------|--------------|
| [Unicycler](https://github.com/rrwick/Unicycler) | Hybrid genome assembly | conda install -c bioconda unicycler |
| [QUAST](https://github.com/ablab/quast) | Assembly QC | conda install -c bioconda quast |
| [BUSCO](https://busco.ezlab.org/) | Completeness assessment | conda install -c bioconda busco |
| [Parsnp](https://harvest.readthedocs.io/) | Core genome alignment | conda install -c bioconda parsnp |
| [Breseq](https://barricklab.org/breseq) | Mutation identification | conda install -c bioconda breseq |
| [Bakta](https://github.com/oschwengers/bakta) | Genome annotation | conda install -c bioconda bakta |
| [MOB-suite](https://github.com/phac-nml/mob-suite) | Plasmid analysis | conda install -c bioconda mob_suite |
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | Sequence alignment | conda install -c bioconda mafft |

---

## Installation Guide

### Step 1: Clone the repository
```bash
git clone https://github.com/Way-Lab/ecoli_natural_immunity.git
cd ecoli_natural_immunity
```

### Step 2: Create Python environment
Using conda (recommended):
```bash
conda create -n ecoli_immunity python=3.12
conda activate ecoli_immunity
pip install pandas numpy matplotlib seaborn scipy statsmodels biopython
```

Or using pip with virtual environment:
```bash
python3 -m venv venv
source venv/bin/activate
pip install pandas numpy matplotlib seaborn scipy statsmodels biopython
```

### Step 3: Install bioinformatics tools (optional, for WGS analysis)
```bash
conda install -c bioconda unicycler quast busco parsnp breseq bakta mob_suite mafft snp-dists
```

### Typical Installation Time
- **Python environment only:** 2-5 minutes
- **Full installation with bioinformatics tools:** 15-30 minutes

---

## Demo

### Demo 1: Statistical Analysis (Primary Demo)

This demo runs the matched case-control analysis of maternal antibodies.

**Data:** Embedded in script (no external files needed)

**Run command:**
```bash
cd statistical_analysis
python logistic_regression_case_control.py --assay OmpA
```

**Expected output:**
```
======================================================================
MATCHED CASE-CONTROL ANTIBODY ANALYSIS SUMMARY - OmpA
======================================================================

SAMPLE SIZE:
  Cases: 100
  Controls: 296
  Matched sets: 100

DESCRIPTIVE STATISTICS:
  Cases - Geometric mean: 265.4
  Controls - Geometric mean: 625.4
  Geometric mean ratio: 2.356

CONDITIONAL LOGISTIC REGRESSION:
  Coefficient: -1.972 (SE: 0.327)
  95% CI: [-2.613, -1.331]
  P-value: 1.6528e-09

ODDS RATIOS:
  Per 10-fold increase: 0.139
  Per 2-fold increase: 0.552

PROTECTIVE THRESHOLDS:
  50% risk reduction: 1346.0
  75% risk reduction: 3320.6
  90% risk reduction: 10090.6
  95% risk reduction: 10327.0
```

**Output files generated:**
  - OmpA_antibody_analysis_figure.pdf
  - OmpA_analysis_summary.txt

**Expected runtime:** < 1 minute on a standard desktop


### Demo 2: OmpA Sequence Analysis

**Run command:**
```bash
cd ompa_analysis
python analyze_selected_ompA_strains.py ../demo_data/selected_ompA_sequences.fasta
```

**Expected output:** Aligned OmpA sequences with peptide loop annotations

**Expected runtime:** < 30 seconds

### Demo 3: Peptide Variability Analysis

**Run command:**
```bash
cd peptide_analysis
python analyze_peptide_variability.py ../demo_data/OmpA_peptide_seqeunces.fasta
```

**Expected output:** Conservation analysis across peptide sequences

**Expected runtime:** < 30 seconds

---

## Instructions for Use

### Running on Your Own Data

#### Statistical Analysis
The statistical analysis script analyzes matched case-control data. To use with your own data:

1. Format your data as tab-separated values with columns: Case, Control1, Control2, Control3
2. Modify the data strings in `logistic_regression_case_control.py` or adapt the script to read from CSV files

#### OmpA/Peptide Analysis
For sequence analysis with your own FASTA files:

```bash
# Extract OmpA sequences from genome assemblies
python ompa_analysis/extract_ompA.py your_assembly.fasta

# Compare OmpA sequences between strains
python ompa_analysis/compare_ompA.py strain1.fasta strain2.fasta

# Analyze peptide variability
python peptide_analysis/analyze_peptide_variability.py your_sequences.fasta
```

#### WGS Analysis Pipeline
For whole genome sequencing analysis:

```bash
# Hybrid assembly (requires Illumina + Nanopore reads)
./wgs_analysis/process_hybrid_assembly.sh sample_prefix num_threads

# Assembly quality control
./wgs_analysis/run_assembly_qc.sh assembly.fasta

# Core genome phylogenetic analysis
./wgs_analysis/run_parsnp_analysis.sh genome_directory reference.fasta

# Mutation identification
./wgs_analysis/run_breseq_analysis.sh reads.fastq reference.gbk
```

---

## Repository Structure

```
ecoli_natural_immunity/
├── README.md                     # This file
├── LICENSE                       # MIT License
├── demo_data/                    # Demo datasets for testing
│   ├── RS218_ompA.fasta         # Reference OmpA sequence
│   ├── selected_ompA_sequences.fasta
│   ├── combined_ompA_all_species.fasta
│   ├── aligned_sequences_with_SCB_authentic.fasta
│   ├── OmpA_peptide_seqeunces.fasta
│   ├── unaligned_for_mafft.fasta
│   ├── ompA_comparison.fasta
│   └── rpoS_comparison.fasta
├── wgs_analysis/                 # Whole genome sequencing pipelines
├── snp_analysis/                 # SNP identification and annotation
├── ompa_analysis/                # OmpA protein analysis
├── peptide_analysis/             # Peptide conservation studies
├── statistical_analysis/         # Case-control statistical analysis
└── utilities/                    # FASTA processing and visualization
```

## Directory Contents

### wgs_analysis/

Shell scripts for hybrid assembly of E coli stock and colonizing strains

| Script | Description |
|--------|-------------|
| `process_hybrid_assembly.sh` | Unicycler hybrid assembly combining Illumina short reads and Nanopore long reads |
| `run_assembly_qc.sh` | Assembly quality control using QUAST and BUSCO |
| `run_parallel_assemblies.sh` | Parallel execution of genome assemblies |
| `run_parsnp_analysis.sh` | Parsnp core genome alignment and phylogenetic tree construction |
| `run_parsnp_with_outgroups.sh` | Parsnp analysis including outgroup reference strains |
| `run_parsnp_without_duplicate.sh` | Parsnp variant excluding duplicate samples |
| `run_comprehensive_parsnp.sh` | Comprehensive parsnp run with full strain collection |
| `run_breseq_analysis.sh` | Breseq analysis for identifying mutations in variant strains |
| `run_all_breseq.sh` | Parallel breseq execution across multiple samples |
| `run_comprehensive_plasmid_analysis.sh` | MOB-suite + Bakta plasmid characterization |
| `run_plasmid_analysis_parallel.sh` | Parallel plasmid analysis processing |
| `run_bakta_plasmids_only.sh` | Bakta annotation for plasmid sequences |
| `prepare_Z6M7Y5_data.sh` | Data preparation for Z6M7Y5 sample |
| `run_Z6M7Y5_assembly_qc.sh` | Assembly QC for Z6M7Y5 |
| `run_Z6M7Y5_complete_workflow.sh` | Complete analysis workflow for Z6M7Y5 |
| `run_Z6M7Y5_hybrid_assemblies.sh` | Hybrid assembly for Z6M7Y5 |
| `run_Z6M7Y5_comparative_analysis.sh` | Comparative analysis of Z6M7Y5 |
| `run_Z6M7Y5_comparative_analysis_with_refs.sh` | Comparative analysis with additional references |

### snp_analysis/

Python scripts for SNP identification, annotation, and visualization between stock and colonizing E. coli isolates

| Script | Description |
|--------|-------------|
| `analyze_comprehensive_parsnp_snps.py` | Analyze comprehensive parsnp SNP data |
| `analyze_filtered_parsnp_snps.py` | Analyze filtered parsnp SNPs with annotations |
| `annotate_snps_with_genes.py` | Annotate SNPs with gene information from Bakta |
| `annotate_filtered_snps_with_genes.py` | Annotate filtered SNPs with functional impact |
| `annotate_snp_effects.py` | Determine amino acid changes from SNPs |
| `create_comprehensive_parsnp_heatmap.py` | Create SNP distance heatmaps |
| `create_ecn_reference_heatmap.py` | Create heatmap using EcN (Nissle) reference |
| `create_filtered_heatmap.py` | Create heatmap from filtered SNP data |
| `create_snp_summary_by_sample.py` | SNP summary statistics by sample |
| `create_filtered_snp_summary.py` | Summary of filtered SNP data |
| `create_wide_snp_table.py` | Wide-format SNP table (nucleotide per strain) |
| `create_nissle_snp_heatmap_from_mapping.py` | Heatmap from Nissle-based mapping |
| `create_nissle_snp_heatmap_from_parsnp_matrix.py` | Heatmap from parsnp distance matrix |
| `extract_nissle_snp_distances.py` | Extract SNP distances from Nissle reference |
| `setup_comprehensive_parsnp_analysis.py` | Configure comprehensive parsnp analysis |
| `setup_parsnp_with_outgroups.py` | Configure parsnp with outgroups |
| `rename_parsnp_files.py` | Rename parsnp output files |
| `analyze_rpoS_gap.py` | Analyze rpoS gene gaps/knockouts |
| `compare_rpoS.py` | Compare rpoS gene sequences |
| `all_samples_summary.py` | Generate summary across all samples |

### ompa_analysis/

Python scripts for OmpA outer membrane protein sequence analysis.

| Script | Description |
|--------|-------------|
| `compare_ompA.py` | Extract and compare ompA gene sequences between strains |
| `extract_ompA.py` | Extract OmpA sequences from scaffold files |
| `analyze_ompA_knockout.py` | Analyze ompA knockout regions and resistance markers |
| `analyze_new_ompA_with_loops.py` | Analyze ompA sequences and compare peptide loop regions |
| `analyze_selected_ompA_strains.py` | Filter OmpA by strain and extract peptide loops |
| `extract_SCB60_SCB61_ompA.py` | Extract ompA from SCB60 and SCB61 strains |
| `clean_ompA_headers.py` | Clean and standardize OmpA FASTA headers |
| `format_ompA_fasta.py` | Format OmpA FASTA files |
| `linearize_SCB60_SCB61.py` | Linearize wrapped OmpA sequences |
| `translate_ompA.py` | Translate ompA nucleotide to protein sequences |
| `download_enterobacteriaceae_ompA.py` | Download OmpA sequences from Enterobacteriaceae via NCBI |
| `download_specific_ompA_strains.py` | Download OmpA from specific bacterial strains |
| `plot_cladogram.py` | Plot phylogenetic cladograms from OmpA alignments |

### peptide_analysis/

Python scripts for peptide sequence analysis and conservation across strains and species.

| Script | Description |
|--------|-------------|
| `analyze_peptide_variability.py` | Analyze variability of peptide sequences within OmpA |
| `generate_peptide_alignments.py` | Generate multi-FASTA alignments for peptide variants |
| `generate_peptide_alignments_with_references.py` | Generate alignments with RS218 reference (dots for matches) |
| `generate_peptide_trees.py` | Generate phylogenetic trees for peptide sequences |
| `generate_full_alignments.py` | Generate full alignments from 130+ sequences |
| `analyze_species_variants.py` | Analyze peptide variants across bacterial species |
| `generate_html_report.py` | Generate HTML reports with embedded graphics |
| `generate_final_summary.py` | Generate final summary report |

### statistical_analysis/

Python scripts for statistical analysis of maternal antibody levels and neonatal sepsis risk.

| Script | Description |
|--------|-------------|
| `logistic_regression_case_control.py` | Matched case-control analysis of maternal antibodies using conditional logistic regression. Analyzes THP-1 opsonization, HL60 opsonization, Mix8 IgG, and OmpA IgG assays. Calculates protective thresholds with bootstrap confidence intervals. |

### utilities/

General-purpose scripts for FASTA processing and visualization.

| Script | Description |
|--------|-------------|
| `clean_headers.py` | Clean FASTA sequence headers |
| `clean_headers_final.py` | Final version of header cleaning |
| `fix_duplicate_headers.py` | Fix duplicate headers in FASTA files |
| `replace_headers.py` | Replace headers with standardized names |
| `format_fasta.py` | Format FASTA files with proper structure |
| `linearize_fasta.py` | Convert wrapped sequences to single-line |
| `make_unique.py` | Remove duplicate sequences |
| `process_fasta.py` | General FASTA processing |
| `remove_sequences.py` | Remove specific sequences from FASTA |
| `remove_sequences_from_processed.py` | Remove sequences from processed files |


## Citation

> Diep, R.E., et al. Natural maternal immunity protects neonates from Escherichia coli sepsis. *Nature* (2025).

## Software Citations

- De Coster, W., D'Hert, S., Schultz, D. T., Cruts, M. & Van Broeckhoven, C. NanoPack: visualizing and processing long-read sequencing data. Bioinformatics 34, 2666–2669 (2018). https://doi.org/10.1093/bioinformatics/bty149
- Bankevich, A. et al. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. J Comput Biol 19, 455–477 (2012). https://doi.org/10.1089/cmb.2012.0021
- Vaser, R., Sovic, I., Nagarajan, N. & Sikic, M. Fast and accurate de novo genome assembly from long uncorrected reads. Genome Res 27, 737–746 (2017). https://doi.org/10.1101/gr.214270.116
- Wick, R. R., Judd, L. M., Gorrie, C. L. & Holt, K. E. Unicycler: Resolving bacterial genome assemblies from short and long sequencing reads. PLoS Comput Biol 13, e1005595 (2017). https://doi.org/10.1371/journal.pcbi.1005595
- Gurevich, A., Saveliev, V., Vyahhi, N. & Tesler, G. QUAST: quality assessment tool for genome assemblies. Bioinformatics 29, 1072–1075 (2013). https://doi.org/10.1093/bioinformatics/btt086
- Katoh, K. & Standley, D. M. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol Biol Evol 30, 772–780 (2013). https://doi.org/10.1093/molbev/mst010
- Deatherage, D.E., Barrick, J.E. Identification of mutations in laboratory-evolved microbes from next-generation sequencing data using breseq. Methods Mol. Biol. 1151: 165–188 (2014).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
