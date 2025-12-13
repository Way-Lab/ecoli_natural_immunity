# E. coli Natural Immunity Analysis Scripts

Scripts supporting the publication: **"Natural maternal immunity protects neonates from Escherichia coli sepsis"**

This repository contains bioinformatics scripts used for whole genome sequencing (WGS) analysis, OmpA protein analysis, and peptide conservation studies of *Escherichia coli* strains.

## Repository Structure

```
ecoli_natural_immunity/
├── wgs_analysis/       # Whole genome sequencing and assembly pipelines
├── ompa_analysis/      # OmpA outer membrane protein analysis
├── peptide_analysis/   # Peptide conservation and variability analysis
├── snp_analysis/       # SNP calling and annotation
└── utilities/          # FASTA processing and visualization tools
```

## Directory Contents

### wgs_analysis/

Shell scripts for E. coli whole genome sequencing analysis pipelines.

| Script | Description |
|--------|-------------|
| `process_hybrid_assembly.sh` | Unicycler hybrid assembly combining Illumina short reads and Nanopore long reads |
| `run_assembly_qc.sh` | Assembly quality control using QUAST and BUSCO |
| `run_parallel_assemblies.sh` | Parallel execution of genome assemblies |
| `run_parsnp_analysis.sh` | Parsnp core genome alignment and phylogenetic tree construction |
| `run_parsnp_with_outgroups.sh` | Parsnp analysis including outgroup reference strains |
| `run_parsnp_without_duplicate.sh` | Parsnp variant excluding duplicate samples |
| `run_comprehensive_parsnp.sh` | Comprehensive parsnp run with full strain collection |
| `run_snippy_core_analysis.sh` | Snippy core genome SNP analysis pipeline |
| `run_snippy_nissle_refseq.sh` | Snippy analysis against Nissle 1917 RefSeq reference |
| `run_breseq_analysis.sh` | Breseq analysis for identifying mutations in variant strains |
| `run_all_breseq.sh` | Parallel breseq execution across multiple samples |
| `run_comprehensive_plasmid_analysis.sh` | MOB-suite + Bakta plasmid characterization |
| `run_plasmid_analysis_parallel.sh` | Parallel plasmid analysis processing |
| `run_bakta_plasmids_only.sh` | Bakta annotation for plasmid sequences |
| `prepare_haslam_illumina.sh` | Data preparation for Haslam Illumina sequencing |
| `prepare_haslam_nanopore.sh` | Data preparation for Haslam Nanopore sequencing |
| `run_haslam_hybrid_assemblies.sh` | Hybrid assemblies for Haslam samples |
| `run_haslam_comparative_analysis.sh` | Comparative analysis of Haslam samples |
| `prepare_Z6M7Y5_data.sh` | Data preparation for Z6M7Y5 sample |
| `run_Z6M7Y5_assembly_qc.sh` | Assembly QC for Z6M7Y5 |
| `run_Z6M7Y5_complete_workflow.sh` | Complete analysis workflow for Z6M7Y5 |
| `run_Z6M7Y5_hybrid_assemblies.sh` | Hybrid assembly for Z6M7Y5 |
| `run_Z6M7Y5_comparative_analysis.sh` | Comparative analysis of Z6M7Y5 |
| `run_Z6M7Y5_comparative_analysis_with_refs.sh` | Comparative analysis with additional references |

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
| `visualize_peptide_conservation.py` | Create conservation heatmaps across Enterobacteriaceae |

### peptide_analysis/

Python scripts for peptide sequence analysis and conservation studies.

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
| `create_lpp_alignment_figure.py` | Create Murein lipoprotein (Lpp) alignment figure |
| `extract_lpp_sequences.py` | Extract Lpp (Braun's lipoprotein) sequences |
| `compare_lpp_multispecies.py` | Compare Lpp across Enterobacteriaceae species |

### snp_analysis/

Python scripts for SNP identification, annotation, and visualization.

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

### utilities/

General-purpose scripts for FASTA processing and visualization.

| Script | Description |
|--------|-------------|
| `align_and_visualize.py` | Multiple sequence alignment and visualization |
| `align_and_visualize_improved.py` | Improved alignment/visualization pipeline |
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
| `visualize_tree.py` | Visualize phylogenetic trees |
| `visualize_final_tree.py` | Visualize final phylogenetic tree |
| `visualize_strains_tree.py` | Visualize strain-specific trees |
| `extract_mlp_comprehensive.py` | Extract Murein lipoprotein via tBLASTn |
| `compare_mlp_all_strains.py` | Compare Murein lipoprotein across strains |

## Software Dependencies

### Assembly and QC
- [Unicycler](https://github.com/rrwick/Unicycler) - Hybrid genome assembly
- [QUAST](https://github.com/ablab/quast) - Assembly quality assessment
- [BUSCO](https://busco.ezlab.org/) - Assembly completeness assessment

### Comparative Genomics
- [Parsnp](https://harvest.readthedocs.io/) - Core genome alignment
- [Snippy](https://github.com/tseemann/snippy) - SNP calling
- [Breseq](https://barricklab.org/breseq) - Mutation identification

### Annotation
- [Bakta](https://github.com/oschwengers/bakta) - Genome annotation
- [MOB-suite](https://github.com/phac-nml/mob-suite) - Plasmid analysis

### Python Libraries
- Biopython
- pandas
- numpy
- matplotlib
- seaborn

## Citation

If you use these scripts, please cite:

> [Authors]. Natural maternal immunity protects neonates from Escherichia coli sepsis. [Journal]. [Year].

## License

These scripts are provided for reproducibility of the published research.
