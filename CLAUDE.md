# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

RNAseq transcriptomics data repository for *Picosynechococcus sp. PCC 11901*, a fast-growing cyanobacterium. Part of the Australian Genome Foundry's research infrastructure.

## Data Structure

```
Data/                           # Raw sequencing data (~131 GB)
  ├── *.fastq.gz               # 240 FASTQ files (60 samples)
  └── AMO17076_md5sum.txt      # MD5 checksums for validation

PCC_11901_annotated genome/    # Reference genome
  ├── GCF_005577135.1_ASM557713v1_genomic.fna   # NCBI RefSeq genome
  ├── PCC11901 genome_fasta.fasta               # Annotated FASTA
  └── genomic.gtf                               # Gene annotations

RNAseq_data_labels.xlsx        # Sample metadata and experimental design
```

## FASTQ File Naming Convention

Pattern: `U[sample]-AMO17076A[id]-22FLL5LT1_S[num]_L[lane]_R[read]_001.fastq.gz`

- **U[sample]**: Sample number (U1-U60)
- **L001/L002**: Sequencing lane
- **R1/R2**: Read pair (paired-end sequencing)
- Each sample has 4 files: 2 lanes × 2 read pairs

## Experimental Design

20 conditions × 3 biological replicates = 60 samples (U1-U60)

### Group 1 - Nutrients (Control: U4,5,6 - MAD media)

| Treatment | Samples |
|-----------|---------|
| Control (MAD media) | U4, U5, U6 |
| 0.75% glycerol | U46, U47, U48 |
| Low Nitrogen | U22, U23, U24 |
| High Nitrogen | U25, U26, U27 |
| Low Phosphate | U28, U29, U30 |
| High Phosphate | U31, U32, U33 |
| Ammonia | U40, U41, U42 |
| Urea | U43, U44, U45 |

### Group 2 - Environmental (Control: U1,2,3 - MAD media)

| Treatment | Samples |
|-----------|---------|
| Control (MAD media) | U1, U2, U3 |
| 9% NaCl | U34, U35, U36 |
| 0.005% H2O2 | U37, U38, U39 |
| Atmospheric CO2 | U7, U8, U9 |
| High CO2 8% | U10, U11, U12 |
| High Temp 38C | U13, U14, U15 |
| Low light 15uE | U16, U17, U18 |
| High light | U19, U20, U21 |

### Group 3 - Circadian Rhythm

| Timepoint | Samples |
|-----------|---------|
| T1 (Light) | U49, U50, U51 |
| T2 (Dark) | U52, U53, U54 |
| T3 (Light) | U55, U56, U57 |
| T4 (Dark) | U58, U59, U60 |

## Sample-to-File Mapping

All FASTQ files are in `Data/`. Each sample has 4 files (2 lanes × 2 reads).

| Group | Treatment | Sample | FASTQ Files (in Data/) |
|-------|-----------|--------|------------------------|
| 1 | Control (MAD) | U4 | `U4-AMO17076A4-22FLL5LT1_S4_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | Control (MAD) | U5 | `U5-AMO17076A5-22FLL5LT1_S5_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | Control (MAD) | U6 | `U6-AMO17076A6-22FLL5LT1_S6_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | 0.75% glycerol | U46 | `U46-AMO17076A46-22FLL5LT1_S46_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | 0.75% glycerol | U47 | `U47-AMO17076A47-22FLL5LT1_S47_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | 0.75% glycerol | U48 | `U48-AMO17076A48-22FLL5LT1_S48_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | Low Nitrogen | U22 | `U22-AMO17076A22-22FLL5LT1_S22_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | Low Nitrogen | U23 | `U23-AMO17076A23-22FLL5LT1_S23_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | Low Nitrogen | U24 | `U24-AMO17076A24-22FLL5LT1_S24_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | High Nitrogen | U25 | `U25-AMO17076A25-22FLL5LT1_S25_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | High Nitrogen | U26 | `U26-AMO17076A26-22FLL5LT1_S26_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | High Nitrogen | U27 | `U27-AMO17076A27-22FLL5LT1_S27_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | Low Phosphate | U28 | `U28-AMO17076A28-22FLL5LT1_S28_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | Low Phosphate | U29 | `U29-AMO17076A29-22FLL5LT1_S29_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | Low Phosphate | U30 | `U30-AMO17076A30-22FLL5LT1_S30_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | High Phosphate | U31 | `U31-AMO17076A31-22FLL5LT1_S31_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | High Phosphate | U32 | `U32-AMO17076A32-22FLL5LT1_S32_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | High Phosphate | U33 | `U33-AMO17076A33-22FLL5LT1_S33_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | Ammonia | U40 | `U40-AMO17076A40-22FLL5LT1_S40_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | Ammonia | U41 | `U41-AMO17076A41-22FLL5LT1_S41_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | Ammonia | U42 | `U42-AMO17076A42-22FLL5LT1_S42_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | Urea | U43 | `U43-AMO17076A43-22FLL5LT1_S43_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | Urea | U44 | `U44-AMO17076A44-22FLL5LT1_S44_L00[1-2]_R[1-2]_001.fastq.gz` |
| 1 | Urea | U45 | `U45-AMO17076A45-22FLL5LT1_S45_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | Control (MAD) | U1 | `U1-AMO17076A1-22FLL5LT1_S1_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | Control (MAD) | U2 | `U2-AMO17076A2-22FLL5LT1_S2_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | Control (MAD) | U3 | `U3-AMO17076A3-22FLL5LT1_S3_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | 9% NaCl | U34 | `U34-AMO17076A34-22FLL5LT1_S34_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | 9% NaCl | U35 | `U35-AMO17076A35-22FLL5LT1_S35_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | 9% NaCl | U36 | `U36-AMO17076A36-22FLL5LT1_S36_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | 0.005% H2O2 | U37 | `U37-AMO17076A37-22FLL5LT1_S37_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | 0.005% H2O2 | U38 | `U38-AMO17076A38-22FLL5LT1_S38_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | 0.005% H2O2 | U39 | `U39-AMO17076A39-22FLL5LT1_S39_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | Atmospheric CO2 | U7 | `U7-AMO17076A7-22FLL5LT1_S7_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | Atmospheric CO2 | U8 | `U8-AMO17076A8-22FLL5LT1_S8_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | Atmospheric CO2 | U9 | `U9-AMO17076A9-22FLL5LT1_S9_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | High CO2 8% | U10 | `U10-AMO17076A10-22FLL5LT1_S10_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | High CO2 8% | U11 | `U11-AMO17076A11-22FLL5LT1_S11_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | High CO2 8% | U12 | `U12-AMO17076A12-22FLL5LT1_S12_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | High Temp 38C | U13 | `U13-AMO17076A13-22FLL5LT1_S13_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | High Temp 38C | U14 | `U14-AMO17076A14-22FLL5LT1_S14_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | High Temp 38C | U15 | `U15-AMO17076A15-22FLL5LT1_S15_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | Low light 15uE | U16 | `U16-AMO17076A16-22FLL5LT1_S16_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | Low light 15uE | U17 | `U17-AMO17076A17-22FLL5LT1_S17_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | Low light 15uE | U18 | `U18-AMO17076A18-22FLL5LT1_S18_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | High light | U19 | `U19-AMO17076A19-22FLL5LT1_S19_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | High light | U20 | `U20-AMO17076A20-22FLL5LT1_S20_L00[1-2]_R[1-2]_001.fastq.gz` |
| 2 | High light | U21 | `U21-AMO17076A21-22FLL5LT1_S21_L00[1-2]_R[1-2]_001.fastq.gz` |
| 3 | T1 (Light) | U49 | `U49-AMO17076A49-22FLL5LT1_S49_L00[1-2]_R[1-2]_001.fastq.gz` |
| 3 | T1 (Light) | U50 | `U50-AMO17076A50-22FLL5LT1_S50_L00[1-2]_R[1-2]_001.fastq.gz` |
| 3 | T1 (Light) | U51 | `U51-AMO17076A51-22FLL5LT1_S51_L00[1-2]_R[1-2]_001.fastq.gz` |
| 3 | T2 (Dark) | U52 | `U52-AMO17076A52-22FLL5LT1_S52_L00[1-2]_R[1-2]_001.fastq.gz` |
| 3 | T2 (Dark) | U53 | `U53-AMO17076A53-22FLL5LT1_S53_L00[1-2]_R[1-2]_001.fastq.gz` |
| 3 | T2 (Dark) | U54 | `U54-AMO17076A54-22FLL5LT1_S54_L00[1-2]_R[1-2]_001.fastq.gz` |
| 3 | T3 (Light) | U55 | `U55-AMO17076A55-22FLL5LT1_S55_L00[1-2]_R[1-2]_001.fastq.gz` |
| 3 | T3 (Light) | U56 | `U56-AMO17076A56-22FLL5LT1_S56_L00[1-2]_R[1-2]_001.fastq.gz` |
| 3 | T3 (Light) | U57 | `U57-AMO17076A57-22FLL5LT1_S57_L00[1-2]_R[1-2]_001.fastq.gz` |
| 3 | T4 (Dark) | U58 | `U58-AMO17076A58-22FLL5LT1_S58_L00[1-2]_R[1-2]_001.fastq.gz` |
| 3 | T4 (Dark) | U59 | `U59-AMO17076A59-22FLL5LT1_S59_L00[1-2]_R[1-2]_001.fastq.gz` |
| 3 | T4 (Dark) | U60 | `U60-AMO17076A60-22FLL5LT1_S60_L00[1-2]_R[1-2]_001.fastq.gz` |

## Reference Genome

- **Organism**: *Picosynechococcus sp. PCC 11901*
- **Assembly**: GCF_005577135.1 (ASM557713v1)
- **Annotation**: GTF format with RefSeq gene IDs, uses translation table 11 (bacterial)

## Data Validation

Verify FASTQ integrity after transfer:
```bash
cd Data/
md5sum -c AMO17076_md5sum.txt
```

## Proposed Folder Structure

```
.
├── Data/                              # Raw FASTQ files (excluded from git)
│   └── AMO17076_md5sum.txt           # Checksums (tracked)
│
├── PCC_11901_annotated genome/        # Reference genome & annotations
│   ├── GCF_005577135.1_ASM557713v1_genomic.fna
│   ├── PCC11901 genome_fasta.fasta
│   └── genomic.gtf
│
├── 01_QC/                             # Quality control outputs
│   ├── fastqc_raw/                   # FastQC on raw reads
│   ├── fastqc_trimmed/               # FastQC after trimming
│   └── multiqc_reports/              # MultiQC summary reports
│
├── 02_trimmed/                        # Trimmed FASTQ files (fastp output)
│   └── *.fastq.gz
│
├── 03_merged/                         # Lane-merged FASTQ files
│   └── {sample}_R[1-2].fastq.gz      # 120 files (60 samples × 2 reads)
│
├── 04_rRNA_filtered/                  # SortMeRNA output
│   ├── non_rRNA/                     # Cleaned reads for downstream
│   └── rRNA/                         # Removed rRNA reads (optional keep)
│
├── 05_salmon/                         # Path A: Quantification
│   ├── index/                        # Salmon transcriptome index
│   └── quants/                       # Per-sample quantification
│       └── {sample}/
│           ├── quant.sf             # Transcript-level counts
│           └── quant.genes.sf       # Gene-level counts
│
├── 06_alignment/                      # Path B: Visualization
│   ├── bowtie2_index/                # Bowtie2 genome index
│   ├── bam/                          # Sorted BAM files (can delete after bigWig)
│   └── bigwig/                       # CPM-normalized coverage tracks
│       └── {sample}.bw
│
├── 07_deseq2/                         # Differential expression results
│   ├── counts_matrix.csv             # Combined count matrix
│   ├── normalized_counts.csv         # DESeq2 normalized counts
│   ├── group1_nutrients/             # Group 1 comparisons
│   ├── group2_environmental/         # Group 2 comparisons
│   └── group3_circadian/             # Time-series analysis
│
├── 08_functional/                     # Functional annotation & enrichment
│   ├── eggnog_annotations/           # eggNOG-mapper output
│   ├── go_enrichment/                # GO term analysis
│   └── kegg_pathways/                # KEGG pathway analysis
│
├── 09_figures/                        # Publication-ready figures
│   ├── qc_plots/                     # PCA, correlation heatmaps
│   ├── volcano_plots/                # Per-comparison volcano plots
│   ├── heatmaps/                     # Expression heatmaps
│   └── interactive/                  # Plotly HTML outputs
│
├── scripts/                           # Analysis scripts
│   ├── 01_qc_trim.sh                 # FastQC + fastp
│   ├── 02_merge_lanes.sh             # Lane merging
│   ├── 03_sortmerna.sh               # rRNA removal
│   ├── 04_salmon_quant.sh            # Salmon quantification
│   ├── 05_bowtie2_align.sh           # Alignment for visualization
│   ├── 06_bigwig.sh                  # deepTools bamCoverage
│   ├── 07_deseq2_analysis.R          # Differential expression
│   ├── 08_functional_analysis.R      # GO/KEGG enrichment
│   └── utils/                        # Helper functions
│
├── logs/                              # Pipeline logs
│   └── {step}_{sample}.log
│
├── environment.yml                    # Conda environment specification
├── RNAseq_data_labels.xlsx           # Sample metadata
├── CLAUDE.md                          # This file
└── README.md                          # Project documentation
```

## Analysis Pipeline

Based on Conesa et al. (2016) "A survey of best practices for RNA-seq data analysis" - Genome Biology

### Phase 1: Pre-Analysis

#### 1.1 Quality Control of Raw Reads
- **Tool**: FastQC + MultiQC
- **Check for**: Sequence quality, GC content, adapter contamination, overrepresented k-mers, duplicates
- **Output**: `01_QC/fastqc_raw/`, `01_QC/multiqc_reports/`

#### 1.2 Trimming
- **Tool**: fastp
- **Action**: Trim adapters and low-quality bases in one step
- **Output**: `02_trimmed/`

#### 1.3 Lane Merging
- Merge L001 and L002 files for each sample
- **Result**: 120 files (60 samples × 2 read pairs R1/R2)
- **Output**: `03_merged/`

#### 1.4 In Silico rRNA Removal (SortMeRNA)
- **Tool**: SortMeRNA
- **Why**: Bacterial samples often retain significant rRNA contamination
- **Databases**: SILVA and Rfam rRNA databases
- **Resources**: 16GB+ RAM per thread; longest step in pipeline
- **Output**: `04_rRNA_filtered/non_rRNA/`

### Phase 2: Core Analysis (Parallel Workflow)

Two parallel paths from rRNA-filtered FASTQs:

```
Path A: FASTQ → Salmon (quasi-mapping) → Counts → DESeq2
Path B: FASTQ → Bowtie2 → BAM → deepTools → bigWig
```

#### 2.1 Path A: Salmon Quantification (Primary)
- **Tool**: Salmon (alignment-free, quasi-mapping mode)
- **Advantages**:
  - Auto-infers library strandedness
  - Fast, lightweight, no large BAM files
  - Built-in GC and positional bias correction
- **Indexing**: Use **decoy-aware mapping** strategy
  - Target: transcriptome sequences (from GTF + genome FASTA)
  - Decoy: whole genome (prevents intergenic reads from false quantification)
  ```bash
  # Generate decoys.txt (chromosome names)
  grep "^>" genome.fna | cut -d " " -f 1 | sed 's/>//g' > decoys.txt
  # Concatenate transcriptome + genome
  cat transcriptome.fa genome.fna > gentrome.fa
  # Build index with decoys
  salmon index -t gentrome.fa -d decoys.txt -i salmon_index -p 8
  ```
- **Output**: `05_salmon/quants/{sample}/`

#### 2.2 Path B: Bowtie2 Alignment (For Visualization Only)
- **Tool**: Bowtie2 (ungapped aligner, appropriate for bacteria)
- **Reference**: PCC 11901 genome
- **Output**: `06_alignment/bam/` → `06_alignment/bigwig/`
- **Note**: BAM files can be deleted after bigWig creation

#### 2.3 Normalization & Differential Expression
- **Tool**: DESeq2 (R)
- **QC plots**: PCA, sample correlation heatmap, dispersion plots
- **Comparisons**:
  - Group 1: Each nutrient condition vs Control (U4,5,6)
  - Group 2: Each environmental condition vs Control (U1,2,3)
  - Group 3: Time-series analysis for circadian data
- **Thresholds**: |log2FC| > 1, adjusted p-value < 0.05
- **Output**: `07_deseq2/`

### Phase 3: Advanced Analysis

#### 3.1 Functional Profiling
- **Problem**: No OrgDb exists for Picosynechococcus PCC 11901
- **Solution**: Custom annotation via eggNOG-mapper
  1. Extract protein sequences from genome
  2. Run eggNOG-mapper (online or local)
  3. Use `clusterProfiler::enricher()` with custom TERM2GENE
- **Output**: `08_functional/`

#### 3.2 Visualization
- **Coverage tracks**: deepTools bamCoverage with CPM normalization
  ```bash
  bamCoverage --normalizeUsing CPM -b sample.bam -o sample.bw
  ```
- **Static plots**: pheatmap, ComplexHeatmap, ggplot2, EnhancedVolcano
- **Interactive plots**: Plotly for exploratory analysis
- **Output**: `09_figures/`

## Tools Summary

| Step | Tool | Notes |
|------|------|-------|
| QC | FastQC, MultiQC | Run before and after trimming |
| Trimming | fastp | Fast, handles adapters + quality in one step |
| rRNA removal | SortMeRNA | Essential for bacterial RNA-seq |
| Alignment | Bowtie2 | Simple for bacteria (no splicing) |
| Quantification | Salmon | Auto-infers strandedness, fast, bias-corrected |
| DE Analysis | DESeq2 | Good for 3 replicates, well-documented |
| bigWig | deepTools | CPM-normalized tracks for IGV |
| Functional | clusterProfiler | With custom eggNOG annotations |
| Visualization | Plotly | Interactive exploratory plots |

## Conda Environment

Activate with:
```bash
conda activate pcc11901_rnaseq
```

See `environment.yml` for full specification.

## Hardware Considerations

Optimized for Apple M2 Max (64GB RAM, 4E + 8P cores):
- **SortMeRNA**: Use 3-4 threads max (~15GB per thread, leaving RAM for OS overhead)
- **Salmon**: 8-12 threads
- **Bowtie2**: 8 threads
- **Python multiprocessing**: Use `joblib` with `n_jobs=8` for parallel sample processing
- **Progress bars**: tqdm for all batch operations

**Apple Silicon Note**: Some bioconda tools lack native ARM64 builds. If installation fails, use Rosetta emulation:
```bash
CONDA_SUBDIR=osx-64 conda env create -f environment.yml
```

## Key Considerations

1. **Bacterial RNA-seq**: rRNA depletion was likely used; in silico removal still needed
2. **3 replicates**: Minimum for DESeq2/edgeR statistical inference
3. **Multiple comparisons**: FDR correction across all tests
4. **Two control groups**: Groups 1 and 2 have separate controls - analyze separately
5. **Circadian data**: Time-series methods (maSigPro or DESeq2 LRT) for Group 3
