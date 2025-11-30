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

## Analysis Tools

Typical RNAseq workflow for this data:
- **QC**: FastQC, MultiQC
- **Alignment**: STAR or Bowtie2 (against PCC 11901 genome)
- **Quantification**: featureCounts or RSEM (using genomic.gtf)
- **Differential expression**: DESeq2 or edgeR (R)
