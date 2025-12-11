# RNA-seq Analysis of *Picosynechococcus sp.* PCC 11901

Transcriptomic analysis of a fast-growing cyanobacterium under various environmental and nutrient stress conditions.

## Study Purpose

This project analyzes RNA-seq data from *Picosynechococcus sp.* PCC 11901 (formerly *Synechococcus sp.*) to identify differentially expressed genes under 20 experimental conditions across three experimental groups:

1. **Nutrient stress** (7 conditions): Glycerol, nitrogen levels, phosphate levels, ammonia, urea
2. **Environmental stress** (7 conditions): CO2 levels, salinity, oxidative stress (H2O2), temperature, light intensity
3. **Circadian rhythm** (4 timepoints): Light/dark cycles

Each condition has 3 biological replicates (60 samples total, 240 FASTQ files).

![Workflow Overview](workflow_overview.png)

---

## Environment Setup

### Prerequisites

- Conda/Mamba package manager
- ~64 GB RAM recommended (SortMeRNA is memory-intensive)
- ~500 GB disk space for intermediate files

### Installation

**Apple Silicon (M1/M2) users must use Rosetta emulation** due to missing ARM64 builds for some bioconda packages:

```bash
# Create environment in x86_64 mode
CONDA_SUBDIR=osx-64 conda env create -f environment.yml

# Activate
conda activate pcc11901_rnaseq

# Install additional R packages (run once)
R -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(c("clusterProfiler", "enrichplot", "GOSemSim", "goseq", "maSigPro"))'
```

### Update Environment

```bash
CONDA_SUBDIR=osx-64 conda env update -f environment.yml --prune
```

### Key Tools

| Tool | Purpose |
|------|---------|
| FastQC, MultiQC | Quality control |
| fastp | Read trimming |
| SortMeRNA | rRNA removal |
| gffread | Transcriptome extraction |
| Salmon | Transcript quantification |
| Bowtie2, samtools | Alignment (for visualization) |
| deepTools | Coverage tracks (bigWig) |
| DESeq2 (R) | Differential expression |
| gseapy (Python) | Pathway enrichment |
| Plotly | Interactive visualization |

---

## Folder Structure

```
.
├── Data/                              # Raw FASTQ files (~131 GB, not in git)
│   ├── *.fastq.gz                    # 240 files (60 samples × 2 lanes × 2 reads)
│   ├── AMO17076_md5sum.txt           # MD5 checksums
│   └── silva_db/                     # SILVA 138.1 rRNA database
│
├── PCC_11901_annotated genome/        # Reference genome
│   ├── GCF_005577135.1_*.fna         # NCBI RefSeq genome
│   ├── genomic.gtf                   # Gene annotations
│   ├── transcriptome.fa              # Extracted transcripts (generated)
│   ├── gentrome.fa                   # Transcriptome + genome (for Salmon)
│   └── decoys.txt                    # Chromosome names (for Salmon)
│
├── CyanoCycDB/                        # CyanoCyc pathway database
│   ├── All-genes-of-*.txt            # Gene annotations with GO terms
│   └── All-pathways-of-*.txt         # Pathway definitions
│
├── 01_QC/                             # Quality control outputs
│   ├── fastqc_raw/                   # FastQC on raw reads
│   ├── fastqc_trimmed/               # FastQC after trimming
│   └── multiqc_reports/              # Aggregated QC reports
│
├── 02_trimmed/                        # Trimmed reads (fastp output)
│
├── 03_merged/                         # Lane-merged reads
│   └── {sample}_R[1-2].fastq.gz      # 120 files (2 per sample)
│
├── 04_rRNA_filtered/                  # SortMeRNA output
│   ├── non_rRNA/                     # Clean reads for analysis
│   └── rRNA/                         # Removed rRNA reads
│
├── 05_salmon/                         # Salmon quantification
│   ├── index/                        # Salmon index
│   └── quants/{sample}/              # Per-sample counts
│       ├── quant.sf                  # Transcript-level counts
│       └── aux_info/meta_info.json   # Mapping statistics
│
├── 06_alignment/                      # Bowtie2 alignment (visualization only)
│   ├── bowtie2_index/                # Genome index
│   ├── bam/                          # Sorted BAM files
│   └── bigwig/                       # Coverage tracks for IGV
│
├── 07_deseq2/                         # Differential expression results
│   ├── counts_matrix.csv             # Combined count matrix
│   ├── group[1-3]_metadata.csv       # Sample metadata per group
│   ├── group[1-3]_normalized_counts.csv
│   ├── group[1-3]_vst_counts.csv     # VST-transformed (for PCA/heatmaps)
│   └── group[1-3]_*_vs_Control.csv   # DE results per comparison
│
├── 08_functional/                     # Functional analysis
│   └── cyanocyc_enrichment/          # Pathway enrichment results
│
├── 09_figures/                        # Publication figures
│   ├── qc_plots/
│   ├── volcano_plots/
│   ├── heatmaps/
│   └── reports/                      # HTML reports per group
│
├── logs/                              # Pipeline logs
│   └── run_sortmerna_all.sh          # Generated SortMeRNA script
│
├── rnaseq_pipeline.ipynb             # Main analysis notebook
├── environment.yml                    # Conda environment
├── CLAUDE.md                          # AI assistant instructions
└── README.md                          # This file
```

---

## Pipeline Phases

The analysis is implemented in `rnaseq_pipeline.ipynb` with 5 phases:

### Phase 1: Pre-processing

#### 1.1 Quality Control (FastQC + MultiQC)

```bash
fastqc -t 8 -o 01_QC/fastqc_raw/ Data/*.fastq.gz
multiqc 01_QC/fastqc_raw/ -o 01_QC/multiqc_reports/
```

| Parameter | Description |
|-----------|-------------|
| `-t 8` | Number of threads |
| `-o DIR` | Output directory |

#### 1.2 Trimming (fastp)

```bash
fastp \
    -i input_R1.fastq.gz \
    -I input_R2.fastq.gz \
    -o output_R1.fastq.gz \
    -O output_R2.fastq.gz \
    --detect_adapter_for_pe \
    --correction \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --thread 8 \
    --json sample.fastp.json \
    --html sample.fastp.html
```

| Parameter | Description |
|-----------|-------------|
| `-i`, `-I` | Input R1 and R2 files |
| `-o`, `-O` | Output R1 and R2 files |
| `--detect_adapter_for_pe` | Auto-detect adapters for paired-end |
| `--correction` | Enable base correction in overlapped regions |
| `--qualified_quality_phred 20` | Quality threshold (Q20 = 99% accuracy) |
| `--length_required 50` | Discard reads shorter than 50 bp |
| `--thread 8` | Number of threads |
| `--json`, `--html` | QC report outputs |

#### 1.3 Lane Merging

Concatenate L001 and L002 files for each sample:

```bash
cat sample_L001_R1.fastq.gz sample_L002_R1.fastq.gz > sample_R1.fastq.gz
cat sample_L001_R2.fastq.gz sample_L002_R2.fastq.gz > sample_R2.fastq.gz
```

#### 1.4 rRNA Removal (SortMeRNA)

```bash
sortmerna \
    --ref silva_db/SILVA_138.1_SSURef_NR99.fasta \
    --ref silva_db/SILVA_138.1_LSURef_NR99.fasta \
    --reads sample_R1.fastq.gz \
    --reads sample_R2.fastq.gz \
    --paired_in \
    --out2 \
    --aligned rRNA/sample \
    --other non_rRNA/sample \
    --fastx \
    --threads 4 \
    --workdir /tmp/sortmerna_sample
```

| Parameter | Description |
|-----------|-------------|
| `--ref FILE` | rRNA reference database (can specify multiple) |
| `--reads FILE` | Input FASTQ files (specify twice for PE) |
| `--paired_in` | Keep both reads if either matches rRNA |
| `--out2` | Output paired reads to separate files |
| `--aligned PREFIX` | Output path for rRNA reads |
| `--other PREFIX` | Output path for non-rRNA reads (clean) |
| `--fastx` | Output in FASTQ format |
| `--threads 4` | Threads (use ≤4, each needs ~15 GB RAM) |
| `--workdir DIR` | Temporary directory (must be unique per sample) |

**Note**: SortMeRNA is the slowest step (~45-60 min/sample). The notebook generates a shell script (`logs/run_sortmerna_all.sh`) for batch execution outside the notebook.

---

### Phase 2: Quantification and Alignment

#### 2.1 Transcriptome Extraction (gffread)

Extract transcript sequences from genome + GTF:

```bash
gffread \
    -w transcriptome.fa \
    -g genome.fna \
    annotation.gtf
```

| Parameter | Description |
|-----------|-------------|
| `-w FILE` | Output transcriptome FASTA |
| `-g FILE` | Input genome FASTA |
| (positional) | Input GTF annotation |

#### 2.2 Salmon Index (decoy-aware)

```bash
# Step 1: Create decoys list (chromosome names)
grep "^>" genome.fna | cut -d " " -f 1 | sed 's/>//g' > decoys.txt

# Step 2: Concatenate transcriptome + genome
cat transcriptome.fa genome.fna > gentrome.fa

# Step 3: Build index
salmon index \
    -t gentrome.fa \
    -d decoys.txt \
    -i salmon_index \
    -k 31 \
    -p 8
```

| Parameter | Description |
|-----------|-------------|
| `-t FILE` | Input sequences (gentrome = transcriptome + genome) |
| `-d FILE` | Decoy sequences file (prevents false quantification) |
| `-i DIR` | Output index directory |
| `-k 31` | k-mer size (31 is appropriate for bacteria) |
| `-p 8` | Number of threads |

#### 2.3 Salmon Quantification

```bash
salmon quant \
    -i salmon_index \
    -l A \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    -o quants/sample \
    --validateMappings \
    --gcBias \
    --seqBias \
    --numBootstraps 100 \
    -p 10
```

| Parameter | Description |
|-----------|-------------|
| `-i DIR` | Salmon index directory |
| `-l A` | Library type: Auto-detect strandedness |
| `-1`, `-2` | Input R1 and R2 files |
| `-o DIR` | Output directory |
| `--validateMappings` | Selective alignment for improved accuracy |
| `--gcBias` | Correct for GC content bias |
| `--seqBias` | Correct for sequence-specific bias |
| `--numBootstraps 100` | Bootstrap samples for uncertainty estimation |
| `-p 10` | Number of threads |

#### 2.4 Bowtie2 Alignment (for visualization only)

```bash
# Build index
bowtie2-build genome.fna bowtie2_index/genome

# Align reads
bowtie2 \
    -x bowtie2_index/genome \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    -p 8 \
    | samtools sort -@ 4 -o sample.bam

# Index BAM
samtools index sample.bam
```

| Parameter | Description |
|-----------|-------------|
| `-x PREFIX` | Index prefix |
| `-1`, `-2` | Input R1 and R2 files |
| `-p 8` | Number of threads |
| `samtools sort -@` | Threads for sorting |
| `samtools sort -o` | Output sorted BAM |

#### 2.5 Coverage Tracks (deepTools)

```bash
bamCoverage \
    -b sample.bam \
    -o sample.bw \
    --normalizeUsing CPM \
    --binSize 10 \
    -p 8
```

| Parameter | Description |
|-----------|-------------|
| `-b FILE` | Input BAM file |
| `-o FILE` | Output bigWig file |
| `--normalizeUsing CPM` | Normalize to Counts Per Million |
| `--binSize 10` | Resolution in bp |
| `-p 8` | Number of threads |

---

### Phase 3: Differential Expression Analysis (DESeq2)

DESeq2 analysis is performed separately for each experimental group via rpy2:

```r
library(DESeq2)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ condition
)

# Filter low counts (≥10 counts in ≥3 samples)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# Run differential expression
dds <- DESeq(dds)

# Extract results for each comparison
res <- results(dds, contrast=c("condition", "Treatment", "Control"))

# VST transformation for visualization
vsd <- vst(dds, blind=TRUE)
vst_counts <- as.data.frame(assay(vsd))
```

**Significance thresholds**: |log2FoldChange| > 1, adjusted p-value < 0.05

**Output files per group**:
- `groupX_normalized_counts.csv` - DESeq2 normalized counts
- `groupX_vst_counts.csv` - VST-transformed counts (for PCA, heatmaps)
- `groupX_Treatment_vs_Control.csv` - DE results per comparison

---

### Phase 4: Visualization

Interactive Plotly visualizations using VST-transformed counts:

- **PCA plots**: Sample clustering (uses VST counts)
- **Correlation heatmaps**: Sample-sample Pearson correlation
- **Volcano plots**: log2FC vs -log10(p-value)
- **Expression heatmaps**: Top DE genes across conditions

All plots use publication-quality styling defined in `PUBLICATION_LAYOUT`.

---

### Phase 5: CyanoCyc Pathway Enrichment

#### 5.1 Over-Representation Analysis (ORA)

Uses gseapy's `enrich()` function with Fisher's exact test:

```python
import gseapy as gp

# gene_list: DEGs (significant genes)
# gene_sets: dict mapping pathway -> gene list
# background: all genes tested

enr = gp.enrich(
    gene_list=deg_list,
    gene_sets=pathway_dict,
    background=background_genes,
    outdir=None
)
```

**Limitation**: ORA requires arbitrary significance cutoffs for DEG selection.

#### 5.2 Gene Set Enrichment Analysis (GSEA Prerank)

Uses all genes ranked by significance, no arbitrary cutoff:

```python
# Rank metric: -log10(pvalue) × sign(log2FoldChange)
ranked_genes = -np.log10(pvalue) * np.sign(log2fc)

gsea_res = gp.prerank(
    rnk=ranked_genes,           # Ranked gene list
    gene_sets=pathway_dict,      # Pathway definitions
    min_size=5,                  # Minimum genes in pathway
    max_size=500,                # Maximum genes in pathway
    permutation_num=1000,        # Permutations for p-value
    outdir=None
)
```

| Parameter | Description |
|-----------|-------------|
| `rnk` | Series/DataFrame with gene IDs as index, rank metric as values |
| `gene_sets` | Dict mapping pathway ID -> list of genes |
| `min_size` | Exclude pathways with fewer genes |
| `max_size` | Exclude pathways with more genes |
| `permutation_num` | Number of permutations (1000 recommended) |

**Output**: Normalized Enrichment Score (NES), FDR q-value, leading edge genes

---

## Running the Analysis

1. **Activate environment**: `conda activate pcc11901_rnaseq`

2. **Open notebook**: `jupyter lab rnaseq_pipeline.ipynb`

3. **Run phases sequentially** (each phase has `RUN_X = False` flags):
   - Phase 1: Set `RUN_FASTP = True`, etc.
   - **SortMeRNA**: Run generated script outside notebook (`bash logs/run_sortmerna_all.sh`)
   - Phase 2-5: Continue in notebook

4. **Outputs**:
   - DESeq2 results: `07_deseq2/`
   - Enrichment results: `08_functional/cyanocyc_enrichment/`
   - HTML reports: `09_figures/reports/`

---

## Hardware Recommendations

| Step | Threads | RAM | Time per sample |
|------|---------|-----|-----------------|
| fastp | 8 | 2 GB | ~2 min |
| SortMeRNA | 4 | 60 GB | 45-60 min |
| Salmon quant | 10 | 8 GB | ~3 min |
| Bowtie2 | 8 | 8 GB | ~10 min |
| DESeq2 | 1 | 16 GB | ~5 min/group |

**Total time estimate**: ~3-4 days for 60 samples (SortMeRNA is the bottleneck)

---

## References

- Conesa et al. (2016). A survey of best practices for RNA-seq data analysis. *Genome Biology*
- Love et al. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*
- Patro et al. (2017). Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods*
- Kopylova et al. (2012). SortMeRNA: fast and accurate filtering of ribosomal RNAs. *Bioinformatics*
