# VIRTUS3: Viral Detection Pipeline for Single-Cell RNA-seq

VIRTUS3 is a specialized bioinformatics pipeline designed to detect **Epstein-Barr Virus (EBV)** sequences in 10x Genomics Chromium single-cell RNA sequencing (scRNA-seq) data. The pipeline combines the robust cell-calling capabilities of Cell Ranger with the sensitive viral transcript quantification of Salmon/Alevin.

## Overview

### Purpose
VIRTUS3 enables the detection and quantification of viral transcripts in scRNA-seq experiments where the primary focus is human transcriptome analysis. It identifies cells infected with EBV and measures viral gene expression at single-cell resolution.

### Key Features
- **Two-step strategy**: Leverages Cell Ranger for human transcriptome analysis and cell identification, then quantifies unmapped reads against viral references
- **Efficient processing**: Only quantifies reads that failed to map to the human transcriptome, reducing computational burden
- **Multi-lane support**: Handles sequencing runs across multiple lanes with proper barcode management
- **Flexible chemistry support**: Compatible with both 3' and 5' 10x Genomics Chromium chemistries
- **Checkpointing**: Resume interrupted runs with the `--skip_exist` flag
- **Multiple output formats**: Generates both scanpy-compatible AnnData objects (.h5ad) and human-readable CSV matrices

## Pipeline Workflow

```
Input FASTQ files (paired-end, 10x Chromium)
         ↓
    [Step 1: Cell Ranger Count]
    - Align to human transcriptome (GRCh38-2020-A)
    - Cell barcode identification and UMI counting
    - Generate BAM file and count matrices
         ↓
    [Step 2: Extract Unmapped Reads]
    - Extract reads that didn't map to human (samtools view -f 4)
    - Convert BAM to FASTQ format
    - Preserve cell barcode and UMI information
         ↓
    [Step 3: Create Barcode Whitelist]
    - Extract valid cell barcodes from Cell Ranger output
    - Create whitelist for viral quantification
         ↓
    [Step 4: Viral Quantification with Alevin]
    - Quantify unmapped reads against EBV reference (per lane)
    - Alignment-free k-mer based quantification
    - Aggregate transcript counts to genes
         ↓
    [Step 5: Result Aggregation]
    - Concatenate results from all lanes
    - Convert to AnnData format (scanpy-compatible)
    - Export as H5AD and CSV matrices
         ↓
Output: Viral UMI counts per cell
```

## Input Requirements

### FASTQ Files
- **Format**: Paired-end FASTQ files (gzip compressed recommended)
- **Expected naming convention**: `{sample}_S{#}_L{lane}_{R1|R2}_{#}.fastq.gz`
  - Example: `MS1713_tonsil_HHT_S1_L001_R1_001.fastq.gz`
- **Chemistry**:
  - R1: Cell barcode + UMI (variable length depending on chemistry)
  - R2: cDNA sequence (typically 91bp for standard protocols)
- **Origin**: Raw output from 10x Genomics Chromium single-cell sequencing

### Reference Data
All reference files are included in the `data/` directory:
- **EBV FASTA**: `NC_007605.1_CDS_EBER12.fa` (96 sequences)
  - 94 protein-coding genes (CDS)
  - 2 non-coding RNAs (EBER1 and EBER2)
  - EBV strain B95-8 (RefSeq: NC_007605.1)
- **Transcript-to-Gene Mapping**: `NC_007605.1_CDS_EBER12.tgMap.tsv`
  - Maps NCBI identifiers to human-readable gene names
- **Salmon Index**: `NC_007605.1_CDS_EBER12_salmon_index/`
  - Pre-built k-31 index for rapid quantification

## Installation & Setup

### System Requirements
- **RAM**: 200+ GB (for Cell Ranger human alignment)
- **CPU**: 8+ cores
- **Disk**: 500+ GB (for intermediate files)
- **Runtime**: 8-12 hours per sample

### Required Tools
1. **Cell Ranger** (v7.x or v8.x)
   - 10x Genomics official tool for scRNA-seq analysis
   - Download: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads

2. **Salmon** (v1.10.0 or later)
   - For alignment-free quantification
   - Installation: `conda install -c bioconda salmon`

3. **SAMtools** (v1.18+)
   - For BAM file manipulation
   - Installation: `conda install -c bioconda samtools`

### Python Environment
```bash
# Create conda environment
conda create -n virtus3 python=3.9

# Activate environment
conda activate virtus3

# Install required packages
conda install -c bioconda salmon samtools
pip install numpy pandas scipy scanpy
```

### Human Reference Genome
The pipeline expects the 10x Genomics formatted human reference:
- **File**: `refdata-gex-GRCh38-2020-A/`
- **Download**: Available from 10x Genomics support website
- **Setup**: Update the `CELLRANGER_REFERENCE` path in SLURM scripts

## Usage

### Basic Command
```bash
python src/virtus3.py \
  --fastq_dir /path/to/fastq/files \
  --output_dir /path/to/output \
  --sample_name my_sample \
  --cellranger_ref /path/to/GRCh38-2020-A \
  --chemistry auto
```

### Command-Line Arguments

| Argument | Type | Required | Description |
|----------|------|----------|-------------|
| `--fastq_dir` | str | Yes | Path to directory containing FASTQ files |
| `--output_dir` | str | Yes | Path where results will be written |
| `--sample_name` | str | Yes | Sample name (used in output file naming) |
| `--cellranger_ref` | str | Yes | Path to 10x Genomics reference genome (GRCh38-2020-A) |
| `--chemistry` | str | No | 10x chemistry (default: "auto"; options: "auto", "SC3Pv3", "SC3Pv2", "SC5P-R2", "ARC-v1") |
| `--lib_alevin` | str | No | Additional Alevin options (e.g., "-l ISF --umiLength 12 --barcodeLength 16 --end 5" for 5' data) |
| `--skip_exist` | flag | No | Skip steps where output already exists (useful for resuming runs) |
| `--n_cores` | int | No | Number of CPU cores for Cell Ranger (default: 8) |

### Example: 5' scRNA-seq Data
```bash
python src/virtus3.py \
  --fastq_dir /path/to/fastq \
  --output_dir /path/to/output \
  --sample_name tonsil_sample \
  --cellranger_ref /gpfs/gibbs/pi/hafler/hc865/references/refdata-gex-GRCh38-2020-A \
  --chemistry ARC-v1 \
  --lib_alevin="-l ISF --umiLength 12 --barcodeLength 16 --end 5" \
  --skip_exist
```

### SLURM Job Submission
The repository includes example SLURM submission scripts. To submit a job:
```bash
sbatch src/2step_cr_al_240627_h10x_Tonsil_MS1713MS1738.sh
```

Edit the script with your specific parameters before submission:
- Update `--fastq_dir` with your input directory
- Update `--output_dir` for desired output location
- Adjust `--chemistry` and `--lib_alevin` based on your sequencing protocol
- Modify resource requests (`--cpus-per-task`, `--mem`) as needed

## Output Files

### Directory Structure
```
output_dir/
├── cellranger_human/          # Cell Ranger human alignment results
│   └── outs/
│       ├── possorted_genome_bam.bam    # Human-aligned reads
│       ├── possorted_genome_bam.bam.bai
│       ├── filtered_feature_bc_matrix/ # Gene expression matrix
│       └── raw_feature_bc_matrix/      # All detected barcodes
├── unmapped_fqs/              # Extracted unmapped reads
├── alevin_virus_lane_001/     # Viral quantification (per lane)
│   └── alevin/
│       ├── quants_mat.mtx.gz
│       ├── quants_mat_rows.txt
│       ├── quants_mat_cols.txt
│       └── quants_mat.log
├── alevin_virus.h5ad          # Final output: AnnData object
├── alevin_virus.csv           # Final output: CSV matrix
└── log.txt                    # Pipeline execution log
```

### Output Files Explained

| File | Format | Description |
|------|--------|-------------|
| `alevin_virus.h5ad` | HDF5/AnnData | Scanpy-compatible single-cell data object containing viral UMI counts, cell barcodes, and gene annotations |
| `alevin_virus.csv` | CSV | Tab-separated matrix: viral genes (rows) × cells (columns) with UMI counts |
| `cellranger_human/outs/filtered_feature_bc_matrix/` | MTX format | Cell Ranger output: human gene expression matrix for quality control |
| `log.txt` | TXT | Execution log with pipeline parameters, step timing, and final viral UMI counts |
| `unmapped_fqs/` | FASTQ.GZ | Extracted unmapped reads for each lane (intermediate files) |

### Loading Results in Python
```python
import scanpy as sc

# Load viral expression data
adata = sc.read_h5ad('output_dir/alevin_virus.h5ad')

# Access viral UMI counts per cell
print(adata.X)  # Gene × Cell sparse matrix

# Cell metadata
print(adata.obs)  # Cell barcodes and metadata

# Gene metadata
print(adata.var)  # EBV gene names
```

## Data Formats

### Input: FASTQ
- Paired-end reads from 10x Chromium
- R1: Cell barcode + UMI (typically 26bp for 3', 28bp for 5')
- R2: cDNA transcript sequence

### Reference: FASTA
EBV coding sequences with structured headers:
```
>lcl|NC_007605.1_cds_YP_401631.1_1 [gene=LMP-2A] [locus_tag=EBV_gp110] [protein=latent_membrane_protein_2A]
ATGAGTCTCGAAGCTCGCCTGATGAATGAAGACCTGGATTACGTAG...
```

### Reference: tgMap (Transcript-to-Gene Mapping)
Two-column TSV file mapping transcript identifiers to gene names:
```
lcl|NC_007605.1_cds_YP_401631.1_1	LMP-2A
lcl|NC_007605.1_cds_YP_401633.1_1	LMP-1
EBER1	EBER1
EBER2	EBER2
```

### Output: AnnData H5AD
HDF5-based single-cell data format used by scanpy:
- `.X`: Sparse matrix of viral UMI counts (genes × cells)
- `.obs`: DataFrame with cell barcodes and metadata
- `.var`: DataFrame with EBV gene annotations
- `.obs_names`: Cell barcodes (10x Chromium format)
- `.var_names`: Gene names (human-readable)

### Output: CSV
Tab-separated matrix for compatibility with standard tools:
- Header row: Cell barcodes
- First column: Gene names
- Values: UMI counts

## Troubleshooting

### Cell Ranger Fails
- Ensure the reference genome path is absolute and correct
- Check that FASTQ files match the expected naming convention
- Verify sufficient disk space for BAM file (~150-200 GB)

### No Viral Reads Detected
- This is normal for many samples (indicates no viral infection)
- Pipeline creates empty matrices with correct dimensions
- Check log.txt for "Viral UMI count: 0"

### Memory Errors
- Increase `--mem` in SLURM script (200GB recommended)
- Reduce number of cores if RAM allocation is limited
- Cell Ranger parallelization requires ~20GB per core

### Alevin Quantification Fails
- Verify Salmon index exists at `data/NC_007605.1_CDS_EBER12_salmon_index/`
- Check that barcode whitelist was generated from Cell Ranger output
- Ensure unmapped FASTQ files contain valid sequences

## Citation

If you use VIRTUS3 in your research, please cite:
- Cell Ranger: Zheng et al., Nature Communications (2017)
- Salmon: Patro et al., Nature Methods (2017)
- AnnData: Wolf et al., Genome Biology (2018)

## Contact & Support

For questions or issues:
1. Check the log.txt file for error messages
2. Verify all input files and reference paths
3. Ensure all required tools (Cell Ranger, Salmon, SAMtools) are installed
4. Review SLURM job logs for resource-related errors

## License

Specify your project's license here.

## Project Status

**Current Branch**: feature/package_creation
- Working on pipeline packaging and distribution improvements

**Recent Updates**:
- Renamed pipeline.py to virtus3.py for clarity
- Consolidated and organized source code structure
- Enhanced documentation and README
