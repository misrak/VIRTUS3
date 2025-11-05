# VIRTUS3: Viral Detection Pipeline for Single-Cell RNA-seq

VIRTUS3 is a specialized bioinformatics pipeline designed to detect viral sequences in 10x Genomics Chromium single-cell RNA sequencing (scRNA-seq) data. 

## Overview

### Purpose
VIRTUS3 enables the detection and quantification of viral transcripts in scRNA-seq experiments. We demonstrate that VIRTUS3 can be used for identification of cells infected with EBV and measures viral gene expression at single-cell resolution. Leverages Cell Ranger for human transcriptome analysis and cell identification, then quantifies unmapped reads against viral references. Only quantifies reads that failed to map to the human transcriptome, reducing the false positive originating from ambiguous reads similar to both human and viral genome. Also cellranger is not the best tool for assesing viral transcripts because they discard reads mapped to multiple genes and viral genomes are compact and has a lot of overlapped genes. VIRTUS3 overcome this by combining quasi-alignment for viral transcripts quantifications.


## Input Requirements

### FASTQ Files
- **Format**: Paired-end FASTQ files (gzip compressed recommended)
- **Expected naming convention**: `{sample}_S{#}_L{lane}_{R1|R2}_{#}.fastq.gz`
  - Example: `SRRxxxx_S1_L001_R1_001.fastq.gz`
- **Chemistry**:
  - R1: Cell barcode + UMI (variable length depending on chemistry)
  - R2: cDNA sequence (typically 91bp for standard protocols)
- **Origin**: Raw output from 10x Genomics Chromium single-cell sequencing

### Reference Data
All reference files for EBV are included in the `data/` directory:
- **EBV FASTA**: `NC_007605.1_CDS_EBER12.fa` (96 sequences)
  - 94 protein-coding genes (CDS)
  - 2 non-coding RNAs (EBER1 and EBER2)
  - EBV strain B95-8 (RefSeq: NC_007605.1)
- **Transcript-to-Gene Mapping**: `NC_007605.1_CDS_EBER12.tgMap.tsv`
  - Maps NCBI identifiers to human-readable gene names
- **Salmon Index**: `NC_007605.1_CDS_EBER12_salmon_index/`
  - Pre-built k-31 index for rapid quantification

Users can also make a custom references.

## Installation & Setup

### Required Tools
1. **Cell Ranger** (tested with v9.0.1)
   - 10x Genomics official tool for scRNA-seq analysis
   - Download: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads

2. **Salmon** (tested with v1.4.0)
   - For alignment-free quantification
   - Installation: `conda install -c bioconda salmon`

3. **SAMtools** (tested with v1.21)
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
- **File**: `refdata-gex-GRCh38-2024-A` etc. Available from 10x Genomics support website.

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
| `--lib_alevin` | str | No | Additional Alevin options (e.g., "-l ISF --umiLength 12 --barcodeLength 16 --end 5" for 5' v3 data) |
| `--skip_exist` | flag | No | Skip steps where output already exists (useful for resuming runs) |
| `--n_cores` | int | No | Number of CPU cores for Cell Ranger (default: 8) |

#### --chemistry: chemistry option for cellranger
Please refer to [the official documents](https://www.10xgenomics.com/support/jp/software/cell-ranger/latest/resources/cr-command-line-arguments).

#### --lib_alevin: chemistry option for alevin
- 3' v3 : `--lib_alevin="-l ISR --chromiumV3"`
- 3' v4 : `to be confirmed, maybe --lib_alevin="-l ISR --umiLength 12 --barcodeLength 16 --end 3"`
- 5' v2 : `--lib_alevin="-l ISF --umiLength 10 --barcodeLength 16 --end 5"`
- 5' v3 : `--lib_alevin="-l ISF --umiLength 12 --barcodeLength 16 --end 5" `

Please refer to [the official documents](https://salmon.readthedocs.io/en/latest/library_type.html).

### Example: 5' scRNA-seq Data
```bash
python ./src/virtus3.py \
    --fastqs /path/to/fastq/palmer_scratch/SRR16976513_data \
    --chemistry_cr ARC-v1 \
    --sample SRR16976513 \
    --lib_alevin="-l ISR --chromiumV3" \
    --output /path/to/output \
    --index_human /path/to/reference/refdata-gex-GRCh38-2024-A \
    --index_virus ./data/NC_007605.1_CDS_EBER12_salmon_index \
    --tgMap ./data/NC_007605.1_CDS_EBER12.tgMap.tsv \
    --cellranger /path/to/CellRanger/9.0.1/cellranger \
    --salmon /path/to/salmon/salmon \
    --cores 20
```


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

## Citation

Manuscript is under preparation.

## Contact & Support

For questions or issues:
1. Check the log.txt file for error messages
2. Verify all input files and reference paths
3. Ensure all required tools (Cell Ranger, Salmon, SAMtools) are installed

## License

Specify your project's license here.