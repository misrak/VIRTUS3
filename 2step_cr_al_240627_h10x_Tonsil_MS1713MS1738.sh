#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20GB
#SBATCH --time=10:00:00
#SBATCH --mail-user=yoshiaki.yasumizu@yale.edu
#SBATCH --mail-type=ALL

set -x

module load CellRanger/8.0.1
module load miniconda
conda activate py3.8
module load SAMtools/1.18-GCC-12.2.0

DIR_PIPELINE=/home/yy693/pi_hafler/EBV_scRNAseq/EBV_scRNAseq_pipeline

ID=240627_h10x_Tonsil_GEX

SAMPLE=MS1713_tonsil_HHT
DIR_DATA=/home/yy693/palmer_scratch/EBV_scRNAseq/original_data/20240627_jm3942_22CFCCLT4_RQ25280_5p
DIR_OUTPUT=/home/yy693/pi_hafler/EBV_scRNAseq/output/$ID/$SAMPLE/2step_cr_Alevin

python3 $DIR_PIPELINE/2step_cr_Alevin/pipeline.py \
    --fastqs $DIR_DATA \
    --chemistry_cr auto \
    --sample $SAMPLE \
    --lib_alevin="-l ISF --umiLength 12 --barcodeLength 16 --end 5" \
    --output $DIR_OUTPUT \
    --index_human /gpfs/gibbs/data/genomes/10xgenomics/refdata-gex-GRCh38-2020-A \
    --index_virus $DIR_PIPELINE/2step_cr_Alevin/data/NC_007605.1_CDS_EBER12_salmon_index \
    --tgMap $DIR_PIPELINE/2step_cr_Alevin/data/NC_007605.1_CDS_EBER12.tgMap.tsv \
    --cellranger /vast/palmer/apps/avx2/software/CellRanger/8.0.1/bin/cellranger \
    --salmon /home/yy693/programs/salmon-latest_linux_x86_64/bin/salmon \
    --cores 40 \
    --skip_exist
