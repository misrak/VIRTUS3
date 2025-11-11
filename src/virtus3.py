#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: yyasumizu
# @Date: 2023-07-20
# @Last Modified time: 2025-11-05


"""
# specify all files in full path
python virtus3.py \
    --fastqs /home/hc865/palmer_scratch/SRR16976513_data \
    --chemistry_cr ARC-v1 \
    --sample SRR16976513 \
    --lib_alevin="-l ISR --chromiumV3" \
    --output /gpfs/gibbs/pi/hafler/hc865/VIRTUS3/test_output \
    --index_human /gpfs/gibbs/pi/hafler/hc865/VIRTUS3/data/refdata-gex-GRCh38-2024-A \
    --index_virus /gpfs/gibbs/pi/hafler/hc865/VIRTUS3/data/NC_007605.1_CDS_EBER12_salmon_index \
    --tgMap /gpfs/gibbs/pi/hafler/hc865/VIRTUS3/data/NC_007605.1_CDS_EBER12.tgMap.tsv \
    --cellranger /vast/palmer/apps/avx2/software/CellRanger/9.0.1/cellranger \
    --salmon /vast/palmer/apps/avx2/software/Salmon/1.4.0-gompi-2020b/bin/salmon \
    --cores 40

# for 5'
    --lib_alevin="-l ISF --umiLength 10 --barcodeLength 16 --end 5" \
"""

import subprocess
import argparse
import os
import datetime
import re
import sys
import glob
from scipy.io import mmread
import pandas as pd
import numpy as np
import scanpy as sc

# Handle imports for different execution modes
try:
    # When run as a package (virtus3 command or python -m virtus3)
    from ._version import __version__
    from .logger_config import setup_logger
except ImportError:
    # When run as a script (python virtus3.py)
    try:
        from _version import __version__
        from logger_config import setup_logger
    except ImportError:
        # Fallback if _version.py is not found
        __version__ = "unknown"
        from logger_config import setup_logger

# Initialize logger
logger = setup_logger(__name__)


def run_command(command):
    logger.info(f"Executing command: {command}")
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Command failed with return code {result.returncode}")
        logger.error(f"stderr: {result.stderr}")
    return result.stdout

def analyze_fastq_name(fastqs):
    list_attributes = []
    for filename in fastqs:
        pattern = r'^(.+)_S\d+_L(\d+)_([RI]\d)_\d+.fastq.gz$'
        match = re.match(pattern, filename)

        if match:
            sample_name = match.group(1)
            lane_number = match.group(2)
            read_number = match.group(3)

            if read_number == "R1":
                logger.info("# samples")
                logger.info(f"Sample name: {sample_name}")
                logger.info(f"Lane number: {lane_number}")
                list_attributes.append([sample_name, lane_number, filename, filename.replace("_R1_", "_R2_")])

        else:
            raise Exception(f"Wired fastq name: {filename}")

    return list_attributes



def pipeline(args):
    log = ""

    # 0. examine input fastqs
    list_fastqs = os.listdir(args.fastqs)

    # retain only fastq files
    list_fastqs = [x for x in list_fastqs if x.endswith(".fastq.gz")]

    logger.info(f'num fastqs: {len(list_fastqs)}')
    logger.info(f'fastqs: {list_fastqs}')
    df_samples = pd.DataFrame(analyze_fastq_name(list_fastqs), columns = ['sample_name', 'lane', 'file1', 'file2'])

    # 0. make output dir and change dir
    os.makedirs(args.output, exist_ok=True)
    os.chdir(args.output)

    # 1. cellranger count for human
    command = "cellranger --version"
    version = run_command(command)
    logger.info(f"CellRanger version: {version.strip()}")

    if int(version.split('-')[-1].split('.')[0]) >= 8:
        opt_createbam = "--create-bam=true"
    else:
        opt_createbam = ""

    command = f"""{args.cellranger} count \
                    --id=cellranger_human \
                    --chemistry={args.chemistry_cr} \
                    --transcriptome={args.index_human} \
                    --fastqs={args.fastqs} \
                    --sample={args.sample} \
                    --localcores={args.cores} --localmem=64 \
                    --include-introns true {opt_createbam}"""
    if (not args.skip_exist) or (not os.path.exists('./cellranger_human/outs/possorted_genome_bam.bam')):
        run_command(command)
    else:
        logger.info("Skipping cellranger count (output already exists)")

    # 2. extract unmapped reads
    os.chdir("cellranger_human/outs")
    command = f"samtools view -f 4 -h -b possorted_genome_bam.bam > unmapped.bam"

    if (not args.skip_exist) or (not os.path.exists('./unmapped.bam')):
        run_command(command)
    else:
        logger.info("Skipping samtools view (unmapped.bam already exists)")

    # Validate that unmapped.bam is not empty
    if os.path.getsize('./unmapped.bam') == 0:
        logger.warning("unmapped.bam is empty. This means no reads were unmapped (all reads mapped to human genome).")
        logger.warning("No viral reads will be detected. This is expected if the sample has no viral content.")

    # os.makedirs("unmapped_fqs", exist_ok=True)
    command = f"{args.cellranger} bamtofastq --nthreads={args.cores} unmapped.bam unmapped_fqs"
    if (not args.skip_exist) or (not os.path.exists('./unmapped_fqs')):
        run_command(command)
    else:
        logger.info("Skipping bamtofastq (unmapped_fqs already exists)")

    # Validate that bamtofastq produced output
    if not os.path.exists('./unmapped_fqs'):
        raise Exception("ERROR: bamtofastq failed to create unmapped_fqs directory. Check that unmapped.bam is a valid BAM file.")

    # 3. make barcode whitelist from cellranger output
    command = f"zcat raw_feature_bc_matrix/barcodes.tsv.gz | sed 's/-1//g' > raw_feature_bc_matrix/barcodes.tsv"
    if (not args.skip_exist) or (not os.path.exists('./raw_feature_bc_matrix/barcodes.tsv')):
        run_command(command)
    else:
        logger.info('Skipping barcode whitelist creation (already exists)')

    # 4. Alevin (per lane)
    f_out_h5ad = f"{args.output}/alevin_virus.h5ad"
    f_out_csv = f"{args.output}/alevin_virus.csv"

    lanes = df_samples.loc[df_samples['sample_name']==args.sample, 'lane'].unique()
    list_unmapped_fqs_R1 = glob.glob('unmapped_fqs/*/*_R1_00*.fastq.gz')
    list_unmapped_fqs_R2 = glob.glob('unmapped_fqs/*/*_R2_00*.fastq.gz')

    # Validate that fastq files were found
    if not list_unmapped_fqs_R1 or not list_unmapped_fqs_R2:
        raise Exception(f"ERROR: No unmapped fastq files found. This could mean:\n"
                       f"  1. The unmapped.bam file was empty (no unmapped reads)\n"
                       f"  2. The glob pattern doesn't match the actual file names\n"
                       f"  3. bamtofastq failed to create fastq files\n"
                       f"Files found - R1: {list_unmapped_fqs_R1}\nFiles found - R2: {list_unmapped_fqs_R2}")

    list_adata = []

    for lane in lanes:
        logger.info(f"Processing lane: {lane}")
        if len(lanes) == 1:
            command = f"""{args.salmon} alevin \
                            {args.lib_alevin} \
                            -1 {' '.join([x for x in list_unmapped_fqs_R1])} \
                            -2 {' '.join([x for x in list_unmapped_fqs_R2])}  \
                            -i {args.index_virus} \
                            -p {args.cores} \
                            -o alevin_virus_lane_{lane} \
                            --tgMap {args.tgMap} \
                            --whitelist raw_feature_bc_matrix/barcodes.tsv \
                            --dumpMtx
                        """
        else:
            command = f"""{args.salmon} alevin \
                            {args.lib_alevin} \
                            -1 {' '.join([x for x in list_unmapped_fqs_R1 if x.split('_')[-3] == 'L'+lane])} \
                            -2 {' '.join([x for x in list_unmapped_fqs_R2 if x.split('_')[-3] == 'L'+lane])}  \
                            -i {args.index_virus} \
                            -p {args.cores} \
                            -o alevin_virus_lane_{lane} \
                            --tgMap {args.tgMap} \
                            --whitelist raw_feature_bc_matrix/barcodes.tsv \
                            --dumpMtx
                        """

        if (not args.skip_exist) or (not os.path.exists(f'alevin_virus_lane_{lane}/logs/salmon_quant.log')):
            run_command(command)
        else:
            logger.info(f"Skipping alevin for lane {lane} (output already exists)")

        # parse salmon log
        f_salmon_log = f'alevin_virus_lane_{lane}/logs/salmon_quant.log'
        with open(f_salmon_log) as f:
            salmon_log = f.read()
        num_reads = re.search(r'Counted ([\d,]+) total reads in the equivalence classes', salmon_log).group(1)
        num_reads = int(num_reads.replace(",", ""))
        is_finish = re.search(r'\[jointLog\] \[info\] finished quantifyLibrary\(\)', salmon_log) != None

        if is_finish & (num_reads > 0):
            logger.info(f"Viral reads detected in lane {lane}, num reads: {num_reads}")
            # convert alevin output to h5ad
            f = f"alevin_virus_lane_{lane}/alevin/quants_mat.mtx.gz"
            f_obs = f"alevin_virus_lane_{lane}/alevin/quants_mat_rows.txt"
            f_var = f"alevin_virus_lane_{lane}/alevin/quants_mat_cols.txt"
            adata = sc.AnnData(X=np.array(mmread(f).todense()), obs=pd.read_csv(f_obs, header=None, index_col=0), var=pd.read_csv(f_var, header=None, index_col=0))
            # adata = sc.AnnData(mmread(f), pd.read_csv(f_obs, header=None, index_col=0), pd.read_csv(f_var, header=None, index_col=0))
            adata.obs.index.name = None
            adata.obs.index  = adata.obs.index + '-1'
            adata.var.index.name = None

            list_adata.append(adata)

        else:
            logger.info(f"No viral reads detected in lane {lane}")

            f_var = f"alevin_virus_lane_{lane}/alevin/quants_mat_cols.txt"
            df_var = pd.read_csv(f_var, header=None, index_col=0)
            adata = sc.AnnData(X=np.empty((0,df_var.shape[0])), var=df_var)
            adata.var.index.name = None
            list_adata.append(adata)

    adata_concat = sc.concat(list_adata)
    adata_concat.to_df().to_csv(f_out_csv)
    adata_concat.write(f_out_h5ad)
    num_viral_reads = adata_concat.to_df().sum().sum()
    logger.info(f"Total num viral reads (UMIs): {num_viral_reads}")
    log += f"Total num viral reads (UMIs): {num_viral_reads}"
    return log


def main(args=None):
    """
    Main entry point for VIRTUS3 pipeline.

    Parameters
    ----------
    args : list, optional
        Command line arguments. If None, uses sys.argv[1:].

    Returns
    -------
    int
        Exit code (0 for success, 1 for error)
    """
    parser = argparse.ArgumentParser(
        prog='virtus3',
        description='VIRTUS3: Detection of viral transcripts in single-cell RNA-seq data',
        epilog='For more information, visit: https://github.com/yyoshiaki/VIRTUS3'
    )
    parser.add_argument("--fastqs", type=str, help="input fastqs file for cellranger", required=True)
    parser.add_argument("--chemistry_cr", "-cc", type=str, help="chemistry for cellranger", required=True)
    parser.add_argument("--sample", "-s", type=str, help="sample name for cellranger", required=True)
    parser.add_argument("--output", "-o", type=str, help="output dir", required=True)
    parser.add_argument('--lib_alevin', '-l', type=str, help="library for alevin (3': ISR, 5': ISF)", required=True)
    parser.add_argument("--index_human", "-ih", type=str, help="index file (human, cellranger)", required=True)
    parser.add_argument("--tgMap", "-tg", type=str, help="tgMap file (virus, alevin)", required=True)
    parser.add_argument("--index_virus", "-iv", type=str, help="index file (virus, salmon)", required=True)
    parser.add_argument("--cellranger", "-c", type=str, help="cellranger path", required=False, default="cellranger")
    parser.add_argument("--salmon", "-sl", type=str, help="salmon path", required=False, default="salmon")
    parser.add_argument("--samtools", "-sam", type=str, help="samtools path", required=False, default="samtools")
    parser.add_argument("--cores", "-p", type=int, help="number of cores", required=False, default=40)
    parser.add_argument("--skip_exist", "-skip", action='store_true', help="skip if output file exists", required=False, default=False)
    parser.add_argument('-v', '--version', action='version', version=__version__, help='Show version and exit')

    try:
        parsed_args = parser.parse_args(args)
    except SystemExit as e:
        # argparse calls sys.exit() on error or --help
        return e.code if e.code else 0

    args_txt = '\n'.join(f'{k}={v}' for k, v in vars(parsed_args).items())
    logger.info(f"VIRTUS3 Arguments:\n{args_txt}")

    log = str(datetime.datetime.now()) + '\n'
    log += '2step_cr_Alevin\n\n'
    log += '*****args*****\n' + args_txt + '\n'

    # Validate that all paths are absolute
    for file in [parsed_args.fastqs, parsed_args.output, parsed_args.index_human,
                 parsed_args.index_virus, parsed_args.tgMap, parsed_args.cellranger,
                 parsed_args.salmon]:
        if not os.path.isabs(file):
            logger.error(f'The path "{file}" is not an absolute path.')
            return 1

    try:
        log += pipeline(parsed_args)

        with open(f"{parsed_args.output}/log.txt", mode='w') as f:
            f.write(log)

        logger.info("Pipeline completed successfully!")
        return 0

    except Exception as e:
        logger.error(f"Pipeline failed with exception: {e}")
        import traceback
        traceback.print_exc()
        return 1


def cli_main():
    """
    Entry point for installed command-line tool (virtus3 command).

    This is the function that setuptools will call when the user runs
    the 'virtus3' command after package installation via pip.
    It simply wraps main() and ensures proper exit code handling.
    """
    sys.exit(main())


if __name__ == "__main__":
    sys.exit(main())
