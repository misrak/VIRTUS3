#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: yyasumizu
# @Date: 2023-07-20
# @Last Modified time: 2023-10-20


"""
# specify all files in full path
python pipeline.py \
    --fastqs /home/yyasumizu/Yale/EBV/EBV_scRNAseq_pipeline/example/SRR16976513/fastq \
    --chemistry_cr ARC-v1 \
    --sample SRR16976513 \
    --lib_alevin="-l ISR --chromiumV3" \
    --output /home/yyasumizu/Yale/EBV/EBV_scRNAseq_pipeline/example/SRR16976513/outputs/2step_cr_Alevin \
    --index_human /home/yyasumizu/reference/refdata-gex-GRCh38-2020-A \
    --index_virus /home/yyasumizu/Yale/EBV/EBV_scRNAseq_pipeline/2step_cr_Alevin/data/NC_007605.1_CDS_EBER12_salmon_index \
    --tgMap /home/yyasumizu/Yale/EBV/EBV_scRNAseq_pipeline/2step_cr_Alevin/data/NC_007605.1_CDS_EBER12.tgMap.tsv \
    --cellranger /home/yyasumizu/Programs/cellranger-7.1.0/cellranger \
    --salmon /home/yyasumizu/Programs/salmon-latest_linux_x86_64/bin/salmon \
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

def run_command(command):
    print(command)
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
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
                print("# samples")
                print("Sample name: ", sample_name)
                print("Lane number: ", lane_number)
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

    print(f'num fastqs: {len(list_fastqs)}')
    print(f'fastqs: {list_fastqs}')
    df_samples = pd.DataFrame(analyze_fastq_name(list_fastqs), columns = ['sample_name', 'lane', 'file1', 'file2'])

    # 0. make output dir and change dir
    os.makedirs(args.output, exist_ok=True)
    os.chdir(args.output)

    # 1. cellranger count for human
    command = "cellranger --version"
    version = run_command(command)
    print(version)

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
        print("skip cellranger count")

    # 2. extract unmapped reads
    os.chdir("cellranger_human/outs")
    command = f"samtools view -f 4 -h -b possorted_genome_bam.bam > unmapped.bam"

    if (not args.skip_exist) or (not os.path.exists('./unmapped.bam')):
        run_command(command)
    else:
        print("skip samtools view (extract unmapped reads)")

    # os.makedirs("unmapped_fqs", exist_ok=True)
    command = f"{args.cellranger} bamtofastq --nthreads={args.cores} unmapped.bam unmapped_fqs"
    if (not args.skip_exist) or (not os.path.exists('./unmapped_fqs')):
        run_command(command)
    else:
        print("skip bamtofastq")

    # 3. make barcode whitelist from cellranger output
    command = f"zcat raw_feature_bc_matrix/barcodes.tsv.gz | sed 's/-1//g' > raw_feature_bc_matrix/barcodes.tsv"
    if (not args.skip_exist) or (not os.path.exists('./raw_feature_bc_matrix/barcodes.tsv')):
        run_command(command)
    else:
        print('skip make barcode whitelist')

    # 4. Alevin (per lane)
    f_out_h5ad = f"{args.output}/alevin_virus.h5ad"
    f_out_csv = f"{args.output}/alevin_virus.csv"

    lanes = df_samples.loc[df_samples['sample_name']==args.sample, 'lane'].unique()
    list_unmapped_fqs_R1 = glob.glob('unmapped_fqs/*/*_R1_00*.fastq.gz')
    list_unmapped_fqs_R2 = glob.glob('unmapped_fqs/*/*_R2_00*.fastq.gz')
    list_adata = []


    for lane in lanes:
        print(f"lane: {lane}")
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
            print(f"skip alevin for lane {lane}")

        # parse salmon log
        f_salmon_log = f'alevin_virus_lane_{lane}/logs/salmon_quant.log'
        with open(f_salmon_log) as f:
            salmon_log = f.read()
        num_reads = re.search(r'Counted ([\d,]+) total reads in the equivalence classes', salmon_log).group(1)
        num_reads = int(num_reads.replace(",", ""))
        is_finish = re.search(r'\[jointLog\] \[info\] finished quantifyLibrary\(\)', salmon_log) != None

        if is_finish & (num_reads > 0):
            print(f"viral read was detected in lane {lane}, num reads: {num_reads}")
            # convert alevin output to h5ad
            f = f"alevin_virus_lane_{lane}/alevin/quants_mat.mtx.gz"
            f_obs = f"alevin_virus_lane_{lane}/alevin/quants_mat_rows.txt"
            f_var = f"alevin_virus_lane_{lane}/alevin/quants_mat_cols.txt"
            adata = sc.AnnData(X=np.array(mmread(f).todense()), obs=pd.read_csv(f_obs, header=None, index_col=0), var=pd.read_csv(f_var, header=None, index_col=0))
            # adata = sc.AnnData(mmread(f), pd.read_csv(f_obs, header=None, index_col=0), pd.read_csv(f_var, header=None, index_col=0))
            adata.obs.index.name = None
            adata.obs.index  = adata.obs.index + '-' + str(int(lane))
            adata.var.index.name = None

            list_adata.append(adata)
        
        else:
            print(f"no viral read was detected in lane {lane}")

            f_var = f"alevin_virus_lane_{lane}/alevin/quants_mat_cols.txt"
            df_var = pd.read_csv(f_var, header=None, index_col=0)
            adata = sc.AnnData(X=np.empty((0,df_var.shape[0])), var=df_var)
            adata.var.index.name = None
            list_adata.append(adata)
    
    adata_concat = sc.concat(list_adata)
    adata_concat.to_df().to_csv(f_out_csv)
    adata_concat.write(f_out_h5ad)
    num_viral_reads = adata_concat.to_df().sum().sum()
    print(f"Total num viral reads (UMIs): {num_viral_reads}")
    log += f"Total num viral reads (UMIs): {num_viral_reads}"
    return log


if __name__ == "__main__": 
    parser = argparse.ArgumentParser(description='2step Cellranger, Alevin')
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

    args = parser.parse_args()
    args_txt = '\n'.join(f'{k}={v}' for k, v in vars(args).items())
    print(args_txt)

    log = str(datetime.datetime.now()) + '\n'
    log += '2step_cr_Alevin\n\n'
    log += '*****args*****\n' + args_txt + '\n'

    for file in [args.fastqs, args.output, args.index_human, args.index_virus, args.tgMap, args.cellranger, args.salmon]:
        if not os.path.isabs(file):
            print(f'Error: The path "{file}" is not an absolute path.')
            sys.exit(1)

    log += pipeline(args)


    with open(f"{args.output}/log.txt", mode='w') as f:
        f.write(log)

    print("finish!")