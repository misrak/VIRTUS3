from scipy.io import mmread
import pandas as pd
import numpy as np
import scanpy as sc

# f = "/home/yyasumizu/NGS_public/GSE189141_GSE159674_EBVPBMCB_tonsil_scRNAseq/SRR16976513_cr_GRCh38/outs/NC_007605.1_CDS_salmon_whitelist/alevin/quants_mat.mtx.gz"
# f_obs = "/home/yyasumizu/NGS_public/GSE189141_GSE159674_EBVPBMCB_tonsil_scRNAseq/SRR16976513_cr_GRCh38/outs/NC_007605.1_CDS_salmon_whitelist/alevin/quants_mat_rows.txt"
# f_var = "/home/yyasumizu/NGS_public/GSE189141_GSE159674_EBVPBMCB_tonsil_scRNAseq/SRR16976513_cr_GRCh38/outs/NC_007605.1_CDS_salmon_whitelist/alevin/quants_mat_cols.txt"
# f_out = "/home/yyasumizu/NGS_public/GSE189141_GSE159674_EBVPBMCB_tonsil_scRNAseq/SRR16976513_cr_GRCh38/outs/NC_007605.1_CDS_salmon_whitelist/alevin.h5ad"

f = "/home/yyasumizu/NGS_public/GSE189141_GSE159674_EBVPBMCB_tonsil_scRNAseq/SRR16976513_cr_GRCh38/outs/NC_007605.1_CDS_salmon_whitelist_filtered/alevin/quants_mat.mtx.gz"
f_obs = "/home/yyasumizu/NGS_public/GSE189141_GSE159674_EBVPBMCB_tonsil_scRNAseq/SRR16976513_cr_GRCh38/outs/NC_007605.1_CDS_salmon_whitelist_filtered/alevin/quants_mat_rows.txt"
f_var = "/home/yyasumizu/NGS_public/GSE189141_GSE159674_EBVPBMCB_tonsil_scRNAseq/SRR16976513_cr_GRCh38/outs/NC_007605.1_CDS_salmon_whitelist_filtered/alevin/quants_mat_cols.txt"
f_out = "/home/yyasumizu/NGS_public/GSE189141_GSE159674_EBVPBMCB_tonsil_scRNAseq/SRR16976513_cr_GRCh38/outs/NC_007605.1_CDS_salmon_whitelist_filtered/alevin.h5ad"

adata = sc.AnnData(X=np.array(mmread(f).todense()), obs=pd.read_csv(f_obs, header=None, index_col=0), var=pd.read_csv(f_var, header=None, index_col=0))
# adata = sc.AnnData(mmread(f), pd.read_csv(f_obs, header=None, index_col=0), pd.read_csv(f_var, header=None, index_col=0))
adata.obs.index.name = None
adata.var.index.name = None

adata.write(f_out)