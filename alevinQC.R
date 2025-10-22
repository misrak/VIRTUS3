library(alevinQC)

alevinQCReport(baseDir = "/home/yyasumizu/NGS_public/GSE189141_GSE159674_EBVPBMCB_tonsil_scRNAseq/SRR16976513_cr_GRCh38/outs/NC_007605.1_CDS_salmon_whitelist",
               sampleId = "SRR16976513", 
               outputFile = "alevinReport.html", 
               outputFormat = "html_document",
               outputDir = tempdir(), forceOverwrite = TRUE)