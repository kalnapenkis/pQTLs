#!/usr/bin/Rscript


library(MatrixEQTL)
library(data.table)
library(dplyr)


protein_files <- c("/path/to/EGCUT_cvd2.ped",
                   "/path/to/EGCUT_cvd3.ped",
                   "/path/to/EGCUT_onc.ped",
                   "/path/to/EGCUT_inf.ped")
protein_arrays <- c("CVD2", "CVD3", "ONC", "INF")

cnv_file <- "/path/to/EGCUT_WGS_GenomeStrip_CNVs_filtered_ASonly_olink_freqfilt0.05_unique.RDS"
cnv <- readRDS(cnv_file) %>% as.matrix()

cov <- fread("/path/to/plink2.eigenvec", data.table=F)
nn <- cov$IID
cov <- cov[-1] %>% as.matrix() %>% t()
colnames(cov) <- nn

excl <- read.table("/path/to/exclude_due_to_relatedness")$V1

all_out <- NULL
for (i in 1:length(protein_files)) {
  prot <- fread(protein_files[i], data.table=F)
  
  samples <- intersect(colnames(cnv), intersect(colnames(cov), prot$IID))
  samples <- samples[!(samples %in% excl)]
  
  cnv_tmp <- cnv[,samples]
  
  cov_tmp <- cov[,match(colnames(cnv_tmp), colnames(cov))]
  
  prot <- prot[match(colnames(cnv_tmp), prot$IID),]
  nn <- prot$IID
  prot <- t(prot[,grepl("resid", colnames(prot))])
  colnames(prot) <- nn
  
  print(dim(prot))
  print(prot[1:10,1:10])
  
  cnvs.data <- SlicedData$new()
  cnvs.data$fileOmitCharacters <- "NA"
  cnvs.data$CreateFromMatrix(cnv_tmp)
  cnvs.data$ResliceCombined(sliceSize = 1000)
  
  cvrt <- SlicedData$new()
  cvrt$CreateFromMatrix(cov_tmp)
  cvrt$ResliceCombined(sliceSize = 1000)
  
  prot.data <- SlicedData$new()
  prot.data$CreateFromMatrix(prot)
  prot.data$ResliceCombined(sliceSize = 1000)
  
  
  me <- Matrix_eQTL_engine(
    snps = cnvs.data,
    gene = prot.data,
    cvrt = cvrt,
    output_file_name = "tmp.txt",
    pvOutputThreshold = 0.05,
    useModel = modelLINEAR,
    errorCovariance = numeric(),
    verbose = F,
    pvalue.hist = F,
    min.pv.by.genesnp = F,
    noFDRsaveMemory = FALSE)
  
  out <- me$all$eqtls
  out$array <- protein_arrays[i]
  all_out <- rbind(all_out, out)
}


all_out <- all_out[order(all_out$pvalue),]
saveRDS(all_out, file = "tmp.RDS")


