#!/usr/bin/Rscript
# Script by Maarja Jõeloo
# CNV-based eQTL analysis

run.eqtl <- function(cnv.file = "/path/to/EXPR_QTL/EGCUT_CN_freqfilter_unique.RDS",
                     expr.file = "/path/to/rnaseq_logcpm.RDS",
                     covariates.file = "/path/to/combined_analysis_data_finalNumeric.RDS",
                     output.dir = "/path/to/EXPR_QTL",
                     pvOutputThreshold = 0.001) {
  require(MatrixEQTL)
  
  # read in data
  cnvs <- readRDS(cnv.file)
  expr <- readRDS(expr.file)
  cov.full <- readRDS(covariates.file)
  
  # format cnvs file
  cnv.meta <- cnvs[1:9] # contains CNV location and other metadata
  cnvs <- cnvs[-(1:9)] # contains genotypes
  names(cnvs) <- gsub("_CN", "", names(cnvs))
  rownames(cnvs) <- paste0(cnv.meta$CHR, ":", cnv.meta$START, "-", cnv.meta$END)
  
  # format covariates and expression files
  w <- which(!(cov.full$Vcode.wgs %in% names(cnvs)))
  expr <- expr[,-w]
  cov.full <- cov.full[-w,]
  cov <- t(as.matrix(cov.full[, c("RBC", "PLT", "Neut", "Mono", "Lymph", "Eo", "Baso", "gender.code",
                                  "age", "bmi", "PC1", "PC2", "PC3", "PC4", "SV1", "SV2")]))
  batch <- cov.full$batch.mixup
  batch.contrasts <- t(model.matrix(~ batch, contrasts = list(batch = "contr.sum"))[,-1])
  
  cov <- rbind(batch.contrasts, cov)
  colnames(cov) <- cov.full$Vcode.rna
  
  # reorder samples in CNV file according to expression file
  m <- match(cov.full$Vcode.wgs, names(cnvs))
  cnvs <- as.matrix(cnvs[,m])
  colnames(cnvs) <- cov.full$Vcode.rna
  
  #print(cnvs[1:10,1:10])
  #print(cov[1:10,1:10])
  #print(expr[1:10,1:10])
  
  
  cnvs.data <- SlicedData$new()
  cnvs.data$fileOmitCharacters <- "NA"
  cnvs.data$CreateFromMatrix(cnvs)
  cnvs.data$ResliceCombined(sliceSize = 1000)
  
  cvrt <- SlicedData$new()
  cvrt$CreateFromMatrix(cov)
  cvrt$ResliceCombined(sliceSize = 1000)
  
  expr.data <- SlicedData$new()
  expr.data$CreateFromMatrix(expr)
  expr.data$ResliceCombined(sliceSize = 1000)
  
  #print(expr.data)
  #print(cnvs.data)
  #print(cvrt)
  
  saveRDS(cnvs, file = "cnvs_all_formatted_for_MatrixEQTL.RDS")
  saveRDS(expr, file = "rnaseq_logcpm_formatted_for_MatrixEQTL.RDS")
  saveRDS(cov, file = "covariates_formatted_for_MatrixEQTL.RDS")
  
  me <- Matrix_eQTL_engine(
    snps = cnvs.data,
    gene = expr.data,
    cvrt = cvrt,
    output_file_name = paste0(output.dir, "/eqtl_pvalues_theshold_1e-3.txt"),
    pvOutputThreshold = pvOutputThreshold,
    useModel = modelLINEAR, 
    errorCovariance = numeric(), 
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = TRUE,
    noFDRsaveMemory = FALSE)
  
  tmp <- me$all$min.pv.gene
  saveRDS(tmp, file = paste0(output.dir, "/eqtl_min_pvalues.RDS"))
}

run.permutations <- function(n = 10000, expr.mat.file = "/path/to/rnaseq_logcpm_formatted_for_MatrixEQTL.RDS", 
                             cnvs.mat.file = "/path/to/cnvs_all_formatted_for_MatrixEQTL.RDS", 
                             cov.mat.file = "/path/to/covariates_formatted_for_MatrixEQTL.RDS", 
                             output.dir = "/path/to/EXPR_QTL", cores = 60) {
  require(MatrixEQTL)
  require(parallel)
  
  #read data
  expr.mat <- readRDS(expr.mat.file)
  cnvs.mat <- readRDS(cnvs.mat.file)
  cov.mat <- readRDS(cov.mat.file)
  
  cnvs <- SlicedData$new()
  cnvs$fileOmitCharacters <- "NA"
  cnvs$CreateFromMatrix(cnvs.mat)
  cnvs$ResliceCombined(sliceSize = 1000)
  
  cvrt <- SlicedData$new()
  cvrt$CreateFromMatrix(cov.mat)
  cvrt$ResliceCombined(sliceSize = 1000)
  
  min.pv <- mclapply(1:n, function(i) {
    # permute expression values
    s <- sample(1:ncol(expr.mat))
    expr.mat <- expr.mat[,s]
    
    expr <- SlicedData$new()
    expr$CreateFromMatrix(expr.mat)
    expr$ResliceCombined(sliceSize = 1000)
    
    me <- Matrix_eQTL_engine(
      snps = cnvs,
      gene = expr,
      cvrt = cvrt,
      output_file_name = paste0(output.dir, "/TMP/perm_", i, "_tmp.txt"),
      pvOutputThreshold = 1e-100,
      useModel = modelLINEAR, 
      errorCovariance = numeric(), 
      verbose = FALSE,
      pvalue.hist = FALSE,
      min.pv.by.genesnp = TRUE,
      noFDRsaveMemory = TRUE)
    min.pv <- me$all$min.pv.gene
    return(min.pv)
  }, mc.cores = cores)
  min.pv <- do.call("rbind", min.pv)
  
  saveRDS(min.pv, file = paste0(output.dir, "/TMP/eqtl_permutations_min_pvalue.RDS"))
}


p.value.correction.main <- function(eqtl.file = "/path/to/EXPR_QTL/eqtl_pvalues_theshold_1e-3.txt",
                                    min.permut.pvalues.file = "/path/to/EXPR_QTL/eqtl_permutations_min_pvalue.RDS",
                                    cores = 40) {
  require(parallel)
  require(dplyr)
  require(ismev)
  
  eqtl <- read.table(eqtl.file, header = TRUE)
  min.permut.pvalues <- readRDS(min.permut.pvalues.file)
  
  results <- mclapply(1:ncol(min.permut.pvalues), function(i) {
    cat(i, "\n")
    
    gene.name <- colnames(min.permut.pvalues)[i]
    min.pv.current.gene <- min.permut.pvalues[,i]
    eqtl.current.gene <- eqtl %>% filter(gene == gene.name)
    if (nrow(eqtl.current.gene) == 0) {
      return(NULL)
    }
    
    pv.perm <- lapply(1:nrow(eqtl.current.gene), function(j) {
      pv.nom <- as.numeric(eqtl.current.gene[j,5])
      
      # ordinary permutations test
      pv.ord <- (sum(pv.nom >= min.pv.current.gene) + 1) / (length(min.pv.current.gene) + 1)
      
      
      # gumbel distribution
      fit <- gum.fit(-log10(min.pv.current.gene), show = FALSE)
      par <- fit$mle
      pv.gumbel <- 1 - exp(-exp(-(-log10(pv.nom) - par[1]) / par[2]))
      if (fit$conv != 0) {
        pv.gumbel = NA
      }
      return(c(pv.ord, par[1], par[2], pv.gumbel))
    })
    pv.perm <- do.call("rbind", pv.perm)
    
    results <- cbind(eqtl.current.gene, pv.perm)
    return(results)
  }, mc.cores = cores)
  results <- do.call("rbind", results)
  names(results)[7:10] <- c("pv.corrected.perm", "gumbel.loc", "gumbel.scale", "pv.corrected.gumbel")
  
  saveRDS(results, file = "eqtl_pvalues_corrected_gumbel.RDS")
}

p.value.correction.min <- function(min.pv.file = "/path/to/EXPR_QTL/eqtl_min_pvalues.RDS",
                                   min.permut.pvalues.file = "/path/to/EXPR_QTL/eqtl_permutations_min_pvalue.RDS",
                                   cores = 40) {
  require(parallel)
  require(ismev)
  
  min.pv <- readRDS(min.pv.file)
  min.permut.pvalues <- readRDS(min.permut.pvalues.file)
  
  results <- mclapply(1:ncol(min.permut.pvalues), function(i) {
    cat(i, "\n")
    
    gene.name <- colnames(min.permut.pvalues)[i]
    min.pv.current.gene <- min.permut.pvalues[,i]
    pv <- min.pv[i]
    
    fit <- gum.fit(-log10(min.pv.current.gene), show = FALSE)
    par <- fit$mle
    pv.gumbel <- 1 - exp(-exp(-(-log10(pv) - par[1]) / par[2]))
    if (fit$conv != 0) {
      pv.gumbel = NA
    }
    
    return(c(pv, pv.gumbel))
  }, mc.cores = cores)
  results <- do.call("rbind", results)
  results <- data.frame(gene = colnames(min.permut.pvalues), pv = results[,1], pv.gumbel = results[,2])
  saveRDS(results, file = "min_pvalues_corrected_gumbel.RDS")
}
