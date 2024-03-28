#!/usr/bin/Rscript

# Author: Maarja Jõeloo
# Date: 2022-09-21

# script used to identify CNV-tagging SNPs


# ==== Libraries

library(data.table)
library(dplyr)
library(parallel)
library(stringr)

# === Variables

query_snps_file <- "/path/to/Sun2022_pQTLs_hg37.txt"
maf_file <- "/path/to/marker_frequencies_EstBB_MAF0.01.frq"
annot_file <- "/path/to/estbb_wgs_cnv_filtered_majAF0.99_VEP104.tsv"
all_snp_dir <- "/path/to/EGV_WGS_genotypes"
cnv_file <- "/path/to/EGCUT_WGS_GenomeStrip_CNVs_filtered_ASonly.RDS"
window <- 500000
r2_cutoff <- 0.2

# === Functions

read_query_snps <- function(query_snps_file) {
	Query_SNPs <- fread(query_snps_file, data.table=F, dec=",") %>%
		filter(chr %in% 1:22) # keep autosomes only

	cat("Input:", nrow(Query_SNPs), "query SNPs read...\n")
	return(Query_SNPs)
}

read_cnvs <- function(cnv_file) {
	CNVs <- readRDS(cnv_file)
	Sample_IDs <- gsub("_CN", "", gsub("\\.", "-", names(CNVs[10:ncol(CNVs)])))

	CNVs <- CNVs %>% filter(unique == 1)
	cat("Input:", nrow(CNVs), "unique CNVs read...\n")
	return(list(Sample_IDs = Sample_IDs, CNVs = CNVs))
}

.read_snp_genotypes <- function(all_snp_dir, chr) {
	bim_file <- paste0(all_snp_dir, "/wgs_genotypes_chrom", chr, ".bim")
	genotype_file <- paste0(all_snp_dir, "/wgs_genotypes_chrom", chr, ".raw.RDS")

	BIM_Chr <- fread(bim_file, data.table=F, header=F)
	GEN_Chr <- readRDS(genotype_file)

	return(list(BIM_Chr = BIM_Chr, GEN_Chr = GEN_Chr))
}

.r2_calculations_per_chr <- function(Query_SNPs, CNVs_Data, all_snp_dir, chr_this, window, r2_cutoff) {
	cat("Analysing: Chromosome", chr_this, "\n")

	Query_SNPs_Chr <- Query_SNPs %>% filter(chr == chr_this) %>%
			arrange(position)
	SNP_genotypes_Chr <- .read_snp_genotypes(all_snp_dir, chr_this)
	cat("=====", sum(Query_SNPs_Chr$position %in% SNP_genotypes_Chr$BIM_Chr$V4), "present in our data and included in calculations...\n")

	
	# query SNP locations in genotype data
	is_query_snp <- which(SNP_genotypes_Chr$BIM_Chr$V4 %in% Query_SNPs_Chr$position)
	BIM_Chr_Query <- SNP_genotypes_Chr$BIM_Chr[is_query_snp,]
	GEN_Chr_Query <- SNP_genotypes_Chr$GEN_Chr[match(CNVs_Data$Sample_IDs, SNP_genotypes_Chr$GEN_Chr$FID), is_query_snp+6]

	# calculate R2
	tagSNP_Chr <- bind_rows(mclapply(1:length(is_query_snp), function(i) {
		CNVs_Flank <- CNVs_Data$CNVs %>% filter(CHR == BIM_Chr_Query$V1[i], 
				START <= BIM_Chr_Query$V4[i] + window,
				END >= BIM_Chr_Query$V4[i] - window)
		if (nrow(CNVs_Flank) == 0) return(NULL)

		CNVs_Flank_Matrix <- t(as.matrix(CNVs_Flank[10:ncol(CNVs_Flank)]))
		R <- cor(GEN_Chr_Query[,i], CNVs_Flank_Matrix, use="pairwise.complete.obs")
		R2 <- R*R

		is_above_cutoff <- which(R2 > r2_cutoff)
		if (length(is_above_cutoff) == 0) return(NULL)

		output <- data.frame(chr = rep(chr_this, length(is_above_cutoff)), 
					position = rep(BIM_Chr_Query$V4[i], length(is_above_cutoff)),
					cnv = paste0(CNVs_Flank$CHR[is_above_cutoff], ":", CNVs_Flank$START[is_above_cutoff], "-", CNVs_Flank$END[is_above_cutoff]),
					r = R[is_above_cutoff],
					r2 = R2[is_above_cutoff])
		return(output)
	}, mc.cores = 40))
	return(tagSNP_Chr)
}

full_snp_cnv_r2_calculations <- function(Query_SNPs, CNVs_Data, all_snp_dir, window, r2_cutoff) {
	tagSNP <- bind_rows(lapply(1:22, function(chr) {
		return(.r2_calculations_per_chr(Query_SNPs, CNVs_Data, all_snp_dir, chr, window, r2_cutoff))
	}))
	return(tagSNP)
}

match_r2_to_query_snps <- function(Query_SNPs, full_r2_table) {
	# get max r2 per SNP
	tagSNP_best_per_SNP <- full_r2_table %>% mutate(snp_id = paste(chr, position, sep="_")) %>%
			group_by(snp_id) %>%
			summarise(nCNV_02 = n(), nCNV_08 = sum(r2 > 0.8), maxR2 = max(r2), maxR2_R = r[which.max(r2)], maxR2_CNV = cnv[which.max(r2)],
					tag_CNV = paste(cnv[which(r2 > 0.8)], collapse=";"))
	print(head(tagSNP_best_per_SNP))

	snp_match <- match(paste(Query_SNPs$chr, Query_SNPs$position, sep="_"), tagSNP_best_per_SNP$snp_id)
	Query_SNPs$nCNV_02 <- tagSNP_best_per_SNP$nCNV_02[snp_match]
	Query_SNPs$nCNV_08 <- tagSNP_best_per_SNP$nCNV_08[snp_match]
	Query_SNPs$maxR2 <- tagSNP_best_per_SNP$maxR2[snp_match]
	Query_SNPs$maxR2_R <- tagSNP_best_per_SNP$maxR2_R[snp_match]
	Query_SNPs$maxR2_CNV <- tagSNP_best_per_SNP$maxR2_CNV[snp_match]
	
	Query_SNPs$tag_CNV <- tagSNP_best_per_SNP$tag_CNV[snp_match]
	Query_SNPs$tag_CNV[Query_SNPs$tag_CNV == ""] <- NA

	return(Query_SNPs)
}

add_est_frequency <- function(Query_SNPs_with_r2, maf_file) {
	maf_table <- fread(maf_file, data.table=F)

	snp_match <- match(Query_SNPs_with_r2$rsID, maf_table$SNP)
	maf_table <- maf_table[snp_match,]

	Query_SNPs_with_r2$A1_Est <- maf_table$A1
	Query_SNPs_with_r2$A2_Est <- maf_table$A2
	Query_SNPs_with_r2$Freq_Est <- ifelse(Query_SNPs_with_r2$A1_Est == Query_SNPs_with_r2$A2, maf_table$MAF, 1-maf_table$MAF)
	return(Query_SNPs_with_r2)
}

.get_annotations_table <- function(annot) {
        annot_tbl <- fread(annot, skip = "#Uploaded_variation", header=T, data.table=F)
        names(annot_tbl)[1] <- "Uploaded_variation"
        annot_tbl$Impact <- str_match(annot_tbl$Extra, "IMPACT=(.*?);")[,2]
        annot_tbl$Impact[is.na(annot_tbl$Impact)] <- gsub("IMPACT=", "", annot_tbl$Extra[is.na(annot_tbl$Impact)])
        annot_tbl$Impact[grepl("coding_sequence_variant", annot_tbl$Consequence)] <- "HIGH"
        return(annot_tbl)
}

.fix_comma_list <- function(comma_list) {
        elements <- unique(sort(unlist(strsplit(as.character(comma_list), ","))))
        new_comma_list <- paste(elements, collapse=",")
        return(new_comma_list)
}

.get_most_severe_impact <- function(impact_list) {
        if ("HIGH" %in% impact_list) return("HIGH")
        if ("MODERATE" %in% impact_list) return("MODERATE")
        if ("LOW" %in% impact_list) return("LOW")
        if ("MODIFIER" %in% impact_list) return("MODIFIER")
        return(NA)
}


add_cnv_annotation <- function(Query_SNPs_with_r2, annot_file) {
	annot_table <- .get_annotations_table(annot_file)

	impact_table <- bind_rows(mclapply(1:nrow(Query_SNPs_with_r2), function(i) {
		if (is.na(Query_SNPs_with_r2$maxR2_CNV[i])) return(data.frame(maxR2_CNV_Impact = NA, maxR2_CNV_Consequence = NA))

		cnv_this <- paste0("CNV_", gsub(":", "_", gsub("-", "_", Query_SNPs_with_r2$maxR2_CNV[i])))
		annot_this <- annot_table %>% filter(Uploaded_variation == cnv_this)
		impact_this <- .get_most_severe_impact(annot_this$Impact)
		cons_this <- .fix_comma_list(paste(annot_this$Consequence, collapse=","))

		return(data.frame(maxR2_CNV_Impact = impact_this, maxR2_CNV_Consequence = cons_this))
	}, mc.cores = 40))

	Query_SNPs_with_r2 <- cbind(Query_SNPs_with_r2, impact_table)
	return(Query_SNPs_with_r2)
}



# === Main

Query_SNPs <- read_query_snps(query_snps_file)
CNVs_Data <- read_cnvs(cnv_file)

full_r2_table <- full_snp_cnv_r2_calculations(Query_SNPs, CNVs_Data, all_snp_dir, window, r2_cutoff)
write.table(full_r2_table, file = "results/Sun2022_pqtl_cnv_correlations_0.2.tsv", row.names=F, quote=F, sep="\t")

# match r2 results with original Query table
Query_SNPs_with_r2 <- match_r2_to_query_snps(Query_SNPs, full_r2_table)
Query_SNPs_with_r2 <- add_est_frequency(Query_SNPs_with_r2, maf_file)
Query_SNPs_with_r2 <- add_cnv_annotation(Query_SNPs_with_r2, annot_file)
write.table(Query_SNPs_with_r2, file = "results/Sun2022_pqtl_cnv_correlations_0.2_full_formatted.tsv", row.names=F, quote=F, sep="\t")

# format final table for R2>0.8
Query_SNPs_with_r2 <- Query_SNPs_with_r2 %>% mutate(A2_freq_Est = Freq_Est) %>% 
	select(chr, position, rsID, A1, A2, target, cis_trans, A2_freq_discovery, A2_freq_replication, A2_freq_Est, maxR2, maxR2_CNV, maxR2_CNV_Impact, maxR2_CNV_Consequence) %>%
	filter(maxR2 > 0.8)
write.table(Query_SNPs_with_r2, file = "results/Sun2022_pqtl_cnv_tagging.tsv", row.names=F, quote=F, sep="\t")


