#######################################################
# 2. QTL Data Harmonisation and Preparation for the Colocalisation Pipeline
# Steps:
#   1. Filter QTL pairs for at least one significant measure
#   2. Harmonise
#   3. Prepare dataset for colocalistion pipeline
#######################################################

#######################################################
# Load Libraries
#######################################################

library(dotenv)
library(data.table)
library(susieR)
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)
library(tidyr)
library(ggplot2)
library(coloc)

#######################################################
# Initialising file paths
#######################################################

load_dot_env("config.env")

interim_data <- Sys.getenv("interimdatadir")

#######################################################
# Read in QTL pairs
#######################################################

joined <- fread(file.path(interim_data, "eQTL_pQTL_matched.csv"))

#######################################################
#Extracting QTL pairs with at least one significant measurement
#######################################################

threshold <- 5*10^-8

filtered <- joined %>% filter(PVAL_pQTLs < threshold | Pvalue_eQTLs < threshold)

# XXX QTL pairs for XXX proteins

#######################################################
#Formatting
#######################################################

pQTL_formatted <- format_data(data.frame(filtered), 
                              type = "exposure",
                              snp_col = "NAME",
                              chr_col = "CHROM_pQTLs",
                              pos_col = "GENPOS_pQTLs",
                              effect_allele_col = "ALLELE1_pQTLs",
                              other_allele_col = "ALLELE0_pQTLs",
                              phenotype_col = "PROTEIN_pQTLs",
                              eaf_col = "A1FREQ_pQTLs",
                              beta_col = "BETA_pQTLs",
                              se_col = "SE_pQTLs",
                              pval_col = "PVAL_pQTLs",
                              samplesize_col = "N_pQTLs",
                              z_col = "ZSCORE_pQTLs",
                              info_col = "INFO_pQTLs"
)
pQTL_formatted$id.exposure <- pQTL_formatted$exposure

# XXX pQTLs for XXX proteins

eQTL_formatted <- format_data(data.frame(filtered), 
                              type = "outcome", 
                              snp_col = "NAME",
                              chr_col = "SNPChr_eQTLs",
                              pos_col = "SNPPos_eQTLs",
                              effect_allele_col = "AssessedAllele_eQTLs",
                              other_allele_col = "OtherAllele_eQTLs",
                              phenotype_col = "GeneSymbol_eQTLs",
                              eaf_col = "AlleleB_all_eQTLs",
                              beta_col = "Beta_eQTLs",
                              se_col = "SE_eQTLs",
                              pval_col = "Pvalue_eQTLs",
                              samplesize_col = "NrSamples_eQTLs",
                              z_col = "Zscore_eQTLs",
                              gene_col = "Gene_eQTLs"
)
eQTL_formatted$id.outcome <- eQTL_formatted$outcome

# XXX eQTLs for XXX proteins

#######################################################
#Harmonise
#######################################################

harmonised <- TwoSampleMR::harmonise_data(pQTL_formatted, eQTL_formatted) %>% 
  filter(outcome == exposure) %>% 
  filter(mr_keep == TRUE)

# XXX QTL pairs for XXX proteins

#######################################################
# Prepare harmonised dataset for input into pipeline
# 1. Extract columns of interest
# 2. drop indels if any (only keep A, C, G, T SNPs)
#######################################################

cols_to_keep <- c("SNP", 
                  "chr.exposure", 
                  "pos.exposure",
                  "effect_allele.exposure", 
                  "other_allele.exposure",
                  "exposure",
                  "eaf.exposure",
                  "beta.exposure", 
                  "se.exposure",
                  "pval.exposure",
                  "samplesize.exposure",
                  "outcome",
                  "eaf.outcome",
                  "beta.outcome", 
                  "se.outcome",
                  "pval.outcome",
                  "samplesize.outcome"
)
harmonised <- harmonised[,cols_to_keep]

alleles <- c("A", "C", "G", "T")
harmonised$keep_alleles <- ifelse(harmonised$effect_allele.exposure %in% alleles & harmonised$other_allele.exposure %in% alleles,
                                  TRUE,
                                  FALSE)
harmonised <- harmonised[harmonised$keep_alleles == TRUE, ]
harmonised$keep_alleles<-NULL

for (i in 1:nrow(harmonised)){
  harmonised$SNP[i] <- unlist(strsplit(harmonised$SNP[i], split = ":"))[2]
}

# XXX QTL pairs for XXX proteins

#######################################################
# Save 
#######################################################

fwrite(harmonised, file.path(interim_data, "QTLs_coloc_ready.csv"))
