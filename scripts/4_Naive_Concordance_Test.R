#######################################################
# 4. LD pruning and Concordance test on Colocalised QTLs
# Steps:
#   1. Extract QTL pairs from harmonised dataset - using pQTL SNPs
#   2. For genes with multiple variants -> LD prune
#   3. Naive Concordance Test
#   4. Repeat using eQTL SNPs

# hit1 = exposure = pQTL
# hit2 = outcome = eQTL
#######################################################

#######################################################
#Load in libraries
#######################################################

library(dotenv)
library(arrow)               # Read parquet files
library(data.table)          # Fast file reading/writing
library(dplyr)               # Data wrangling
library(tidyr)               # Tidying data
library(ggplot2)             # Plotting
library(patchwork)           # Combining plots
library(ieugwasr)            # GWAS summary utilities
library(gt)                  # Pretty tables
library(grid)                # For unit()
library(extrafont)           # Font handling for publication-quality plots
library(stringr)             # For str_replace

# Load fonts for postscript output

#font_import()               # Run once if fonts are not installed
loadfonts(device = "postscript")

#######################################################
# Initialising file paths
#######################################################

load_dot_env("config.env")

interim_data <- Sys.getenv("interimdatadir")
results_data <- Sys.getenv("resultsdir")
docs_data <- Sys.getenv("docsdir")

#######################################################
# Read in colocalised QTL pairs (H4 PP > 0.7)
# (Output from SNP_coloc_all_genes.R)
#######################################################

QTLs_coloc <- fread(file.path(interim_data, "coloc_output.csv")) # Colocalised pQTLs and eQTLs

harmonised <- fread(file.path(interim_data, "QTLs_coloc_ready.csv")) # Harmonised pQTLs and eQTLs

#######################################################
# FUNCTIONS
#######################################################
#######################################################
# Naive Concordance Test
#######################################################

check_concordance <- function(df, threshold){
  df$type <- ifelse(
    (df$beta.pQTL >= 0 & df$beta.eQTL >= 0) |
      (df$beta.pQTL < 0 & df$beta.eQTL < 0),
    "concordant", "discordant"
  )
  df$type[df$pval.pQTL > threshold] <- "pQTL_dropped"
  df$type[df$pval.eQTL > threshold] <- "eQTL_dropped"
  df
}

#######################################################
# Concordance Table
#######################################################

make_concordance_table <- function(df){
  npairs <- nrow(df)
  nconc <- sum(df$type == "concordant")
  ndis <- sum(df$type == "discordant")
  npdrop <- sum(df$type == "pQTL_dropped")
  nedrop <- sum(df$type == "eQTL_dropped")
  
  data.frame(
    Category = c("Number", "Proportion %"),
    Pairs = c(npairs, "100"),
    Concordant = c(nconc, nconc/npairs*100),
    Discordant = c(ndis, ndis/npairs*100),
    pQTL_dropped = c(npdrop, npdrop/npairs*100),
    eQTL_dropped = c(nedrop, nedrop/npairs*100)
  ) 
}

#######################################################
# Taking hit1 SNPs (pQTLs) 
#######################################################

pQTLs <- tidyr::separate(QTLs_coloc, col = hit1, sep="_", into=c("SNP", "A1", "A2"), remove=F)

names(pQTLs)[names(pQTLs) == "gene"] <- "outcome"

pQTLs_harmonised_coloc <- merge(pQTLs, harmonised, by = c("SNP", "outcome")) %>% 
  arrange(outcome)

#######################################################
# Clumping - pQTLs
#######################################################

plink <- genetics.binaRies::get_plink_binary() #For LD clumping

pQTLs_clumped <- ld_clump(
  dplyr::tibble(rsid = pQTLs_harmonised_coloc$SNP, 
                pval = pQTLs_harmonised_coloc$pval.exposure, 
                id = pQTLs_harmonised_coloc$outcome, 
                chr = pQTLs_harmonised_coloc$chr.exposure, 
                pos = pQTLs_harmonised_coloc$pos.exposure,
                effect_allele = pQTLs_harmonised_coloc$effect_allele.exposure, 
                other_allele = pQTLs_harmonised_coloc$other_allele.exposure, 
                eaf.pQTL = pQTLs_harmonised_coloc$eaf.exposure, 
                beta.pQTL = pQTLs_harmonised_coloc$beta.exposure, 
                se.pQTL = pQTLs_harmonised_coloc$se.exposure,
                pval.pQTL = pQTLs_harmonised_coloc$pval.exposure,
                eaf.eQTL = pQTLs_harmonised_coloc$eaf.outcome,
                beta.eQTL = pQTLs_harmonised_coloc$beta.outcome,
                se.eQTL = pQTLs_harmonised_coloc$se.outcome,
                pval.eQTL = pQTLs_harmonised_coloc$pval.outcome),
  clump_r2 = 0.001,
  plink_bin = plink,
  bfile = "/Volumes/MRC-IEU-research/projects/ieu3/p3/023/working/data/1000gp/EUR" #Clumping
)

#######################################################
# Naive Concordance Test - pQTLs
#######################################################

pQTL_concordance <- check_concordance(pQTLs_clumped, threshold = 5e-8)

pQTL_concordance_table <- make_concordance_table(pQTL_concordance) %>%
  gt() %>%
  tab_header(title = md("Table of QTL Concordance"),
             subtitle = md("Selected on pQTL SNPs")) %>%
  fmt_number(rows = 1, decimals = 0) %>%
  fmt_number(rows = 2, decimals = 2)

#######################################################
# Save - pQTLs
#######################################################

fwrite(pQTL_concordance, file.path(results_data, "pQTL_concordance.csv"))

gtsave(pQTL_concordance_table, file.path(docs_data, "data/Coloc/pQTL_concordance.png"))

#######################################################
# Taking hit2 SNPs (eQTLs) 
#######################################################

eQTLs <- tidyr::separate(QTLs, col = hit2, sep="_", into=c("SNP", "A1", "A2"), remove=F)

names(eQTLs)[names(eQTLs) == "gene"] <- "outcome"

eQTLs_harmonised_coloc <- merge(eQTLs, harmonised, by = c("SNP", "outcome")) %>% 
  arrange(outcome)

#######################################################
# Clumping - eQTLs
#######################################################

eQTLs_clumped <- ld_clump(
  dplyr::tibble(rsid = eQTLs_harmonised_coloc$SNP, 
                pval = eQTLs_harmonised_coloc$pval.outcome, 
                id = eQTLs_harmonised_coloc$outcome, 
                chr = eQTLs_harmonised_coloc$chr.exposure, 
                pos = eQTLs_harmonised_coloc$pos.exposure,
                effect_allele = eQTLs_harmonised_coloc$effect_allele.exposure, 
                other_allele = eQTLs_harmonised_coloc$other_allele.exposure, 
                eaf.pQTL = eQTLs_harmonised_coloc$eaf.exposure, 
                beta.pQTL = eQTLs_harmonised_coloc$beta.exposure, 
                se.pQTL = eQTLs_harmonised_coloc$se.exposure,
                pval.pQTL = eQTLs_harmonised_coloc$pval.exposure,
                eaf.eQTL = eQTLs_harmonised_coloc$eaf.outcome,
                beta.eQTL = eQTLs_harmonised_coloc$beta.outcome,
                se.eQTL = eQTLs_harmonised_coloc$se.outcome,
                pval.eQTL = eQTLs_harmonised_coloc$pval.outcome),
  clump_r2 = 0.001,
  plink_bin = plink,
  bfile = "/Volumes/MRC-IEU-research/projects/ieu3/p3/023/working/data/1000gp/EUR" #Clumping
)

#######################################################
# Naive Concordance Test - eQTLs
#######################################################

eQTL_concordance <- check_concordance(eQTLs_clumped, threshold = 5e-8)

eQTL_concordance_table <- make_concordance_table(eQTL_concordance) %>%
  gt() %>%
  tab_header(title = md("Table of QTL Concordance"),
             subtitle = md("Selected on eQTL SNPs")) %>%
  fmt_number(rows = 1, decimals = 0) %>%
  fmt_number(rows = 2, decimals = 2)

#######################################################
# Save - eQTLs
#######################################################

fwrite(eQTL_concordance, "data/Coloc/eQTL_concordance.csv")

gtsave(eQTL_concordance_table, "data/Coloc/eQTL_concordance.png")
