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
library(webshot)             # To save gt images

# Load fonts for postscript output

#font_import()               # Run once if fonts are not installed
loadfonts(device = "postscript")

#######################################################
# Initialising file paths
#######################################################

load_dot_env("config.env")

interim_data <- Sys.getenv("interimdatadir")
processed_data <- Sys.getenv("processeddatadir")
docs_data <- Sys.getenv("docsdir")

#######################################################
# Read in colocalised QTL pairs (H4 PP > 0.7)
#######################################################

QTLs_coloc <- fread(file.path(interim_data, "coloc_output.csv")) # Colocalised pQTLs and eQTLs

harmonised <- fread(file.path(interim_data, "QTLs_coloc_ready.csv")) # Harmonised pQTLs and eQTLs

#######################################################
# FUNCTIONS
#######################################################
#######################################################
# ld_clump, but catching the errors resulting from null results 
#######################################################

safe_ld_clump <- function(dat) {
  tryCatch(
    {
      ld_clump(dat,
               clump_r2 = 0.001,
               plink_bin = plink,
               bfile = "/Volumes/MRC-IEU-research/projects/ieu3/p3/023/working/data/raw/1000gp/EUR")
    },
    error = function(e) {
      message("No clumped SNPs = dropping outcome")
      return(NULL)
    }
  )
}

#######################################################
# Naive Concordance Test
#######################################################

check_concordance <- function(df, threshold){
  df$type <- ifelse(
    (df$beta_pQTL >= 0 & df$beta_eQTL >= 0) |
      (df$beta_pQTL < 0 & df$beta_eQTL < 0),
    "concordant", "discordant"
  )
  df$type[df$pval_pQTL > threshold] <- "pQTL_dropped"
  df$type[df$pval_eQTL > threshold] <- "eQTL_dropped"
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

names(pQTLs)[names(pQTLs) == "gene"] <- "exposure"

pQTLs_harmonised_coloc <- merge(pQTLs, harmonised, by = c("SNP", "exposure")) %>% 
  arrange(exposure)

#######################################################
# Clumping - pQTLs
#######################################################

plink <- genetics.binaRies::get_plink_binary() #For LD clumping
pQTL_clumped <- data.frame()

pQTLS_to_clump <- dplyr::tibble(rsid = pQTLs_harmonised_coloc$SNP,
                                pval = pQTLs_harmonised_coloc$pval.exposure, 
                                id = pQTLs_harmonised_coloc$exposure, 
                                chr = pQTLs_harmonised_coloc$chr.exposure, 
                                pos = pQTLs_harmonised_coloc$pos.exposure,
                                effect_allele = pQTLs_harmonised_coloc$effect_allele.exposure, 
                                other_allele = pQTLs_harmonised_coloc$other_allele.exposure, 
                                eaf_pQTL = pQTLs_harmonised_coloc$eaf.exposure, 
                                beta_pQTL = pQTLs_harmonised_coloc$beta.exposure, 
                                se_pQTL = pQTLs_harmonised_coloc$se.exposure,
                                pval_pQTL = pQTLs_harmonised_coloc$pval.exposure,
                                eaf_eQTL = pQTLs_harmonised_coloc$eaf.outcome,
                                beta_eQTL = pQTLs_harmonised_coloc$beta.outcome,
                                se_eQTL = pQTLs_harmonised_coloc$se.outcome,
                                pval_eQTL = pQTLs_harmonised_coloc$pval.outcome)

for (i in unique(pQTLS_to_clump$id)){
  dat <- pQTLS_to_clump %>% filter(id == i)
  
  dat_clumped <- safe_ld_clump(dat)

  pQTL_clumped <- rbind(pQTL_clumped, dat_clumped)
}

pQTL_clumped <- unique(pQTL_clumped)

#######################################################
# Naive Concordance Test - pQTLs
#######################################################

pQTL_concordance <- check_concordance(pQTL_clumped, threshold = 5e-8)

pQTL_concordance_table <- make_concordance_table(pQTL_concordance) %>%
  gt() %>%
  tab_header(title = md("Table of QTL Concordance"),
             subtitle = md("Selected on pQTL SNPs")) %>%
  fmt_number(rows = 1, decimals = 0) %>%
  fmt_number(rows = 2, decimals = 2)

#######################################################
# Naive Concordance Test - pQTLs - 1 SNP per gene
#######################################################

pQTL_concordance_1_SNP <- data.table()

for (i in unique(pQTL_concordance$id)){
  dat <- pQTL_concordance %>% 
    filter(id == i)
  
  # If a gene has multiple associated SNPs, take the lowest pQTL pvalue
  if (nrow(dat) > 1){
    dat <- dat %>% filter(pval == min(dat$pval)) 
  }
  # If there are still multiple SNPs, take the lowest eQTL pvalue
  if (nrow(dat) > 1){
    dat <- dat %>% filter(pval_eQTL == min(dat$pval_eQTL))
    # If there are still multiple SNPs, take first row
    dat <- dat[1,] 
  }
  
  pQTL_concordance_1_SNP <- rbind(pQTL_concordance_1_SNP, dat)
  
}

pQTL_concordance_table_1_SNP <- make_concordance_table(pQTL_concordance_1_SNP) %>%
  gt() %>%
  tab_header(title = md("Table of QTL Concordance"),
             subtitle = md("Selected on pQTL SNPs")) %>%
  fmt_number(rows = 1, decimals = 0) %>%
  fmt_number(rows = 2, decimals = 2)

#######################################################
# Save - pQTLs
#######################################################

fwrite(pQTL_concordance, file.path(processed_data, "Concordance/pQTL_concordance.csv"))
gtsave(pQTL_concordance_table, file.path(docs_data, "Concordance/pQTL_concordance.png"))

fwrite(pQTL_concordance_1_SNP, file.path(processed_data, "Concordance/pQTL_concordance_1_SNP.csv"))
gtsave(pQTL_concordance_table_1_SNP, file.path(docs_data, "Concordance/pQTL_concordance_1_SNP.png"))

#######################################################
# Taking hit2 SNPs (eQTLs) 
#######################################################

eQTLs <- tidyr::separate(QTLs_coloc, col = hit2, sep="_", into=c("SNP", "A1", "A2"), remove=F)

names(eQTLs)[names(eQTLs) == "gene"] <- "outcome"

eQTLs_harmonised_coloc <- merge(eQTLs, harmonised, by = c("SNP", "outcome")) %>% 
  arrange(outcome)

#######################################################
# Clumping - eQTLs
#######################################################

eQTL_clumped <- data.frame()

eQTLS_to_clump <- dplyr::tibble(rsid = eQTLs_harmonised_coloc$SNP,
                                pval = eQTLs_harmonised_coloc$pval.outcome, 
                                id = eQTLs_harmonised_coloc$outcome, 
                                chr = eQTLs_harmonised_coloc$chr.exposure, 
                                pos = eQTLs_harmonised_coloc$pos.exposure,
                                effect_allele = eQTLs_harmonised_coloc$effect_allele.exposure, 
                                other_allele = eQTLs_harmonised_coloc$other_allele.exposure, 
                                eaf_pQTL = eQTLs_harmonised_coloc$eaf.exposure, 
                                beta_pQTL = eQTLs_harmonised_coloc$beta.exposure, 
                                se_pQTL = eQTLs_harmonised_coloc$se.exposure,
                                pval_pQTL = eQTLs_harmonised_coloc$pval.exposure,
                                eaf_eQTL = eQTLs_harmonised_coloc$eaf.outcome,
                                beta_eQTL = eQTLs_harmonised_coloc$beta.outcome,
                                se_eQTL = eQTLs_harmonised_coloc$se.outcome,
                                pval_eQTL = eQTLs_harmonised_coloc$pval.outcome)

for (i in unique(eQTLS_to_clump$id)){
  dat <- eQTLS_to_clump %>% filter(id == i)
  
  dat_clumped <- safe_ld_clump(dat)
  
  eQTL_clumped <- rbind(eQTL_clumped, dat_clumped)
}

eQTL_clumped <- unique(eQTL_clumped)

#######################################################
# Naive Concordance Test - eQTLs
#######################################################

eQTL_concordance <- check_concordance(eQTL_clumped, threshold = 5e-8)

eQTL_concordance_table <- make_concordance_table(eQTL_concordance) %>%
  gt() %>%
  tab_header(title = md("Table of QTL Concordance"),
             subtitle = md("Selected on eQTL SNPs")) %>%
  fmt_number(rows = 1, decimals = 0) %>%
  fmt_number(rows = 2, decimals = 2)

#######################################################
# Naive Concordance Test - eQTLs - 1 SNP per gene
#######################################################

eQTL_concordance_1_SNP <- data.table()

for (i in unique(eQTL_concordance$id)){
  dat <- eQTL_concordance %>% 
    filter(id == i)
  
  # If a gene has multiple associated SNPs, take the lowest eQTL pvalue
  if (nrow(dat) > 1){
    dat <- dat %>% filter(pval == min(dat$pval)) 
  }
  # If there are still multiple SNPs, take the lowest pQTL pvalue
  if (nrow(dat) > 1){
    dat <- dat %>% filter(pval_eQTL == min(dat$pval_pQTL))
    # If there are still multiple SNPs, take first row
    dat <- dat[1,] 
  }
  
  eQTL_concordance_1_SNP <- rbind(eQTL_concordance_1_SNP, dat)
  
}

eQTL_concordance_table_1_SNP <- make_concordance_table(eQTL_concordance_1_SNP) %>%
  gt() %>%
  tab_header(title = md("Table of QTL Concordance"),
             subtitle = md("Selected on eQTL SNPs")) %>%
  fmt_number(rows = 1, decimals = 0) %>%
  fmt_number(rows = 2, decimals = 2)

#######################################################
# Save - eQTLs
#######################################################

fwrite(eQTL_concordance, file.path(processed_data, "Concordance/eQTL_concordance.csv"))
gtsave(eQTL_concordance_table, file.path(docs_data, "Concordance/eQTL_concordance.png"))

fwrite(eQTL_concordance_1_SNP, file.path(processed_data, "Concordance/eQTL_concordance_1_SNP.csv"))
gtsave(eQTL_concordance_table_1_SNP, file.path(docs_data, "Concordance/eQTL_concordance_1_SNP.png"))
