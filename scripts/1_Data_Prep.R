#######################################################
# 1. QTL Data Processing
# Steps:
#   1. Load and format eQTL and pQTL datasets
#   2. Match SNPs between datasets (1-to-1 dbSNP rsID matching)
#######################################################

#######################################################
# Load in libraries
#######################################################

library(dotenv)
library(data.table)          # Fast file reading/writing
library(R.utils)
#library(TwoSampleMR)         # Mendelian randomisation tools
library(dplyr)               # Data wrangling
library(magrittr)            # Piping
library(gt)                  # Table formatting
library(webshot)

#######################################################
# Initialising file paths
#######################################################

load_dot_env("config.env")

raw_data <- Sys.getenv("rawdatadir")
interim_data <- Sys.getenv("interimdatadir")

#######################################################
# Load pQTL data - UKBBPPP
#######################################################

pQTL_files <- list.files(
  file.path(raw_data, "filtered_pQTLs"), 
  full.names = TRUE)

# 2,921 files (proteins)

# dbSNP rsID maps

olink <- list.files(
  file.path(raw_data, "SNP RSID maps"), 
  full.names = TRUE)[1:24]

#######################################################
# Load eQTL data - eQTLgen Phase I
#######################################################

eQTL_raw <- fread(
  file.path(raw_data, "Raw_eQTL_data/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt"),
  data.table = F)            # eQTL data

# 127,341,798 eQTL measures for 19,127 genes

# Effect allele frequencies (EAFs)

eQTL_eaf <- fread(
  file.path(raw_data, "Raw_eQTL_data/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt"),
  data.table = F)

# 10,823,015 measures

#######################################################
# Format pQTL data
#######################################################

header <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", 
            "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA")

all_pQTLs <- data.frame()
pQTL_rsIDs <- data.frame() 

# Read, clean, and combine all pQT files

for (f in pQTL_files){
  if (file.size(f) > 0){
    pQTL_data <- fread(f)
    colnames(pQTL_data) <- header
    
    # Extract protein name from filename
    
    protein <- unlist(strsplit(basename(f), "_"))[2] %>%
      strsplit("[.]") %>% 
      unlist() %>% 
      .[1]
    pQTL_data$PROTEIN <- protein
    
    #Calculating p value from log10p values given
    
    pQTL_data$PVAL <- 10 ^ (-(pQTL_data$LOG10P))
    
    #Calculate Z score
    
    pQTL_data$ZSCORE <- abs(pQTL_data$BETA / pQTL_data$SE)
    
    #Combine all pQTLs
    
    all_pQTLs <- rbind(all_pQTLs, pQTL_data)
    
    print(f)
  } else{
    message(f, " has no cis-SNPs")
  }
}

# 104,867,846 pQTL measures for 2,831 proteins

# Assigning dbSNP rsIDS

for(f in olink){
  rsIDs <- fread(f)
  print(f)
  IDs <- rsIDs[rsIDs$ID %in% all_pQTLs$ID, ]
  pQTL_rsIDs <- rbind(pQTL_rsIDs, IDs)
}

pQTL_mapped <- merge(all_pQTLs, pQTL_rsIDs, by = "ID")
pQTL_mapped <- pQTL_mapped[grep("rs", pQTL_mapped$rsid), ] # Only keep SNPs with rsIDs

# 99,797,456 pQTL measures for 2,831 proteins

#######################################################
# Format eQTL data
#######################################################

# Merge with allele frequency data

eQTL_merged <- merge(eQTL_raw, eQTL_eaf, by = "SNP")

# 127,341,798 eQTL measures for 19,127 proteins

#Identify palindromic SNPs

palindromic <- with(eQTL_merged,
                    (AssessedAllele == "A" & OtherAllele == "T") |
                      (AssessedAllele == "T" & OtherAllele == "A") |
                      (AssessedAllele == "G" & OtherAllele == "C") |
                      (AssessedAllele == "C" & OtherAllele == "G"))
eQTL_merged$palindromic <- ifelse(palindromic, "TRUE", "FALSE")

# Identifying ambiguous SNPs (0.42 < EAF < 0.58)

eQTL_merged$ambiguous <- ifelse(
  eQTL_merged$AlleleB_all > 0.42 & eQTL_merged$AlleleB_all < 0.58,
  "TRUE", "FALSE"
)

# Harmonising allele frequencies

effect_diff <- which(eQTL_merged$AssessedAllele != eQTL_merged$AlleleB)
eQTL_merged$AlleleB_all[effect_diff] <- 1 - eQTL_merged$AlleleB_all[effect_diff]

# Back-calculating the Beta and SE from Z-scores.

eQTL_eaf <- eQTL_merged %>%
  dplyr::select(SNP, Pvalue, SNPChr, SNPPos, AssessedAllele, OtherAllele, Zscore, Gene, 
                GeneSymbol, NrSamples, AlleleB_all, palindromic, ambiguous)

eQTL_eaf$Beta <- eQTL_eaf$Zscore / sqrt(2 * eQTL_eaf$AlleleB_all * (1 - eQTL_eaf$AlleleB_all) * 
                                          (eQTL_eaf$NrSamples + eQTL_eaf$Zscore^2))

eQTL_eaf$SE <- 1 / sqrt(2 * eQTL_eaf$AlleleB_all * (1 - eQTL_eaf$AlleleB_all) * 
                          (eQTL_eaf$NrSamples + eQTL_eaf$Zscore^2))

#######################################################
# 1-to-1 dbSNP rsID matching between pQTLs and eQTLs
#######################################################

pQTLs <- pQTL_mapped
eQTLs <- eQTL_eaf

# Create unique identifiers for SNP-protein/gene pairs

pQTLs$NAME <- paste(pQTLs$PROTEIN, pQTLs$rsid, sep = ":")
eQTLs$NAME <- paste(eQTLs$GeneSymbol, eQTLs$SNP, sep = ":")

# Align column names before joining

colnames(pQTLs)[1:22] <- paste0(colnames(pQTLs)[1:22], "_pQTLs")
colnames(eQTLs)[1:15] <- paste0(colnames(eQTLs)[1:15], "_eQTLs")

#Join two datasets based on unique identifiers: keep only matched SNPs
joined <- pQTLs[eQTLs, on = "NAME", nomatch = 0]

# 13,314,372 QTL pairs for 2,113 proteins

#######################################################
# Save
#######################################################

fwrite(joined, file.path(interim_data, "eQTL_pQTL_matched.csv"))