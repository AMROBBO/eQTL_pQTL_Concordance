#Formatting eQTL and pQTL data
#1-1 SNP matching
#Primary categorisation
#Naive concordance test
#Sensitivity Analysis

#######################################################
#Load in libraries
#######################################################

library(data.table) #For fread
library(TwoSampleMR)
library(dplyr)
library(biomaRt)
library(ieugwasr)
library(gwasvcf)
library(VariantAnnotation)
library(GenomicRanges)
library(magrittr)
library(ggforestplot)
library(ggplot2) #For plotting aesthetics
library(R.utils)
library(forestplot)
library(patchwork) #For plot layout
library(miamiplot)
library(gt)
library(webshot)

#######################################################
#Set Working Directory
#######################################################

setwd("")

#######################################################
#Read in eQTLgen Phase I data
#######################################################

eQTL_raw <- fread("data/Raw_eQTL_data/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",
                  data.table = F) #eQTL data
eQTL_eaf <- fread("data/Raw_eQTL_data/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt",
                  data.table = F) #EAF data

#######################################################
#Reading in pQTL data files
#######################################################

pQTL_files <- list.files("data/filtered_pQTLs/", full.names = TRUE) #pQTL files

olink <- list.files("data/SNP RSID maps/", full.names = TRUE) #rsIDs for SNPs
olink <- olink[1:24] #Removing last file

#######################################################
#Formatting pQTL data
#######################################################

header <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", 
            "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA")

all_pQTLs <- data.frame()
pQTL_rsIDs <- data.frame() 

#In a loop
for (f in pQTL_files){
  if ((file.size(f) > 0) == TRUE){
    pQTL_data <- fread(f) #Read in file
    colnames(pQTL_data) <- header #Name columns
    
    #Assign protein name
    protein <- unlist(strsplit(f, split = "/"))[4]
    protein <- unlist(strsplit(protein, split = "_"))[2]
    protein <- unlist(strsplit(protein, split = "[.]"))[1]
    pQTL_data$PROTEIN <- protein
    
    #Calculating p value from log10p values given
    pQTL_data$PVAL <- 10 ^ (-(pQTL_data$LOG10P))
    
    if (nrow(pQTL_data) > 0){
      #Calculate Z score
      pQTL_data$ZSCORE <- abs(pQTL_data$BETA / pQTL_data$SE)
      
      #Combine all pQTLs
      all_pQTLs <- rbind(all_pQTLs, pQTL_data)
      
      print(f)
    } else{
      print(paste(f, "has no significant SNPs"))
    }
  } else{
    print(paste(f, "has no cis-SNPs"))
  }
}

## Assigning rsIDS
for(f in olink){
  rsIDs <- fread(f)
  print(f)
  IDs <- rsIDs[which(rsIDs$ID %in% all_pQTLs$ID)]
  pQTL_rsIDs <- rbind(pQTL_rsIDs, IDs)
}

pQTL_mapped <- merge(all_pQTLs, pQTL_rsIDs, by = "ID")
pQTL_mapped <- pQTL_mapped[grep("rs", pQTL_mapped$rsid), ] #Only using SNPs with rsIDs

#######################################################
#Formatting eQTL data
#######################################################

#Adding EAF to eQTL data
eQTL_merged <- merge(eQTL_raw, eQTL_eaf, by = "SNP")

#Identifying palindromic SNPs
palindromic_at <- which(eQTL_merged$AssessedAllele == "A" & eQTL_merged$OtherAllele == "T")
palindromic_ta <- which(eQTL_merged$AssessedAllele == "T" & eQTL_merged$OtherAllele == "A")
palindromic_gc <- which(eQTL_merged$AssessedAllele == "G" & eQTL_merged$OtherAllele == "C")
palindromic_cg <- which(eQTL_merged$AssessedAllele == "C" & eQTL_merged$OtherAllele == "G")

eQTL_merged$palindromic <- "FALSE"
eQTL_merged$palindromic[palindromic_at] <- "TRUE"
eQTL_merged$palindromic[palindromic_ta] <- "TRUE"
eQTL_merged$palindromic[palindromic_gc] <- "TRUE"
eQTL_merged$palindromic[palindromic_cg] <- "TRUE"

#Identifying ambiguous SNPs (0.42 < EAF < 0.58)
ambiguous_lower <- which(eQTL_merged$AlleleB_all < 0.42)
ambiguous_higher <- which(eQTL_merged$AlleleB_all > 0.58)

eQTL_merged$ambiguous <- "TRUE"
eQTL_merged$ambiguous[ambiguous_lower] <- "FALSE"
eQTL_merged$ambiguous[ambiguous_higher] <- "FALSE"

#Harmonising data
effect_diff <- which(eQTL_merged$AssessedAllele != eQTL_merged$AlleleB) #Position of SNPs with different effect alleles
eQTL_merged$AlleleB_all[effect_diff] <- 1 - eQTL_merged$AlleleB_all[effect_diff]

#Back Calculating the Beta and SE.

eQTL_eaf <- eQTL_merged %>% select(c("SNP", "Pvalue", "SNPChr", "SNPPos", "AssessedAllele", "OtherAllele", "Zscore", "Gene", "GeneSymbol", "NrSamples", "AlleleB_all", "palindromic", "ambiguous"))

eQTL_eaf$Beta <- (eQTL_eaf$Zscore) / sqrt(2*eQTL_eaf$AlleleB_all*(1 - eQTL_eaf$AlleleB_all)*(eQTL_eaf$NrSamples + eQTL_eaf$Zscore^2))
eQTL_eaf$SE <- 1 / sqrt(2 * eQTL_eaf$AlleleB_all * (1 - eQTL_eaf$AlleleB_all) * (eQTL_eaf$NrSamples + eQTL_eaf$Zscore^2))

#######################################################
#1-to-1 SNP match
#######################################################

pQTLs <- pQTL_mapped
eQTLs <- eQTL_eaf

#Combining rsIDs and genes to create unique names
pQTLs$NAME <- paste(pQTLs$PROTEIN, pQTLs$rsid, sep = ":")
eQTLs$NAME <- paste(eQTLs$GeneSymbol, eQTLs$SNP, sep = ":")

#Join two datasets based on unique names
colnames(pQTLs)[1:22] <- paste(colnames(pQTLs)[1:22], "pQTLs", sep = "_")
colnames(eQTLs)[1:15] <- paste(colnames(eQTLs)[1:15], "eQTLs", sep = "_")

joined <- pQTLs[eQTLs, on = "NAME", nomatch = 0]

#######################################################
#Extracting QTL pairs with at least one significant measurement
#######################################################

threshold <- 5*10^-8

filtered <- joined %>% filter(PVAL_pQTLs < threshold | Pvalue_eQTLs < threshold)

#######################################################
#Formatting
#######################################################

pQTL_formatted <- format_data(data.frame(filtered), 
                              type = "exposure", 
                              snp_col = "NAME", 
                              phenotype_col = "PROTEIN_pQTLs",
                              beta_col = "BETA_pQTLs",
                              se_col = "SE_pQTLs",
                              eaf_col = "A1FREQ_pQTLs",
                              effect_allele_col = "ALLELE1_pQTLs",
                              other_allele_col = "ALLELE0_pQTLs",
                              pval_col = "PVAL_pQTLs",
                              z_col = "ZSCORE_pQTLs",
                              chr_col = "CHROM_pQTLs",
                              pos_col = "GENPOS_pQTLs",
                              samplesize_col = "N_pQTLs",
                              info_col = "INFO_pQTLs"
)
pQTL_formatted$id.exposure <- pQTL_formatted$exposure

eQTL_formatted <- format_data(data.frame(filtered), 
                              type = "outcome", 
                              snp_col = "NAME", 
                              phenotype_col = "GeneSymbol_eQTLs",
                              beta_col = "Beta_eQTLs",
                              se_col = "SE_eQTLs",
                              eaf_col = "AlleleB_all_eQTLs",
                              effect_allele_col = "AssessedAllele_eQTLs",
                              other_allele_col = "OtherAllele_eQTLs",
                              pval_col = "Pvalue_eQTLs",
                              z_col = "Zscore_eQTLs",
                              samplesize_col = "NrSamples_eQTLs",
                              chr_col = "SNPChr_eQTLs",
                              pos_col = "SNPPos_eQTLs",
                              gene_col = "Gene_eQTLs"
)
eQTL_formatted$id.outcome <- eQTL_formatted$outcome

#######################################################
#Harmonise
#######################################################

harmonised <- TwoSampleMR::harmonise_data(pQTL_formatted, eQTL_formatted) %>% filter(outcome == exposure)
harmonised <- harmonised[which(harmonised$mr_keep == TRUE),] 

#######################################################
#Categorizing into primary and non-primary QTL pairs
#######################################################

primary_QTLs <- data.frame()
non_primary_QTLs <- data.frame()

for (f in unique(harmonised$exposure)){
  QTL_tmp <- harmonised[which(harmonised$exposure == f),]
  
  primary <- QTL_tmp[which(QTL_tmp$z.exposure == max(QTL_tmp$z.exposure)),]
  non_primary <- QTL_tmp[which(QTL_tmp$z.exposure != max(QTL_tmp$z.exposure)),]
  
  if (nrow(primary) != 1){ #for the QTL pairs with multiple maximum pQTL Zscores -> taking highest eQTL score
    primary_2 <- primary[which(abs(primary$z.outcome) == max(abs(primary$z.outcome))),] #Used the absolute here as the eQTL Zscores are not absolute, whereas the pQTL Z scores were
    non_primary_2 <- primary[which(abs(primary$z.outcome) != max(abs(primary$z.outcome))),]
    
    primary <- primary_2
    non_primary <- rbind(non_primary, non_primary_2)
    
    if (nrow(primary) != 1){ #For the QTL pairs with multiple maximum pQTL and eQTL Zscores -> taking first instance
      primary_2 <- primary[1,]
      non_primary_2 <- primary[2:nrow(primary),]
      
      primary <- primary_2
      non_primary <- rbind(non_primary, non_primary_2)
    }
  }
  primary_QTLs <- rbind(primary_QTLs, primary)
  non_primary_QTLs <- rbind(non_primary_QTLs, non_primary)
}

#######################################################
#Naive Concordance Test - Primary QTLs
#######################################################

primary_QTLs$type <- ifelse(
  (primary_QTLs$beta.exposure >= 0 & primary_QTLs$beta.outcome >= 0) | (primary_QTLs$beta.exposure < 0 & primary_QTLs$beta.outcome < 0),
  "concordant",
  "discordant"
)

pQTL_dropped <- which(primary_QTLs$pval.exposure > 5*10^-8)
eQTL_dropped <- which(primary_QTLs$pval.outcome > 5*10^-8)

primary_QTLs$type[pQTL_dropped] <- "pQTL_dropped"
primary_QTLs$type[eQTL_dropped] <- "eQTL_dropped"

fwrite(primary_QTLs, "1-1_concordance_primary.csv")

#######################################################
#Naive Concordance Test - Non-Primary QTLs
#######################################################

non_primary_QTLs$type <- ifelse(
  (non_primary_QTLs$beta.exposure >= 0 & non_primary_QTLs$beta.outcome >= 0) | (non_primary_QTLs$beta.exposure < 0 & non_primary_QTLs$beta.outcome < 0),
  "concordant",
  "discordant"
)

pQTL_dropped <- which(non_primary_QTLs$pval.exposure > 5*10^-8)
eQTL_dropped <- which(non_primary_QTLs$pval.outcome > 5*10^-8)

non_primary_QTLs$type[pQTL_dropped] <- "pQTL_dropped"
non_primary_QTLs$type[eQTL_dropped] <- "eQTL_dropped"

#######################################################
#Concordance Table - Primary
#######################################################

npairs <- nrow(primary_QTLs)
nconc <- length(which(primary_QTLs$type == "concordant"))
ndis <- length(which(primary_QTLs$type == "discordant"))
npdropped <- length(which(primary_QTLs$type == "pQTL_dropped"))
nedropped <- length(which(primary_QTLs$type == "eQTL_dropped"))

concordance_primary = data.frame(
  Primary = c("Number", "Proportion %", "Sub Proportion %"),
  Pairs = c(npairs, "100", "100"),
  Concordant = c(nconc,(nconc/npairs * 100), (nconc/(nconc+ndis) * 100)),
  Discordant = c(ndis, (ndis/npairs * 100), (ndis/(nconc+ndis) * 100)),
  pQTL_dropped = c(npdropped, (npdropped/npairs * 100), (npdropped/(npdropped+nedropped) * 100)),
  eQTL_dropped = c(nedropped, (nedropped/npairs * 100), (nedropped/(npdropped+nedropped) * 100))
)

concordance_table <- concordance_primary %>%
  gt() %>%
  tab_header(title = md("Table of QTL Concordance"),
             subtitle = md("Primary QTL pairs")) %>%
  fmt_number(
    rows = 1,
    decimals = 0
  ) %>%
  fmt_number(
    rows = c(2, 3),
    decimals = 2
  )

gtsave(concordance_table, "")

#######################################################
#Concordance Table - Sensitivity
#######################################################

npairs <- nrow(non_primary_QTLs)
nconc <- length(which(non_primary_QTLs$type == "concordant"))
ndis <- length(which(non_primary_QTLs$type == "discordant"))
npdropped <- length(which(non_primary_QTLs$type == "pQTL_dropped"))
nedropped <- length(which(non_primary_QTLs$type == "eQTL_dropped"))

concordance_sens = data.frame(
  Primary = c("Number", "Proportion %", "Sub Proportion %"),
  Pairs = c(npairs, "100", "100"),
  Concordant = c(nconc, (nconc/npairs * 100), (nconc/(nconc+ndis) * 100)),
  Discordant = c(ndis, (ndis/npairs * 100), (ndis/(nconc+ndis) * 100)),
  pQTL_dropped = c(npdropped, (npdropped/npairs * 100), (npdropped/(npdropped+nedropped) * 100)),
  eQTL_dropped = c(nedropped, (nedropped/npairs * 100), (nedropped/(npdropped+nedropped) * 100))
)

concordance_table_sens <- concordance_sens %>%
  gt() %>%
  tab_header(title = md("Table of QTL Concordance"),
             subtitle = md("Non-Primary QTL pairs")) %>%
  fmt_number(
    rows = 1,
    decimals = 0
  ) %>%
  fmt_number(
    rows = c(2, 3),
    decimals = 2
  )

gtsave(concordance_table_sens, "")

#######################################################
#Chi-Square test
#######################################################

primary_chi <- concordance_primary[1,]
non_primary_chi <- concordance_sens[1,]

chi <- rbind(primary_chi, non_primary_chi)
rownames(chi) <- c("Primary", "Non_Primary")
chi <- chi[,-(1:2)]

# Perform chi-square test
chi_square_test <- chisq.test(chi)

#Save

chi_tab <- c(chi_square_test$statistic, chi_square_test$parameter, chi_square_test$p.value)
chi_tab <- as.data.table(chi_tab) %>% transpose()
colnames(chi_tab) <- c("X-squared", "df", "p-value")

chi_tab <- chi_tab%>%
  gt() %>%
  tab_header(title = md("Chi-Square Test"),
             subtitle = md("Primary VS Non-Primary QTL Pairs")) %>%
  fmt_number(
    columns = 1,
    n_sigfig = 5)

gtsave(chi_tab, "")
