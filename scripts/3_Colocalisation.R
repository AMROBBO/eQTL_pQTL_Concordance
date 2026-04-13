#######################################################
# 3. QTL Colocalisation
# Steps:
#   1. Get LD matrix for the region
#   2. Harmonise LD reference panel against the QTL pairs
#   3. Run the finemapping diagnostics (kriging)
#   4. Use the kriging diagnostics to perform QC
#   5. Compute lambda and regularise LD matrix
#   6. Perform susie finemapping on the cleaned data
#   7. Perform susie coloc on the finemap dataset
#   8. Extract SNPs with H4 > 0.7
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

raw_data <- Sys.getenv("rawdatadir")
interim_data <- Sys.getenv("interimdatadir")

#######################################################
# Read in harmonised QTL pairs
#######################################################

harmonised <- fread(file.path(interim_data, "QTLs_coloc_ready.csv"))

#######################################################
# FUNCTIONS
#######################################################
#######################################################
# Step 1: Get LD matrix for the region
# 1. Take unique rsids from the prepared harmonised file in step 1
# 2. Compute the LD matrix (signed r values) using 1,000G EUR LD reference panel
#######################################################

# harmon.df = gene by gene chunk of harmonised

get_ld_matrix <- function(harmon.df, path_to_bfile){
  
  snp_list <- unique(harmon.df$SNP)
  
  myld<-NULL
  
  try(myld <- ieugwasr::ld_matrix(
    unique(snp_list),
    with_alleles = T,
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = path_to_bfile
  ))
  
  return(myld)
}

#######################################################
# Step 2: Harmonise LD reference panel against the GWAS
# 1. Merges the harmonised summary statistics against the 1000G reference panel
# 2. Check and resolve strand flip issues  (does not account for palindromic SNPs at this step, assumes these are on the same strand as ref panel)
# 3. Drop any alleles where the GWAS and LD panel do not match
# 4. Flip the beta on the GWAS to match the A2 plink assigned allele in the LD reference panel
#######################################################

LDref_harmon <- function(Gene_1, myld){
  
  #Get A1 and A2 from the matrix
  df1 <- data.frame(row.names(myld))
  names(df1) <- c("LDpanel.ID")
  
  df1 <- tidyr::separate(df1, col=LDpanel.ID, sep="_", into=c("SNP", "A1", "A2"), remove=F)
  df1$A1 <- ifelse(df1$A1 == "TRUE", 
                   dplyr::recode(df1$A1, 
                                 "TRUE" = "T"),
                   df1$A1) # For some reason one of the panels has TRUE as an allele, I presume this is T? CACNA1C
  
  df1$A2 <- ifelse(df1$A2 == "TRUE", 
                   dplyr::recode(df1$A2, 
                                 "TRUE" = "T"),
                   df1$A2) # RNF31 has TRUE as an allele
  
  #Merge the harmon data.
  df2 <- merge(x=df1, y=Gene_1, by="SNP")
  
  #Check alleles in the GWAS match the alleles in the LD reference panel (either on forward or reverse strand).
  #----------------------
  #Forward strand.
  #Get the alleles for the SNP from the LD reference panel.
  
  df2$G11 <- paste0(df2$A1, df2$A2)
  df2$G12 <- paste0(df2$A2, df2$A1)
  
  #Get the alleles for the SNP from GWAS.
  df2$G2 <- paste0(df2$other_allele.exposure, df2$effect_allele.exposure)
  
  #Check if match
  
  df2$geno.match <- ifelse(df2$G11 == df2$G2 | df2$G12 == df2$G2,
                           TRUE,
                           FALSE)
  
  #Reverse strand.
  #Now check the strand flip.
  #-----------------------
  
  #Find the palindromic SNPs as do not want to flip these.
  
  df2$palindromic <- ifelse(df2$G2 == "AT" | df2$G2 == "TA" | df2$G2 == "GC" | df2$G2 == "CG",
                            TRUE,
                            FALSE)
  
  
  
  
  df2$A1.flip <- ifelse(df2$palindromic==TRUE, df2$A1, dplyr::recode(df2$A1, 
                                                                     A = "T",
                                                                     C = "G", 
                                                                     G = "C",
                                                                     T = "A"))
  
  df2$A2.flip <- ifelse(df2$palindromic==TRUE, df2$A2, dplyr::recode(df2$A2, 
                                                                     A = "T",
                                                                     C = "G", 
                                                                     G = "C",
                                                                     T = "A"))
  
  
  df2$G21 <- paste0(df2$A1.flip, df2$A2.flip)
  df2$G22 <- paste0(df2$A2.flip, df2$A1.flip)
  
  #Check if match
  
  df2$geno.flip.match <- ifelse(df2$G21 == df2$G2 | df2$G22 == df2$G2,
                                TRUE,
                                FALSE)
  
  #Exclude mis-matching alleles and repair any strand flips.
  #------------
  
  df2$geno.keep <- ifelse(df2$geno.match == TRUE | df2$geno.flip.match == TRUE,
                          TRUE,
                          FALSE)
  
  #Remove the alleles that do not match with the LD reference panel.
  
  print("Number of alleles that match between the GWAS and LD reference.  Those that do not will be dropped.")
  print(table(df2$geno.keep))
  
  df2 <- df2[df2$geno.keep==TRUE,]
  
  #Fix the A1 and A2 to be on the same strand as the GWAS.
  
  df2$A1 <- ifelse(df2$geno.flip.match == TRUE, df2$A1.flip, df2$A1)
  df2$A2 <- ifelse(df2$geno.flip.match == TRUE, df2$A2.flip, df2$A2)
  
  #Tidy up data frame.
  
  drop <- c("geno.match", "G11", "G12", "geno.flip.match", "G21", "G22", "G2", "geno.keep", "A1.flip", "A2.flip")
  df2 <- df2[,!(names(df2) %in% drop)]
  #-----------
  
  #Align effect allele and effect size in GWAS to A2 in the LD reference panel
  #---------------------------
  
  #Check if need to flip the beta to be on the A2 allele for the reference panel.
  
  df2$need.A2flip <- ifelse(df2$A2 == df2$effect_allele.exposure, FALSE, TRUE)
  
  print("Number of betas to flip in the GWAS to match the A2 allele in the LD reference.")
  print(table(df2$need.A2flip))
  
  #Flip the alleles and beta.
  
  df2$LDflip.EA <- ifelse(df2$need.A2flip==TRUE, df2$other_allele.exposure, df2$effect_allele.exposure)
  df2$LDflip.nonEA <- ifelse(df2$need.A2flip==TRUE, df2$effect_allele.exposure, df2$other_allele.exposure)
  df2$LDflip.beta.exposure <- ifelse(df2$need.A2flip==TRUE, -1*df2$beta.exposure, df2$beta.exposure)
  df2$LDflip.beta.outcome <- ifelse(df2$need.A2flip==TRUE, -1*df2$beta.outcome, df2$beta.outcome)
  
  #Return the dataset.
  #-------------------
  #Fix the ordering of the LD panel to be the same as the GWAS file.
  
  ### Adding this ifelse, as this doesn't work if only one SNP in myld
  #Original: LD.result <- myld[df2$LDpanel.ID, df2$LDpanel.ID]
  
  if(nrow(myld) == 1){
    LD.result <- myld
  } else{
    LD.result <-myld[df2$LDpanel.ID, df2$LDpanel.ID]
  }
  
  #Sanity check, do not return if items mis-ordered.
  if(identical(rownames(LD.result), df2$LDpanel.ID)==FALSE){
    print("LD matrix not ordered correctly.")
    return( list("harmon_data" = NULL,
                 "LD_data" = NULL))
  }
  
  if(identical(rownames(LD.result), df2$LDpanel.ID)==TRUE){
    return( list("harmon_data" = df2,
                 "LD_data" = LD.result))
  }
  
}

#######################################################
# Step 3: Run the finemapping diagnostics (kriging).
# Compares the observed to expected z scores based on the LD from the reference panel.
# Kriging plots highlights SNP outliers in red that need to be QCed.
#######################################################

kriging_diagnostic <- function(harmon.df, ld.mat){
  
  #Exposure (quant)
  
  diagnostic <- susieR::kriging_rss(z=harmon.df$LDflip.beta.exposure/harmon.df$se.exposure, 
                                    R= ld.mat, 
                                    n= median(harmon.df$samplesize.exposure))
  
  dist1.df <- diagnostic$conditional_dist
  dist1.df$trait <- unique(harmon.df$exposure)
  dist1.df$SNP <- harmon.df$SNP
  dist1.df$palindromic <- harmon.df$palindromic
  
  #Outcome (discrete)
  diagnostic <- susieR::kriging_rss(z=harmon.df$LDflip.beta.outcome/harmon.df$se.outcome, 
                                    R=ld.mat, 
                                    n=median(harmon.df$samplesize.outcome))
  
  dist2.df <- diagnostic$conditional_dist 
  dist2.df$trait <- unique(harmon.df$outcome)
  dist2.df$SNP <- harmon.df$SNP
  dist2.df$palindromic <- harmon.df$palindromic
  
  #Merge and add in additional flags.
  
  dist.df <- rbind(dist1.df, dist2.df)
  dist.df$outlier <- ifelse(dist.df$logLR > 2 & abs(dist.df$z) > 2, TRUE, FALSE)
  dist.df$samesign <- (dist.df$z * dist.df$condmean) > 0
  
  return(dist.df)
}

#Scatterplot for the kriging diagnostic

kriging_diagnostic_plot <- function(kriging.df){
  
  cols <- c("TRUE" = "red", "FALSE" = "black")
  p <- ggplot(kriging.df, aes(x= condmean, y= z, col=outlier, shape=palindromic)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    scale_colour_manual(values = cols) +
    xlab("Expected value") + ylab("Observed z score") +
    facet_wrap(~trait, ncol = 1, scales = "free") +
    theme_bw()
  
  return(p)
}

#######################################################
# Step 4: Use the kriging diagnostics to perform QC
# This function returns always a list with two elements: harmon_data and LD_data.
#######################################################

#QC the kriging dataset.
kriging_QC <- function(harmon_gwas_LD, harmon_LD_matrix, kriging.df){
  
  print("Applying QC to the outliers from the kriging diagnostic")
  print(table(kriging.df$outlier))
  
  QC_flag <- sum(kriging.df$outlier)
  
  harmon.df <- harmon_gwas_LD
  ld.mat <- harmon_LD_matrix
  
  #No QC required.
  if(QC_flag == 0){
    message("No outliers. QC not required")
    
    krig_list <- list("harmon_data" = harmon.df,
                      "LD_data" = ld.mat)
  }
  
  #Do the QC and return the harmon data and matrix.
  if(QC_flag > 0){
    message("QC required")
    #Label the SNPs to QC on the harmon df.
    snps_to_QC <- unique(kriging.df[kriging.df$outlier==TRUE,]$SNP)
    harmon.df$kriging.outlier <- ifelse(harmon.df$SNP %in% snps_to_QC, TRUE, FALSE)
    
    #Label the SNPS that are outliers for both.
    trait.names <- unique(kriging.df$trait)
    trait1.outlier <- kriging.df %>% filter(trait==trait.names[1],outlier==TRUE) %>% select(SNP)
    trait2.outlier <- kriging.df %>% filter(trait==trait.names[2],outlier==TRUE) %>% select(SNP)
    bothtrait.outlier <- intersect(trait1.outlier, trait2.outlier)
    harmon.df$kriging.outlierboth <- ifelse(harmon.df$SNP %in% bothtrait.outlier, TRUE, FALSE)
    
    #print(table(harmon.df$kriging.outlierboth))
    
    #If the outlier is palindromic switch the beta otherwise keep the alleles in the ordering.
    if(any(harmon.df$palindromic)==TRUE){
      message("Palindromic SNPs found in the data.")
      harmon.df$LDflip.EA <- ifelse(harmon.df$palindromic & harmon.df$kriging.outlierboth , harmon.df$LDflip.nonEA,  harmon.df$LDflip.EA)
      harmon.df$LDflip.nonEA <- ifelse(harmon.df$LDflip.EA==harmon.df$other_allele.exposure, harmon.df$effect_allele.exposure, harmon.df$other_allele.exposure)
      harmon.df$LDflip.beta.exposure <- ifelse(harmon.df$palindromic & harmon.df$kriging.outlierboth , -1*harmon.df$LDflip.beta.exposure, harmon.df$LDflip.beta.exposure)
      harmon.df$LDflip.beta.outcome <- ifelse(harmon.df$palindromic & harmon.df$kriging.outlierboth , -1*harmon.df$LDflip.beta.outcome, harmon.df$LDflip.beta.outcome)
    } else{
      message("No palindromic SNPs found in the data.")
    }
    
    #Sort out the rules for dropping - if outlier in both and palindromic flip the sign, drop all other SNPs.
    palindromic.keep <- harmon.df[harmon.df$palindromic & harmon.df$kriging.outlierboth,]$SNP
    snps_to_remove <- setdiff(snps_to_QC, palindromic.keep)    
    harmon.df$kriging.remove <- ifelse(harmon.df$SNP %in% snps_to_remove, TRUE, FALSE)
    
    #Remove the SNPs from the LD matrix if needed. 
    if(length(snps_to_remove) == 0){  
      message("No SNPs to remove")
      krig_list <- list("harmon_data" = harmon.df,
                        "LD_data" = ld.mat)
    }else{
      message("Number of SNPs to remove: ", length(snps_to_remove))
      harmon.df <- harmon.df[harmon.df$kriging.remove == FALSE, ]
      ld.mat <- ld.mat[harmon.df$LDpanel.ID, harmon.df$LDpanel.ID]
      
      krig_list <-list("harmon_data" = harmon.df,
                       "LD_data" = ld.mat)    
    }
  }
  return(krig_list)
} 

#######################################################
# Step 5: Compute lambda and regularise LD matrix
#######################################################

regularize_LDmatrix <- function(harmon.df, ld.mat){
  #Exposure
  exp_lambda <- susieR::estimate_s_rss(z=harmon.df$LDflip.beta.exposure/harmon.df$se.exposure, 
                                       R= ld.mat, 
                                       n= median(harmon.df$samplesize.exposure))
  
  #Adjust the LD matrix
  exp_R <- ld.mat*(1-exp_lambda)
  diag(exp_R) <- diag(exp_R) + exp_lambda
  
  #Outcome
  out_lambda <- susieR::estimate_s_rss(z=harmon.df$LDflip.beta.outcome/harmon.df$se.outcome, 
                                       R= ld.mat, 
                                       n= median(harmon.df$samplesize.outcome))
  
  #Adjust the LD matrix
  out_R <- ld.mat*(1-out_lambda)
  diag(out_R) <- diag(out_R) + out_lambda
  
  return( list("exposure_lambda" = exp_lambda,
               "outcome_lambda" = out_lambda,
               "exposure_LD_data" = exp_R,
               "outcome_LD_data" = out_R))
}

#######################################################
# Step 6: Perform susie finemapping on the cleaned data.
#######################################################

do_susie_finemap <- function(harmon.df, ld.mat){
  
  #Do not proceed if LD matrices dims are not correct.
  
  #Perfect square, rownames == colnames
  stopifnot(identical(rownames(ld.mat$exposure_LD_data), colnames(ld.mat$exposure_LD_data)))
  stopifnot(identical(rownames(ld.mat$outcome_LD_data), colnames(ld.mat$outcome_LD_data)))
  
  #Exposure and outcome regularized matrix should be ordered the same.
  stopifnot(identical(rownames(ld.mat$exposure_LD_data), rownames(ld.mat$outcome_LD_data)))
  
  #SNPs in the harmon GWAS file should be in the same order as the LD ref.
  stopifnot(identical(rownames(ld.mat$exposure_LD_data), harmon.df$LDpanel.ID))
  
  #Exposure (quant).
  
  exposureforcoloc <- list(
    beta = harmon.df$LDflip.beta.exposure,
    varbeta = harmon.df$se.exposure^2,
    snp = harmon.df$LDpanel.ID,
    position = harmon.df$pos.exposure,
    MAF = ifelse(harmon.df$eaf.exposure > 0.5, 1 - harmon.df$eaf.exposure, harmon.df$eaf.exposure), 
    type = "quant",
    N = median(harmon.df$samplesize.exposure), 
    LD = ld.mat$exposure_LD_data
  )
  
  exposuresusie <- coloc::runsusie(exposureforcoloc)    
  print(summary(exposuresusie))
  
  #Outcome (quant).
  
  outcomeforcoloc <- list(
    beta = harmon.df$LDflip.beta.outcome,
    varbeta = harmon.df$se.outcome^2,
    snp = harmon.df$LDpanel.ID,
    position = harmon.df$pos.exposure,
    MAF =ifelse(harmon.df$eaf.outcome > 0.5, 1 - harmon.df$eaf.outcome, harmon.df$eaf.outcome), 
    type = "quant",
    N = median(harmon.df$samplesize.outcome),
    LD = ld.mat$outcome_LD_data
  )
  
  outcomesusie <- coloc::runsusie(outcomeforcoloc)    
  print(summary(outcomesusie))
  
  #Return the susie finemapping results.
  
  return(list("exposure_pre_finemap" = exposureforcoloc,
              "exposure_susie_finemap" = exposuresusie,
              "outcome_pre_finemap" = outcomeforcoloc,
              "outcome_susie_finemap" = outcomesusie))
  
}

#######################################################
# Step 7: Perform susie coloc on the finemap dataset
#######################################################

do_susie_coloc <- function(susie.finemap.result){
  
  exposure.finemap <- susie.finemap.result$exposure_susie_finemap
  outcome.finemap <- susie.finemap.result$outcome_susie_finemap
  
  coloc_out <- NULL
  
  try(coloc_out <- coloc.susie(exposure.finemap, outcome.finemap))
  
  if(!is.null(coloc_out)){
    print(coloc_out$summary)
  }
  return(coloc_out)
  
}

#######################################################
# Running through pipeline
#######################################################

path_to_bfile <- file.path(raw_data, "1000gp/EUR")

top_coloc_snps <- data.frame(matrix(ncol = 11, nrow = 0))
col_names <- c("nsnps", "hit1", "hit2", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf",
               "PP.H4.abf","idx1", "idx2", "gene")
colnames(top_coloc_snps) <- col_names

old_warn <- getOption("warn")
options(warn = 1)

for (i in unique(harmonised$exposure)){

  print(paste0("Processing gene: ", i))

  Gene_data <- harmonised[harmonised$exposure == i,]
  
  # Step 2: Get LD matrix
  
  snp_list <- unique(Gene_data$SNP)
  
  myld<-NULL
  
  myld <- try(ieugwasr::ld_matrix(
    unique(snp_list),
    with_alleles = T,
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = path_to_bfile
  ))
  
  if (inherits(myld, "try-error")) {
    message("Error for element ", i, ", skipping.")
    next   # continue to next iteration
  }
  
  # Step 3: Harmonising reference panel against GWAS
  
  LD.harm <- LDref_harmon(Gene_data, myld)
  
  # Step 4: Run the Finemap diagnostic (kriging)
  
  harmon.df <- LD.harm$harmon_data
  ld.mat <- LD.harm$LD_data
  
  kriging.df <- kriging_diagnostic(harmon.df, ld.mat)
  
  #Plot
  
  p <- kriging_diagnostic_plot(kriging.df)
  
  # Step 5: Use kriging diagnostics to perform QC
  
  krig_list <- kriging_QC(harmon.df, ld.mat, kriging.df)
  
  # Step 6: Compute lambda and regularise LD matrix
  
  harmon.df <- krig_list$harmon_data
  ld.mat <- krig_list$LD_data
  
  ld_reg <- regularize_LDmatrix(harmon.df, ld.mat)
  
  # Step 7: Perform finemapping on the cleaned data
  
  ld.mat <- ld_reg
  
  ## Susie can only be run on more than 1 SNP
  
  if(nrow(harmon.df) > 1){
    
    susie.finemap.result <- do_susie_finemap(harmon.df, ld.mat)
  
    # Step 8: Perform susie coloc on finemap dataset
  
    coloc_out <- do_susie_coloc(susie.finemap.result)
    
  } else{
    
    exposureforcoloc <- list(
      beta = harmon.df$LDflip.beta.exposure,
      varbeta = harmon.df$se.exposure^2,
      snp = harmon.df$LDpanel.ID,
      position = harmon.df$pos.exposure,
      MAF = ifelse(harmon.df$eaf.exposure > 0.5, 1 - harmon.df$eaf.exposure, harmon.df$eaf.exposure), 
      type = "quant",
      N = median(harmon.df$samplesize.exposure)
    )
    
    outcomeforcoloc <- list(
      beta = harmon.df$LDflip.beta.outcome,
      varbeta = harmon.df$se.outcome^2,
      snp = harmon.df$LDpanel.ID,
      position = harmon.df$pos.exposure,
      MAF =ifelse(harmon.df$eaf.outcome > 0.5, 1 - harmon.df$eaf.outcome, harmon.df$eaf.outcome), 
      type = "quant",
      N = median(harmon.df$samplesize.outcome)
    )
    
    coloc_out <- coloc::coloc.abf(exposureforcoloc, outcomeforcoloc)
    
    if(!is.null(coloc_out)){
      print(coloc_out$summary)
      
      coloc_out$summary <- as.data.frame(t(coloc_out$summary))
      
    }
  }
  
  # Step 9: Extracting SNPs with H4 > 0.7
  
  if (!is.null(coloc_out$summary)){
    
    strong_pairs <- subset(coloc_out$summary, PP.H4.abf > 0.7)
    
    if (nrow(strong_pairs) > 0){
      print("Credible set pairs with strong colocalization:")
      print(strong_pairs)
    
      strong_pairs$gene <- i
    
      top_coloc_snps <- rbind(top_coloc_snps, strong_pairs)
      
    } else {
      print(paste0(i, " has no common causal variants with H4 > 0.7"))
    }
  } else {
    print(paste0(i, " has no common causal variants"))
  }
}

#######################################################
# Save
#######################################################

fwrite(top_coloc_snps, file.path(interim_data, "coloc_output.csv"))
