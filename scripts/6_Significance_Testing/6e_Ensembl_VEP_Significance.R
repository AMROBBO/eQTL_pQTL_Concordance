#######################################################
# 6e. Significance Testing for Ensembl VEP Enrichment Analysis
# Creating 2x2 contingency tables and running the Fisher's exact test to 
# see significance in the enrichment analyses
#######################################################

#######################################################
#Load in libraries
#######################################################

library(dotenv)
library(data.table)
library(dplyr)
library(gt)

#######################################################
# Initialising file paths
#######################################################

load_dot_env("config.env")

processed_data <- Sys.getenv("processeddatadir")
docs_data <- Sys.getenv("docsdir")

#######################################################
# Reading in Ensembl labelled QTLs
#######################################################

con_QTLs <- fread(file.path(processed_data, "Ensembl_VEP/VEP_labelled_QTLs_con.csv"))
dis_QTLs <- fread(file.path(processed_data, "Ensembl_VEP/VEP_labelled_QTLs_dis.csv"))

#######################################################
# Creating 2X2 Contingency Tables
#######################################################

VEP_summary <- colSums(con_QTLs[,7:ncol(con_QTLs)])
VEP_summary_2 <- colSums(dis_QTLs[,7:ncol(dis_QTLs)])

#Combine

VEP_summary <- merge(VEP_summary, VEP_summary_2, by = 0, all = T)

VEP_summary[is.na(VEP_summary)] <- 0
colnames(VEP_summary) <- c("Variant Effect", "Concordant", "Discordant")

# Unique Protein Classes
subgroups <- unique(VEP_summary$`Variant Effect`)

# Prepare results list
results <- list()

for (sub in subgroups) {
  # Create logical vectors
  in_sub <- VEP_summary[VEP_summary$`Variant Effect` == sub,]
  out_sub <- VEP_summary[VEP_summary$`Variant Effect` != sub,]
  
  # Build the 2x2 table
  a <- sum(in_sub$Concordant)
  b <- sum(in_sub$Discordant)
  c <- sum(out_sub$Concordant)
  d <- sum(out_sub$Discordant)
  
  table_2x2 <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  
#######################################################
# Fisher's exact test
#######################################################
  
  fisher_result <- fisher.test(table_2x2)
  
  # Append results
  results[[sub]] <- data.frame(
    Subgroup = sub,
    A_B1 = a,
    A_B2 = b,
    NotA_B1 = c,
    NotA_B2 = d,
    Odds_Ratio = fisher_result$estimate,
    P_Value = fisher_result$p.value
  )
}

# Combine into one data frame
results_df <- bind_rows(results) 

# Apply Bonferroni correction
results_df <- results_df %>%
  mutate(Adj_P_Bonferroni = p.adjust(P_Value, method = "bonferroni")) %>%
  dplyr::select("Subgroup", "P_Value", "Adj_P_Bonferroni")

#######################################################
# Create gt table
#######################################################

results_gt <- results_df %>%
  gt() %>%
  tab_header(
    title = md("Fisher's Exact Test Results"),
    subtitle = md("Predicted Variant Effects")
  ) %>%
  cols_label(
    Subgroup = "Predicted Variant Effect",
    P_Value = "P-Value",
    Adj_P_Bonferroni = "Adjusted P-Value"
  )

#######################################################
# Save
#######################################################

gtsave(results_gt, file.path(docs_data, "Ensembl_VEP/Ensembl_Fisher_Test.png"))
