#######################################################
# 6b. Significance Testing for Human Protein Atlas Enrichment Analysis
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
# Reading in HPA labelled QTLs
#######################################################

QTLs_HPA <- fread(file.path(processed_data, "Human_Protein_Atlas/HPA_labelled_QTLs.csv"))

#######################################################
# Summing number of QTLs in each protein class for each concordance group
#######################################################

QTLs_HPA_con <- QTLs_HPA[which(QTLs_HPA$type == "concordant"),][,5:ncol(QTLs_HPA)]
QTLs_HPA_dis <- QTLs_HPA[which(QTLs_HPA$type == "discordant"),][,5:ncol(QTLs_HPA)]

HPA_summary <- colSums(QTLs_HPA_con)
HPA_summary <- cbind(HPA_summary, colSums(QTLs_HPA_dis))

HPA_summary <- as.data.frame(HPA_summary)
HPA_summary$Protein_Class <- rownames(HPA_summary)

colnames(HPA_summary) <- c("Concordant", "Discordant", "Protein Class")

#######################################################
# Creating 2X2 Contingency Tables
#######################################################

# Unique Protein Classes
subgroups <- unique(HPA_summary$`Protein Class`)

# Prepare results list
results <- list()

for (sub in subgroups) {
  # Create logical vectors
  in_sub <- HPA_summary[HPA_summary$`Protein Class` == sub,]
  out_sub <- HPA_summary[HPA_summary$`Protein Class` != sub,]
  
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
    subtitle = md("Protein Class")
  ) %>%
  fmt_number(
    columns = c(P_Value),
    decimals = 4
  ) %>%
  cols_label(
    Subgroup = "Protein Class",
    P_Value = "P-Value",
    Adj_P_Bonferroni = "Adjusted P-Value"
  )

#######################################################
# Save
#######################################################

gtsave(results_gt, file.path(docs_data, "Human_Protein_Atlas/HPA_Fisher_Test.png"))
