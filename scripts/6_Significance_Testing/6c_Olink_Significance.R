#######################################################
# 6c. Significance Testing for Olink Enrichment Analysis
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
# Reading in Olink labelled QTLs
#######################################################

olink_QTLs <- fread(file.path(processed_data, "Olink/Olink_labelled_QTLs.csv"))

#######################################################
# Creating 2X2 Contingency Tables
#######################################################

olink_QTLs <- olink_QTLs %>%
  dplyr::select("type", "Class")

olink_QTLs <- olink_QTLs[which(olink_QTLs$type == "concordant" | olink_QTLs$type == "discordant"),]

# Unique Blood Groups
subgroups <- unique(olink_QTLs$Class)

# Prepare results list
results <- list()

for (sub in subgroups) {
  # Create logical vectors
  in_sub <- olink_QTLs$Class == sub
  in_Con <- olink_QTLs$type == 'concordant'
  in_Dis <- olink_QTLs$type == 'discordant'
  
  # Build the 2x2 table
  a <- sum(in_sub & in_Con)
  b <- sum(in_sub & in_Dis)
  c <- sum(!in_sub & in_Con)
  d <- sum(!in_sub & in_Dis)
  
  table_2x2 <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)

#######################################################
# Fisher's exact test
#######################################################
  
  # Fisher's exact test
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

# Create gt table
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

gtsave(results_gt, file.path(docs_data, "Olink/Olink_Fisher_Test.png"))
