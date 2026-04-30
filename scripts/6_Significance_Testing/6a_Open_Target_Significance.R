#######################################################
# 6a. Significance Testing for Open Target Enrichment Analysis
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
# Reading in OT labelled QTLs
#######################################################

QTLs_OT <- fread(file.path(processed_data, "Open_Target/OT_labelled_QTLs.csv"))

#######################################################
# Creating 2X2 Contingency Tables
#######################################################

QTLs_OT <- QTLs_OT %>%
  dplyr::select("Concordance_Group", "Threshold")

QTLs_OT <- QTLs_OT[which(QTLs_OT$Concordance_Group == "concordant" | QTLs_OT$Concordance_Group == "discordant"),]

# Unique Blood Groups
subgroups <- unique(QTLs_OT$Threshold)

# Prepare results list
results <- list()

for (sub in subgroups) {
  # Create logical vectors
  in_sub <- QTLs_OT$Threshold == sub
  in_Con <- QTLs_OT$Concordance_Group == 'concordant'
  in_Dis <- QTLs_OT$Concordance_Group == 'discordant'
  
  # Build the 2x2 table
  a <- sum(in_sub & in_Con)
  b <- sum(in_sub & in_Dis)
  c <- sum(!in_sub & in_Con)
  d <- sum(!in_sub & in_Dis)
  
  table_2x2 <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  
#######################################################
# Fisher's exact test
#######################################################
  
  fisher_result <- fisher.test(table_2x2)
  
  sub <- paste0(">", sub, "TMP")
  
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

results_df$Subgroup <- factor(results_df$Subgroup, levels = c(">0TMP", ">10TMP", ">100TMP", ">1000TMP", ">10000TMP", ">1e+05TMP"))
results_df <- results_df[order(results_df$Subgroup), ]

#######################################################
# Create gt table
#######################################################

results_gt <- results_df %>%
  gt() %>%
  tab_header(
    title = md("Fisher's Exact Test Results"),
    subtitle = md("Blood Expression")
  ) %>%
  fmt_number(
    columns = c(P_Value),
    decimals = 4
  ) %>%
  cols_label(
    Subgroup = "Blood Group",
    P_Value = "P-Value",
    Adj_P_Bonferroni = "Adjusted P-Value"
  )

#######################################################
# Save
#######################################################

gtsave(results_gt, file.path(docs_data, "Open_Target/OT_Fisher_Test.png"))
