#######################################################
# 6d. Significance Testing for Reactome Enrichment Analysis
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
# Reading in Reactome labelled QTLs
#######################################################

# Reading in Reactome labelled QTLs
QTLs_reactome <- fread(file.path(processed_data, "Reactome/Reactome_labelled_QTLs.csv"))

#######################################################
# Creating 2X2 Contingency Tables
#######################################################

QTLs_reactome <- QTLs_reactome %>%
  dplyr::select("type", "Event (Pathway or Reaction) Name")
QTLs_reactome <- QTLs_reactome[which(QTLs_reactome$type == "concordant" | QTLs_reactome$type == "discordant"),]

# Unique Blood Groups
subgroups <- unique(QTLs_reactome$`Event (Pathway or Reaction) Name`)

# Prepare results list
results <- list()

for (sub in subgroups) {
  # Create logical vectors
  in_sub <- QTLs_reactome$`Event (Pathway or Reaction) Name` == sub
  in_Con <- QTLs_reactome$type == 'concordant'
  in_Dis <- QTLs_reactome$type == 'discordant'
  
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

results_gt <- results_df %>%
  gt() %>%
  tab_header(
    title = md("Fisher's Exact Test Results"),
    subtitle = md("Pathway Involvement")
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

gtsave(results_gt, file.path(docs_data, "Reactome/Reactome_Fisher_Test.png"))

#######################################################
# For significant pathways
#######################################################

reactome_summary <- fread(file.path(processed_data, "Reactome/Reactome_summary.csv"))

# Summary table for pathways with >10 Genes involved

reactome_summary <- reactome_summary[rowSums(reactome_summary[,2:5]) > 10,]

results_df_sig <- results_df %>% 
  filter(Subgroup %in% reactome_summary$`Event (Pathway or Reaction) Name`) %>% 
  arrange(Subgroup)

#######################################################
# Create gt table
#######################################################

results_gt_sig <- results_df_sig %>%
  gt() %>%
  tab_header(
    title = md("Fisher's Exact Test Results"),
    subtitle = md("Pathway Involvement")
  ) %>%
  fmt_number(
    columns = c(P_Value),
    decimals = 4
  ) %>%
  cols_label(
    Subgroup = "Protein Class",
    P_Value = "P-Value",
    Adj_P_Bonferroni = "Adjusted P-Value"
  ) %>%
  tab_options(
    data_row.padding = px(2),
    table.font.size = px(11)
  )

#######################################################
# Save
#######################################################

gtsave(results_gt_sig, file.path(docs_data, "Reactome/Reactome_Fisher_Test_Sig.png"))
