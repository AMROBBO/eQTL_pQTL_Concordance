#######################################################
# 6d. Significance Testing for Reactome Enrichment Analysis
# Concordant vs Discordant:
# Creating 2x2 contingency tables and running the Fisher's exact test to 
# see significance in the enrichment analyses
# All categories:
# Comparing distribution of 4 categories within each class against distribution
# within all other classes - Fisher test
#######################################################
# To note:

# Some of the genes are members of multiple group. This is against the assumptions
# of the fisher test. It is worth noting this and interpreting with caution

#######################################################
#Load in libraries
#######################################################

library(dotenv)
library(data.table)
library(dplyr)
library(gt)
library(webshot)

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
# Concordance vs Discordant
#######################################################
#######################################################
# Creating 2X2 Contingency Tables
#######################################################

QTLs_reactome_discon <- QTLs_reactome %>%
  dplyr::select("type", "Event (Pathway or Reaction) Name") %>% 
  filter(type == "concordant" | type == "discordant")

# Unique Blood Groups
subgroups <- unique(QTLs_reactome_discon$`Event (Pathway or Reaction) Name`)

# Prepare results list
results <- list()

for (sub in subgroups) {
  # Create logical vectors
  in_sub <- QTLs_reactome_discon$`Event (Pathway or Reaction) Name` == sub
  in_Con <- QTLs_reactome_discon$type == 'concordant'

  # Build the 2x2 table
  a <- sum(in_sub & in_Con)
  b <- sum(in_sub & !in_Con)
  c <- sum(!in_sub & in_Con)
  d <- sum(!in_sub & !in_Con)
  
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
  mutate(Adj_P_BH_FDR = p.adjust(P_Value, method = "BH")) %>% 
  dplyr::select("Subgroup", "P_Value", "Adj_P_Bonferroni", "Adj_P_BH_FDR")

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
    columns = c(P_Value, Adj_P_Bonferroni, Adj_P_BH_FDR),
    decimals = 4
  ) %>%
  cols_label(
    Subgroup = "Pathway",
    P_Value = "P-Value",
    Adj_P_Bonferroni = "Bonferroni Adjusted P-Value",
    Adj_P_BH_FDR = "FDR Adjusted P-Value"
  ) %>%
  tab_options(
    data_row.padding = px(2),
    table.font.size = px(11)
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

reactome_summary_sig <- reactome_summary[rowSums(reactome_summary[,2:5]) > 10,]

results_df_sig <- results_df %>% 
  filter(Subgroup %in% reactome_summary_sig$`Event (Pathway or Reaction) Name`) %>% 
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
    columns = c(P_Value, Adj_P_Bonferroni, Adj_P_BH_FDR),
    decimals = 4
  ) %>%
  cols_label(
    Subgroup = "Pathway",
    P_Value = "P-Value",
    Adj_P_Bonferroni = "Bonferroni Adjusted P-Value",
    Adj_P_BH_FDR = "FDR Adjusted P-Value"
  ) %>%
  tab_options(
    data_row.padding = px(2),
    table.font.size = px(11)
  )

#######################################################
# Save
#######################################################

gtsave(results_gt_sig, file.path(docs_data, "Reactome/Reactome_Fisher_Test_Sig.png"))

#######################################################
# All Categories
#######################################################

#######################################################
# Creating 4X2 Contingency Tables
#######################################################

pvals <- c()

for (i in 1:nrow(reactome_summary)){
  group <- as.numeric(reactome_summary[i, -1])
  
  others <- as.numeric(colSums(reactome_summary[-i, -1]))
  
  tbl <- rbind(group, others)
  
  colnames(tbl) <- colnames(reactome_summary[ ,-1])
  
  rownames(tbl) <- c(reactome_summary$`Event (Pathway or Reaction) Name`[i], "Others")
  
#######################################################
# Fisher's Exact Test
#######################################################
  
  test <- fisher.test(tbl)
  
  pvals[i] <- test$p.value
}

# Combine
results <- data.frame(
  Group = reactome_summary$`Event (Pathway or Reaction) Name`,
  P_Value = pvals
)

# Apply Multiple Testing Correction
results <- results %>%
  mutate(Adj_P_Bonferroni = p.adjust(P_Value, method = "bonferroni")) %>%
  mutate(Adj_P_BH_FDR = p.adjust(P_Value, method = "BH"))

results <- results[order(results$Adj_P_BH_FDR), ]

#######################################################
# Create gt table
#######################################################

# Create gt table
results_gt_all <- results %>%
  gt() %>%
  tab_header(
    title = md("Fisher's Exact Test Results"),
    subtitle = md("Pathway Involvement")
  ) %>%
  fmt_number(
    columns = c(P_Value, Adj_P_Bonferroni, Adj_P_BH_FDR),
    decimals = 4
  ) %>%
  cols_label(
    Group = "Protein Class",
    P_Value = "P-Value",
    Adj_P_Bonferroni = "Bonferroni Adjusted P-Value",
    Adj_P_BH_FDR = "FDR Adjusted P-Value"
  )

#######################################################
# Save
#######################################################

gtsave(results_gt_all, file.path(docs_data, "Reactome/Reactome_Fisher_Test_All.png"))

#######################################################
# For significant pathways
#######################################################

results_sig <- results %>% 
  filter(Group %in% reactome_summary_sig$`Event (Pathway or Reaction) Name`) %>% 
  arrange(Group)

#######################################################
# Create gt table
#######################################################

results_gt_all_sig <- results_sig %>%
  gt() %>%
  tab_header(
    title = md("Fisher's Exact Test Results"),
    subtitle = md("Pathway Involvement")
  ) %>%
  fmt_number(
    columns = c(P_Value, Adj_P_Bonferroni, Adj_P_BH_FDR),
    decimals = 4
  ) %>%
  cols_label(
    Group = "Pathway",
    P_Value = "P-Value",
    Adj_P_Bonferroni = "Bonferroni Adjusted P-Value",
    Adj_P_BH_FDR = "FDR Adjusted P-Value"
  ) %>%
  tab_options(
    data_row.padding = px(2),
    table.font.size = px(11)
  )

#######################################################
# Save
#######################################################

gtsave(results_gt_all_sig, file.path(docs_data, "Reactome/Reactome_Fisher_Test_All_Sig.png"))
