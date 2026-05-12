#######################################################
# 6b. Significance Testing for Human Protein Atlas Enrichment Analysis
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
# Reading in HPA labelled QTLs
#######################################################

QTLs_HPA <- fread(file.path(processed_data, "Human_Protein_Atlas/HPA_labelled_QTLs.csv"))

#######################################################
# Summing number of QTLs in each protein class for each concordance group
#######################################################

QTLs_HPA_con <- QTLs_HPA[which(QTLs_HPA$type == "concordant"),][,5:ncol(QTLs_HPA)] %>% 
  colSums()
QTLs_HPA_dis <- QTLs_HPA[which(QTLs_HPA$type == "discordant"),][,5:ncol(QTLs_HPA)] %>% 
  colSums()
QTLs_HPA_pQTL_dropped <- QTLs_HPA[which(QTLs_HPA$type == "pQTL_dropped"),][,5:ncol(QTLs_HPA)] %>% 
  colSums()
QTLs_HPA_eQTL_dropped <- QTLs_HPA[which(QTLs_HPA$type == "eQTL_dropped"),][,5:ncol(QTLs_HPA)] %>% 
  colSums()

HPA_summary <- data.frame(
  Protein_Class = colnames(QTLs_HPA)[5:ncol(QTLs_HPA)],
  Concordant = QTLs_HPA_con,
  Discordant = QTLs_HPA_dis,
  pQTL_Dropped = QTLs_HPA_pQTL_dropped,
  eQTL_Dropped = QTLs_HPA_eQTL_dropped
)

#######################################################
# Concordance vs Discordant
#######################################################
#######################################################
# Creating 2X2 Contingency Tables
#######################################################

# Unique Protein Classes
subgroups <- unique(HPA_summary$Protein_Class)

# Prepare results list
results <- list()

for (sub in subgroups) {
  # Create logical vectors
  in_sub <- HPA_summary %>% 
    filter(Protein_Class == sub)
  out_sub <- HPA_summary %>% 
    filter_out(Protein_Class == sub)
  
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
  mutate(Adj_P_BH_FDR = p.adjust(P_Value, method = "BH")) %>% 
  dplyr::select("Subgroup", "P_Value", "Adj_P_Bonferroni", "Adj_P_BH_FDR")

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
    columns = c(P_Value, Adj_P_Bonferroni, Adj_P_BH_FDR),
    decimals = 4
  ) %>%
  cols_label(
    Subgroup = "Protein Class",
    P_Value = "P-Value",
    Adj_P_Bonferroni = "Bonferroni Adjusted P-Value",
    Adj_P_BH_FDR = "FDR Adjusted P-Value"
  )

#######################################################
# Save
#######################################################

gtsave(results_gt, file.path(docs_data, "Human_Protein_Atlas/HPA_Fisher_Test.png"))

#######################################################
# All Categories
#######################################################

#######################################################
# Creating 4X2 Contingency Tables
#######################################################

pvals <- c()

for (i in 1:nrow(HPA_summary)){
  group <- as.numeric(HPA_summary[i, -1])
  
  others <- as.numeric(colSums(HPA_summary[-i, -1]))
  
  tbl <- rbind(group, others)
  
  colnames(tbl) <- colnames(HPA_summary[ ,-1])
  
  rownames(tbl) <- c(HPA_summary$Protein_Class[i], "Others")
  
#######################################################
# Fisher's Exact Test
#######################################################
  
  test <- fisher.test(tbl)
  
  pvals[i] <- test$p.value
}

# Combine
results <- data.frame(
  Group = HPA_summary$Protein_Class,
  P_Value = pvals
)

# Apply Multiple Testing Correction
results <- results %>%
  mutate(Adj_P_Bonferroni = p.adjust(P_Value, method = "bonferroni")) %>%
  mutate(Adj_P_BH_FDR = p.adjust(P_Value, method = "BH"))

#######################################################
# Create gt table
#######################################################

# Create gt table
results_gt_all <- results %>%
  gt() %>%
  tab_header(
    title = md("Fisher's Exact Test Results"),
    subtitle = md("Protein Class")
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

gtsave(results_gt_all, file.path(docs_data, "Human_Protein_Atlas/HPA_Fisher_Test_All.png"))
