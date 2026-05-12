#######################################################
# 6e. Significance Testing for Ensembl VEP Enrichment Analysis
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
# Reading in Ensembl labelled QTLs
#######################################################

con_QTLs <- fread(file.path(processed_data, "Ensembl_VEP/VEP_labelled_QTLs_con.csv"))
dis_QTLs <- fread(file.path(processed_data, "Ensembl_VEP/VEP_labelled_QTLs_dis.csv"))
pQTL_dropped_QTLs <- fread(file.path(processed_data, "Ensembl_VEP/VEP_labelled_QTLs_pQTL_dropped.csv"))
eQTL_dropped_QTLs <- fread(file.path(processed_data, "Ensembl_VEP/VEP_labelled_QTLs_eQTL_dropped.csv"))

#######################################################
# Creating 2X2 Contingency Tables
#######################################################

VEP_summary <- colSums(con_QTLs[,7:ncol(con_QTLs)])
VEP_summary_2 <- colSums(dis_QTLs[,7:ncol(dis_QTLs)])
VEP_summary_3 <- colSums(pQTL_dropped_QTLs[,7:ncol(pQTL_dropped_QTLs)])
VEP_summary_4 <- colSums(eQTL_dropped_QTLs[,7:ncol(eQTL_dropped_QTLs)])

#Combine

VEP_summary <- merge(VEP_summary, VEP_summary_2, by = 0, all = T)
VEP_summary_5 <- merge(VEP_summary_3, VEP_summary_4, by = 0, all = T)

VEP_summary <- merge(VEP_summary, VEP_summary_5, by = "Row.names", all = T)

VEP_summary[is.na(VEP_summary)] <- 0
colnames(VEP_summary) <- c("Variant Effect", "Concordant", "Discordant", "pQTL_Dropped", "eQTL_Dropped")

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
  mutate(Adj_P_BH_FDR = p.adjust(P_Value, method = "BH")) %>% 
  dplyr::select("Subgroup", "Odds_Ratio", "P_Value", "Adj_P_Bonferroni", "Adj_P_BH_FDR")

#######################################################
# Create gt table
#######################################################

results_gt <- results_df %>%
  gt() %>%
  tab_header(
    title = md("Fisher's Exact Test Results"),
    subtitle = md("Predicted Variant Effects")
  ) %>%
  fmt_number(
    columns = c(P_Value, Adj_P_Bonferroni, Adj_P_BH_FDR),
    decimals = 4
  ) %>%
  cols_label(
    Subgroup = "Predicted Variant Effect",
    P_Value = "P-Value",
    Adj_P_Bonferroni = "Bonferroni Adjusted P-Value",
    Adj_P_BH_FDR = "FDR Adjusted P-Value"
  )

#######################################################
# Save
#######################################################

gtsave(results_gt, file.path(docs_data, "Ensembl_VEP/Ensembl_Fisher_Test.png"))

#######################################################
# All Categories
#######################################################

#######################################################
# Creating 4X2 Contingency Tables
#######################################################

pvals <- c()

for (i in 1:nrow(VEP_summary)){
  group <- as.numeric(VEP_summary[i, -1])
  
  others <- as.numeric(colSums(VEP_summary[-i, -1]))
  
  tbl <- rbind(group, others)
  
  colnames(tbl) <- colnames(VEP_summary[ ,-1])
  
  rownames(tbl) <- c(VEP_summary$`Variant Effect`[i], "Others")
  
#######################################################
# Fisher's Exact Test
#######################################################
  
  test <- fisher.test(tbl)
  
  pvals[i] <- test$p.value
}

# Combine
results <- data.frame(
  Group = VEP_summary$`Variant Effect`,
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
    subtitle = md("Predicted Variant Effect")
  ) %>%
  fmt_number(
    columns = c(P_Value, Adj_P_Bonferroni, Adj_P_BH_FDR),
    decimals = 4
  ) %>%
  cols_label(
    Group = "Predicted Variant Effect",
    P_Value = "P-Value",
    Adj_P_Bonferroni = "Bonferroni Adjusted P-Value",
    Adj_P_BH_FDR = "FDR Adjusted P-Value"
  )

#######################################################
# Save
#######################################################

gtsave(results_gt_all, file.path(docs_data, "Ensembl_VEP/Ensembl_Fisher_Test_All.png"))
