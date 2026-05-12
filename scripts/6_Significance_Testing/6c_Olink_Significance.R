#######################################################
# 6c. Significance Testing for Olink Enrichment Analysis
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
# Reading in Olink labelled QTLs
#######################################################

olink_QTLs <- fread(file.path(processed_data, "Olink/Olink_labelled_QTLs.csv"))

#######################################################
# Concordance vs Discordant
#######################################################
#######################################################
# Creating 2X2 Contingency Tables
#######################################################

olink_QTLs_discon <- olink_QTLs %>%
  dplyr::select("type", "Class") %>% 
  filter(type == "concordant" | type == "discordant")

# Unique Protein Groups
subgroups <- unique(olink_QTLs_discon$Class)

# Prepare results list
results <- list()

for (sub in subgroups) {
  # Create logical vectors
  in_sub <- olink_QTLs_discon$Class == sub
  in_Con <- olink_QTLs_discon$type == 'concordant'

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

# Apply multiple testing correction
results_df <- results_df %>%
  mutate(Adj_P_Bonferroni = p.adjust(P_Value, method = "bonferroni")) %>%
  mutate(Adj_P_BH_FDR = p.adjust(P_Value, method = "BH")) %>% 
  dplyr::select("Subgroup", "P_Value", "Adj_P_Bonferroni", "Adj_P_BH_FDR")

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

gtsave(results_gt, file.path(docs_data, "Olink/Olink_Fisher_Test.png"))

#######################################################
# All Categories
#######################################################

#######################################################
# Creating Summary Table
#######################################################

olink_QTLs_all <- dcast(
  olink_QTLs[, .N, by = .(Class, type)],
  Class ~ type,
  value.var = "N",
  fill = 0
)

pvals <- c()

#######################################################
# Creating 4X2 Contingency Tables
#######################################################

for (i in 1:nrow(olink_QTLs_all)){
  group <- as.numeric(olink_QTLs_all[i, -1])
  
  others <- as.numeric(colSums(olink_QTLs_all[-i, -1]))
  
  tbl <- rbind(group, others)
  
  colnames(tbl) <- colnames(olink_QTLs_all[ ,-1])
  
  rownames(tbl) <- c(olink_QTLs_all$Class[i], "Others")
  
#######################################################
# Fisher's Exact Test
#######################################################
  
  test <- fisher.test(tbl)
  
  pvals[i] <- test$p.value
}

# Combine
results <- data.frame(
  Group = olink_QTLs_all$Class,
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

gtsave(results_gt_all, file.path(docs_data, "Olink/Olink_Fisher_Test_All.png"))
