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
# Concordance vs Discordant
#######################################################
#######################################################
# Creating 2X2 Contingency Tables
#######################################################

QTLs_OT_discon <- QTLs_OT %>%
  dplyr::select("Concordance_Group", "top_organ") %>% 
  filter(Concordance_Group == "concordant" | Concordance_Group == "discordant")

# Unique Blood Groups
subgroups <- unique(QTLs_OT_discon$top_organ)

# Prepare results list
results <- list()

for (sub in subgroups) {
  # Create logical vectors
  in_sub <- QTLs_OT_discon$top_organ == sub
  in_Con <- QTLs_OT_discon$Concordance_Group == 'concordant'
  
  # Build the 2x2 table
  a <- sum(in_sub & in_Con)
  b <- sum(in_sub & !in_Con)
  c <- sum(!in_sub & in_Con)
  d <- sum(!in_sub & !in_Con)
  
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
    subtitle = md("Organ Expression")
  ) %>%
  fmt_number(
    columns = c(P_Value, Adj_P_Bonferroni, Adj_P_BH_FDR),
    decimals = 4
  ) %>%
  cols_label(
    Subgroup = "Blood Expression",
    P_Value = "P-Value",
    Adj_P_Bonferroni = "Bonferroni Adjusted P-Value",
    Adj_P_BH_FDR = "FDR Adjusted P-Value"
  )

#######################################################
# Save
#######################################################

gtsave(results_gt, file.path(docs_data, "Open_Target/OT_Organ_Fisher_Test.png"))

#######################################################
# All Categories
#######################################################

#######################################################
# Creating Summary Table
#######################################################

QTLs_OT_all <- dcast(
  QTLs_OT[, .N, by = .(top_organ, Concordance_Group)],
  top_organ ~ Concordance_Group,
  value.var = "N",
  fill = 0
)

pvals <- c()

#######################################################
# Creating 4X2 Contingency Tables
#######################################################

for (i in 1:nrow(QTLs_OT_all)){
  group <- as.numeric(QTLs_OT_all[i, -1])
  
  others <- as.numeric(colSums(QTLs_OT_all[-i, -1]))
  
  tbl <- rbind(group, others)
  
  colnames(tbl) <- colnames(QTLs_OT_all[ ,-1])
  
  rownames(tbl) <- c(QTLs_OT_all$top_organ[i], "Others")
  
  #######################################################
  # Fisher's Exact Test
  #######################################################
  
  test <- fisher.test(tbl)
  
  pvals[i] <- test$p.value
}

# Combine
results <- data.frame(
  Group = QTLs_OT_all$top_organ,
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
    subtitle = md("Organ Expression")
  ) %>%
  fmt_number(
    columns = c(P_Value, Adj_P_Bonferroni, Adj_P_BH_FDR),
    decimals = 4
  ) %>%
  cols_label(
    Group = "Blood Expression",
    P_Value = "P-Value",
    Adj_P_Bonferroni = "Bonferroni Adjusted P-Value",
    Adj_P_BH_FDR = "FDR Adjusted P-Value"
  )

#######################################################
# Save
#######################################################

gtsave(results_gt_all, file.path(docs_data, "Open_Target/OT_Organ_Fisher_Test_All.png"))
