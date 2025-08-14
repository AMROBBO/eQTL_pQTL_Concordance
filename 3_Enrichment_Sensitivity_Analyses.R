##Creating 2x2 contigency tables and running the Fisher's exact test to see significance in the enrichment analyses

#######################################################
#Load in libraries
#######################################################

library(data.table)
library(dplyr)
library(gt)

#######################################################
#Set Working Directory
#######################################################

setwd("")

#######################################################
#Open Target Enrichment
#######################################################

# Reading in OT labelled QTLs
QTLs_OT <- fread("OT_labelled_QTLs.csv")

QTLs_OT <- QTLs_OT %>%
  dplyr::select("Concordance_Group", "Blood_Group")
QTLs_OT <- QTLs_OT[which(QTLs_OT$Concordance_Group == "concordant" | QTLs_OT$Concordance_Group == "discordant"),]

# Unique Blood Groups
subgroups <- unique(QTLs_OT$Blood_Group)

# Prepare results list
results <- list()

for (sub in subgroups) {
  # Create logical vectors
  in_sub <- QTLs_OT$Blood_Group == sub
  in_Con <- QTLs_OT$Concordance_Group == 'concordant'
  in_Dis <- QTLs_OT$Concordance_Group == 'discordant'
  
  # Build the 2x2 table
  a <- sum(in_sub & in_Con)
  b <- sum(in_sub & in_Dis)
  c <- sum(!in_sub & in_Con)
  d <- sum(!in_sub & in_Dis)
  
  table_2x2 <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  
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

# Create gt table
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

# Save
gtsave(results_gt, "")

#######################################################
#Human Protein Atlas Enrichment
#######################################################

# Reading in HPA labelled QTLs
QTLs_HPA <- fread("HPA_labelled_QTLs.csv")

#Separating QTLs into concordance type

QTLs_HPA_con <- QTLs_HPA[which(QTLs_HPA$type == "concordant"),][,6:29]
QTLs_HPA_dis <- QTLs_HPA[which(QTLs_HPA$type == "discordant"),][,6:29]

#Summing number of QTLs in each protein class for each concordance group

HPA_summary <- colSums(QTLs_HPA_con)
HPA_summary <- cbind(HPA_summary, colSums(QTLs_HPA_dis))

HPA_summary <- as.data.frame(HPA_summary)
HPA_summary$Protein_Class <- rownames(HPA_summary)

colnames(HPA_summary) <- c("Concordant", "Discordant", "Protein Class")

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

# Save
gtsave(results_gt, "")

#######################################################
#Olink Enrichment
#######################################################

# Reading in Olink labelled QTLs
olink_QTLs <- fread("Olink_labelled_QTLs.csv")

olink_QTLs <- olink_QTLs %>%
  dplyr::select("type", "class")
olink_QTLs <- olink_QTLs[which(olink_QTLs$type == "concordant" | olink_QTLs$type == "discordant"),]

# Unique Blood Groups
subgroups <- unique(olink_QTLs$class)

# Prepare results list
results <- list()

for (sub in subgroups) {
  # Create logical vectors
  in_sub <- olink_QTLs$class == sub
  in_Con <- olink_QTLs$type == 'concordant'
  in_Dis <- olink_QTLs$type == 'discordant'
  
  # Build the 2x2 table
  a <- sum(in_sub & in_Con)
  b <- sum(in_sub & in_Dis)
  c <- sum(!in_sub & in_Con)
  d <- sum(!in_sub & in_Dis)
  
  table_2x2 <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  
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

# Save
gtsave(results_gt, "")

#######################################################
#Reactome Enrichment
#######################################################

# Reading in Reactome labelled QTLs
QTLs_reactome <- fread("Reactome_labelled_QTLs.csv")

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

#Save CSV
fwrite(results_df, "")

# Create gt table
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

# Save
gtsave(results_gt, "")

##For significant pathways

sig_pathways <- fread("sig_pathways.csv") %>%
  dplyr::select("Event (Pathway or Reaction) Name")

results_df_sig <- results_df[results_df$Subgroup %in% sig_pathways$`Event (Pathway or Reaction) Name`,]

# Create gt table
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
  )

# Save
gtsave(results_gt_sig, "")

#######################################################
#Ensembl VEP Enrichment
#######################################################

# Reading in Ensembl labelled QTLs
con_QTLs <- fread("VEP_labelled_QTLs_con.csv")
dis_QTLs <- fread("VEP_labelled_QTLs_dis.csv")

VEP_summary <- colSums(con_QTLs[,11:ncol(con_QTLs)])
VEP_summary_2 <- colSums(dis_QTLs[,11:ncol(dis_QTLs)])

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

# Create gt table
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

# Save
gtsave(results_gt, "")
