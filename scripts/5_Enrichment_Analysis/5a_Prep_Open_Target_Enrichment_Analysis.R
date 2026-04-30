#######################################################
# 5a pre-step - 
# Extracting Ensemble gene names in order to extract data from Open Target
#######################################################

#######################################################
# Load Libraries
#######################################################

library(dotenv)
library(data.table)
library(dplyr)
library(tidyr)

#######################################################
# Initialising file paths
#######################################################

load_dot_env("config.env")

interim_data <- Sys.getenv("interimdatadir")

#######################################################
# Read in QTL pairs
#######################################################

joined <- fread(file.path(interim_data, "eQTL_pQTL_matched.csv"))

#######################################################
# Extract Gene names
#######################################################

gene_names <- joined %>% select(c("GeneSymbol_eQTLs", "Gene_eQTLs"))
gene_names <- unique(gene_names)

#######################################################
# Save
#######################################################

fwrite(gene_names, file.path(interim_data, "gene_names.csv"))