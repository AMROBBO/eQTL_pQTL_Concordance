#######################################################
# 5a. Enrichment Analysis of Colocalised QTL Pairs
# Datasets:
#   a. Open Targets: Baseline Expression (Blood vs non blood)
##
# As there is little variation between concordance of datasets selected on pQTLs
# and eQTLs and those with all available independent SNPs and those with only 1
# SNP per gene, the rest of the enrichment analysis is carried out on all available
# independent SNPs selected from pQTL data:
# pQTL_concordance.csv from 4_Naive_Concordance_Test.R
##
#######################################################

#######################################################
#Load in libraries
#######################################################

library(dotenv)
library(arrow)               # Read parquet files
library(data.table)          # Fast file reading/writing
library(dplyr)               # Data wrangling
library(tidyr)               # Tidying data
library(ggplot2)             # Plotting
library(patchwork)           # Combining plots
library(gt)                  # Pretty tables
library(grid)                # For unit()
library(extrafont)           # Font handling for publication-quality plots

# Load fonts for postscript output

#font_import()               # Run once if fonts are not installed
loadfonts(device = "postscript")

#######################################################
# Initialising file paths
#######################################################

load_dot_env("config.env")

raw_data <- Sys.getenv("rawdatadir")
interim_data <- Sys.getenv("interimdatadir")
processed_data <- Sys.getenv("processeddatadir")
docs_data <- Sys.getenv("docsdir")

#######################################################
# Read in colocalised QTL pairs 
#######################################################

QTLs <- fread(file.path(processed_data, "Concordance/pQTL_concordance.csv"))

# Gene names
gene_names <- fread(file.path(interim_data, "gene_names.csv"))

#######################################################
# Open Targets Data: Baseline expression (blood vs non-blood)
#######################################################

OT_files <- list.files(file.path(raw_data, "open_targets_expression"), full.names = TRUE) #Expression files from Open Target
OT_files <- OT_files[-1] #First file is empty

# Read and combine parquet files

OT_data <- read_parquet(OT_files[1], as_tibble = TRUE)

for (f in OT_files[-1]){
  tmp <- read_parquet(f, as_tibble = TRUE)
  OT_data <- bind_rows(OT_data, tmp)
}

#######################################################
# Labeling QTL data with Ensembl gene names
#######################################################

colnames(gene_names) <- c("id", "Gene_name")

QTLs <- QTLs %>% 
  left_join(gene_names, by = "id")

#######################################################
# Labeling QTL data for Blood Expression
#######################################################

QTLs_OT <- QTLs %>%
  dplyr::select("rsid", "id", "Gene_name", "type") %>% 
  filter(Gene_name %in% OT_data$id) # Filter for QTLs in OT data

OT_data <- OT_data %>% 
  filter(id %in% QTLs_OT$Gene_name)

QTLs_OT$blood_expression <- NA
QTLs_OT$nonblood_expression <- NA
QTLs_OT$top_organ <- NA

#Extract top blood and non-blood expression and most expressed organ for each gene

for (i in unique(QTLs_OT$Gene_name)){  
  
  QTL_row <- which(QTLs_OT$Gene_name == i)
  OT_row <- which(OT_data$id == i)
  
  OT_table <- as.data.frame(OT_data[[2]][[OT_row]])
  OT_rna <- as.data.frame(OT_table[[5]])
  
  bexp <- character()
  bnonexp <- character()
  
  for (j in 1:nrow(OT_table)){
    if (length(grep("blood", OT_table$organs[[j]])) > 0){ 
      bexp <- append(bexp, OT_rna[j,1])        # All blood associated expression
    }else{ 
      bnonexp <- append(bnonexp, OT_rna[j, 1]) # All non blood associated expression
    }
  }
  
  QTLs_OT$blood_expression[QTL_row] <- max(as.numeric(bexp))                                          # Highest blood associated expression
  QTLs_OT$nonblood_expression[QTL_row] <- max(as.numeric(bnonexp))                                    # Highest non blood associated expression
  QTLs_OT$top_organ[QTL_row] <- OT_table$organs[grep(max(OT_table$rna$value), OT_table$rna$value)][1] # Top expressed organ
  
  print(i)
}

#######################################################
# Assign Expression Thresholds
#######################################################

thresholds <- c(0, 10, 100, 1000, 10000, 100000)

QTLs_OT$expression_threshold <- 0
for (i in thresholds){
  QTLs_OT[which(QTLs_OT$blood_expression > i),]$expression_threshold <- i
}

# Save Labelled QTL

QTLs_OT <- QTLs_OT %>%
  rename(
    Protein = id,
    Concordance_Group = type,
    `Blood Expression/TPM` = blood_expression,
    Threshold = expression_threshold
  )

fwrite(QTLs_OT, file.path(processed_data, "Open_Target/OT_labelled_QTLs.csv"))

#######################################################
#Creating OT Summary Table
#######################################################

OT_summary <- data.frame(Blood_Expression_Threshold = thresholds)

for (i in 1:length(thresholds)){
  limit <- thresholds[i]
  
  OT_summary$Concordant[i] <- sum(QTLs_OT$Concordance_Group == "concordant" & QTLs_OT$`Blood Expression/TPM` >= limit)
  OT_summary$Discordant[i] <- sum(QTLs_OT$Concordance_Group == "discordant" & QTLs_OT$`Blood Expression/TPM` >= limit)
  OT_summary$eQTL_dropped[i] <- sum(QTLs_OT$Concordance_Group == "eQTL_dropped" & QTLs_OT$`Blood Expression/TPM` >= limit)
  OT_summary$pQTL_dropped[i] <- sum(QTLs_OT$Concordance_Group == "pQTL_dropped" & QTLs_OT$`Blood Expression/TPM` >= limit)
}

#Save summary as a table

OT_summary_tab <- OT_summary %>%
  gt() %>%
  tab_header(title = md("Blood Expression Analysis"),
             subtitle = md("Open Target Data")) %>%
  fmt_number(
    columns = 1,
    decimals = 0,
    use_seps = TRUE
  )

gtsave(OT_summary_tab, file.path(docs_data, "Open_Target/OT_summary.png"))

#######################################################
# Plotting Blood Expression
#######################################################

# Proportions per concordance group

total_counts <- table(QTLs_OT$Concordance_Group)

OT_plot <- OT_summary %>%
  mutate(
    Concordant = Concordant / total_counts["concordant"] * 100,
    Discordant = Discordant / total_counts["discordant"] * 100,
    eQTL_dropped = eQTL_dropped / total_counts["eQTL_dropped"] * 100,
    pQTL_dropped = pQTL_dropped / total_counts["pQTL_dropped"] * 100
  )

OT_plot <- OT_plot %>%
  pivot_longer(cols = c(Concordant, Discordant, eQTL_dropped, pQTL_dropped),
               names_to = "type_plot", values_to = "prop_plot") %>%
  mutate(blood_plot = c(rep(">0", 4), rep(">10", 4), rep(">100", 4), rep(">1000", 4), rep(">10000", 4), rep(">100000", 4))) 

# Open Target Plot 1a - all Concordance groups

OT_prop_plot <- ggplot(OT_plot, aes(x = blood_plot, y = prop_plot, fill = type_plot)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(text = element_text(size = 8, family = "Times"),
        axis.text.x = element_text(size = 8)) +
  labs(y = "Proportion of Concordance Category", x = "Blood Expression/TPM", fill = "Concordance Category")

# Open Target Plot 1b - Con/Dis groups only

OT_plot_discon <- OT_plot %>% 
  filter(type_plot == "Concordant" | type_plot == "Discordant")

OT_prop_plot_discon <- ggplot(OT_plot_discon, aes(x=blood_plot, y=prop_plot, fill=type_plot)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(text = element_text(size=8, family = "Times"),
        axis.text.x = element_text(size=8)) +
  labs(y = "Proportion of Concordance Category", x = "Blood Expression/TPM", fill = "Concordance Category")

#Plotting concordance categories as a proportion of each blood group
#Proportion = N QTLs in blood group / total number of QTLs in blood group

OT_plot_stacked <- OT_summary

for (i in 1:nrow(OT_plot_stacked)){
  blood_sum <- sum(OT_plot_stacked[i,2:5]) #total number of QTLs in blood group
  OT_plot_stacked$Concordant[i] <- (OT_plot_stacked$Concordant[i]/blood_sum) * 100
  OT_plot_stacked$Discordant[i] <- (OT_plot_stacked$Discordant[i]/blood_sum) * 100
  OT_plot_stacked$eQTL_dropped[i] <- (OT_plot_stacked$eQTL_dropped[i]/blood_sum) * 100
  OT_plot_stacked$pQTL_dropped[i] <- (OT_plot_stacked$pQTL_dropped[i]/blood_sum) * 100
}

#Rearrange to plot

OT_plot_stacked <- OT_plot_stacked %>%
  pivot_longer(cols = c(Concordant, Discordant, eQTL_dropped, pQTL_dropped),
               names_to = "type_plot", values_to = "prop_plot") %>%
  mutate(blood_plot = c(rep(">0", 4), rep(">10", 4), rep(">100", 4), rep(">1000", 4), rep(">10000", 4), rep(">100000", 4)),
         num_plot = c(OT_summary[1,2:5], OT_summary[2,2:5], OT_summary[3,2:5], OT_summary[4,2:5], OT_summary[5,2:5], OT_summary[6,2:5])) 

OT_plot_stacked$num_plot <- gsub("^0$", "", OT_plot_stacked$num_plot)

# Open Target Plot 2 - Stacked plot

OT_prop_plot_stacked <- ggplot(OT_plot_stacked, aes(fill=type_plot, y=prop_plot, x=blood_plot, label = num_plot)) + 
  geom_bar(position="stack", stat="identity") +
  theme(text = element_text(size=8, family = "Times"),
        axis.text.x = element_text(size=8)) +
  geom_text(size = 2, position = position_stack(vjust = 0.5), check_overlap = TRUE, family = "Times") +
  labs(y = "Proportion of Blood Expression Group", x = "Blood Expression/TPM", fill = "Concordance Category")

#######################################################
#Saving OT Figures
#######################################################

###FIGURE 2

fig2a <- OT_prop_plot +
  theme(legend.position = "none") +
  coord_flip()

fig2b <- OT_prop_plot_stacked +
  theme(axis.title.y = element_blank()) +
  coord_flip()

fig2 <- fig2a + fig2b

fig2 <- fig2 & theme(plot.margin = unit(c(2, 2, 2, 2), "pt"))

ggsave(file = file.path(docs_data, "Open_Target/Fig2.eps"), plot = fig2, width = 5.2, height = 4, 
       units = "in", dpi = 600, family = "Times", device = postscript, paper = "special", horizontal = FALSE)

###FIGURE S1

figs1 <- OT_prop_plot_discon +
  coord_flip()

ggsave(file = file.path(docs_data, "Open_Target/FigS1.eps"), plot = figs1, width = 3.5, height = 4, 
       units = "in", dpi = 600, family = "Times", device = postscript, paper = "special", horizontal = FALSE)
