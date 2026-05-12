#######################################################
# 5b. Enrichment Analysis of Colocalised QTL Pairs
# Datasets:
#   b. Human Protein Atlas (HPA): Protein Class annotation
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

#######################################################
#Read in HPA Data
#######################################################

HPA <- fread(file.path(raw_data, "proteinatlas.tsv")) %>% 
  dplyr::select("Gene", "Ensembl", "Protein class")

#######################################################
#Labeling QTL data for Protein Class
#######################################################

QTLs_HPA <- QTLs %>%
  dplyr::select("rsid", "id", "type")
colnames(QTLs_HPA)[2] <- "Gene"

#Assigning protein classes to genes based on Ensembl ID

QTLs_HPA <- merge(QTLs_HPA, HPA, by = "Gene")

#Making a count for which protein classes appear for each QTL

QTLs_HPA$`Protein class` <- strsplit(QTLs_HPA$'Protein class', split = ", ")

QTLs_HPA <- QTLs_HPA %>%
  mutate(row = row_number()) %>%  # Add a row identifier
  unnest('Protein class') %>%     # Unnest the list column
  mutate(present = 1) %>%         # Mark presence
  pivot_wider(
    names_from = 'Protein class',
    values_from = present,
    values_fill = 0               # Fill missing combinations with 0
  ) %>%
  dplyr::select(-row)             # Drop helper column

# Save

fwrite(QTLs_HPA, file.path(processed_data, "Human_Protein_Atlas/HPA_labelled_QTLs.csv"))

#######################################################
# Creating HPA Summary Table
#######################################################

#Separating QTLs into concordance type

HPA_summary <- colSums(QTLs_HPA[QTLs_HPA$type == "concordant", 5:ncol(QTLs_HPA)])
HPA_summary <- cbind(HPA_summary, colSums(QTLs_HPA[QTLs_HPA$type == "discordant", 5:ncol(QTLs_HPA)]))
HPA_summary <- cbind(HPA_summary, colSums(QTLs_HPA[QTLs_HPA$type == "eQTL_dropped", 5:ncol(QTLs_HPA)]))
HPA_summary <- cbind(HPA_summary, colSums(QTLs_HPA[QTLs_HPA$type == "pQTL_dropped", 5:ncol(QTLs_HPA)]))

#Summing number of QTLs in each protein class for each concordance group

HPA_summary <- as.data.frame(HPA_summary)
HPA_summary$Protein_Class <- rownames(HPA_summary)
HPA_summary <- HPA_summary %>%
  relocate(Protein_Class, .before = HPA_summary)
colnames(HPA_summary) <- c("Protein Class", "Concordant", "Discordant", "eQTL_dropped", "pQTL_dropped")

#Save

HPA_summary_tab <- HPA_summary%>%
  gt() %>%
  tab_header(title = md("Protein Class Analysis"),
             subtitle = md("Human Protein Atlas"))


gtsave(HPA_summary_tab, file.path(docs_data, "Human_Protein_Atlas/HPA_summary.png"))

#######################################################
#Plotting Protein Classes
#######################################################

#Total number of QTLs in each concordance group

total_counts <- table(QTLs_HPA$type)

HPA_plot <- HPA_summary %>%
  mutate(
    Concordant = Concordant / total_counts["concordant"] * 100,
    Discordant = Discordant / total_counts["discordant"] * 100,
    eQTL_dropped = eQTL_dropped / total_counts["eQTL_dropped"] * 100,
    pQTL_dropped = pQTL_dropped / total_counts["pQTL_dropped"] * 100
  )

#Plotting Protein Class involvement as a proportion of each concordance category:
#Proportion = N QTLs in protein class / total number of QTLs in concordance group

#Rearrange to plot
HPA_plot <- HPA_plot %>%
  pivot_longer(cols = c(Concordant, Discordant, eQTL_dropped, pQTL_dropped),
               names_to = "type_plot", values_to = "prop_plot")

#Human Protein Atlas Plot 1a - all Concordance groups

HPA_prop_plot <- ggplot(HPA_plot, aes(fill=type_plot, y=prop_plot, x=`Protein Class`)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(text = element_text(size=8, family = "Times"),
        axis.text.x = element_text(size=8)) +
  labs(y = "Proportion of Concordance Category", x = "Protein Class", fill = "Concordance Category")

#Human Protein Atlas Plot 1b - Con/Dis Groups only

HPA_plot_discon <- HPA_plot[which(HPA_plot$type_plot == "Concordant" | HPA_plot$type_plot == "Discordant"),]

HPA_prop_plot_discon <- ggplot(HPA_plot_discon, aes(fill=type_plot, y=prop_plot, x=`Protein Class`)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(text = element_text(size=8, family = "Times"),
        axis.text.x = element_text(size=8)) +
  labs(y = "Proportion of Concordance Category", x = "Protein Class", fill = "Concordance Category")

#Plotting concordance categories as a proportion of each protein class
#Proportion = N QTLs in protein class / total number of QTLs in protein class

HPA_plot_stacked <- HPA_summary

for (i in 1:nrow(HPA_plot_stacked)){
  blood_sum <- sum(HPA_plot_stacked[i,2:5])
  HPA_plot_stacked$Concordant[i] <- (HPA_plot_stacked$Concordant[i]/blood_sum) * 100
  HPA_plot_stacked$Discordant[i] <- (HPA_plot_stacked$Discordant[i]/blood_sum) * 100
  HPA_plot_stacked$eQTL_dropped[i] <- (HPA_plot_stacked$eQTL_dropped[i]/blood_sum) * 100
  HPA_plot_stacked$pQTL_dropped[i] <- (HPA_plot_stacked$pQTL_dropped[i]/blood_sum) * 100
}

HPA_num <- HPA_summary %>% 
  pivot_longer(cols = c(Concordant, Discordant, eQTL_dropped, pQTL_dropped),
               names_to = "type_plot", values_to = "num_plot")

#Rearrange to plot
HPA_plot_stacked <- HPA_plot_stacked %>%
  pivot_longer(cols = c(Concordant, Discordant, eQTL_dropped, pQTL_dropped),
               names_to = "type_plot", values_to = "prop_plot") %>% 
  mutate(num_plot = HPA_num$num_plot)

HPA_plot_stacked$num_plot <- gsub("^0$", "", HPA_plot_stacked$num_plot)

#Human Protein Atlas Plot 3 - Stacked plot

HPA_prop_plot_stacked <- ggplot(HPA_plot_stacked, aes(fill=type_plot, y=prop_plot, x=`Protein Class`, label=num_plot)) + 
  geom_bar(position="stack", stat="identity") +
  theme(text = element_text(size=8, family = "Times"),
        axis.text.x = element_text(size=8)) +
  geom_text(size = 2, position = position_stack(vjust = 0.5), family = "Times") +
  labs(y = "Proportion of Protein Class", x = "Protein Class", fill = "Concordance Category")

#######################################################
#Saving HPA Figures
#######################################################

###FIGURE S2

figs2a <- HPA_prop_plot +
  theme(legend.position = "none") +
  coord_flip()

figs2b <- HPA_prop_plot_stacked +
  theme(axis.title.y = element_blank()) +
  coord_flip()

figs2 <- figs2a + figs2b

figs2 <- figs2 & theme(plot.margin = unit(c(2, 2, 2, 2), "pt"))

ggsave(file = file.path(docs_data, "Human_Protein_Atlas/FigS2.eps"), plot = figs2, width = 7.5, height = 5, 
       units = "in", dpi = 600, family = "Times", device = postscript, paper = "special", horizontal = FALSE)

###FIGURE S3

figs3 <- HPA_prop_plot_discon +
  coord_flip()

ggsave(file = file.path(docs_data, "Human_Protein_Atlas/FigS3.eps"), plot = figs3, width = 5, height = 5, 
       units = "in", dpi = 600, family = "Times", device = postscript, paper = "special", horizontal = FALSE)
