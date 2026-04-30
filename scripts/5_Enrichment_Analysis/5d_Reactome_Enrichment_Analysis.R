#######################################################
# 5d. Enrichment Analysis of Colocalised QTL Pairs
# Datasets:
#   d. Reactome: Pathway annotation
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
#Read in Reactome Data
#######################################################

reactome <- fread(file.path(raw_data, "Ensembl2Reactome.tsv"), header = FALSE)

#######################################################
# Labeling QTL data with Ensembl gene names
#######################################################

colnames(gene_names) <- c("id", "Gene_name")

QTLs <- QTLs %>% 
  left_join(gene_names, by = "id")

#######################################################
#Labeling QTL data for pathways
#######################################################

QTLs_reactome <- QTLs %>%
  dplyr::select("rsid", "id", "Gene_name", "type")

#Formatting pathway data - Homo sapiens only

header <- c("Source database identifier", "Reactome Pathway Stable identifier", 
            "URL", "Event (Pathway or Reaction) Name", "Evidence Code", "Species")
colnames(reactome) <- header

reactome <- reactome %>% 
  filter(Species == "Homo sapiens") %>%
  dplyr::select("Source database identifier", "Event (Pathway or Reaction) Name")

#Assigning pathways to genes based on Ensembl ID

colnames(reactome)[1] <- "Gene_name"

QTLs_reactome <- merge(QTLs_reactome, reactome, by = "Gene_name")
QTLs_reactome <- unique(QTLs_reactome)

fwrite(QTLs_reactome, file.path(processed_data, "Reactome/Reactome_labelled_QTLs.csv"))

#######################################################
#Creating Reactome Summary Table
#######################################################

#Categorising QTLs into concordance groups and summing pathways present in each

QTLs_reactome_con <- QTLs_reactome[which(QTLs_reactome$type == "concordant"),]
QTLs_reactome_con <- QTLs_reactome_con[, .N, by = 'Event (Pathway or Reaction) Name']
colnames(QTLs_reactome_con)[2] <- "Concordant"

QTLs_reactome_dis <- QTLs_reactome[which(QTLs_reactome$type == "discordant"),]
QTLs_reactome_dis <- QTLs_reactome_dis[, .N, by = 'Event (Pathway or Reaction) Name']
colnames(QTLs_reactome_dis)[2] <- "Discordant"

QTLs_reactome_eQTL_dropped <- QTLs_reactome[which(QTLs_reactome$type == "eQTL_dropped"),]
QTLs_reactome_eQTL_dropped <- QTLs_reactome_eQTL_dropped[, .N, by = 'Event (Pathway or Reaction) Name']
colnames(QTLs_reactome_eQTL_dropped)[2] <- "eQTL_dropped"

QTLs_reactome_pQTL_dropped <- QTLs_reactome[which(QTLs_reactome$type == "pQTL_dropped"),]
QTLs_reactome_pQTL_dropped <- QTLs_reactome_pQTL_dropped[, .N, by = 'Event (Pathway or Reaction) Name']
colnames(QTLs_reactome_pQTL_dropped)[2] <- "pQTL_dropped"

#Join into one table

reactome_summary <- merge(QTLs_reactome_con, QTLs_reactome_dis, by = 'Event (Pathway or Reaction) Name', all = T)
reactome_summary <- merge(reactome_summary, QTLs_reactome_eQTL_dropped, by = 'Event (Pathway or Reaction) Name', all = T)
reactome_summary <- merge(reactome_summary, QTLs_reactome_pQTL_dropped, by = 'Event (Pathway or Reaction) Name', all = T)

reactome_summary[is.na(reactome_summary)] <- 0 #Removing NAs

#Save

fwrite(reactome_summary, file.path(processed_data, "Reactome/Reactome_summary.csv"))

#######################################################
# Summary table for pathways with >10 Genes involved
#######################################################

reactome_sig <- reactome_summary[rowSums(reactome_summary[,2:5]) > 10,]

reactome_summary_tab <- reactome_sig %>%
  gt() %>%
  tab_header(title = md("Pathway Analysis"),
             subtitle = md("Reactome"))


gtsave(reactome_summary_tab, file.path(docs_data, "Reactome/Reactome_summary.png"))

#######################################################
#Plotting Pathways
#######################################################

#Number of QTLs in each concordance type
total_con <- QTLs_reactome %>% 
  filter(type == "concordant") %>% 
  select(rsid) %>% 
  unique() %>% 
  nrow()
total_dis <- QTLs_reactome %>% 
  filter(type == "discordant") %>% 
  select(rsid) %>% 
  unique() %>% 
  nrow()
total_eQTL_dropped <- QTLs_reactome %>% 
  filter(type == "eQTL_dropped") %>% 
  select(rsid) %>% 
  unique() %>% 
  nrow()
total_pQTL_dropped <- QTLs_reactome %>% 
  filter(type == "pQTL_dropped") %>% 
  select(rsid) %>% 
  unique() %>% 
  nrow()

#Plotting Pathway involvement as a proportion of each concordance category:
#Proportion = N QTLs in pathway / total number of QTLs in concordance group

reactome_plot <- reactome_summary %>%
  mutate(
    Concordant = Concordant / total_con * 100,
    Discordant = Discordant / total_dis * 100,
    eQTL_dropped = eQTL_dropped / total_eQTL_dropped * 100,
    pQTL_dropped = pQTL_dropped / total_pQTL_dropped * 100
  )

reactome_plot_sig <- reactome_plot[rowSums(reactome_plot[,2:5]) > 10,]

#Shortening some of the names for plotting
#reactome_summary$`Event (Pathway or Reaction) Name`[105] <- "Regulation of Insulin-like Growth Factor (IGF) **"
#reactome_summary$`Event (Pathway or Reaction) Name`[49] <- "Immunoregulatory interactions *"

#Rearrange to plot

reactome_plot_sig <- reactome_plot_sig %>%
  pivot_longer(cols = c(Concordant, Discordant, eQTL_dropped, pQTL_dropped),
               names_to = "type_plot", values_to = "prop_plot")

#Reactome Plot 1a - all Concordance groups

reactome_prop_plot <- ggplot(reactome_plot_sig, aes(fill=type_plot, y=prop_plot, x=`Event (Pathway or Reaction) Name`)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(text = element_text(size=10, family = "Times"),
        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5, family = "Times")) +
  labs(y = "Proportion of Concordance Category", x = "Pathway", fill = "Concordance Category")

#Reactome Plot 1b - For pathways with over 2% of at least one con group in it

#reactome_summary_sig <- reactome_summary %>% 
#  filter(Concordant_prop > 2 | Discordant_prop > 2 | eQTL_dropped_prop > 2 | pQTL_dropped_prop > 2)

#Rearrange to plot

#pathway_plot <- rep(reactome_summary_sig$`Event (Pathway or Reaction) Name`, 4)
#type_plot <- c(rep("Concordant", 25), rep("Discordant", 25), rep("eQTL_dropped", 25), rep("pQTL_dropped", 25))
#prop_plot <- c(reactome_summary_sig$Concordant_prop, reactome_summary_sig$Discordant_prop, reactome_summary_sig$eQTL_dropped_prop, reactome_summary_sig$pQTL_dropped_prop)

#reactome_plot_sig <- data.frame(pathway_plot, type_plot, prop_plot)

#Reactome Plot 1b - all Concordance groups, significant pathways

#reactome_prop_plot_sig <- ggplot(reactome_plot_sig, aes(fill=type_plot, y=prop_plot, x=pathway_plot)) + 
#  geom_bar(position="dodge", stat="identity") +
#  theme(text = element_text(size=10, family = "Times"),
#        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5)) +
#  labs(y = "Proportion of Concordance Category", x = "Pathway", fill = "Concordance Category")

#Save pathways

#fwrite(reactome_summary_sig, file.path(results_data, "Reactome/sig_pathways.csv"))

#Reactome Plot 1c - Con/Dis groups only

reactome_plot_discon <- reactome_plot %>% 
  select(`Event (Pathway or Reaction) Name`, Concordant, Discordant) %>% 
  filter(Concordant > 0 | Discordant > 0)

reactome_plot_discon_sig <- reactome_plot_discon[rowSums(reactome_plot_discon[,2:3]) > 5,]

#Rearrange to plot

reactome_plot_discon_sig <- reactome_plot_discon_sig %>%
  pivot_longer(cols = c(Concordant, Discordant),
               names_to = "type_plot", values_to = "prop_plot")


reactome_prop_plot_discon <- ggplot(reactome_plot_discon_sig, aes(fill=type_plot, y=prop_plot, x=`Event (Pathway or Reaction) Name`)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(text = element_text(size=10, family = "Times"),
        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5)) +
  labs(y = "Proportion of Concordance Category", x = "Pathway", fill = "Concordance Category")

#Reactome Plot 1d - Pathways with multiple QTLs involved, Discon only

#reactome_plot_discon_sig <- reactome_plot_sig[which(reactome_plot_sig$type_plot == "Concordant" | reactome_plot_sig$type_plot == "Discordant"),]

#reactome_prop_plot_discon_sig <- ggplot(reactome_plot_discon_sig, aes(fill=type_plot, y=prop_plot, x=pathway_plot)) + 
#  geom_bar(position="dodge", stat="identity") +
#  theme(text = element_text(size=10, family = "Times"),
#        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5)) +
#  labs(y = "Proportion of Concordance Category", x = "Pathway", fill = "Concordance Category")

#Plotting concordance categories as a proportion of each pathway
#Proportion = N QTLs in pathway / total number of QTLs in pathway

for (i in 1:nrow(reactome_sig)){
  pathway_sum <- sum(reactome_sig[i,2:5]) #Number of QTLs in pathway
  reactome_sig$Concordant[i] <- (reactome_sig$Concordant[i]/pathway_sum) * 100
  reactome_sig$Discordant[i] <- (reactome_sig$Discordant[i]/pathway_sum) * 100
  reactome_sig$eQTL_dropped[i] <- (reactome_sig$eQTL_dropped[i]/pathway_sum) * 100
  reactome_sig$pQTL_dropped[i] <- (reactome_sig$pQTL_dropped[i]/pathway_sum) * 100
}

#Rearrange to plot

reactome_plot_stacked <- reactome_sig %>%
  pivot_longer(cols = c(Concordant, Discordant, eQTL_dropped, pQTL_dropped),
               names_to = "type_plot", values_to = "prop_plot")

reactome_num <- reactome_summary %>% 
  filter(`Event (Pathway or Reaction) Name` %in% reactome_plot_stacked$`Event (Pathway or Reaction) Name`) %>% 
  pivot_longer(cols = c(Concordant, Discordant, eQTL_dropped, pQTL_dropped),
               names_to = "type_plot", values_to = "num_plot")

reactome_plot_stacked <- reactome_plot_stacked %>% 
  left_join(reactome_num, by = c("Event (Pathway or Reaction) Name", "type_plot"))

#Reactome Plot 2 - Stacked Plot

reactome_prop_plot_stacked <- ggplot(reactome_plot_stacked, aes(fill=type_plot, y=prop_plot, x=`Event (Pathway or Reaction) Name`, label=num_plot)) + 
  geom_bar(position="stack", stat="identity") +
  theme(text = element_text(size=10, family = "Times"),
        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5)) +
  geom_text(size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE, angle = 90, family = "Times") +
  labs(y = "Proportion of Pathway", x = "Pathway", fill = "Concordance Category")

#Reactome Plot 2b - Pathways with multiple QTLs involved

#reactome_plot_stacked_sig <- reactome_plot[which(reactome_plot$pathway_plot %in% reactome_plot_sig$pathway_plot),]

#reactome_plot_stacked_sig <- ggplot(reactome_plot_stacked_sig, aes(fill=type_plot, y=prop_plot, x=pathway_plot, label=num_plot)) + 
#  geom_bar(position="stack", stat="identity") +
#  theme(text = element_text(size=10, family = "Times"),
#        axis.text.x = element_text(size=8, angle = 90, hjust = 1, vjust = 0.5)) +
#  geom_text(size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE, angle = 90, family = "Times") +
#  labs(y = "Proportion of Pathway", x = "Pathway", fill = "Concordance Category")

#######################################################
#Saving Reactome Figures
#######################################################

###FIGURE 4

fig4a <- reactome_prop_plot +
  theme(legend.position = "none")

fig4b <- reactome_prop_plot_stacked

fig4 <- fig4a + fig4b

fig4 <- fig4 & theme(plot.margin = unit(c(2, 2, 2, 2), "pt"))

ggsave(file = file.path(docs_data, "Reactome/Fig4.eps"), plot = fig4, width = 7.5, height = 5.5, 
       units = "in", dpi = 600, family = "Times", device = postscript, paper = "special", horizontal = FALSE)

###FIGURE S5

#figs5a <- reactome_prop_plot +
#  theme(legend.position = "none")

#figs5b <- reactome_prop_plot_stacked

#figs5 <- figs5a + figs5b

#ggsave(file = file.path(docs_data, "Reactome/FigS5.eps"), plot = figs5, width = 27, height = 8.75, 
#       units = "in", dpi = 600, family = "Times", device = postscript, paper = "special", horizontal = FALSE)

###FIGURE S6

figs6a <- reactome_prop_plot_discon #+
#  theme(legend.position = "none")

#figs6b <- reactome_prop_plot_discon_sig

figs6 <- figs6a #+ figs6b

ggsave(file = file.path(docs_data, "Reactome/FigS6.eps"), plot = figs6, width = 27, height = 8.75, 
       units = "in", dpi = 600, family = "Times", device = postscript, paper = "special", horizontal = FALSE)
