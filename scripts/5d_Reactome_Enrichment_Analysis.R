#######################################################
# 5d. Enrichment Analysis of Colocalised QTL Pairs
# Datasets:
#   d. Reactome: Pathway annotation

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
results_data <- Sys.getenv("resultsdir")
docs_data <- Sys.getenv("docsdir")

#######################################################
# Read in colocalised QTL pairs 
#######################################################

QTLs <- fread(file.path(results_data, "pQTL_concordance.csv"))

#######################################################
#Read in Reactome Data
#######################################################

reactome <- fread(file.path(raw_data, "Ensembl2Reactome.tsv", header = FALSE))

#######################################################
#Labeling QTL data for pathways
#######################################################

QTLs_reactome <- QTLs %>%
  dplyr::select("SNP", "id.outcome", "gene.outcome", "type")

#Formatting pathway data - Homo sapiens only

header <- c("Source database identifier", "Reactome Pathway Stable identifier", 
            "URL", "Event (Pathway or Reaction) Name", "Evidence Code", "Species")
colnames(reactome) <- header

reactome_hs <- reactome[which(reactome$Species == "Homo sapiens"),] %>%
  dplyr::select("Source database identifier", "Event (Pathway or Reaction) Name")

#Assigning pathways to genes based on Ensembl ID

colnames(QTLs_reactome)[3] <- "Ensembl"
colnames(reactome_hs)[1] <- "Ensembl"

QTLs_reactome <- merge(QTLs_reactome, reactome_hs, by = "Ensembl")
QTLs_reactome <- unique(QTLs_reactome)

fwrite(QTLs_reactome, file.path(results_data, "Reactome/Reactome_labelled_QTLs.csv"))

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

reactome_summary <- merge(QTLs_reactome_con, QTLs_reactome_dis, by = 'Event (Pathway or Reaction) Name')
reactome_summary <- merge(reactome_summary, QTLs_reactome_eQTL_dropped, by = 'Event (Pathway or Reaction) Name')
reactome_summary <- merge(reactome_summary, QTLs_reactome_pQTL_dropped, by = 'Event (Pathway or Reaction) Name')

#Save

fwrite(reactome_summary, file.path(results_data, "Reactome/Reactome_summary.csv"))

reactome_summary_tab <- reactome_summary%>%
  gt() %>%
  tab_header(title = md("Pathway Analysis"),
             subtitle = md("Reactome"))


gtsave(reactome_summary_tab, file.path(docs_data, "Reactome/Reactome_summary.png"))

#######################################################
#Plotting Pathways
#######################################################

#Number of QTLs in each concordance type

total_con <- nrow(unique(QTLs_reactome[which(QTLs_reactome$type == "concordant"), 1]))
total_dis <- nrow(unique(QTLs_reactome[which(QTLs_reactome$type == "discordant"), 1]))
total_eQTL_dropped <- nrow(unique(QTLs_reactome[which(QTLs_reactome$type == "eQTL_dropped"), 1]))
total_pQTL_dropped <- nrow(unique(QTLs_reactome[which(QTLs_reactome$type == "pQTL_dropped"), 1]))

#Plotting Pathway involvement as a proportion of each concordance category:
#Proportion = N QTLs in pathway / total number of QTLs in concordance group

reactome_summary$Concordant_prop <- (reactome_summary$Concordant/total_con) * 100
reactome_summary$Discordant_prop <- (reactome_summary$Discordant/total_dis) * 100
reactome_summary$eQTL_dropped_prop <- (reactome_summary$eQTL_dropped/total_eQTL_dropped) * 100
reactome_summary$pQTL_dropped_prop <- (reactome_summary$pQTL_dropped/total_pQTL_dropped) * 100

#Shortening some of the names for plotting
reactome_summary$`Event (Pathway or Reaction) Name`[105] <- "Regulation of Insulin-like Growth Factor (IGF) **"
reactome_summary$`Event (Pathway or Reaction) Name`[49] <- "Immunoregulatory interactions *"

#Rearrange to plot

pathway_plot <- rep(reactome_summary$`Event (Pathway or Reaction) Name`, 4)
type_plot <- c(rep("Concordant", 138), rep("Discordant", 138), rep("eQTL_dropped", 138), rep("pQTL_dropped", 138))
prop_plot <- c(reactome_summary$Concordant_prop, reactome_summary$Discordant_prop, reactome_summary$eQTL_dropped_prop, reactome_summary$pQTL_dropped_prop)

reactome_plot <- data.frame(pathway_plot, type_plot, prop_plot)

#Reactome Plot 1a - all Concordance groups

reactome_prop_plot <- ggplot(reactome_plot, aes(fill=type_plot, y=prop_plot, x=pathway_plot)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(text = element_text(size=10, family = "Times"),
        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5, family = "Times")) +
  labs(y = "Proportion of Concordance Category", x = "Pathway", fill = "Concordance Category")

#Reactome Plot 1b - For pathways with over 2% of at least one con group in it

reactome_summary_sig <- reactome_summary %>% 
  filter(Concordant_prop > 2 | Discordant_prop > 2 | eQTL_dropped_prop > 2 | pQTL_dropped_prop > 2)

#Rearrange to plot

pathway_plot <- rep(reactome_summary_sig$`Event (Pathway or Reaction) Name`, 4)
type_plot <- c(rep("Concordant", 25), rep("Discordant", 25), rep("eQTL_dropped", 25), rep("pQTL_dropped", 25))
prop_plot <- c(reactome_summary_sig$Concordant_prop, reactome_summary_sig$Discordant_prop, reactome_summary_sig$eQTL_dropped_prop, reactome_summary_sig$pQTL_dropped_prop)

reactome_plot_sig <- data.frame(pathway_plot, type_plot, prop_plot)

#Reactome Plot 1b - all Concordance groups, significant pathways

reactome_prop_plot_sig <- ggplot(reactome_plot_sig, aes(fill=type_plot, y=prop_plot, x=pathway_plot)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(text = element_text(size=10, family = "Times"),
        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5)) +
  labs(y = "Proportion of Concordance Category", x = "Pathway", fill = "Concordance Category")

#Save pathways

fwrite(reactome_summary_sig, file.path(results_data, "Reactome/sig_pathways.csv"))

#Reactome Plot 1c - Con/Dis groups only

reactome_plot_discon <- reactome_plot[which(reactome_plot$type_plot == "Concordant" | reactome_plot$type_plot == "Discordant"),]

reactome_prop_plot_discon <- ggplot(reactome_plot_discon, aes(fill=type_plot, y=prop_plot, x=pathway_plot)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(text = element_text(size=10, family = "Times"),
        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5)) +
  labs(y = "Proportion of Concordance Category", x = "Pathway", fill = "Concordance Category")

#Reactome Plot 1d - Pathways with multiple QTLs involved, Discon only

reactome_plot_discon_sig <- reactome_plot_sig[which(reactome_plot_sig$type_plot == "Concordant" | reactome_plot_sig$type_plot == "Discordant"),]

reactome_prop_plot_discon_sig <- ggplot(reactome_plot_discon_sig, aes(fill=type_plot, y=prop_plot, x=pathway_plot)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(text = element_text(size=10, family = "Times"),
        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5)) +
  labs(y = "Proportion of Concordance Category", x = "Pathway", fill = "Concordance Category")

#Plotting concordance categories as a proportion of each pathway
#Proportion = N QTLs in pathway / total number of QTLs in pathway

for (i in 1:nrow(reactome_summary)){
  pathway_sum <- sum(reactome_summary[i,2:5]) #Number of QTLs in pathway
  reactome_summary$Concordant_prop[i] <- (reactome_summary$Concordant[i]/pathway_sum) * 100
  reactome_summary$Discordant_prop[i] <- (reactome_summary$Discordant[i]/pathway_sum) * 100
  reactome_summary$eQTL_dropped_prop[i] <- (reactome_summary$eQTL_dropped[i]/pathway_sum) * 100
  reactome_summary$pQTL_dropped_prop[i] <- (reactome_summary$pQTL_dropped[i]/pathway_sum) * 100
}

#Rearrange to plot

pathway_plot <- rep(reactome_summary$`Event (Pathway or Reaction) Name`, 4)
type_plot <- c(rep("Concordant", 138), rep("Discordant", 138), rep("eQTL_dropped", 138), rep("pQTL_dropped", 138))
prop_plot <- c(reactome_summary$Concordant_prop, reactome_summary$Discordant_prop, reactome_summary$eQTL_dropped_prop, reactome_summary$pQTL_dropped_prop)
num_plot <- c(reactome_summary$Concordant, reactome_summary$Discordant, reactome_summary$eQTL_dropped, reactome_summary$pQTL_dropped)

reactome_plot <- data.frame(pathway_plot, type_plot, prop_plot, num_plot)

#Reactome Plot 2 - Stacked Plot

reactome_prop_plot_stacked <- ggplot(reactome_plot, aes(fill=type_plot, y=prop_plot, x=pathway_plot, label=num_plot)) + 
  geom_bar(position="stack", stat="identity") +
  theme(text = element_text(size=10, family = "Times"),
        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5)) +
  geom_text(size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE, angle = 90, family = "Times") +
  labs(y = "Proportion of Pathway", x = "Pathway", fill = "Concordance Category")

#Reactome Plot 2b - Pathways with multiple QTLs involved

reactome_plot_stacked_sig <- reactome_plot[which(reactome_plot$pathway_plot %in% reactome_plot_sig$pathway_plot),]

reactome_plot_stacked_sig <- ggplot(reactome_plot_stacked_sig, aes(fill=type_plot, y=prop_plot, x=pathway_plot, label=num_plot)) + 
  geom_bar(position="stack", stat="identity") +
  theme(text = element_text(size=10, family = "Times"),
        axis.text.x = element_text(size=8, angle = 90, hjust = 1, vjust = 0.5)) +
  geom_text(size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE, angle = 90, family = "Times") +
  labs(y = "Proportion of Pathway", x = "Pathway", fill = "Concordance Category")

#######################################################
#Saving Reactome Figures
#######################################################

###FIGURE 4

fig4a <- reactome_prop_plot_sig +
  theme(legend.position = "none")

fig4b <- reactome_plot_stacked_sig

fig4 <- fig4a + fig4b

fig4 <- fig4 & theme(plot.margin = unit(c(2, 2, 2, 2), "pt"))

ggsave(file = file.path(docs_data, "Reactome/Fig4.eps", plot = fig4, width = 7.5, height = 5.5, 
       units = "in", dpi = 600, family = "Times", device = postscript, paper = "special", horizontal = FALSE)

###FIGURE S5

figs5a <- reactome_prop_plot +
  theme(legend.position = "none")

figs5b <- reactome_prop_plot_stacked

figs5 <- figs5a + figs5b

ggsave(file = file.path(docs_data, "Reactome/FigS5.eps"), plot = figs5, width = 27, height = 8.75, 
       units = "in", dpi = 600, family = "Times", device = postscript, paper = "special", horizontal = FALSE)

###FIGURE S6

figs6a <- reactome_prop_plot_discon +
  theme(legend.position = "none")

figs6b <- reactome_prop_plot_discon_sig

figs6 <- figs6a + figs6b

ggsave(file = file.path(docs_data, "Reactome/FigS6.eps"), plot = figs6, width = 27, height = 8.75, 
       units = "in", dpi = 600, family = "Times", device = postscript, paper = "special", horizontal = FALSE)


