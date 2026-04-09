#######################################################
# 5c. Enrichment Analysis of Colocalised QTL Pairs
# Datasets:
#   c. Olink: Protein Class annotation

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
#Read in Olink Data
#######################################################

olink_files <- list.files(file.path(raw_data, "olink"), pattern = "Olink", full.names = T)

#######################################################
#Labeling QTL data for Protein Function
#######################################################

QTLs_olink <- QTLs %>%
  dplyr::select("SNP", "id.outcome", "gene.outcome", "type")

#Reading in each olink panel and assigning function to those SNPs present in panel

QTLs_protein_function <- data.frame()

for (f in olink_files){
  
  tmp <- fread(f) #Read in SNPs in each file
  
  class <- unlist(strsplit(f, split = " ")) #Extract protein function from file name
  class <- class[4:(length(class) - 2)]
  
  QTL_class <- QTLs_olink[QTLs_olink$id.outcome %in% tmp$Gene,] #Assigning protein class to QTLs present in panel
  QTL_class$class <- paste(class, collapse = "_")
  
  QTLs_protein_function <- rbind(QTLs_protein_function, QTL_class) #Combine to one dataset
}

#Combining subcategories into one

QTLs_protein_function[QTLs_protein_function$class == "Cardiovascular_II" | QTLs_protein_function$class == "Cardiovascular_III",]$class <- "Cardiovascular"
QTLs_protein_function[QTLs_protein_function$class == "Oncology_II" | QTLs_protein_function$class == "Oncology_III",]$class <- "Oncology"

fwrite(QTLs_protein_function, file.path(results_data, "Olink/Olink_labelled_QTLs.csv"))

#######################################################
#Creating Olink Summary Table
#######################################################

#Grouping by concordance type and summing number of protein functions present

QTLs_protein_function_con <- QTLs_protein_function[QTLs_protein_function$type == "concordant",]
QTLs_protein_function_con <- QTLs_protein_function_con[, .N, by = 'class']
colnames(QTLs_protein_function_con)[2] <- "Concordant"

QTLs_protein_function_dis <- QTLs_protein_function[which(QTLs_protein_function$type == "discordant"),]
QTLs_protein_function_dis <- QTLs_protein_function_dis[, .N, by = 'class']
colnames(QTLs_protein_function_dis)[2] <- "Discordant"

QTLs_protein_function_eQTL_dropped <- QTLs_protein_function[which(QTLs_protein_function$type == "eQTL_dropped"),]
QTLs_protein_function_eQTL_dropped <- QTLs_protein_function_eQTL_dropped[, .N, by = 'class']
colnames(QTLs_protein_function_eQTL_dropped)[2] <- "eQTL_dropped"

QTLs_protein_function_pQTL_dropped <- QTLs_protein_function[which(QTLs_protein_function$type == "pQTL_dropped"),]
QTLs_protein_function_pQTL_dropped <- QTLs_protein_function_pQTL_dropped[, .N, by = 'class']
colnames(QTLs_protein_function_pQTL_dropped)[2] <- "pQTL_dropped"

#Combining into one table

olink_summary <- merge(QTLs_protein_function_con, QTLs_protein_function_dis, by = 'class')
olink_summary <- merge(olink_summary, QTLs_protein_function_eQTL_dropped, by = 'class')
olink_summary <- merge(olink_summary, QTLs_protein_function_pQTL_dropped, by = 'class', all = T)
olink_summary[is.na(olink_summary)] <- 0 #Removing NAs

#Save

olink_summary_tab <- olink_summary%>%
  gt() %>%
  tab_header(title = md("Protein Class Analysis"),
             subtitle = md("Olink"))


gtsave(olink_summary_tab, file.path(docs_data, "Olink/Olink_summary.png"))

#######################################################
#Plotting Protein Functions
#######################################################

#Total number of QTLs in each concordance category

total_con <- nrow(unique(QTLs_protein_function[which(QTLs_protein_function$type == "concordant"), 1]))
total_dis <- nrow(unique(QTLs_protein_function[which(QTLs_protein_function$type == "discordant"), 1]))
total_eQTL_dropped <- nrow(unique(QTLs_protein_function[which(QTLs_protein_function$type == "eQTL_dropped"), 1]))
total_pQTL_dropped <- nrow(unique(QTLs_protein_function[which(QTLs_protein_function$type == "pQTL_dropped"), 1]))

#Plotting Protein Class involvement as a proportion of each concordance category:
#Proportion = N QTLs in protein class / total number of QTLs in concordance group

olink_summary$Concordant_prop <- (olink_summary$Concordant/total_con) * 100
olink_summary$Discordant_prop <- (olink_summary$Discordant/total_dis) * 100
olink_summary$eQTL_dropped_prop <- (olink_summary$eQTL_dropped/total_eQTL_dropped) * 100
olink_summary$pQTL_dropped_prop <- (olink_summary$pQTL_dropped/total_pQTL_dropped) * 100

#Rearrange to plot

olink_plot <- rep(olink_summary$class, 4)
type_plot <- c(rep("Concordant", 12), rep("Discordant", 12), rep("eQTL_dropped", 12), rep("pQTL_dropped", 12))
prop_plot <- c(olink_summary$Concordant_prop, olink_summary$Discordant_prop, olink_summary$eQTL_dropped_prop, olink_summary$pQTL_dropped_prop)

olink_plot <- data.frame(olink_plot, type_plot, prop_plot)

#Olink Plot 1a - all Concordance groups

olink_prop_plot <- ggplot(olink_plot, aes(fill=type_plot, y=prop_plot, x=olink_plot)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(text = element_text(size=10, family = "Times"),
        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5)) +
  labs(y = "Proportion of Concordance Category", x = "Protein Class", fill = "Concordance Category")

#Olink Plot 1b - Con/Dis Groups only

olink_plot_discon <- olink_plot[which(olink_plot$type_plot == "Concordant" | olink_plot$type_plot == "Discordant"),]

olink_prop_plot_discon <- ggplot(olink_plot_discon, aes(fill=type_plot, y=prop_plot, x=olink_plot)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(text = element_text(size=10, family = "Times"),
        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5)) +
  labs(y = "Proportion of Concordance Category", x = "Protein Class", fill = "Concordance Category")


#Plotting concordance categories as a proportion of each protein class
#Proportion = N QTLs in protein class / total number of QTLs in protein class

for (i in 1:nrow(olink_summary)){
  olink_sum <- sum(olink_summary[i,2:5])
  olink_summary$Concordant_prop[i] <- (olink_summary$Concordant[i]/olink_sum) * 100
  olink_summary$Discordant_prop[i] <- (olink_summary$Discordant[i]/olink_sum) * 100
  olink_summary$eQTL_dropped_prop[i] <- (olink_summary$eQTL_dropped[i]/olink_sum) * 100
  olink_summary$pQTL_dropped_prop[i] <- (olink_summary$pQTL_dropped[i]/olink_sum) * 100
}

#Rearrange to plot

olink_plot <- rep(olink_summary$class, 4)
type_plot <- c(rep("Concordant", 12), rep("Discordant", 12), rep("eQTL_dropped", 12), rep("pQTL_dropped", 12))
prop_plot <- c(olink_summary$Concordant_prop, olink_summary$Discordant_prop, olink_summary$eQTL_dropped_prop, olink_summary$pQTL_dropped_prop)
num_plot <- c(olink_summary$Concordant, olink_summary$Discordant, olink_summary$eQTL_dropped, olink_summary$pQTL_dropped)

olink_plot <- data.frame(olink_plot, type_plot, prop_plot, num_plot)

#Olink Plot 2 - Stacked plot

olink_prop_plot_stacked <- ggplot(olink_plot, aes(fill=type_plot, y=prop_plot, x=olink_plot, label=num_plot)) + 
  geom_bar(position="stack", stat="identity") +
  theme(text = element_text(size=10, family = "Times"),
        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5)) +
  geom_text(size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE, angle = 90, family = "Times") +
  labs(y = "Proportion of Protein Class", x = "Protein Class", fill = "Concordance Category")

#######################################################
#Saving Olink Figures
#######################################################

###FIGURE 3

fig3a <- olink_prop_plot +
  theme(legend.position = "none")

fig3b <- olink_prop_plot_stacked

fig3 <- fig3a + fig3b

fig3 <- fig3 & theme(plot.margin = unit(c(2, 2, 2, 2), "pt"))

ggsave(file = file.path(docs_data, "Olink/Fig3.eps"), plot = fig3, width = 5.2, height = 4, 
       units = "in", dpi = 600, family = "Times", device = postscript, paper = "special", horizontal = FALSE)

###FIGURE S4

figs4 <- olink_prop_plot_discon

ggsave(file = file.path(docs_data, "Olink/FigS4.eps"), plot = figs4, width = 4, height = 4, 
       units = "in", dpi = 600, family = "Times", device = postscript, paper = "special", horizontal = FALSE)
