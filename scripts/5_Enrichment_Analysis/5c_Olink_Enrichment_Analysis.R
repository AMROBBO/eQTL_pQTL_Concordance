#######################################################
# 5c. Enrichment Analysis of Colocalised QTL Pairs
# Datasets:
#   c. Olink: Protein Class annotation
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
#Read in Olink Data
#######################################################

olink_files <- list.files(file.path(raw_data, "olink"), pattern = "Olink", full.names = T)

#######################################################
#Labeling QTL data for Protein Function
#######################################################

QTLs_olink <- QTLs %>%
  dplyr::select("rsid", "id", "type")

#Reading in each olink panel and assigning function to those SNPs present in panel

QTLs_protein_function <- data.frame()

for (f in olink_files){
  
  tmp <- fread(f) #Read in SNPs in each file
  
  class <- unlist(strsplit(f, split = " ")) #Extract protein function from file name
  class <- class[4:(length(class) - 2)]
  
  QTL_class <- QTLs_olink[QTLs_olink$id %in% tmp$Gene,] #Assigning protein class to QTLs present in panel
  QTL_class$Class <- paste(class, collapse = "_")
  
  QTLs_protein_function <- rbind(QTLs_protein_function, QTL_class) #Combine to one dataset
}

#Combining subcategories into one

QTLs_protein_function[QTLs_protein_function$Class == "Cardiovascular_II" | QTLs_protein_function$Class == "Cardiovascular_III",]$Class <- "Cardiovascular"
QTLs_protein_function[QTLs_protein_function$Class == "Oncology_II" | QTLs_protein_function$Class == "Oncology_III",]$Class <- "Oncology"

fwrite(QTLs_protein_function, file.path(processed_data, "Olink/Olink_labelled_QTLs.csv"))

#######################################################
#Creating Olink Summary Table
#######################################################

#Grouping by concordance type and summing number of protein functions present

QTLs_protein_function_con <- QTLs_protein_function[QTLs_protein_function$type == "concordant",]
QTLs_protein_function_con <- QTLs_protein_function_con[, .N, by = 'Class']
colnames(QTLs_protein_function_con)[2] <- "Concordant"

QTLs_protein_function_dis <- QTLs_protein_function[which(QTLs_protein_function$type == "discordant"),]
QTLs_protein_function_dis <- QTLs_protein_function_dis[, .N, by = 'Class']
colnames(QTLs_protein_function_dis)[2] <- "Discordant"

QTLs_protein_function_eQTL_dropped <- QTLs_protein_function[which(QTLs_protein_function$type == "eQTL_dropped"),]
QTLs_protein_function_eQTL_dropped <- QTLs_protein_function_eQTL_dropped[, .N, by = 'Class']
colnames(QTLs_protein_function_eQTL_dropped)[2] <- "eQTL_dropped"

QTLs_protein_function_pQTL_dropped <- QTLs_protein_function[which(QTLs_protein_function$type == "pQTL_dropped"),]
QTLs_protein_function_pQTL_dropped <- QTLs_protein_function_pQTL_dropped[, .N, by = 'Class']
colnames(QTLs_protein_function_pQTL_dropped)[2] <- "pQTL_dropped"

#Combining into one table

olink_summary <- merge(QTLs_protein_function_con, QTLs_protein_function_dis, by = 'Class', all = T)
olink_summary <- merge(olink_summary, QTLs_protein_function_eQTL_dropped, by = 'Class', all = T)
olink_summary <- merge(olink_summary, QTLs_protein_function_pQTL_dropped, by = 'Class', all = T)
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

total_con <- QTLs_protein_function %>% 
  filter(type == "concordant") %>% 
  select(rsid) %>% 
  unique() %>% 
  nrow()
total_dis <- QTLs_protein_function %>% 
  filter(type == "discordant") %>% 
  select(rsid) %>% 
  unique() %>% 
  nrow()
total_eQTL_dropped <- QTLs_protein_function %>% 
  filter(type == "eQTL_dropped") %>% 
  select(rsid) %>% 
  unique() %>% 
  nrow()
total_pQTL_dropped <- QTLs_protein_function %>% 
  filter(type == "pQTL_dropped") %>% 
  select(rsid) %>% 
  unique() %>% 
  nrow()

#Plotting Protein Class involvement as a proportion of each concordance category:
#Proportion = N QTLs in protein class / total number of QTLs in concordance group

olink_plot <- olink_summary %>%
  mutate(
    Concordant = Concordant / total_con * 100,
    Discordant = Discordant / total_dis * 100,
    eQTL_dropped = eQTL_dropped / total_eQTL_dropped * 100,
    pQTL_dropped = pQTL_dropped / total_pQTL_dropped * 100
  )

#Rearrange to plot

olink_plot <- olink_plot %>%
  pivot_longer(cols = c(Concordant, Discordant, eQTL_dropped, pQTL_dropped),
               names_to = "type_plot", values_to = "prop_plot")

#Olink Plot 1a - all Concordance groups

olink_prop_plot <- ggplot(olink_plot, aes(fill=type_plot, y=prop_plot, x=Class)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(text = element_text(size=10, family = "Times"),
        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5)) +
  labs(y = "Proportion of Concordance Category", x = "Protein Class", fill = "Concordance Category")

#Olink Plot 1b - Con/Dis Groups only

olink_plot_discon <- olink_plot[which(olink_plot$type_plot == "Concordant" | olink_plot$type_plot == "Discordant"),]

olink_prop_plot_discon <- ggplot(olink_plot_discon, aes(fill=type_plot, y=prop_plot, x=Class)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(text = element_text(size=10, family = "Times"),
        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5)) +
  labs(y = "Proportion of Concordance Category", x = "Protein Class", fill = "Concordance Category")

#Plotting concordance categories as a proportion of each protein class
#Proportion = N QTLs in protein class / total number of QTLs in protein class

olink_plot_stacked <- olink_summary

for (i in 1:nrow(olink_plot_stacked)){
  olink_sum <- sum(olink_plot_stacked[i,2:5])
  olink_plot_stacked$Concordant[i] <- (olink_plot_stacked$Concordant[i]/olink_sum) * 100
  olink_plot_stacked$Discordant[i] <- (olink_plot_stacked$Discordant[i]/olink_sum) * 100
  olink_plot_stacked$eQTL_dropped[i] <- (olink_plot_stacked$eQTL_dropped[i]/olink_sum) * 100
  olink_plot_stacked$pQTL_dropped[i] <- (olink_plot_stacked$pQTL_dropped[i]/olink_sum) * 100
}

#Rearrange to plot

olink_num <- olink_summary %>% 
  pivot_longer(cols = c(Concordant, Discordant, eQTL_dropped, pQTL_dropped),
               names_to = "type_plot", values_to = "num_plot")


olink_plot_stacked <- olink_plot_stacked %>%
  pivot_longer(cols = c(Concordant, Discordant, eQTL_dropped, pQTL_dropped),
               names_to = "type_plot", values_to = "prop_plot") %>% 
  mutate(num_plot = olink_num$num_plot)


#Olink Plot 2 - Stacked plot

olink_prop_plot_stacked <- ggplot(olink_plot_stacked, aes(fill=type_plot, y=prop_plot, x=Class, label=num_plot)) + 
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
