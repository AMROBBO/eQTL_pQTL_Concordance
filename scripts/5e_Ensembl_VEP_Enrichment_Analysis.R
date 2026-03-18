#######################################################
# 5e. Enrichment Analysis of Colocalised QTL Pairs
# Datasets:
#   e. Ensembl VEP: Variant effect annotation pipeline (inputs and outputs)

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
# Extracting Concordant and Discordant SNPs
#######################################################

QTLs_VEP <- QTLs %>%
  dplyr::select("SNP", "effect_allele.exposure", "other_allele.exposure", "id.outcome", "gene.outcome", "type")

#Saving list of SNPs to extract proxies using plink - Concordant and discordant only

QTLs_VEP_con <- QTLs_VEP[which(QTLs_VEP$type == "concordant"),]
QTLs_VEP_dis <- QTLs_VEP[which(QTLs_VEP$type == "discordant"),]

for (i in 1:nrow(QTLs_VEP_con)){
  QTLs_VEP_con$rsID[i] <- unlist(strsplit(QTLs_VEP_con$SNP[i], split=":"))[2]
}

for (i in 1:nrow(QTLs_VEP_dis)){
  QTLs_VEP_dis$rsID[i] <- unlist(strsplit(QTLs_VEP_dis$SNP[i], split=":"))[2]
}

QTLs_VEP_con_SNPs <- QTLs_VEP_con$rsID
QTLs_VEP_dis_SNPs <- QTLs_VEP_dis$rsID

#Save SNPs

write.table(QTLs_VEP_con_SNPs, 
            file = file.path(interim_data, "SNPs_con.txt"), 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

write.table(QTLs_VEP_dis_SNPs, 
            file = file.path(interim_data, "SNPs_dis.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

#SNP lists were read into Ensembl VEP
##############################################################################################################
# Will really need to check this as method has changed

#######################################################
#Reading back in the VEP SNPs
#######################################################

VEP_con <- fread(file.path(interim_data, "VEP_con_res.txt"))
VEP_dis <- fread(file.path(interim_data, "VEP_dis_res.txt")

#######################################################
#Labeling QTLs with Predicted Variant Effects
#######################################################

#Linking proxies with VEPs back to original SNPs

VEP_con <- VEP_con %>%
  select("#Uploaded_variation", "Allele", "Consequence", "SYMBOL")
VEP_dis <- VEP_dis %>%
  select("#Uploaded_variation", "Allele", "Consequence", "SYMBOL")

#Merge QTLs and VEP by rsID names, filtering for matching Gene names and effect alleles

VEP_con_QTLs <- merge(QTLs_VEP_con, VEP_con, by = "SNP_A") %>% 
  filter(id.outcome == SYMBOL) %>%
  filter(effect_allele.exposure == Allele)

VEP_dis_QTLs <- merge(QTLs_VEP_dis, VEP_dis, by = "SNP_A") %>% 
  filter(id.outcome == SYMBOL) %>%
  filter(effect_allele.exposure == Allele)

#Assigning VEPs to QTLs

VEP_con_QTLs$Consequence <- strsplit(VEP_con_QTLs$Consequence, split = ",")
VEP_dis_QTLs$Consequence <- strsplit(VEP_dis_QTLs$Consequence, split = ",")

con_QTLs <- VEP_con_QTLs %>%
  mutate(row = row_number()) %>%  # Add a row identifier
  unnest('Consequence') %>%     # Unnest the list column
  mutate(present = 1) %>%         # Mark presence
  pivot_wider(
    names_from = 'Consequence',
    values_from = present,
    values_fill = 0               # Fill missing combinations with 0
  ) %>%
  dplyr::select(-row)                    # Drop helper column

dis_QTLs <- VEP_dis_QTLs %>%
  mutate(row = row_number()) %>%  # Add a row identifier
  unnest('Consequence') %>%     # Unnest the list column
  mutate(present = 1) %>%         # Mark presence
  pivot_wider(
    names_from = 'Consequence',
    values_from = present,
    values_fill = 0               # Fill missing combinations with 0
  ) %>%
  dplyr::select(-row)                    # Drop helper column

#Save

fwrite(con_QTLs, file.path(results_data, "Ensembl_VEP/VEP_labelled_QTLs_con.csv"))
fwrite(dis_QTLs, file.path(results_data, "Ensembl_VEP/VEP_labelled_QTLs_dis.csv")

#######################################################
#Creating VEP Summary Table
#######################################################

#Sum number of effects in each Concordance Group

VEP_summary <- colSums(con_QTLs[,11:ncol(con_QTLs)])
VEP_summary_2 <- colSums(dis_QTLs[,11:ncol(dis_QTLs)])

#Combine

VEP_summary <- merge(VEP_summary, VEP_summary_2, by = 0, all = T)

VEP_summary[is.na(VEP_summary)] <- 0
colnames(VEP_summary) <- c("Variant Effect", "Concordant", "Discordant")

#Save

VEP_summary_tab <- VEP_summary%>%
  gt() %>%
  tab_header(title = md("Variant Effect Prediction Analysis"),
             subtitle = md("Ensembl VEP"))


gtsave(VEP_summary_tab, file.path(docs_data, "Ensembl_VEP/Ensembl_summary.png")

#######################################################
#Plotting Variant Effects
#######################################################

#Total number of QTLs in each Concordance group

total_con <- length(unique(VEP_con_QTLs$SNP.x))
total_dis <- length(unique(VEP_dis_QTLs$SNP.x))

#Plotting VEP as a proportion of each concordance category:
#Proportion = N QTLs in Variant Effect / total number of QTLs in concordance group

VEP_summary$Concordant_prop <- (VEP_summary$Concordant/total_con) * 100
VEP_summary$Discordant_prop <- (VEP_summary$Discordant/total_dis) * 100

#Rearrange to Plot

variant_plot <- rep(VEP_summary$`Variant Effect`, 2)
type_plot <- c(rep("Concordant", 18), rep("Discordant", 18))
prop_plot <- c(VEP_summary$Concordant_prop, VEP_summary$Discordant_prop)

VEP_plot <- data.frame(variant_plot, type_plot, prop_plot)

#VEP Plot 1

VEP_prop_plot <- ggplot(VEP_plot, aes(fill=type_plot, y=prop_plot, x=variant_plot)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(text = element_text(size=10, family = "Times"),
        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5)) +
  labs(y = "Proportion of Concordance Category", x = "Predicted Variant Effect", fill = "Concordance Category")

#Plotting concordance categories as a proportion of each variant effect
#Proportion = N QTLs with variant effect / total number of QTLs with variant effect

for (i in 1:nrow(VEP_summary)){
  variant_sum <- sum(VEP_summary[i,2:3]) #Total QTLs with variant effect
  VEP_summary$Concordant_prop[i] <- (VEP_summary$Concordant[i]/variant_sum) * 100
  VEP_summary$Discordant_prop[i] <- (VEP_summary$Discordant[i]/variant_sum) * 100
}

#Rearrange for plotting

variant_plot <- rep(VEP_summary$`Variant Effect`, 2)
type_plot <- c(rep("Concordant", 18), rep("Discordant", 18))
prop_plot <- c(VEP_summary$Concordant_prop, VEP_summary$Discordant_prop)
num_plot <- c(VEP_summary$Concordant, VEP_summary$Discordant)

VEP_plot <- data.frame(variant_plot, type_plot, prop_plot, num_plot)

#VEP plot 2 - Stacked

VEP_prop_plot_stacked <- ggplot(VEP_plot, aes(fill=type_plot, y=prop_plot, x=variant_plot, label=num_plot)) + 
  geom_bar(position="stack", stat="identity") +
  theme(text = element_text(size=10, family = "Times"),
        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5, family = "Times")) +
  geom_text(size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE, angle = 90, family = "Times") +
  labs(y = "Proportion of Predicted Variant Effect", x = "Predicted Variant Effect", fill = "Concordance Category")

#######################################################
#Saving Ensembl VEP Figures
#######################################################

###FIGURE 5

fig5a <- VEP_prop_plot +
  theme(legend.position = "none")

fig5b <- VEP_prop_plot_stacked

fig5 <- fig5a + fig5b

fig5 <- fig5 & theme(plot.margin = unit(c(2, 2, 2, 2), "pt"))

ggsave(file = file.path(docs_data, "Ensembl_VEP/Fig5.eps"), plot = fig5, width = 5.2, height = 5, 
       units = "in", dpi = 600, family = "Times", device = postscript, paper = "special", horizontal = FALSE)
