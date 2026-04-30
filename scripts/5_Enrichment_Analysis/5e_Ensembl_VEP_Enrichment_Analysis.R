#######################################################
# 5e. Enrichment Analysis of Colocalised QTL Pairs
# Datasets:
#   e. Ensembl VEP: Variant effect annotation pipeline (inputs and outputs)
##
# For VEP enrichment analysis, only variants with H4 = 1 were used 
# ie. same rsid for both pQTL and eQTL
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
library(ieugwasr)            # GWAS summary utilities

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

QTLs_coloc <- fread(file.path(interim_data, "coloc_output.csv")) # Colocalised pQTLs and eQTLs

QTLs <- fread(file.path(processed_data, "Concordance/pQTL_concordance.csv"))

harmonised <- fread(file.path(interim_data, "QTLs_coloc_ready.csv")) # Harmonised pQTLs and eQTLs

#######################################################
# FUNCTIONS
#######################################################
#######################################################
# ld_clump, but catching the errors resulting from null results 
#######################################################

safe_ld_clump <- function(dat) {
  tryCatch(
    {
      ld_clump(dat,
               clump_r2 = 0.001,
               plink_bin = plink,
               bfile = "/Volumes/MRC-IEU-research/projects/ieu3/p3/023/working/data/raw/1000gp/EUR")
    },
    error = function(e) {
      message("No clumped SNPs = dropping outcome")
      return(NULL)
    }
  )
}

#######################################################
# Naive Concordance Test
#######################################################

check_concordance <- function(df, threshold){
  df$type <- ifelse(
    (df$beta_pQTL >= 0 & df$beta_eQTL >= 0) |
      (df$beta_pQTL < 0 & df$beta_eQTL < 0),
    "concordant", "discordant"
  )
  df$type[df$pval_pQTL > threshold] <- "pQTL_dropped"
  df$type[df$pval_eQTL > threshold] <- "eQTL_dropped"
  df
}

#######################################################
# Extracting QTLs pairs with coloc of H4 = 1 (Same rsid)
#######################################################

QTLs_h4 <- QTLs_coloc %>% 
  filter(hit1 == hit2)

#######################################################
# Merging
#######################################################

QTLs_h4 <- tidyr::separate(QTLs_h4, col = hit1, sep = "_", into = c("SNP", "A1", "A2"), remove = F)

names(QTLs_h4)[names(QTLs_h4) == "gene"] <- "exposure"

QTLs_h4_harmonised <- merge(QTLs_h4, harmonised, by = c("SNP", "exposure")) %>% 
  arrange(exposure)

#######################################################
# Clumping
#######################################################

plink <- genetics.binaRies::get_plink_binary() #For LD clumping
QTL_clumped <- data.frame()

QTLS_to_clump <- dplyr::tibble(rsid = QTLs_h4_harmonised$SNP,
                                pval = QTLs_h4_harmonised$pval.exposure, 
                                id = QTLs_h4_harmonised$exposure, 
                                chr = QTLs_h4_harmonised$chr.exposure, 
                                pos = QTLs_h4_harmonised$pos.exposure,
                                effect_allele = QTLs_h4_harmonised$effect_allele.exposure, 
                                other_allele = QTLs_h4_harmonised$other_allele.exposure, 
                                eaf_pQTL = QTLs_h4_harmonised$eaf.exposure, 
                                beta_pQTL = QTLs_h4_harmonised$beta.exposure, 
                                se_pQTL = QTLs_h4_harmonised$se.exposure,
                                pval_pQTL = QTLs_h4_harmonised$pval.exposure,
                                eaf_eQTL = QTLs_h4_harmonised$eaf.outcome,
                                beta_eQTL = QTLs_h4_harmonised$beta.outcome,
                                se_eQTL = QTLs_h4_harmonised$se.outcome,
                                pval_eQTL = QTLs_h4_harmonised$pval.outcome)

for (i in unique(QTLS_to_clump$id)){
  dat <- QTLS_to_clump %>% filter(id == i)
  
  dat_clumped <- safe_ld_clump(dat)
  
  QTL_clumped <- rbind(QTL_clumped, dat_clumped)
}

QTL_clumped <- unique(QTL_clumped)

#######################################################
# Naive Concordance Test - pQTLs
#######################################################

QTL_concordance <- check_concordance(QTL_clumped, threshold = 5e-8)

#######################################################
# Extracting Concordant and Discordant SNPs
#######################################################

QTLs_VEP <- QTL_concordance %>%
  dplyr::select("rsid", "effect_allele", "other_allele", "id", "type")

#Saving list of SNPs for Ensembl VEP

QTLs_VEP_con <- QTLs_VEP %>% 
  filter(type == "concordant")
QTLs_VEP_dis <- QTLs_VEP %>% 
  filter(type == "discordant")
QTLs_VEP_pQTL_dropped <- QTLs_VEP %>% 
  filter(type == "pQTL_dropped")
QTLs_VEP_eQTL_dropped <- QTLs_VEP %>% 
  filter(type == "eQTL_dropped")

#Save SNPs

write.table(QTLs_VEP_con$rsid, 
            file = file.path(interim_data, "Ensembl_VEP/SNPs_con.txt"), 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

write.table(QTLs_VEP_dis$rsid, 
            file = file.path(interim_data, "Ensembl_VEP/SNPs_dis.txt"), 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

write.table(QTLs_VEP_pQTL_dropped$rsid, 
            file = file.path(interim_data, "Ensembl_VEP/SNPs_pQTL_dropped.txt"), 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

write.table(QTLs_VEP_eQTL_dropped$rsid, 
            file = file.path(interim_data, "Ensembl_VEP/SNPs_eQTL_dropped.txt"), 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

#SNP lists were read into Ensembl VEP

#######################################################
#Reading back in the VEP SNPs
#######################################################

VEP_con <- fread(file.path(interim_data, "Ensembl_VEP/VEP_result_con.txt"))
VEP_dis <- fread(file.path(interim_data, "Ensembl_VEP/VEP_result_dis.txt"))
VEP_pQTL_dropped <- fread(file.path(interim_data, "Ensembl_VEP/VEP_result_pQTL_dropped.txt"))
VEP_eQTL_dropped <- fread(file.path(interim_data, "Ensembl_VEP/VEP_result_eQTL_dropped.txt"))

#######################################################
# Labeling QTLs with Predicted Variant Effects
#######################################################

#Linking proxies with VEPs back to original SNPs

VEP_con <- VEP_con %>%
  select("#Uploaded_variation", "Allele", "Consequence", "SYMBOL") %>% 
  unique()
VEP_dis <- VEP_dis %>%
  select("#Uploaded_variation", "Allele", "Consequence", "SYMBOL") %>% 
  unique()
VEP_pQTL_dropped <- VEP_pQTL_dropped %>%
  select("#Uploaded_variation", "Allele", "Consequence", "SYMBOL") %>% 
  unique()
VEP_eQTL_dropped <- VEP_eQTL_dropped %>%
  select("#Uploaded_variation", "Allele", "Consequence", "SYMBOL") %>% 
  unique()

#Merge QTLs and VEP by rsID names, filtering for matching Gene names and effect alleles

VEP_colnames <- c("rsid", "effect_allele", "Consequence", "id")

names(VEP_con) <- VEP_colnames
names(VEP_dis) <- VEP_colnames
names(VEP_pQTL_dropped) <- VEP_colnames
names(VEP_eQTL_dropped) <- VEP_colnames

VEP_con <- VEP_con %>% 
  merge(QTLs_VEP_con, by = c("rsid", "effect_allele", "id"))

VEP_dis <- VEP_dis %>% 
  merge(QTLs_VEP_dis, by = c("rsid", "effect_allele", "id"))

VEP_pQTL_dropped <- VEP_pQTL_dropped %>% 
  merge(QTLs_VEP_pQTL_dropped, by = c("rsid", "effect_allele", "id"))

VEP_eQTL_dropped <- VEP_eQTL_dropped %>% 
  merge(QTLs_VEP_eQTL_dropped, by = c("rsid", "effect_allele", "id"))

#Assigning VEPs to QTLs

VEP_con$Consequence <- strsplit(VEP_con$Consequence, split = ",")
VEP_dis$Consequence <- strsplit(VEP_dis$Consequence, split = ",")
VEP_pQTL_dropped$Consequence <- strsplit(VEP_pQTL_dropped$Consequence, split = ",")
VEP_eQTL_dropped$Consequence <- strsplit(VEP_eQTL_dropped$Consequence, split = ",")

VEP_con <- VEP_con %>%
  mutate(row = row_number()) %>%  # Add a row identifier
  unnest('Consequence') %>%     # Unnest the list column
  mutate(present = 1) %>%         # Mark presence
  pivot_wider(
    names_from = 'Consequence',
    values_from = present,
    values_fill = 0               # Fill missing combinations with 0
  ) %>%
  dplyr::select(-row)                    # Drop helper column

VEP_dis <- VEP_dis %>%
  mutate(row = row_number()) %>%  # Add a row identifier
  unnest('Consequence') %>%     # Unnest the list column
  mutate(present = 1) %>%         # Mark presence
  pivot_wider(
    names_from = 'Consequence',
    values_from = present,
    values_fill = 0               # Fill missing combinations with 0
  ) %>%
  dplyr::select(-row)                    # Drop helper column

VEP_pQTL_dropped <- VEP_pQTL_dropped %>%
  mutate(row = row_number()) %>%  # Add a row identifier
  unnest('Consequence') %>%     # Unnest the list column
  mutate(present = 1) %>%         # Mark presence
  pivot_wider(
    names_from = 'Consequence',
    values_from = present,
    values_fill = 0               # Fill missing combinations with 0
  ) %>%
  dplyr::select(-row)                    # Drop helper column

VEP_eQTL_dropped <- VEP_eQTL_dropped %>%
  mutate(row = row_number()) %>%  # Add a row identifier
  unnest('Consequence') %>%     # Unnest the list column
  mutate(present = 1) %>%         # Mark presence
  pivot_wider(
    names_from = 'Consequence',
    values_from = present,
    values_fill = 0               # Fill missing combinations with 0
  ) %>%
  dplyr::select(-row)                    # Drop helper column

VEP_con$rsid_id <- paste(VEP_con$rsid, VEP_con$id, sep = ":")
VEP_con <- VEP_con %>%
  group_by(rsid_id) %>%
  summarise(across(everything(), max))

VEP_dis$rsid_id <- paste(VEP_dis$rsid, VEP_dis$id, sep = ":")
VEP_dis <- VEP_dis %>%
  group_by(rsid_id) %>%
  summarise(across(everything(), max))

VEP_pQTL_dropped$rsid_id <- paste(VEP_pQTL_dropped$rsid, VEP_pQTL_dropped$id, sep = ":")
VEP_pQTL_dropped <- VEP_pQTL_dropped %>%
  group_by(rsid_id) %>%
  summarise(across(everything(), max))

VEP_eQTL_dropped$rsid_id <- paste(VEP_eQTL_dropped$rsid, VEP_eQTL_dropped$id, sep = ":")
VEP_eQTL_dropped <- VEP_eQTL_dropped %>%
  group_by(rsid_id) %>%
  summarise(across(everything(), max))

#Save

fwrite(VEP_con, file.path(processed_data, "Ensembl_VEP/VEP_labelled_QTLs_con.csv"))
fwrite(VEP_dis, file.path(processed_data, "Ensembl_VEP/VEP_labelled_QTLs_dis.csv"))
fwrite(VEP_pQTL_dropped, file.path(processed_data, "Ensembl_VEP/VEP_labelled_QTLs_pQTL_dropped.csv"))
fwrite(VEP_eQTL_dropped, file.path(processed_data, "Ensembl_VEP/VEP_labelled_QTLs_eQTL_dropped.csv"))

#######################################################
#Creating VEP Summary Table
#######################################################

#Sum number of effects in each Concordance Group

VEP_summary <- colSums(VEP_con[,7:ncol(VEP_con)])
VEP_summary_2 <- colSums(VEP_dis[,7:ncol(VEP_dis)])
VEP_summary_3 <- colSums(VEP_pQTL_dropped[,7:ncol(VEP_pQTL_dropped)])
VEP_summary_4 <- colSums(VEP_eQTL_dropped[,7:ncol(VEP_eQTL_dropped)])

#Combine

VEP_summary <- merge(VEP_summary, VEP_summary_2, by = 0, all = T)
VEP_summary_3 <- merge(VEP_summary_3, VEP_summary_4, by = 0, all = T)
VEP_summary <- merge(VEP_summary, VEP_summary_3, by = "Row.names", all = T)


VEP_summary[is.na(VEP_summary)] <- 0
colnames(VEP_summary) <- c("Variant Effect", "Concordant", "Discordant", "pQTL_Dropped", "eQTL_Dropped")

#Save

VEP_summary_tab <- VEP_summary%>%
  gt() %>%
  tab_header(title = md("Variant Effect Prediction Analysis"),
             subtitle = md("Ensembl VEP"))


gtsave(VEP_summary_tab, file.path(docs_data, "Ensembl_VEP/Ensembl_summary.png"))

#######################################################
#Plotting Variant Effects
#######################################################

#Total number of QTLs in each Concordance group

total_con <- length(unique(VEP_con$rsid_id))
total_dis <- length(unique(VEP_dis$rsid_id))
total_pQTL_dropped <- length(unique(VEP_pQTL_dropped$rsid_id))
total_eQTL_dropped <- length(unique(VEP_eQTL_dropped$rsid_id))

#Plotting VEP as a proportion of each concordance category:
#Proportion = N QTLs in Variant Effect / total number of QTLs in concordance group

VEP_plot <- VEP_summary %>%
  mutate(
    Concordant = Concordant / total_con * 100,
    Discordant = Discordant / total_dis * 100,
    eQTL_Dropped = eQTL_Dropped / total_eQTL_dropped * 100,
    pQTL_Dropped = pQTL_Dropped / total_pQTL_dropped * 100
  )

#Rearrange to Plot

VEP_plot <- VEP_plot %>%
  pivot_longer(cols = c(Concordant, Discordant, eQTL_Dropped, pQTL_Dropped),
               names_to = "type_plot", values_to = "prop_plot")


#VEP Plot 1a

VEP_prop_plot <- ggplot(VEP_plot, aes(fill=type_plot, y=prop_plot, x=`Variant Effect`)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(text = element_text(size=10, family = "Times"),
        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5)) +
  labs(y = "Proportion of Concordance Category", x = "Predicted Variant Effect", fill = "Concordance Category")

# VEP Plot 1b - Con/Dis Groups only

VEP_plot_discon <- VEP_plot[which(VEP_plot$type_plot == "Concordant" | VEP_plot$type_plot == "Discordant"),]

VEP_prop_plot_discon <- ggplot(VEP_plot_discon, aes(fill=type_plot, y=prop_plot, x=`Variant Effect`)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(text = element_text(size=10, family = "Times"),
        axis.text.x = element_text(size=8, angle=90, hjust = 1, vjust = 0.5)) +
  labs(y = "Proportion of Concordance Category", x = "Protein Class", fill = "Concordance Category")

#Plotting concordance categories as a proportion of each variant effect
#Proportion = N QTLs with variant effect / total number of QTLs with variant effect

VEP_plot_stacked <- VEP_summary

for (i in 1:nrow(VEP_plot_stacked)){
  VEP_sum <- sum(VEP_plot_stacked[i,2:5])
  VEP_plot_stacked$Concordant[i] <- (VEP_plot_stacked$Concordant[i]/VEP_sum) * 100
  VEP_plot_stacked$Discordant[i] <- (VEP_plot_stacked$Discordant[i]/VEP_sum) * 100
  VEP_plot_stacked$eQTL_Dropped[i] <- (VEP_plot_stacked$eQTL_Dropped[i]/VEP_sum) * 100
  VEP_plot_stacked$pQTL_Dropped[i] <- (VEP_plot_stacked$pQTL_Dropped[i]/VEP_sum) * 100
}

#Rearrange for plotting

VEP_num <- VEP_summary %>% 
  pivot_longer(cols = c(Concordant, Discordant, eQTL_Dropped, pQTL_Dropped),
               names_to = "type_plot", values_to = "num_plot")


VEP_plot_stacked <- VEP_plot_stacked %>%
  pivot_longer(cols = c(Concordant, Discordant, eQTL_Dropped, pQTL_Dropped),
               names_to = "type_plot", values_to = "prop_plot") %>% 
  mutate(num_plot = VEP_num$num_plot)

#VEP plot 2 - Stacked

VEP_prop_plot_stacked <- ggplot(VEP_plot_stacked, aes(fill=type_plot, y=prop_plot, x=`Variant Effect`, label=num_plot)) + 
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

###FIGURE S7

figs7 <- VEP_prop_plot_discon

ggsave(file = file.path(docs_data, "Ensembl_VEP/FigS7.eps"), plot = figs7, width = 4, height = 4, 
       units = "in", dpi = 600, family = "Times", device = postscript, paper = "special", horizontal = FALSE)
