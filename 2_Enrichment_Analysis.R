#Enrichment Analysis of 1-1 Matched Primary QTL Pairs

#Open Target Data - Baseline Expression (Blood vs non blood)
#HPA Data - Protein Function annotation
#Olink - Alternative Protein Function annotation
#Reactome - Pathway annotation
#Ensembl VEP - Variant effects

#######################################################
#Load in libraries
#######################################################

library(arrow)
library(data.table)
library(biomaRt)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(rlist)
library(TwoSampleMR)
library(ieugwasr)
library(gt)
library(pheatmap)

#######################################################
#Set Working Directory
#######################################################

setwd("")

#######################################################
#Read in harmonised QTL pairs
#######################################################

QTLs <- fread("1-1_concordance_primary.csv")

#######################################################
#Read in Open Target Data
#######################################################

OT_files <- list.files("data/open_targets_expression/", full.names = TRUE) #Expression files from Open Target
OT_files <- OT_files[-1] #First file is empty

#Read in each file and combine into one table

OT_data <- read_parquet(OT_files[1], col_select = NULL, as_tibble = TRUE)

for (f in 2:length(OT_files)){
  tmp <- read_parquet(OT_files[f])
  OT_data <- bind_rows(OT_data, tmp)
}

#######################################################
#Labeling QTL data for Blood Expression
#######################################################

QTLs_OT <- QTLs %>%
  dplyr::select("SNP", "id.outcome", "gene.outcome", "type") 

#Taking just the genes in both datasets

OT <- OT_data[which(OT_data$id %in% QTLs_OT$gene.outcome),]
QTLs_OT <- QTLs_OT[which(QTLs_OT$gene.outcome %in% OT$id),]

#Ordering genes in both datasets

OT <- OT[order(OT$id),]
QTLs_OT <- QTLs_OT[order(gene.outcome),]

#Extracting top blood expression, top non-blood expression and most expressed organ for each gene

QTLs_OT$blood.expression <- NA
QTLs_OT$nonblood.expression <- NA
QTLs_OT$top.organ <- NA

for (i in 1:nrow(OT)){
  bexp <- character()
  bnonexp <- character()
  
  OT_table <- as.data.frame(OT[[2]][[i]])
  OT_rna <- as.data.frame(OT_table[[5]])
  
  if (OT$id[i] != QTLs_OT$gene.outcome[i]){ #Making sure the genes match(if they don't, STOP)
    print("Genes don't match")
  }
  
  for (j in 1:nrow(OT_table)){
    if (length(grep("blood", OT_table$organs[[j]])) > 0){ 
      bexp <- append(bexp, OT_rna[j,1]) #All blood associated expression
    }else{ 
      bnonexp <- append(bnonexp, OT_rna[j, 1]) #All non blood associated expression
    }
  }
  
  QTLs_OT$blood.expression[i] <- max(as.numeric(bexp)) #Highest blood associated expression
  QTLs_OT$nonblood.expression[i] <- max(as.numeric(bnonexp)) #Highest non blood associated expression
  QTLs_OT$top.organ[i] <- OT_table$organs[grep(max(OT_table$rna$value), OT_table$rna$value)][1] #Top expressed organ
  
  print(i)
}

#Calculating and catagorising blood expression for each gene

##prop_blood = top blood associated expression / top non blood associated expression
#Group 1 = <25%
#Group 2 = 25 - 50%
#Group 3 = 50 - 75%
#Group 4 = 75 = 100%
#Group 5 = >100% (blood is top organ)

QTLs_OT$prop_blood <- (QTLs_OT$blood.expression/QTLs_OT$nonblood.expression) * 100

QTLs_OT$blood_group <- NA

for (i in 1:nrow(QTLs_OT)){
  if(QTLs_OT$prop_blood[i] < 25){
    QTLs_OT$blood_group[i] <- 1
  }
  if (QTLs_OT$prop_blood[i] > 25 & QTLs_OT$prop_blood[i] < 50){
    QTLs_OT$blood_group[i] <- 2
  } 
  if (QTLs_OT$prop_blood[i] > 50 & QTLs_OT$prop_blood[i] < 75){
    QTLs_OT$blood_group[i] <- 3
  } 
  if (QTLs_OT$prop_blood[i] > 75 & QTLs_OT$prop_blood[i] < 100){
    QTLs_OT$blood_group[i] <- 4
  } 
  if (QTLs_OT$prop_blood[i] > 100){
    QTLs_OT$blood_group[i] <- 5
  }
  print(i)
}

#Save

fwrite(QTLs_OT, "QTLs_OT.csv")

#######################################################
#Saving OT Labelled QTL
#######################################################

QTLs_OT_tab <- QTLs_OT %>%
  dplyr::select("SNP", "id.outcome", "type", "top.organ", "prop_blood", "blood_group")

colnames(QTLs_OT_tab) <- c("Protein:rsID", "Protein", "Concordance_Group", "Top_organ", "Proportion_Blood_Expression", "Blood_Group")

fwrite(QTLs_OT_tab, "OT_labelled_QTLs.csv")

#######################################################
#Creating OT Summary Table
#######################################################

#Separating concordance groups and calculating number of QTLs in each blood group

QTLs_OT_con <- QTLs_OT[which(QTLs_OT$type == "concordant"),]
QTLs_OT_con <- QTLs_OT_con[, .N, by = 'blood_group'] 
colnames(QTLs_OT_con)[2] <- "Concordant"

QTLs_OT_dis <- QTLs_OT[which(QTLs_OT$type == "discordant"),]
QTLs_OT_dis <- QTLs_OT_dis[, .N, by = 'blood_group']
colnames(QTLs_OT_dis)[2] <- "Discordant"

QTLs_OT_eQTL_dropped <- QTLs_OT[which(QTLs_OT$type == "eQTL_dropped"),]
QTLs_OT_eQTL_dropped <- QTLs_OT_eQTL_dropped[, .N, by = 'blood_group']
colnames(QTLs_OT_eQTL_dropped)[2] <- "eQTL_dropped"

QTLs_OT_pQTL_dropped <- QTLs_OT[which(QTLs_OT$type == "pQTL_dropped"),]
QTLs_OT_pQTL_dropped <- QTLs_OT_pQTL_dropped[, .N, by = 'blood_group']
colnames(QTLs_OT_pQTL_dropped)[2] <- "pQTL_dropped"

#Combine into one summary table

OT_summary <- merge(QTLs_OT_con, QTLs_OT_dis, by = 'blood_group')
OT_summary <- merge(OT_summary, QTLs_OT_eQTL_dropped, by = 'blood_group')
OT_summary <- merge(OT_summary, QTLs_OT_pQTL_dropped, by = 'blood_group')

#Save

OT_summary_tab <- OT_summary%>%
  gt() %>%
  tab_header(title = md("Blood Expression Analysis"),
             subtitle = md("Open Target Data")) %>%
  tab_footnote(
    footnote = "Where blood groups are catagorised based on protein expression in blood:"
  ) %>%
  tab_footnote(
    footnote = "Group 1 = <25%"
  ) %>%
  tab_footnote(
    footnote = "Group 2 = 25 - 50%"
  ) %>%
  tab_footnote(
    footnote = "Group 3 = 50 - 75%"
  ) %>%
  tab_footnote(
    footnote = "Group 4 = 75 - 100%"
  ) %>%
  tab_footnote(
    footnote = "Group 5 = >100% (blood is top organ)"
  )


gtsave(OT_summary_tab, "")

#######################################################
#Plotting Blood Expression
#######################################################

#Total number of QTLs in each concordance group

total_con <- nrow(QTLs_OT[which(QTLs_OT$type == "concordant"),])
total_dis <- nrow(QTLs_OT[which(QTLs_OT$type == "discordant"),])
total_eQTL_dropped <- nrow(QTLs_OT[which(QTLs_OT$type == "eQTL_dropped"),])
total_pQTL_dropped <- nrow(QTLs_OT[which(QTLs_OT$type == "pQTL_dropped"),])

#Plotting blood expression as a proportion of each concordance category
#Proportion = N QTLs in blood group / total number of QTLs in concordance group

OT_summary$Concordant_prop <- (OT_summary$Concordant/total_con) * 100
OT_summary$Discordant_prop <- (OT_summary$Discordant/total_dis) * 100
OT_summary$eQTL_dropped_prop <- (OT_summary$eQTL_dropped/total_eQTL_dropped) * 100
OT_summary$pQTL_dropped_prop <- (OT_summary$pQTL_dropped/total_pQTL_dropped) * 100

#Rearrange to plot

blood_plot <- rep(c("1", "2", "3", "4", "5"), 4)
type_plot <- c(rep("Concordant", 5), rep("Discordant", 5), rep("eQTL_dropped", 5), rep("pQTL_dropped", 5))
prop_plot <- c(OT_summary$Concordant_prop, OT_summary$Discordant_prop, OT_summary$eQTL_dropped_prop, OT_summary$pQTL_dropped_prop)

OT_plot <- data.frame(blood_plot, type_plot, prop_plot)

#Open Target Plot 1a - all Concordance groups

OT_prop_plot <- ggplot(OT_plot, aes(fill=type_plot, y=prop_plot, x=blood_plot)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=18),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))

png("", height = 1200, width = 882)

OT_prop_plot + 
  plot_annotation(title = '1. Proportion of Concordance Type QTLs pairs in each Blood Group - As a proportion of each concordance category')

dev.off()

#Open Target Plot 1b - Con/Dis groups only

OT_plot_discon <- OT_plot[which(OT_plot$type_plot == "Concordant" | OT_plot$type_plot == "Discordant"),]

OT_prop_plot_discon <- ggplot(OT_plot_discon, aes(fill=type_plot, y=prop_plot, x=blood_plot)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=18),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))

png("", height = 1200, width = 882)

OT_prop_plot_discon + 
  plot_annotation(title = '1. Proportion of Concordance Type QTLs pairs in each Blood Group - As a proportion of each concordance category')

dev.off()

#Plotting concordance categories as a proportion of each blood group
#Proportion = N QTLs in blood group / total number of QTLs in blood group

for (i in 1:nrow(OT_summary)){
  blood_sum <- sum(OT_summary[i,2:5]) #total number of QTLs in blood group
  OT_summary$Concordant_prop[i] <- (OT_summary$Concordant[i]/blood_sum) * 100
  OT_summary$Discordant_prop[i] <- (OT_summary$Discordant[i]/blood_sum) * 100
  OT_summary$eQTL_dropped_prop[i] <- (OT_summary$eQTL_dropped[i]/blood_sum) * 100
  OT_summary$pQTL_dropped_prop[i] <- (OT_summary$pQTL_dropped[i]/blood_sum) * 100
}

#Rearrange to plot

blood_plot <- rep(c("1", "2", "3", "4", "5"), 4)
type_plot <- c(rep("Concordant", 5), rep("Discordant", 5), rep("eQTL_dropped", 5), rep("pQTL_dropped", 5))
prop_plot <- c(OT_summary$Concordant_prop, OT_summary$Discordant_prop, OT_summary$eQTL_dropped_prop, OT_summary$pQTL_dropped_prop)
num_plot <- c(OT_summary$Concordant, OT_summary$Discordant, OT_summary$eQTL_dropped, OT_summary$pQTL_dropped)

OT_plot <- data.frame(blood_plot, type_plot, prop_plot, num_plot)

#Open Target Plot 2 - Stacked plot

OT_prop_plot_stacked <- ggplot(OT_plot, aes(fill=type_plot, y=prop_plot, x=blood_plot, label = num_plot)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=18),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18)) +
  geom_text(size = 5, position = position_stack(vjust = 0.5))


png("", height = 1200, width = 882)

OT_prop_plot_stacked + 
  plot_annotation(title = '2. Proportion of Concordance Type QTLs pairs in each Blood Group - As a proportion of the blood group')

dev.off()


#######################################################
#Read in HPA Data
#######################################################

HPA <- fread("data/human_protein_atlas/proteinatlas.tsv")

#######################################################
#Labeling QTL data for Protein Class
#######################################################

HPA_protein_class <- HPA %>% 
  dplyr::select("Gene", "Ensembl", "Protein class")

QTLs_HPA <- QTLs %>%
  dplyr::select("SNP", "id.outcome", "gene.outcome", "type")
colnames(QTLs_HPA)[3] <- "Ensembl"

#Assigning protein classes to genes based on Ensembl ID

QTLs_HPA <- merge(QTLs_HPA, HPA_protein_class, by = "Ensembl")

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
  dplyr::select(-row)                    # Drop helper column

fwrite(QTLs_HPA, "HPA_labelled_QTLs.csv")

#######################################################
#Creating HPA Summary Table
#######################################################

#Separating QTLs into concordance type

QTLs_HPA_con <- QTLs_HPA[which(QTLs_HPA$type == "concordant"),][,6:29]
QTLs_HPA_dis <- QTLs_HPA[which(QTLs_HPA$type == "discordant"),][,6:29]
QTLs_HPA_eQTL_dropped <- QTLs_HPA[which(QTLs_HPA$type == "eQTL_dropped"),][,6:29]
QTLs_HPA_pQTL_dropped <- QTLs_HPA[which(QTLs_HPA$type == "pQTL_dropped"),][,6:29]

#Summing number of QTLs in each protein class for each concordance group

HPA_summary <- colSums(QTLs_HPA_con)
HPA_summary <- cbind(HPA_summary, colSums(QTLs_HPA_dis))
HPA_summary <- cbind(HPA_summary, colSums(QTLs_HPA_eQTL_dropped))
HPA_summary <- cbind(HPA_summary, colSums(QTLs_HPA_pQTL_dropped))
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


gtsave(HPA_summary_tab, "")

#######################################################
#Plotting Protein Classes
#######################################################

#Total number of QTLs in each concordance group

total_con <- nrow(QTLs_HPA[which(QTLs_HPA$type == "concordant"),])
total_dis <- nrow(QTLs_HPA[which(QTLs_HPA$type == "discordant"),])
total_eQTL_dropped <- nrow(QTLs_HPA[which(QTLs_HPA$type == "eQTL_dropped"),])
total_pQTL_dropped <- nrow(QTLs_HPA[which(QTLs_HPA$type == "pQTL_dropped"),])

#Plotting Protein Class involvement as a proportion of each concordance category:
#Proportion = N QTLs in protein class / total number of QTLs in concordance group

HPA_summary$Concordant_prop <- (HPA_summary$Concordant/total_con) * 100
HPA_summary$Discordant_prop <- (HPA_summary$Discordant/total_dis) * 100
HPA_summary$eQTL_dropped_prop <- (HPA_summary$eQTL_dropped/total_eQTL_dropped) * 100
HPA_summary$pQTL_dropped_prop <- (HPA_summary$pQTL_dropped/total_pQTL_dropped) * 100

#Rearrange to plot

protein_plot <- rep(HPA_summary$`Protein Class`, 4)
type_plot <- c(rep("Concordant", 24), rep("Discordant", 24), rep("eQTL_dropped", 24), rep("pQTL_dropped", 24))
prop_plot <- c(HPA_summary$Concordant_prop, HPA_summary$Discordant_prop, HPA_summary$eQTL_dropped_prop, HPA_summary$pQTL_dropped_prop)

HPA_plot <- data.frame(protein_plot, type_plot, prop_plot)

#Human Protein Atlas Plot 1a - all Concordance groups

HPA_prop_plot <- ggplot(HPA_plot, aes(fill=type_plot, y=prop_plot, x=protein_plot)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(size=13, angle=90),
        axis.text.y = element_text(size=13),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10))

png("", height = 1200, width = 882)

HPA_prop_plot + 
  plot_annotation(title = '1. Proportion of Concordance Type QTLs pairs in each Protein Class - As a proportion of each concordance category')

dev.off()

#Human Protein Atlas Plot 1b - Con/Dis Groups only

HPA_plot_discon <- HPA_plot[which(HPA_plot$type_plot == "Concordant" | HPA_plot$type_plot == "Discordant"),]

HPA_prop_plot_discon <- ggplot(HPA_plot_discon, aes(fill=type_plot, y=prop_plot, x=protein_plot)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text = element_text(size=15, angle=90),
        axis.title = element_text(size=18),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))

png("", height = 1200, width = 882)

HPA_prop_plot_discon + 
  plot_annotation(title = '1. Proportion of Concordance Type QTLs pairs in each Blood Group - As a proportion of each concordance category')

dev.off()


#Plotting concordance categories as a proportion of each protein class
#Proportion = N QTLs in protein class / total number of QTLs in protein class

for (i in 1:nrow(HPA_summary)){
  protein_sum <- sum(HPA_summary[i,2:5])
  HPA_summary$Concordant_prop[i] <- (HPA_summary$Concordant[i]/protein_sum) * 100
  HPA_summary$Discordant_prop[i] <- (HPA_summary$Discordant[i]/protein_sum) * 100
  HPA_summary$eQTL_dropped_prop[i] <- (HPA_summary$eQTL_dropped[i]/protein_sum) * 100
  HPA_summary$pQTL_dropped_prop[i] <- (HPA_summary$pQTL_dropped[i]/protein_sum) * 100
}

#Rearrange to plot
protein_plot <- rep(HPA_summary$`Protein Class`, 4)
type_plot <- c(rep("Concordant", 24), rep("Discordant", 24), rep("eQTL_dropped", 24), rep("pQTL_dropped", 24))
prop_plot <- c(HPA_summary$Concordant_prop, HPA_summary$Discordant_prop, HPA_summary$eQTL_dropped_prop, HPA_summary$pQTL_dropped_prop)
num_plot <- c(HPA_summary$Concordant, HPA_summary$Discordant, HPA_summary$eQTL_dropped, HPA_summary$pQTL_dropped)

HPA_plot <- data.frame(protein_plot, type_plot, prop_plot, num_plot)

#Human Protein Atlas Plot 3 - Stacked plot

HPA_prop_plot_stacked <- ggplot(HPA_plot, aes(fill=type_plot, y=prop_plot, x=protein_plot, label=num_plot)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(size=10, angle=90)) +
  geom_text(size = 5, position = position_stack(vjust = 0.5))

png("", height = 1200, width = 882)

HPA_prop_plot_stacked + 
  plot_annotation(title = '2. Proportion of Concordance Type QTLs pairs in each Protein Class - As a proportion of the protein class')

dev.off()

#######################################################
#Read in Olink Data
#######################################################

olink_files <- list.files("data/olink/input/", full.names = T) 

#######################################################
#Labeling QTL data for Protein Function
#######################################################

QTLs_olink <- QTLs %>%
  dplyr::select("SNP", "id.outcome", "gene.outcome", "type")

#Reading in each olink panel and assigning function to those SNPs present in panel

QTLs_protein_function <- data.frame()

for (f in olink_files){
  
  print(f)
  tmp <- fread(f) #Read in SNPs in each file
  
  class <- unlist(strsplit(f, split = " ")) #Extract protein function from file name
  class <- class[4:(length(class) - 2)]
  print(class)
  
  QTL_class <- QTLs_olink[QTLs_olink$id.outcome %in% tmp$Gene,] #Assigning protein class to QTLs present in panel
  QTL_class$class <- paste(class, collapse = "_")
  
  QTLs_protein_function <- rbind(QTLs_protein_function, QTL_class) #Combine to one dataset
}

#Combining subcategories into one

QTLs_protein_function[which(QTLs_protein_function$class == "Cardiovascular_II" | QTLs_protein_function$class == "Cardiovascular_III"),]$class <- "Cardiovascular"
QTLs_protein_function[which(QTLs_protein_function$class == "Oncology_II" | QTLs_protein_function$class == "Oncology_III"),]$class <- "Oncology"

fwrite(QTLs_protein_function, "Olink_labelled_QTLs.csv")

#######################################################
#Creating Olink Summary Table
#######################################################

#Grouping by concordance type and summing number of protein functions present

QTLs_protein_function_con <- QTLs_protein_function[which(QTLs_protein_function$type == "concordant"),]
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


gtsave(olink_summary_tab, "")

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
  theme(axis.text.x = element_text(size=15, angle=90),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=18),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))

png("", height = 1200, width = 882)

olink_prop_plot + 
  plot_annotation(title = '1. Proportion of Concordance Type QTLs pairs in each Protein Function - As a proportion of each concordance category')

dev.off()

#Olink Plot 1b - Con/Dis Groups only

olink_plot_discon <- olink_plot[which(olink_plot$type_plot == "Concordant" | olink_plot$type_plot == "Discordant"),]

olink_prop_plot_discon <- ggplot(olink_plot_discon, aes(fill=type_plot, y=prop_plot, x=olink_plot)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(size=15, angle=90),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=18),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))

png("", height = 1200, width = 882)

olink_prop_plot_discon + 
  plot_annotation(title = '1. Proportion of Concordance Type QTLs pairs in each Blood Group - As a proportion of each concordance category')

dev.off()


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
  theme(axis.text.x = element_text(size=15, angle=90),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=18),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18)) +
  geom_text(size = 5, position = position_stack(vjust = 0.5))

png("", height = 1200, width = 882)

olink_prop_plot_stacked + 
  plot_annotation(title = '2. Proportion of Concordance Type QTLs pairs in each Protein Function - As a proportion of the function')

dev.off()

#######################################################
#Read in Reactome Data
#######################################################

reactome <- fread("data/reactome/Ensembl2Reactome.tsv", header = FALSE)

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

fwrite(QTLs_reactome, "Reactome_labelled_QTLs.csv")

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

fwrite(reactome_summary, "")

reactome_summary_tab <- reactome_summary%>%
  gt() %>%
  tab_header(title = md("Pathway Analysis"),
             subtitle = md("Reactome"))


gtsave(reactome_summary_tab, "")

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

#Rearrange to plot

pathway_plot <- rep(reactome_summary$`Event (Pathway or Reaction) Name`, 4)
type_plot <- c(rep("Concordant", 138), rep("Discordant", 138), rep("eQTL_dropped", 138), rep("pQTL_dropped", 138))
prop_plot <- c(reactome_summary$Concordant_prop, reactome_summary$Discordant_prop, reactome_summary$eQTL_dropped_prop, reactome_summary$pQTL_dropped_prop)

reactome_plot <- data.frame(pathway_plot, type_plot, prop_plot)

#Reactome Plot 1a - all Concordance groups

reactome_prop_plot <- ggplot(reactome_plot, aes(fill=type_plot, y=prop_plot, x=pathway_plot)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(size=6, angle=90))

png("", height = 1200, width = 882)

reactome_prop_plot + 
  plot_annotation(title = '1. Proportion of Concordance Type QTLs pairs in each Pathway - As a proportion of each concordance category')

dev.off()

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
  theme(axis.text.x = element_text(size=10, angle=90))

png("", height = 1200, width = 882)

reactome_prop_plot_sig + 
  plot_annotation(title = '1. Proportion of Concordance Type QTLs pairs in each Pathway - As a proportion of each concordance category')

dev.off()

#Save pathways

fwrite(reactome_summary_sig, "sig_pathways.csv")

#Reactome Plot 1c - Con/Dis groups only

reactome_plot_discon <- reactome_plot[which(reactome_plot$type_plot == "Concordant" | reactome_plot$type_plot == "Discordant"),]

reactome_prop_plot_discon <- ggplot(reactome_plot_discon, aes(fill=type_plot, y=prop_plot, x=pathway_plot)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(size=6, angle=90),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=18),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))

png("", height = 1200, width = 882)

reactome_prop_plot_discon + 
  plot_annotation(title = '1. Proportion of Concordance Type QTLs pairs in each Blood Group - As a proportion of each concordance category')

dev.off()

#Reactome Plot 1d - Pathways with multiple QTLs involved, Discon only

reactome_plot_discon_sig <- reactome_plot_sig[which(reactome_plot_sig$type_plot == "Concordant" | reactome_plot_sig$type_plot == "Discordant"),]

reactome_prop_plot_discon_sig <- ggplot(reactome_plot_discon_sig, aes(fill=type_plot, y=prop_plot, x=pathway_plot)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(size=10, angle=90),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=18),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))

png("/", height = 1200, width = 882)

reactome_prop_plot_discon_sig + 
  plot_annotation(title = '1. Proportion of Concordance Type QTLs pairs in each Blood Group - As a proportion of each concordance category')

dev.off()

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
  theme(axis.text.x = element_text(size=8, angle=90),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=18),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18)) +
  geom_text(size = 4, position = position_stack(vjust = 0.5))

png("", height = 2200, width = 1200)

reactome_prop_plot_stacked + 
  plot_annotation(title = '2. Proportion of Concordance Type QTLs pairs in each Protein Class - As a proportion of the protein class')

dev.off()

#Reactome Plot 2b - Pathways with multiple QTLs involved

reactome_plot_stacked_sig <- reactome_plot[which(reactome_plot$pathway_plot %in% reactome_plot_sig$pathway_plot),]

reactome_plot_stacked_sig <- ggplot(reactome_plot_stacked_sig, aes(fill=type_plot, y=prop_plot, x=pathway_plot, label=num_plot)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(size=12, angle=90),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=18),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18)) +
  geom_text(size = 4, position = position_stack(vjust = 0.5))

png("", height = 1500, width = 882)

reactome_plot_stacked_sig + 
  plot_annotation(title = '2. Proportion of Concordance Type QTLs pairs in each Protein Class - As a proportion of the protein class')

dev.off()

#######################################################
#Variant Functional Annotations
#######################################################
#######################################################
#Extracting Concordant and Discordant SNPs for Proxy Search
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

write.table(QTLs_VEP_con_SNPs, file = "data/SNPs_for_proxies/SNPs_con.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

write.table(QTLs_VEP_dis_SNPs, file = "data/SNPs_for_proxies/SNPs_dis.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

#plink --bfile /data/references/EUR --ld-snp-list SNPs_con.txt --ld-window-kb 500 --ld-window 99999 --ld-window-r2 0.8 --out SNPs_con_ld
#plink --bfile /data/references/EUR --ld-snp-list SNPs_dis.txt --ld-window-kb 500 --ld-window 99999 --ld-window-r2 0.8 --out SNPs_dis_ld

#######################################################
#Reading back in Con/Dis SNPs with proxies - Saving for VEP
#######################################################

VEP_con_prox <- fread("data/SNPs_for_proxies/SNPs_con_ld.ld")
VEP_dis_prox <- fread("data/SNPs_for_proxies/SNPs_dis_ld.ld")

#Extracting all proxy SNPs

VEP_con_prox_SNPs <- unique(VEP_con_prox$SNP_B)
VEP_dis_prox_SNPs <- unique(VEP_dis_prox$SNP_B)

#Save

write.table(VEP_con_prox_SNPs, file = "data/EnsemblVEP/VEP_con_prox_SNPs.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

write.table(VEP_dis_prox_SNPs, file = "data/EnsemblVEP/VEP_dis_prox_SNPs.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

#SNP lists were read into Ensembl VEP

#######################################################
#Reading back in the VEP SNPs
#######################################################

VEP_con <- fread("data/EnsemblVEP/VEP_con_res.txt")
VEP_dis <- fread("data/EnsemblVEP/VEP_dis_res.txt")

#######################################################
#Labeling QTLs with Predicted Variant Effects
#######################################################

#Linking proxies with VEPs back to original SNPs

VEP_con <- VEP_con %>%
  select("#Uploaded_variation", "Allele", "Consequence", "SYMBOL")
VEP_dis <- VEP_dis %>%
  select("#Uploaded_variation", "Allele", "Consequence", "SYMBOL")

VEP_con_prox <- VEP_con_prox %>%
  select("SNP_A", "SNP_B")
VEP_dis_prox <- VEP_dis_prox %>%
  select("SNP_A", "SNP_B")

colnames(VEP_con)[1] <- "SNP"
colnames(VEP_dis)[1] <- "SNP"

colnames(VEP_con_prox)[2] <- "SNP"
colnames(VEP_dis_prox)[2] <- "SNP"

VEP_con <- merge(VEP_con, VEP_con_prox, by = "SNP") #Merge based on rsIDs
VEP_dis <- merge(VEP_dis, VEP_dis_prox, by = "SNP") #Seemed to increase rows, not sure why...

#Linking SNPs with VEPs back to the Gene they were associated with

colnames(QTLs_VEP_con)[7] <- "SNP_A"
colnames(QTLs_VEP_dis)[7] <- "SNP_A"

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

fwrite(con_QTLs, "VEP_labelled_QTLs_con.csv")
fwrite(dis_QTLs, "VEP_labelled_QTLs_dis.csv")

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


gtsave(VEP_summary_tab, "")

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
  theme(axis.text.x = element_text(size=13, angle=90),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))

png("", height = 1200, width = 882)

VEP_prop_plot + 
  plot_annotation(title = '1. Proportion of Concordance Type QTLs pairs in each Protein Class - As a proportion of each concordance category')

dev.off()

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
  theme(axis.text.x = element_text(size=15, angle=90),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=18),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18)) +
  geom_text(size = 4, position = position_stack(vjust = 0.5))

png("", height = 1200, width = 882)

VEP_prop_plot_stacked + 
  plot_annotation(title = '2. Proportion of Concordance Type QTLs pairs in each Protein Class - As a proportion of the protein class')

dev.off()
