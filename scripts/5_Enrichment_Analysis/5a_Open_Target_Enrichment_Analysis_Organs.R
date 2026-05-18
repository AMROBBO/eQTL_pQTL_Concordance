#######################################################
# 5a. Enrichment Analysis of Colocalised QTL Pairs
# Datasets:
#   a. Open Targets: Top Organ Expression
##
# As there is little variation between concordance of datasets selected on pQTLs
# and eQTLs and those with all available independent SNPs and those with only 1
# SNP per gene, the rest of the enrichment analysis is carried out on all available
# independent SNPs selected from pQTL data:
# pQTL_concordance.csv from 4_Naive_Concordance_Test.R
##
#######################################################


#######################################################
# Initialising file paths
#######################################################

load_dot_env("config.env")

raw_data <- Sys.getenv("rawdatadir")
interim_data <- Sys.getenv("interimdatadir")
processed_data <- Sys.getenv("processeddatadir")
docs_data <- Sys.getenv("docsdir")

#######################################################
# Read in OT labelled QTL pairs 
#######################################################

QTLs_OT <- fread(file.path(processed_data, "Open_Target/OT_labelled_QTLs.csv"))

#######################################################
# Creating summary table for top organ
#######################################################

organs <- unique(QTLs_OT$top_organ)

OT_summary <- data.frame(Top_Organ = organs)

for (i in 1:length(organs)){
  top_organ <- organs[i]
  
  OT_summary$Concordant[i] <- sum(QTLs_OT$Concordance_Group == "concordant" & QTLs_OT$top_organ == top_organ)
  OT_summary$Discordant[i] <- sum(QTLs_OT$Concordance_Group == "discordant" & QTLs_OT$top_organ == top_organ)
  OT_summary$eQTL_dropped[i] <- sum(QTLs_OT$Concordance_Group == "eQTL_dropped" & QTLs_OT$top_organ == top_organ)
  OT_summary$pQTL_dropped[i] <- sum(QTLs_OT$Concordance_Group == "pQTL_dropped" & QTLs_OT$top_organ == top_organ)
}

#Save summary as a table

OT_summary_tab <- OT_summary %>%
  gt() %>%
  tab_header(title = md("Organ Expression Analysis"),
             subtitle = md("Open Target Data")) %>%
  fmt_number(
    columns = 1,
    decimals = 0,
    use_seps = TRUE
  )

gtsave(OT_summary_tab, file.path(docs_data, "Open_Target/OT_summary_organs.png"))

#######################################################
# Plotting Organ Expression
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
               names_to = "type_plot", values_to = "prop_plot") 

# Open Target Plot 1a - all Concordance groups

OT_prop_plot <- ggplot(OT_plot, aes(x = Top_Organ, y = prop_plot, fill = type_plot)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(text = element_text(size = 8, family = "Times"),
        axis.text.x = element_text(size = 8)) +
  labs(y = "Proportion of Concordance Category", x = "Blood Expression/TPM", fill = "Concordance Category")

# Open Target Plot 1b - Con/Dis groups only

OT_plot_discon <- OT_plot %>% 
  filter(type_plot == "Concordant" | type_plot == "Discordant")

OT_prop_plot_discon <- ggplot(OT_plot_discon, aes(x=Top_Organ, y=prop_plot, fill=type_plot)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(text = element_text(size=8, family = "Times"),
        axis.text.x = element_text(size=8)) +
  labs(y = "Proportion of Concordance Category", x = "Blood Expression/TPM", fill = "Concordance Category")

#Plotting concordance categories as a proportion of each blood group
#Proportion = N QTLs in blood group / total number of QTLs in blood group

OT_plot_stacked <- OT_summary

for (i in 1:nrow(OT_plot_stacked)){
  organ_sum <- sum(OT_plot_stacked[i,2:5]) #total number of QTLs in blood group
  OT_plot_stacked$Concordant[i] <- (OT_plot_stacked$Concordant[i]/organ_sum) * 100
  OT_plot_stacked$Discordant[i] <- (OT_plot_stacked$Discordant[i]/organ_sum) * 100
  OT_plot_stacked$eQTL_dropped[i] <- (OT_plot_stacked$eQTL_dropped[i]/organ_sum) * 100
  OT_plot_stacked$pQTL_dropped[i] <- (OT_plot_stacked$pQTL_dropped[i]/organ_sum) * 100
}

#Rearrange to plot

OT_num <- OT_summary %>% 
  pivot_longer(cols = c(Concordant, Discordant, eQTL_dropped, pQTL_dropped),
               names_to = "type_plot", values_to = "num_plot")

OT_plot_stacked <- OT_plot_stacked %>%
  pivot_longer(cols = c(Concordant, Discordant, eQTL_dropped, pQTL_dropped),
               names_to = "type_plot", values_to = "prop_plot") %>%
  mutate(num_plot = OT_num$num_plot) 

OT_plot_stacked$num_plot <- gsub("^0$", "", OT_plot_stacked$num_plot)

# Open Target Plot 2 - Stacked plot

OT_prop_plot_stacked <- ggplot(OT_plot_stacked, aes(fill=type_plot, y=prop_plot, x=Top_Organ, label = num_plot)) + 
  geom_bar(position="stack", stat="identity") +
  theme(text = element_text(size=8, family = "Times"),
        axis.text.x = element_text(size=8)) +
  geom_text(size = 2, position = position_stack(vjust = 0.5), family = "Times") +
  labs(y = "Proportion of Blood Expression Group", x = "Blood Expression/TPM", fill = "Concordance Category")

#######################################################
#Saving OT Figures
#######################################################

###FIGURE S7

figs7a <- OT_prop_plot +
  theme(legend.position = "none") +
  coord_flip()

figs7b <- OT_prop_plot_stacked +
  theme(axis.title.y = element_blank()) +
  coord_flip()

figs7 <- figs7a + figs7b

figs7 <- fig2 & theme(plot.margin = unit(c(2, 2, 2, 2), "pt"))

ggsave(file = file.path(docs_data, "Open_Target/FigS7.eps"), plot = figs7, width = 10, height = 4, 
       units = "in", dpi = 600, family = "Times", device = postscript, paper = "special", horizontal = FALSE)

###FIGURE S8

figs8 <- OT_prop_plot_discon +
  coord_flip()

ggsave(file = file.path(docs_data, "Open_Target/FigS8.eps"), plot = figs8, width = 6, height = 4, 
       units = "in", dpi = 600, family = "Times", device = postscript, paper = "special", horizontal = FALSE)
