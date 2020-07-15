library(tidyverse)

# load connectivity map results

setwd("~/OneDrive - King's College London/ConnectivityMap/CLUE")


clue_results <- read_csv("clue_results.csv")

# rename cols
colnames(clue_results)[11:17] <- c("Post_Mortem","Progerin","HPC_6hr",
                                   "HPC_72hr", "HPC_144hr",	"HPC_1hr",	
                                   "HPC_24hr")

# Create correlation matrix to show 6hr similarity (and 24hr dissimilarity) to microarray and progerin signatures


clue_df <- clue_results[, 11:17]

clue_df$HPC_6hr <- as.numeric(clue_df$HPC_6hr)
clue_df$Post_Mortem <- as.numeric(clue_df$Human_Microarray)
clue_df$Progerin <- as.numeric(clue_df$Progerin)
clue_df$HPC_72hr <- as.numeric(clue_df$HPC_72hr)
clue_df$HPC_144hr <- as.numeric(clue_df$HPC_144hr)
clue_df$HPC_1hr <- as.numeric(clue_df$HPC_1hr)
clue_df$HPC_24hr <- as.numeric(clue_df$HPC_24hr)

clue_mat <- as.matrix(clue_df)

clue_corr <- cor(clue_mat, method = "pearson", use = "pairwise.complete.obs")

library(pheatmap)

CMap_scores_heatmap <- pheatmap(clue_corr,
                             fontsize = 12,
                             display_numbers = T,
                             fontsize_number = 16,
                             fontsize_col = 18,
                             fontsize_row = 18, 
                             angle_col = 45,
                             cluster_rows = T,
                             scale = "none",
                             border_color = "black")


svg("Correlations_CMap_Scores.svg", width=9, height=8)
CMap_scores_heatmap # plot should show decreasing disp with increasing mean
dev.off()


# Get top and bottom N drug hits sorted by high rowSums and low rowVars between 6hr, microarray and progerin

# add unique identifier to summarise on
uid <- 1:nrow(clue_results)
scores_df <-  cbind(clue_results[, -c(14:17)], uid) %>% 
  as_tibble()
scores_df <- scores_df[, -10] # get rid of this column as we will calculate ourselves
# convert tau scores to numeric class
scores_df$Human_Microarray <- as.numeric(scores_df$Human_Microarray)
scores_df$Progerin <- as.numeric(scores_df$Progerin)
scores_df$HPC_6hr <- as.numeric(scores_df$HPC_6hr)

# calculated sum of tau scores and variance of tau scores for each drug  
tau_sum_calc <- scores_df %>% 
  pivot_longer(cols = c("HPC_6hr", "Human_Microarray", "Progerin"), 
               names_to = "dataset", values_to = "tau_score") %>% 
  group_by(uid,id, name) %>% 
  summarise(sum_tau = round(sum(tau_score), digits = 1), variance = round(var(tau_score), digits = 1))

# Add sum and variance to df and filter to get drugs with consistent changes
consistent_scores_tbl <- as_tibble(cbind(scores_df, tau_sum_calc[, 4:5])) %>% 
  filter(HPC_6hr > 0 & Human_Microarray > 0 & Progerin > 0 |
           HPC_6hr < 0 & Human_Microarray < 0 & Progerin < 0)
# There are 1,296 drugs with consistent directional tau scores

write_tsv(consistent_scores_tbl, "consistent_direction_Tau_scores.txt")

# create top and bottom 100 corr/anti-corr compounds

bottom_100 <- consistent_scores_tbl %>% 
  arrange(sum_tau) %>% 
  dplyr::slice(1:100)

top_100 <- consistent_scores_tbl %>% 
  arrange(desc(sum_tau)) %>% 
  dplyr::slice(1:100)
  
top_bottom_100 <- rbind(bottom_100,top_100)

write_tsv(top_bottom_100, "TopBottom100_compounds_sortedBy_sumTau.txt")

# Save compounds with tau scores <|> +-90. These are TOP HITS
top_consistent_CMap_hits <- consistent_scores_tbl %>% 
  filter(HPC_6hr>90 & Human_Microarray>90 & Progerin>90|
                 HPC_6hr< -90 & Human_Microarray< -90 & Progerin< -90|
                 HPC_6hr>90 & Human_Microarray>90 & Progerin< -90|
                 HPC_6hr>90 & Human_Microarray< -90 & Progerin>90|
                 HPC_6hr< -90 & Human_Microarray>90 & Progerin>90|
                 HPC_6hr>90 & Human_Microarray< -90 & Progerin< -90|
                 HPC_6hr< -90 & Human_Microarray< -90 & Progerin>90|
                 HPC_6hr< -90 & Human_Microarray>90 & Progerin< -90)

write_tsv(top_consistent_CMap_hits, "top_compounds_over_90corr_90anticorr.txt")

tau75_consistent_CMap_hits <- consistent_scores_tbl %>% 
  filter(HPC_6hr>75 & Human_Microarray>75 & Progerin>75|
           HPC_6hr< -75 & Human_Microarray< -75 & Progerin< -75|
           HPC_6hr>75 & Human_Microarray>75 & Progerin< -75|
           HPC_6hr>75 & Human_Microarray< -75 & Progerin>75|
           HPC_6hr< -75 & Human_Microarray>75 & Progerin>75|
           HPC_6hr>75 & Human_Microarray< -75 & Progerin< -75|
           HPC_6hr< -75 & Human_Microarray< -75 & Progerin>75|
           HPC_6hr< -75 & Human_Microarray>75 & Progerin< -75)

write_tsv(tau75_consistent_CMap_hits, "top_compounds_over_75corr_75anticorr.txt")




# some compounds might be inverse scoring for different datasets. Check!

scores_df %>% 
  filter(HPC_6hr>90 & Human_Microarray>90 & Progerin>90|
           HPC_6hr< -90 & Human_Microarray< -90 & Progerin< -90|
           HPC_6hr>90 & Human_Microarray>90 & Progerin< -90|
           HPC_6hr>90 & Human_Microarray< -90 & Progerin>90|
           HPC_6hr< -90 & Human_Microarray>90 & Progerin>90|
           HPC_6hr>90 & Human_Microarray< -90 & Progerin< -90|
           HPC_6hr< -90 & Human_Microarray< -90 & Progerin>90|
           HPC_6hr< -90 & Human_Microarray>90 & Progerin< -90) %>% 
  view()

# not the case, still 16 compounds



-------### PLOTTING DATASET TAU SCORE AGAINST EACH OTHER ###-------------

# LINCS pert information

setwd("~/OneDrive - King's College London/ConnectivityMap/LINCS Drug Information")

drug_metadata <- read_csv("Drugs_metadata.csv")


# Load drugAge drugs and corresponding pubchemCIDs


setwd("~/OneDrive - King's College London/DrugAge")

drugage <- read_tsv("drugAge_pubchem_exchange.txt")
colnames(drugage) <- c("name", "pubchem_CID")



# Map pubchem ids for all BRD pert_id in clue results

common <- intersect(clue_results$id, drug_metadata$pert_id)

drug_metadata <- drug_metadata[match(common, drug_metadata$pert_id),]
clue_results <- clue_results[match(common, clue_results$id),]

# Now both drug metadata and cmap drugs are in the same order and same length
# we cbind the pubchem_id column onto CMap

cmap_pubchem <- cbind(drug_metadata$pubchem_cid, clue_results)

#rename column
colnames(cmap_pubchem)[1] <- "PubChem_CID"

# Check for overlap between CMap hits and DrugAge compounds using common PubChem_CIDs

common <- intersect(cmap_pubchem$PubChem_CID, drugage$pubchem_CID)

cmap_drugage_overlap <- cmap_pubchem[match(common, cmap_pubchem$PubChem_CID),]


# label cmap summary table with column identifying DrugAge Matches
# and label consistent Tau scores > 90 and > 75

idx <- match(cmap_drugage_overlap$PubChem_CID, cmap_pubchem$PubChem_CID)

cmap_pubchem2 <- cmap_pubchem %>% 
  mutate(in_drugage = ifelse(row_number() %in% idx, "In DrugAge", "")) %>%
  mutate(consistent_top_hits = ifelse(HPC_6hr>90 & Post_Mortem >90 & Progerin>90|
                             HPC_6hr< -90 & Post_Mortem< -90 & Progerin< -90|
                             HPC_6hr>90 & Post_Mortem>90 & Progerin< -90|
                             HPC_6hr>90 & Post_Mortem< -90 & Progerin>90|
                             HPC_6hr< -90 & Post_Mortem>90 & Progerin>90|
                             HPC_6hr>90 & Post_Mortem< -90 & Progerin< -90|
                             HPC_6hr< -90 & Post_Mortem< -90 & Progerin>90|
                             HPC_6hr< -90 & Post_Mortem>90 & Progerin< -90,
                           "High (Tau > 90)", ifelse(HPC_6hr>75 & Post_Mortem>75 & Progerin>75|
                                                     HPC_6hr< -75 & Post_Mortem< -75 & Progerin< -75|
                                                     HPC_6hr>75 & Post_Mortem>75 & Progerin< -75|
                                                     HPC_6hr>75 & Post_Mortem< -75 & Progerin>75|
                                                     HPC_6hr< -75 & Post_Mortem>75 & Progerin>75|
                                                     HPC_6hr>75 & Post_Mortem< -75 & Progerin< -75|
                                                     HPC_6hr< -75 & Post_Mortem< -75 & Progerin>75|
                                                     HPC_6hr< -75 & Post_Mortem>75 & Progerin< -75,
                                                   "Moderate (Tau > 75 )", "Low (Tau < 75)")))
  


# calcuate an integrated tau score for each drug falling between 0 and 100
# squaring tau values does two things; makes positive and assigns more weight to higher Tau scores
cmap_results_scaled <- cmap_pubchem2 %>% 
  select(name,description,target, Post_Mortem, Progerin, HPC_6hr) %>% 
  pivot_longer(cols = 4:6, names_to = "dataset", values_to = "tau") %>% 
  group_by(name, description, target) %>% 
  mutate(squared = tau^2, mean_tau_sqrd = mean(abs(squared))) 

# rescale values between 0 and 100  
cmap_results_scaled$mean_tau_sqrd <- rescale(cmap_results_scaled$mean_tau_sqrd, to = c(0,100))

# collapse drugs to one row

cmap_results_scaled <- cmap_results_scaled %>% 
  group_by(name,description,target) %>% 
  summarise(mean_tau_sqrd = min(mean_tau_sqrd)) %>% 
  arrange(desc(mean_tau_sqrd)) %>% 
  mutate(target_fill = ifelse(str_detect(target, "^CDK"), "CDK", ifelse(
    str_detect(target, "^TOP"), "TOP", ifelse(str_detect(target, "^HDAC"), "HDAC",
                                               ifelse(str_detect(target, "^SIRT"), "SIRT", ifelse(
                                                 str_detect(target, "^MAP"), "MAPK", ifelse(
                                                   str_detect(target, "^GSK"), "GSK", "Other")))))))

setwd("~/OneDrive - King's College London/ConnectivityMap/CLUE")

write_tsv(cmap_results_scaled, "absolute_standardised_connectivity_scores.txt")


# create bar graph of top50 drugs sorted by absolute scaled mean Tau score

library(ggplot2)

top50_drugs <- ggplot(cmap_results_scaled[1:50,], aes(reorder(name,mean_tau_sqrd), mean_tau_sqrd)) +
  geom_col(aes(fill = target_fill), show.legend = T, col = "black",size=0.3) +
  coord_flip() +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = NULL,y= "Absolute Standardised Connectivity Score (Tau)") +
  scale_fill_brewer(palette = "Paired", name = "Target", na.value = "darkgrey") +
  theme_minimal() +
  theme(text = element_text(color = "black"),
        axis.title.x = element_text(size = 15, vjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.background = element_rect(color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.ticks.y = element_line(color="black"))

setwd("~/OneDrive - King's College London/ConnectivityMap/CLUE")

ggsave("top50_highest_Tau_drug_hits.svg", plot = top50_drugs, 
       device = "svg", dpi = 600, width = 9, height = 10)  

# create bar graph of 13 top drugs hits (><90 tau)

top_hits_over_90tau <- cmap_pubchem2 %>% 
  filter(HPC_6hr > 90 & Post_Mortem > 90 & Progerin > 90 |
           HPC_6hr < -90 & Post_Mortem < -90 & Progerin < -90) %>% 
  mutate(target_fill = ifelse(str_detect(target, "^CDK"), "CDK", ifelse(
    str_detect(target, "^TOP"), "TOP", ifelse(str_detect(target, "^HDAC"), "HDAC",
                                              ifelse(str_detect(target, "^SIRT"), "SIRT", ifelse(
                                                str_detect(target, "^MAP"), "MAPK", ifelse(
                                                  str_detect(target, "^GSK"), "GSK", "Other"))))))) 
  

top_hits_over_75tau <- cmap_pubchem2 %>% 
  filter(HPC_6hr > 75 & Post_Mortem > 75 & Progerin > 75 |
           HPC_6hr < -75 & Post_Mortem < -75 & Progerin < -75) %>% 
  mutate(target_fill = ifelse(str_detect(target, "^CDK"), "CDK", ifelse(
    str_detect(target, "^TOP"), "TOP", ifelse(str_detect(target, "^HDAC"), "HDAC",
                                              ifelse(str_detect(target, "^SIRT"), "SIRT", ifelse(
                                                str_detect(target, "^MAP"), "MAPK", ifelse(
                                                  str_detect(target, "^GSK"), "GSK", "Other")))))))



write_tsv(top_hits_over_90tau, "top_compounds_over_90corr_90anticorr.txt")
write_tsv(top_hits_over_75tau, "top_compounds_over_75corr_75anticorr.txt")

# calculate mean 

over90_plot_data <- top_hits_over_90tau %>% 
  select(name,description,target, Post_Mortem, Progerin, HPC_6hr) %>% 
  pivot_longer(cols = 4:6, names_to = "dataset", values_to = "tau") %>% 
  group_by(name, description, target) %>% 
  mutate(mean_tau = mean(tau)) %>% 
  pivot_wider(names_from = "dataset", values_from = "tau") %>% 
  arrange(desc(mean_tau))

over75_plot_data <- top_hits_over_75tau %>% 
  select(name,description,target, Post_Mortem, Progerin, HPC_6hr, target_fill) %>% 
  pivot_longer(cols = 4:6, names_to = "dataset", values_to = "tau") %>% 
  group_by(name, description, target, target_fill) %>% 
  mutate(mean_tau = mean(tau)) %>% 
  group_by(name, description, target, target_fill) %>% 
  summarise(mean_tau = min(mean_tau))

# barplots

top_hits_90_barplot <- ggplot(over90_plot_data, aes(reorder(name,abs(mean_tau)) , abs(mean_tau))) +
  geom_col(aes(fill = target), show.legend = T, col = "black",size=0.3) +
  coord_flip() +
  scale_y_continuous(expand = c(0,0), breaks = c(0,25,50,75,100),
                     limits = c(0,100)) +
  labs(x = NULL,y= "Mean Absolute Tau Score") +
  scale_fill_brewer(palette = "Paired", name = "Target", na.value = "darkgrey") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 15, vjust = 0.5, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.background = element_rect(color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.ticks.y = element_line(color="black"),
        axis.ticks.x = element_line(color="black"))

ggsave("top_hits_over90_tau_barplot.svg", plot = top_hits_90_barplot,
       device = "svg", dpi = 600, width = 15, height = 5)

top_hits_75_barplot <- ggplot(over75_plot_data, aes(reorder(name,abs(mean_tau)) , abs(mean_tau))) +
  geom_col(aes(fill = target_fill), show.legend = T, col = "black",size=0.3) +
  coord_flip() +
  scale_y_continuous(expand = c(0,0), breaks = c(0,25,50,75,100),
                     limits = c(0,100)) +
  labs(x = NULL,y= "Mean Absolute Tau Score") +
  scale_fill_brewer(palette = "Paired", name = "Target", na.value = "darkgrey") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 15, vjust = 0.5, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.background = element_rect(color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.ticks.y = element_line(color="black"),
        axis.ticks.x = element_line(color="black"))

ggsave("top_hits_over75_tau_barplot.svg", plot = top_hits_75_barplot,
       device = "svg", dpi = 600, width = 10, height = 7)

# top 50 drugs absolute mean tau

# get more colors for r color brewer

colourCount <- length(unique(plot_data$target_fill))
getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# make plotting data
plot_data <- cmap_pubchem2 %>% 
  select(name, description, target, Post_Mortem,HPC_6hr,Progerin, in_drugage) %>% 
  pivot_longer(cols = 4:6, names_to = "dataset", values_to = "tau") %>%
  group_by(name) %>% 
  summarise(abs_mean_tau = abs(mean(tau)), description = min(description),
            target = min(target), in_drugage = min(in_drugage)) %>% 
  ungroup() %>% 
  arrange(desc(abs_mean_tau)) %>% 
  mutate(target_fill = ifelse(str_detect(target, "^CDK"), "CDK", ifelse(
    str_detect(target, "^TOP2A"), "TOP2A", ifelse(str_detect(target, "^HDAC"), "HDAC",
                                              ifelse(str_detect(target, "^SIRT"), "SIRT", ifelse(
                                                str_detect(target, "^MAP"), "MAPK", ifelse(
                                                  str_detect(target, "^GSK"), "GSK3B", ifelse(
                                                    str_detect(target, "^AKT"), "AKT", ifelse(
                                                      str_detect(target, "^MTOR"), "MTOR", ifelse(
                                                        str_detect(target, "^PPAR"), "PPAR", ifelse(
                                                          str_detect(target, "^PDGF"), "PDGFR", ifelse(
                                                            str_detect(target, "TUB"), "TUBB", ifelse(
                                                              str_detect(target, "^TOP1"), "TOP1", "Other"))))))))))))) %>%  
  dplyr::slice(1:50)

# construct top 50 barplot

top50_meantau_barplot <- ggplot(plot_data, aes(reorder(name,abs_mean_tau), abs_mean_tau)) +
  geom_col(aes(fill = target_fill), show.legend = T, col = "black",size=0.3) +
  coord_flip() +
  scale_y_continuous(expand = c(0,0), breaks = c(0,25,50,75,100),
                     limits = c(0,100)) +
  labs(x = NULL,y= "Mean Absolute Tau Score") +
  scale_fill_manual(values = getPalette(colourCount), name = "Target") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 26, vjust = 0.5, colour = "black"),
        axis.text.x = element_text(size = 20, colour = "black"),
        axis.text.y = element_text(size = 24, colour = "black"),
        legend.background = element_rect(color = "black"),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        axis.ticks.y = element_line(color="black"),
        axis.ticks.x = element_line(color="black"))
  
ggsave("top_50_meanTau_barplot.svg", plot = top50_meantau_barplot,
       device = "svg", dpi = 600, width = 12, height = 15)



# Pairwise Scatter Plots of Samples'  Tau scores

# Create column to map size onto points
size_mapping <- cmap_pubchem2 %>% 
  pivot_longer(cols = 12:14, names_to = "dataset", values_to = "tau") %>% 
  group_by(id,name,pc,ts_pc, median_tau_score) %>% 
  summarise(squared = tau^2, mean_tau_sqrd = mean(abs(squared))) %>% 
  ungroup() %>% 
  group_by(id,name,pc,ts_pc, median_tau_score) %>% 
  summarise(mean_tau_sqrd = min(mean_tau_sqrd))

common <- intersect(size_mapping$id, cmap_pubchem2$id)
size_mapping <- size_mapping[match(common, size_mapping$id),]
cmap_pubchem2 <- cmap_pubchem2[match(common, cmap_pubchem2$id),]

# bind absolute mean Tau column to cmap results
gg_cmap_pubchem2 <- cbind(cmap_pubchem2, size_mapping[, "mean_tau_sqrd"])

# override drugAge column with actual drug names for ggplot labels

idx2 <- str_which(gg_cmap_pubchem2$in_drugage, "In DrugAge")
gg_cmap_pubchem2$in_drugage[idx2] <- gg_cmap_pubchem2$name[idx2]

# Add extra column for drugAge drugs with high absolute mean Tau score > 3000

temp <- gg_cmap_pubchem2 %>% 
  filter(str_detect(in_drugage, "[[:alpha:]]") & 
           mean_tau_sqrd > 3000)

                           
                                                  
# reorder camp hits to same as labelling idex
idx4 <- match(temp$id, gg_cmap_pubchem2$id)

gg_cmap_pubchem2[idx4,]

# create blank column then fill with labels

gg_cmap_pubchem2$cmap_high_drugage <- ""
gg_cmap_pubchem2$cmap_high_drugage[idx4] <- gg_cmap_pubchem2$name[idx4]



# for human microarry, hpc6hr and Progerin pairwise comparisons
library(ggrepel)

# Post Mortem Microarry vs Progerin

PMvProgerin <- ggplot(gg_cmap_pubchem2, aes(Progerin, Human_Microarray)) +
  geom_point(aes(size = mean_tau_sqrd,col = consistent_top_hits),alpha=0.4) +
  guides(alpha = FALSE) +
  labs(x="Progerin CMap Score (Tau)",y="Post-Mortem CMap Score (Tau)",
       title = "Progerin vs Post-Mortem Hippocampus") +
  scale_size_area(breaks = c(1,5000,7000),max_size = 1.8,
                  limits = c(0,9600)) +
  scale_color_manual(values = c("red", "black","blue"),
                     name = "Consistency in Absolute Scores Between 
Post-Mortem, Progerin and HPC Signatures") +
  geom_smooth(method = "lm", col = "black", size = 1, alpha=0.4) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_label_repel(aes(label = cmap_high_drugage), size = 3.5, label.size = 0.2,
                   point.padding = 0.1, label.padding = 0.15,force = 2,
                   nudge_x = 0.3,nudge_y = 0.3) +
  guides(col = guide_legend(override.aes = list(size=8)),
         size = FALSE) +
  theme_light() +
  theme(legend.text = element_text(size=12),
        plot.title = element_text(size=30, hjust = 0.5,vjust = 1.5),
        axis.title = element_text(size=22),
        axis.text = element_text(size=16, color = "black"),
        axis.ticks = element_line(color = "black", size = 1.5),
        legend.title = element_text(size = 13),
        panel.grid = element_blank(),
        panel.border = element_rect(size=1.5,color="black")) +
  annotate(geom="label", x=10, y=-10, label="r = 0.50, P < 2.2e-16",
           color="red", size = 5)

# Progerin vs HPC

ProgerinvHPC6hr <- ggplot(gg_cmap_pubchem2, aes(Progerin, HPC_6hr)) +
  geom_point(aes(size = mean_tau_sqrd,col = consistent_top_hits),alpha=0.4) +
  guides(alpha = FALSE) +
  labs(x="Progerin CMap Score (Tau)",y="HPC 6hr CMap Score (Tau)",
       title = "Progerin vs HPC 6hr") +
  scale_size_area(breaks = c(1,5000,7000),max_size = 1.8,
                  limits = c(0,9600)) +
  scale_color_manual(values = c("red", "black","blue"),
                     name = "Consistency in Absolute Scores Between 
Post-Mortem, Progerin and HPC Signatures") +
  geom_smooth(method = "lm", col = "black", size = 1, alpha=0.4) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_label_repel(aes(label = cmap_high_drugage), size = 3.5, label.size = 0.2,
                   point.padding = 0.2, label.padding = 0.25,force = 4,
                   nudge_x = 0.3,nudge_y = 0.3) +
  guides(col = guide_legend(override.aes = list(size=8)),
         size = FALSE) +
  theme_light() +
  theme(legend.text = element_text(size=12),
        plot.title = element_text(size=30, hjust = 0.5),
        axis.title = element_text(size=22),
        axis.text = element_text(size=16, color = "black"),
        axis.ticks = element_line(color = "black", size = 1.5),
        legend.title = element_text(size = 13),
        panel.grid = element_blank(),
        panel.border = element_rect(size=1.5,color="black")) +
  annotate(geom="label", x=10, y=-10, label="r = 0.35, P < 2.2e-16",
           color="red", size = 5)

# HPC 6hr vs Post Mort Microarray

PMvHPC6hr <- ggplot(gg_cmap_pubchem2, aes(Human_Microarray, HPC_6hr)) +
  geom_point(aes(size = mean_tau_sqrd,col = consistent_top_hits),alpha=0.4) +
  guides(alpha = FALSE) +
  labs(x="Post-Mortem CMap Score (Tau)",y="HPC 6hr CMap Score (Tau)",
       title = "Post-Mortem Hippocampus vs HPC 6hr") +
  scale_size_area(breaks = c(1,5000,7000),max_size = 1.8,
                  limits = c(0,9600)) +
  scale_color_manual(values = c("red", "black","blue"),
                     name = "Consistency in Absolute Scores Between 
Post-Mortem, Progerin and HPC Signatures") +
  geom_smooth(method = "lm", col = "black", size = 1, alpha=0.4) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_label_repel(aes(label = cmap_high_drugage), size = 3.5, label.size = 0.2,
                   point.padding = 0.2, label.padding = 0.25,force = 4,
                   nudge_x = 0.1,nudge_y = 0.3) +
  guides(col = guide_legend(override.aes = list(size=8)),
         size = FALSE) +
  theme_light() +
  theme(legend.text = element_text(size=12),
        plot.title = element_text(size=30, hjust = 0.5),
        axis.title = element_text(size=22),
        axis.text = element_text(size=16, color = "black"),
        axis.ticks = element_line(color = "black", size = 1.5),
        legend.title = element_text(size = 13),
        panel.grid = element_blank(),
        panel.border = element_rect(size=1.5,color="black")) +
  annotate(geom="label", x=10, y=-10, label="r = 0.28, P < 2.2e-16",
           color="red", size = 5)

# Add pearson coefficient and p value to graph
cor.test(gg_cmap_pubchem2$HPC_6hr, gg_cmap_pubchem2$Human_Microarray, method = "pearson")

setwd("~/OneDrive - King's College London/ConnectivityMap/CLUE")  
ggsave("scatter_postMort_vs_Progerin_Tau.svg", plot = PMvProgerin,
       device = svg,dpi = 600,width = 15,height = 8)
ggsave("scatter_Progerin_vs_HPC6hr_Tau.svg", plot = ProgerinvHPC6hr,
       device = svg,dpi = 600,width = 15,height = 8)
ggsave("scatter_postMort_vs_HPC6hr_Tau.svg", plot = PMvHPC6hr,
       device = svg,dpi = 600,width = 15,height = 8)



# for other comparisons
ggplot(gg_cmap_pubchem2, aes(HPC_24hr, Human_Microarray)) +
  geom_point(size=1, alpha=0.3) +
  geom_smooth(method = "lm") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

cor.test(gg_cmap_pubchem2$HPC_24hr, gg_cmap_pubchem2$Human_Microarray, method = "pearson")
