
----------# FILTER POST-MORTEM DATA TO LINCS L1000 GENES ----------


library(tidyverse)

setwd("~/OneDrive - King's College London/Data/Bioinformatic/IPA/Data")

# Load post_mortem microarray combined signature


ber_kang <- read_tsv("Berch_Kang_Combined_Signature.txt")


# filter gene lists to only include BING genes from CMap 2.0

setwd("~/OneDrive - King's College London/ConnectivityMap")

clue_genes <- read_tsv("CLUE/clue.io_all_genes.txt") %>% 
  filter(Type %in% c("best inferred", "landmark"))

# find matching entrez between clue genes and geneset entrez symbols
# match() function finds the indices of those string matches. Using that index vector we can than subset the
# geneset up set to only include rows with common entrez

common <- intersect(clue_genes$`Entrez ID`, ber_kang$entrez)
ber_kang <- ber_kang[match(common, ber_kang$entrez), ]

# There were 13893 common entrez before filtering. After, there are 8809 entrez


----------# FILTER PROGERIN SIGNATURE TO LINCS L1000 GENES ----------


setwd("~/OneDrive - King's College London/Data/Bioinformatic/External/RNA-seq/Progerin_signatures/Miller_GSE524a31/GSE524a31_Results")

progerin <- read_tsv("progerin_signature.txt")

common <- intersect(clue_genes$`Entrez ID`, progerin$entrez)
progerin <- progerin[match(common, progerin$entrez), ]


----------# FILTER SERUM-TREATED HIPPOCAMPAL PROGENITOR SIGNATURE TO LINCS L1000 GENES ----------


## Read in data ##

setwd("~/OneDrive - King's College London/Dan Array Data/Data for Limma/TYTUS/TopDE genes/TopDEGs/bycoeff_contrast/probeSet_summarised_before_LmFit")

OldvYoungTP1 <- read_tsv("hippocampal_progenitors_1hr.txt", na = "NA")
OldvYoungTP6 <- read_tsv("hippocampal_progenitors_6hr.txt", na = "NA")
OldvYoungTP24 <- read_tsv("hippocampal_progenitors_24hr.txt", na = "NA")
OldvYoungTP72 <- read_tsv("hippocampal_progenitors_72hr.txt", na = "NA")
OldvYoungTP144 <- read_tsv("hippocampal_progenitors_144hr.txt", na = "NA")

common1 <- intersect(OldvYoungTP1$ENTREZ_GENE_ID, clue_genes$`Entrez ID`)
OldvYoungTP1 <- OldvYoungTP1[match(common1,OldvYoungTP1$ENTREZ_GENE_ID),]

common6 <- intersect(OldvYoungTP6$ENTREZ_GENE_ID, clue_genes$`Entrez ID`)
OldvYoungTP6 <- OldvYoungTP6[match(common6,OldvYoungTP6$ENTREZ_GENE_ID),]

common24 <- intersect(OldvYoungTP24$ENTREZ_GENE_ID, clue_genes$`Entrez ID`)
OldvYoungTP24 <- OldvYoungTP24[match(common24,OldvYoungTP24$ENTREZ_GENE_ID),]

common72 <- intersect(OldvYoungTP72$ENTREZ_GENE_ID, clue_genes$`Entrez ID`)
OldvYoungTP72 <- OldvYoungTP72[match(common72,OldvYoungTP72$ENTREZ_GENE_ID),]

common144 <- intersect(OldvYoungTP144$ENTREZ_GENE_ID, clue_genes$`Entrez ID`)
OldvYoungTP144 <- OldvYoungTP144[match(common144,OldvYoungTP144$ENTREZ_GENE_ID),]


-----------# SELECTING GENES FOR CONNECTIVITY MAPPING #-----------

# CLUE online tool (https://clue.io/query). Only accepts max 150 genes as entrez ids.

post_mort_top150_UP <- ber_kang %>% 
  filter(meanlogFC > 0) %>%
  arrange(desc(meanlogFC)) %>% 
  dplyr::slice(1:150) %>% 
  dplyr::rename("berch_kang_UP" = entrezgene_id) %>% 
  dplyr::select(berch_kang_UP)

post_mort_top150_DN <- ber_kang %>% 
  filter(meanlogFC < 0) %>%
  arrange(meanlogFC) %>% 
  dplyr::slice(1:150) %>% 
  dplyr::rename("berch_kang_DN" = entrezgene_id) %>% 
  dplyr::select(berch_kang_DN)

progerin_top150_UP <- progerin %>% 
  filter(logFC > 0) %>%
  arrange(desc(logFC)) %>% 
  slice(1:150) %>% 
  dplyr::select(entrez) %>% 
  rename("progerin_UP" = entrez)

progerin_top150_DN <- progerin %>% 
  filter(logFC < 0) %>%
  arrange(logFC) %>% 
  slice(1:150) %>% 
  dplyr::select(entrez) %>% 
  rename("progerin_UP" = entrez)

serum_progenitor_1hr_UP <- OldvYoungTP1 %>% 
  filter(P.Value < 0.05, logFC > 0) %>% 
  slice(1:150) %>% 
  dplyr::select(ENTREZ_GENE_ID) %>% 
  rename("entrez" = ENTREZ_GENE_ID)

serum_progenitor_1hr_DN <- OldvYoungTP1 %>% 
  filter(P.Value < 0.05,logFC < 0) %>% 
  slice(1:150) %>% 
  dplyr::select(ENTREZ_GENE_ID) %>% 
  rename("entrez" = ENTREZ_GENE_ID)

serum_progenitor_6hr_UP <- OldvYoungTP6 %>% 
  filter(P.Value < 0.05,logFC > 0) %>% 
  slice(1:150) %>% 
  dplyr::select(ENTREZ_GENE_ID) %>% 
  rename("entrez" = ENTREZ_GENE_ID)

serum_progenitor_6hr_DN <- OldvYoungTP6 %>% 
  filter(P.Value < 0.05,logFC < 0) %>% 
  slice(1:150) %>% 
  dplyr::select(ENTREZ_GENE_ID) %>% 
  rename("entrez" = ENTREZ_GENE_ID)

serum_progenitor_24hr_UP <- OldvYoungTP24 %>% 
  filter(P.Value < 0.05, logFC > 0) %>% 
  slice(1:150) %>% 
  dplyr::select(ENTREZ_GENE_ID) %>% 
  rename("entrez" = ENTREZ_GENE_ID)

serum_progenitor_24hr_DN <- OldvYoungTP24 %>% 
  filter(P.Value < 0.05, logFC < 0) %>% 
  slice(1:150) %>% 
  dplyr::select(ENTREZ_GENE_ID) %>% 
  rename("entrez" = ENTREZ_GENE_ID)

serum_progenitor_72hr_UP <- OldvYoungTP72 %>% 
  filter(P.Value < 0.05, logFC > 0) %>% 
  slice(1:150) %>% 
  dplyr::select(ENTREZ_GENE_ID) %>% 
  rename("entrez" = ENTREZ_GENE_ID)

serum_progenitor_72hr_DN <- OldvYoungTP72 %>% 
  filter(P.Value < 0.05,logFC < 0) %>% 
  slice(1:150) %>% 
  dplyr::select(ENTREZ_GENE_ID) %>% 
  rename("entrez" = ENTREZ_GENE_ID)

serum_progenitor_144hr_UP <- OldvYoungTP144 %>% 
  filter(P.Value < 0.05,logFC > 0) %>% 
  slice(1:150) %>% 
  dplyr::select(ENTREZ_GENE_ID) %>% 
  rename("entrez" = ENTREZ_GENE_ID)

serum_progenitor_144hr_DN <- OldvYoungTP144 %>% 
  filter(P.Value < 0.05, logFC < 0) %>% 
  slice(1:150) %>% 
  dplyr::select(ENTREZ_GENE_ID) %>% 
  rename("entrez" = ENTREZ_GENE_ID)


# iLINCS online tool (no gene limit & accepts logFCs, http://www.ilincs.org/ilincs/signatures/main/upload)

post_mortem_ilincs <- ber_kang %>% 
  dplyr::select(entrezgene_id, meanlogFC)

serum1hr_ilincs <- OldvYoungTP1[, c("ENTREZ_GENE_ID", "logFC", "P.Value")]
serum6hr_ilincs <- OldvYoungTP6[, c("ENTREZ_GENE_ID", "logFC", "P.Value")]
serum24hr_ilincs <- OldvYoungTP24[, c("ENTREZ_GENE_ID", "logFC", "P.Value")]
serum72hr_ilincs <- OldvYoungTP72[, c("ENTREZ_GENE_ID", "logFC", "P.Value")]
serum144hr_ilincs <- OldvYoungTP144[, c("ENTREZ_GENE_ID", "logFC", "P.Value")]

progerin_ilincs <- progerin %>% 
  dplyr::select(entrez, logFC, padj)





