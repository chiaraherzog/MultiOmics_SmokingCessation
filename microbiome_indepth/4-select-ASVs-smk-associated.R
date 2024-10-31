
## Prepare sequence, relabund. and tax data for plot relevant ASVs
# 1. filter based on prevalence (in smokers control/all smoker participants)
## Threshold saliva 90%
## Threshold saliva 50%
# 2. filter prevalent ASVs that also are significantly associated with smoking 
# test relative abundance using unpaired two-sample Wilcoxon test never smokers vs smokers baseline
# write away sequences for phylogenetic tree
# gather also taxonomy info, including species assignment
# gather abundance data
## heatmap 1: relative abundance, never smoker control (n=24), smoker baseline (n=24)
## heatmap 2: relative abundance, change from baseline gather M6 smoker no cessation (n=12), M6 smoker cessation (n=12)
# Indicate which ASVs also significantly modulated by smoking cessation (or not) at M6
## high compliance linear model (time:compliance in interaction model generally higher p-values)
# Author: Charlotte Vavourakis

library(tidyverse)
library(ampvis2)
library(Biostrings)
library(MultiAssayExperiment)

dir.create("./4-output")

# Saliva #-----

## get relevant samples, tax info, filter prevalence #-----

# select samples
pheno_saliva_extend <- readRDS("./1-output/pheno_saliva_S_ctrl.Rds")

# tax info, select ASVs
ASVtable_saliva_extend <- readRDS("./1-output/ASVtable_saliva_S_ctrl.Rds")
taxinfo_saliva <- ASVtable_saliva_extend %>%
  dplyr::select(OTU,Phylum,Family,Genus,Species)

# abundance
exp <- "Saliva microbiome: ASVs"
load("../../data/data_raw.Rdata")

abund_saliva <- as.data.frame(longFormat(data[,,exp])) %>%
  dplyr::mutate(sampleId = paste0(primary,"SM")) %>%
  dplyr::filter(sampleId %in% unique(c(pheno_saliva_extend$sampleId))) %>%
  pivot_wider(id_cols = rowname, names_from = sampleId, values_from = value) %>%
  column_to_rownames(var="rowname")

# filter prevalence
prevalence <- rowSums(abund_saliva != 0) / ncol(abund_saliva)
abund_saliva <- abund_saliva[prevalence >= 0.9, ] #61 prevalent ASVs

## smoking Wilcoxon-tests #-----

# Group1: Never smokers, group2: smokers baseline
group1 <- pheno_saliva_extend %>%
  filter(group == "Never smoker (control)") %>%
  pull(sampleId)
group1 <- abund_saliva %>% dplyr::select(group1)

group2 <- pheno_saliva_extend %>%
  filter(group == "Smoker baseline") %>%
  pull(sampleId)
group2 <- abund_saliva %>% dplyr::select(group2)

identical(rownames(group1),rownames(group2))

result_Wilcox <- data.frame(
  feature = rep(NA, nrow(group1)),
  pval = rep(NA, nrow(group1)),
  median_group1 = rep(NA, nrow(group1)),
  median_group2 = rep(NA, nrow(group1))
)

for (i in 1:nrow(group1)) {
  result_Wilcox[i,1] <- rownames(group1)[i]
  result_Wilcox[i,2] <- wilcox.test(as.numeric(group1[i,]), as.numeric(group2[i,]), alternative = "two.sided")$p.value
  result_Wilcox[i,3] <- median(as.numeric(group1[i,]))
  result_Wilcox[i,4] <- median(as.numeric(group2[i,]))
}

result_Wilcox <- result_Wilcox %>%
  mutate(association_smoking = case_when(median_group1 > median_group2 ~ "-",
                                         median_group2 > median_group1 ~ "+",
                                         TRUE ~ ""))
smoking_feat <- result_Wilcox %>%
  dplyr::filter(pval <= 0.05) #19 features
  
## update taxinfo #-----
taxinfo_saliva <- taxinfo_saliva %>%
  dplyr::filter(OTU %in% c(smoking_feat$feature))

## gather and write away all seqs for phylogenetic tree #----
fasta <- readDNAStringSet("./2-output/saliva_S_ctrl_ASV.fasta")
fasta <- fasta[c(taxinfo_saliva$OTU)]
writeXStringSet(fasta, file = "./4-output/saliva_ASV.fasta", format = "fasta")

## ASV table and metadata relative abundance heatmap 1 #-----
metadat1 <- pheno_saliva_extend %>%
  dplyr::filter(group %in% c("Never smoker (control)","Smoker baseline")) %>%
  dplyr::select(sampleId,group, age_at_consent, bmi_at_consent, smoking_py, nic_replacement)

abund_saliva <- abund_saliva %>% 
  dplyr::filter(rownames(abund_saliva) %in% c(smoking_feat$feature)) %>%
  dplyr::select(all_of(c(metadat1$sampleId)))

## ASV table change from baseline heatmap 2 #-----

metadat2 <- pheno_saliva_extend %>%
  dplyr::filter(group %in% c("M6 no smoking cessation","M6 smoking cessation")) %>%
  dplyr::select(sampleId,group, age_at_consent, bmi_at_consent, smoking_py, nic_replacement)

load("../../data/data_baseline_change.Rdata")

abund_baseline_saliva <- as.data.frame(longFormat(data[,,exp])) %>%
  dplyr::mutate(sampleId = paste0(primary,"SM")) %>%
  dplyr::filter(sampleId %in% c(metadat2$sampleId)) %>%
  pivot_wider(id_cols = rowname, names_from = sampleId, values_from = value) %>%
  dplyr::filter(rowname %in% c(smoking_feat$feature)) %>%
  column_to_rownames(var="rowname")

## p-values from lmm models #-----

# # interaciton model
# exp <- "Saliva microbiome: ASVs"
# load("../out_lmm_factor.Rdata")
# dat1 <- as.data.frame(out_lmm$`Interaction model with packyears`) %>%
#   filter(str_detect(x, exp))
# dat1$x <- gsub(paste0(exp,"_"),"",dat1$x)
# dat1 <- dat1 %>% dplyr::filter(x %in% c(smoking_feat$feature))
# 
# lmm_pval <- data.frame(
#   ASV = dat1$x,
#   `time_compliance_interaction` = dat1$`p.value_visitIdM6:compliancehigher compliance`)

# high compliacne model
exp <- "Saliva microbiome: ASVs"
load("../out_lmm_factor.Rdata")
dat1 <- as.data.frame(out_lmm$`Higher compliance only with packyears`) %>%
  filter(str_detect(x, exp))
dat1$x = gsub(paste0(exp,"_"),"",dat1$x)
dat1 <- dat1 %>% dplyr::filter(x %in% c(smoking_feat$feature))

lmm_pval <- data.frame(
  ASV = dat1$x,
  `association_cessation` = gtools::stars.pval(dat1$`p.value_visitIdM6`)) 

## save #-----
list_saliva <- list(taxinfo_saliva,abund_saliva,metadat1,abund_baseline_saliva,metadat2,lmm_pval, smoking_feat)
saveRDS(list_saliva,file="./4-output/list_saliva.Rds")

# Gut #-----

## get relevant samples, tax info, filter prevalence #-----

# select samples
pheno_stool_extend <- readRDS("./1-output/pheno_stool_S_ctrl.Rds")

# tax info, select ASVs
ASVtable_stool_extend <- readRDS("./1-output/ASVtable_stool_S_ctrl.Rds")
taxinfo_stool <- ASVtable_stool_extend %>%
  dplyr::select(OTU,Phylum,Family,Genus,Species)

# abundance
exp <- "Stool microbiome: ASVs"
load("../../data/data_raw.Rdata")

abund_stool <- as.data.frame(longFormat(data[,,exp])) %>%
  dplyr::mutate(sampleId = paste0(primary,"FM")) %>%
  dplyr::filter(sampleId %in% unique(c(pheno_stool_extend$sampleId))) %>%
  pivot_wider(id_cols = rowname, names_from = sampleId, values_from = value) %>%
  column_to_rownames(var="rowname")

# filter prevalence
prevalence <- rowSums(abund_stool != 0) / ncol(abund_stool)
abund_stool <- abund_stool[prevalence >= 0.5, ] #149 prevalent ASVs

## smoking Wilcoxon-tests #-----

# Group1: Never smokers, group2: smokers baseline
group1 <- pheno_stool_extend %>%
  filter(group == "Never smoker (control)") %>%
  pull(sampleId)
group1 <- abund_stool %>% dplyr::select(group1)

group2 <- pheno_stool_extend %>%
  filter(group == "Smoker baseline") %>%
  pull(sampleId)
group2 <- abund_stool %>% dplyr::select(group2)

identical(rownames(group1),rownames(group2))

result_Wilcox <- data.frame(
  feature = rep(NA, nrow(group1)),
  pval = rep(NA, nrow(group1)),
  median_group1 = rep(NA, nrow(group1)),
  median_group2 = rep(NA, nrow(group1))
)

for (i in 1:nrow(group1)) {
  result_Wilcox[i,1] <- rownames(group1)[i]
  result_Wilcox[i,2] <- wilcox.test(as.numeric(group1[i,]), as.numeric(group2[i,]), alternative = "two.sided")$p.value
  result_Wilcox[i,3] <- median(as.numeric(group1[i,]))
  result_Wilcox[i,4] <- median(as.numeric(group2[i,]))
}

result_Wilcox <- result_Wilcox %>%
  mutate(association_smoking = case_when(median_group1 > median_group2 ~ "-",
                                         median_group2 > median_group1 ~ "+",
                                         TRUE ~ ""))
smoking_feat <- result_Wilcox %>%
  dplyr::filter(pval <= 0.05) # 16 features

## update taxinfo #-----
taxinfo_stool <- taxinfo_stool %>%
  dplyr::filter(OTU %in% c(smoking_feat$feature))

## gather and write away all seqs for phylogenetic tree #----
fasta <- readDNAStringSet("./2-output/stool_S_ctrl_ASV.fasta")
fasta <- fasta[c(taxinfo_stool$OTU)]
writeXStringSet(fasta, file = "./4-output/stool_ASV.fasta", format = "fasta")

## ASV table and metadata relative abundance heatmap 1 #-----
metadat1 <- pheno_stool_extend %>%
  dplyr::filter(group %in% c("Never smoker (control)","Smoker baseline")) %>%
  dplyr::select(sampleId,group, age_at_consent, bmi_at_consent, smoking_py, nic_replacement)

abund_stool <- abund_stool %>% 
  dplyr::filter(rownames(abund_stool) %in% c(smoking_feat$feature)) %>%
  dplyr::select(all_of(c(metadat1$sampleId)))

## ASV table change from baseline heatmap 2 #-----

metadat2 <- pheno_stool_extend %>%
  dplyr::filter(group %in% c("M6 no smoking cessation","M6 smoking cessation")) %>%
  dplyr::select(sampleId,group, age_at_consent, bmi_at_consent, smoking_py, nic_replacement)

load("../../data/data_baseline_change.Rdata")

abund_baseline_stool <- as.data.frame(longFormat(data[,,exp])) %>%
  dplyr::mutate(sampleId = paste0(primary,"FM")) %>%
  dplyr::filter(sampleId %in% c(metadat2$sampleId)) %>%
  pivot_wider(id_cols = rowname, names_from = sampleId, values_from = value) %>%
  dplyr::filter(rowname %in% c(smoking_feat$feature)) %>%
  column_to_rownames(var="rowname")

## p-values from lmm models #-----

# high compliacne model
exp <- "Stool microbiome: ASVs"
load("../out_lmm_factor.Rdata")
dat1 <- as.data.frame(out_lmm$`Higher compliance only with packyears`) %>%
  filter(str_detect(x, exp))
dat1$x = gsub(paste0(exp,"_"),"",dat1$x)
dat1 <- dat1 %>% dplyr::filter(x %in% c(smoking_feat$feature))

lmm_pval <- data.frame(
  ASV = dat1$x,
  `association_cessation` = gtools::stars.pval(dat1$`p.value_visitIdM6`)) 

## save #-----
list_stool <- list(taxinfo_stool,abund_stool,metadat1,abund_baseline_stool,metadat2,lmm_pval, smoking_feat)
saveRDS(list_stool,file="./4-output/list_stool.Rds")
