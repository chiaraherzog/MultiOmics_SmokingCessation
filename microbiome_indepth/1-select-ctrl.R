
# select 24 M0 stool samples of non-smokers, with the healthiest bmi
# preprocessing was done on the full microbiome dataset, including both Intermittent Fasting and Smoking Cessation study arms
# prepare ASVtables (needed for full taxonomy) and covariates for datasets with/without control group
## spurious ASVs removed
## only complete cases M0-M6 kept
# prepare corresponding pheno

# Author: Charlotte Vavourakis

library(tidyverse)
library(MultiAssayExperiment)

dir.create('./1-output')

# select samples #----
load("../../data/data_raw.Rdata")
ctrl <- colData(data) %>% 
  as.data.frame() %>%
  filter(interventionId != 'S') %>%
  filter(visitId== "M0" & smoking_ever=="no") %>%
  arrange(bmi_at_consent) %>%
  dplyr::slice(1:24) %>%
  rownames_to_column(var='primary') %>%
  pull(primary)
saveRDS(ctrl, file='./1-output/cltr_primary.Rds')

# make count tables including control group #----
## saliva ##----
ASVtable_saliva <- readRDS("<path-preprocessed-data>/ASVtable_saliva_S.Rds")
smokers <- colnames(ASVtable_saliva)[grep("SM", colnames(ASVtable_saliva))]
cls_get <- c('OTU',smokers,paste0(ctrl,'SM'),'Kingdom','Phylum','Class','Order','Family','Genus','Species')
ASVtable_saliva_extend <- readRDS("<path-preprocessed-data>/ASVtable_saliva_filtered.Rds") %>%
  select(all_of(cls_get)) %>%
  mutate(sum_count = rowSums(select(., ends_with("SM")))) %>% # remove 0 abundance
  filter(sum_count > 0) %>%
  select(-sum_count)
saveRDS(ASVtable_saliva_extend, file = "./1-output/ASVtable_saliva_S_ctrl.Rds")

## stool ##----
ASVtable_stool <- readRDS("<path-preprocessed-data>/ASVtable_stool_S.Rds")
smokers <- colnames(ASVtable_stool)[grep("FM", colnames(ASVtable_stool))]
cls_get <- c('OTU',smokers,paste0(ctrl,'FM'),'Kingdom','Phylum','Class','Order','Family','Genus','Species')
ASVtable_stool_extend <- readRDS("<path-preprocessed-data>/ASVtable_stool_filtered.Rds") %>%
  select(all_of(cls_get)) %>%
  mutate(sum_count = rowSums(select(., ends_with("FM")))) %>% # remove 0 abundance
  filter(sum_count > 0) %>%
  select(-sum_count)
saveRDS(ASVtable_stool_extend, file = "./1-output/ASVtable_stool_S_ctrl.Rds")

# corresponding pheno and count tables without control group #----
features <- c("interventionId", "subjectId", "visitId", "time", "compliance", "smoking_py","nic_replacement","smkstop","time", "mpstatrs", "age_at_consent", "bmi_at_consent", "cig_before", "etohu_curr", "diet", "preg_ever", "ocp_curr", "hrt_curr", "intactcurr")
load("../../data/data_raw.Rdata")
metadat <- colData(data) %>% 
  as.data.frame() %>%
  dplyr::select(all_of(features)) %>%
  rownames_to_column(var="sampleId")

## Smoking only  ##----
ASVtable_saliva <- readRDS("<path-preprocessed-data>/ASVtable_saliva_S.Rds")
pheno_saliva <- metadat %>%
  mutate(sampleId = paste0(sampleId,"SM")) %>%
  filter(sampleId %in% colnames(ASVtable_saliva)[str_detect(colnames(ASVtable_saliva), "SM")]) %>%
  droplevels()
saveRDS(ASVtable_saliva, file = "./1-output/ASVtable_saliva_S.Rds")
saveRDS(pheno_saliva, file = "./1-output/pheno_saliva_S.Rds")

ASVtable_stool <- readRDS("<path-preprocessed-data>/ASVtable_stool_S.Rds")
pheno_stool <- metadat %>%
  mutate(sampleId = paste0(sampleId,"FM")) %>%
  filter(sampleId %in% colnames(ASVtable_stool)[str_detect(colnames(ASVtable_stool), "FM")]) %>%
  droplevels()
saveRDS(ASVtable_stool, file = "./1-output/ASVtable_stool_S.Rds")
saveRDS(pheno_stool, file = "./1-output/pheno_stool_S.Rds")

## Pheno control group  ##----
#ASVtable_saliva_extend <- readRDS("./1-output/ASVtable_saliva_S_ctrl.Rds")

pheno_saliva_extend <- metadat %>%
  mutate(sampleId = paste0(sampleId,"SM")) %>%
  filter(sampleId %in% colnames(ASVtable_saliva_extend)[str_detect(colnames(ASVtable_saliva_extend), "SM")]) %>%
  droplevels() %>%
  mutate(visitId = case_when(interventionId == 'S' ~ visitId, TRUE ~ 'Never smoker (control)')) %>% # set visitId for ctrl group to control
  mutate(compliance = case_when(interventionId == 'S' ~ compliance, TRUE ~ 'Never smoker (control)')) %>% # set compliance for ctrl group to control
  mutate(group=case_when(
    visitId == "Never smoker (control)" ~ visitId,
    visitId == "M0" ~ "Smoker baseline",
    smkstop=="yes" & !visitId %in% c("Never smoker (control)", "M0")  ~ paste0(visitId," smoking cessation"),
    TRUE ~ paste0(visitId," no smoking cessation")
  ))# add group variable for testing
saveRDS(pheno_saliva_extend, file = "./1-output/pheno_saliva_S_ctrl.Rds")

#ASVtable_stool_extend <- readRDS("./1-output/ASVtable_stool_S_ctrl.Rds")

pheno_stool_extend <- metadat %>%
  mutate(sampleId = paste0(sampleId,"FM")) %>%
  filter(sampleId %in% colnames(ASVtable_stool_extend)[str_detect(colnames(ASVtable_stool_extend), "FM")]) %>%
  droplevels() %>%
  mutate(visitId = case_when(interventionId == 'S' ~ visitId, TRUE ~ 'Never smoker (control)')) %>% # set visitId for ctrl group to control
  mutate(compliance = case_when(interventionId == 'S' ~ compliance, TRUE ~ 'Never smoker (control)')) %>% # remove compliance data for control group
  mutate(group=case_when(
    visitId == "Never smoker (control)" ~ visitId,
    visitId == "M0" ~ "Smoker baseline",
    smkstop=="yes" & !visitId %in% c("Never smoker (control)", "M0")  ~ paste0(visitId," smoking cessation"),
    TRUE ~ paste0(visitId," no smoking cessation")
  ))# add group variable for testing
saveRDS(pheno_stool_extend, file = "./1-output/pheno_stool_S_ctrl.Rds")


