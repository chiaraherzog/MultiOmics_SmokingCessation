
# get sequences for Smoking arm ASV tables extended with ctrl group for calculating unifrac distances
# Author: Charlotte Vavourakis

library(tidyverse)
library(MultiAssayExperiment)
library(Biostrings)

dir.create("./2-output")

fasta_path <- "~/Dropbox/tg-data-prep/prep-scripts/11-preprocessed/2-microbiome/1-extended-taxonomy/1-output"
fasta <- readDNAStringSet(file.path(fasta_path,"tg_all_microbiome_asvs.fasta"))
names(fasta) <- gsub(" .*", "", names(fasta))

# stool ##----
ASVtab <- readRDS("./1-output/ASVtable_stool_S_ctrl.Rds")
fasta2 <- fasta[c(ASVtab$OTU)]
out <- ASVtab %>% filter(Kingdom=="Archaea") %>% pull(OTU) # archaea for root
write.table(out, file="./2-output/archaeal_ASV.txt", row.names=F,col.names=F, quote=F)
writeXStringSet(fasta2, file = "./2-output/stool_S_ctrl_ASV.fasta", format = "fasta")

# saliva ##----
ASVtab <- readRDS("./1-output/ASVtable_saliva_S_ctrl.Rds")
fasta2 <- fasta[c(ASVtab$OTU)]
out <- ASVtab %>% filter(Phylum=="Patescibacteria") %>% pull(OTU) # cpr for root
write.table(out, file="./2-output/cpr_ASV.txt", row.names=F,col.names=F, quote=F)
writeXStringSet(fasta2, file = "./2-output/saliva_S_ctrl_ASV.fasta", format = "fasta")