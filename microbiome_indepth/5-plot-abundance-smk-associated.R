
# Author: Charlotte Vavourakis

library(tidyverse)
library(ape)
library(ggtreeExtra)
library(ggtree)
library(ggnewscale)
library(rcartocolor)
library(ggpubr)
library(patchwork)
library(viridis)

dir.create("./5-output")

### Saliva ###----

#### relabel tree and plot #----
taxlab1 <- readRDS("./4-output/list_saliva.Rds")[[1]]
taxlab1 <- taxlab1 %>%
  mutate(label = case_when(
    is.na(Species) ~ paste0(Genus," sp."),
    grepl("/",Species) ~ paste0(Genus," sp."),
    TRUE ~ paste0(Genus," ",Species)
    )) %>%
  dplyr::mutate(label = case_when(
    Genus %in% c("Streptococcus", "Veillonella") ~ paste0(label," ",OTU),
    TRUE ~ label
  ))

# mark recovery yes/no
recovery <- readRDS("./4-output/list_saliva.Rds")[[6]] %>%
  mutate(recovery = case_when(
    association_cessation %in% c(".","*","**","***") ~ "yes", TRUE ~ "no"
  ))
cols_recovery <- c("black","white")
names(cols_recovery) <- c("yes","no")
shapes_recovery <- c(21,25)
names(shapes_recovery) <- c("yes","no")
taxlab1 <- left_join(taxlab1,recovery, by = c("OTU" = "ASV")) 

phy_tree1 <- read.tree("./5-output/saliva_ASV.nwk")
taxlab1 <- taxlab1[match(phy_tree1$tip.label, taxlab1$OTU),]
identical(phy_tree1$tip.label,c(taxlab1$OTU))
phy_tree1$tip.label <- c(taxlab1$label)
taxlab1 <- taxlab1 %>% relocate(label)

p1 <- ggtree(phy_tree1,
             branch.length="none",
             layout = "rectangular", 
             size = 0.2) %<+%
  taxlab1 +
  geom_tippoint(aes(shape = recovery,
                    fill = recovery)) +
  scale_fill_manual(values=cols_recovery) +
  scale_shape_manual(values=shapes_recovery) +
  theme(plot.margin=margin(0,0,0,0),
         legend.position="bottom") +
  guides(fill=guide_legend(ncol=3),
         shape=guide_legend(ncol=1))+
  geom_tiplab(size=5,offset = 2.3) + 
  coord_cartesian(clip = "off") + # allow plotting outside of plot panel (see: https://github.com/tidyverse/ggplot2/issues/2536)
  theme(plot.margin = margin(0,250,20,0, unit = "pt")) # increase margins

#### Relative abundance heatmaps #----
# Never smokers
metadat1 <- readRDS("./4-output/list_saliva.Rds")[[3]] %>%
  dplyr::filter(group == "Never smoker (control)") 
dat1 <- readRDS("./4-output/list_saliva.Rds")[[2]] %>%
  dplyr::select(c(metadat1$sampleId)) %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxlab1[,1:2], by = c("ASV" = "OTU")) %>%
  dplyr::select(-ASV,) %>%
  pivot_longer(cols = -label, names_to = "sampleId", values_to = "abund") %>%
  group_by(label) %>%
  summarise(median_value = median(abund, na.rm = TRUE))

# Smoker baseline
metadat2 <- readRDS("./4-output/list_saliva.Rds")[[3]] %>%
  dplyr::filter(group == "Smoker baseline") 
dat2 <- readRDS("./4-output/list_saliva.Rds")[[2]] %>%
  dplyr::select(c(metadat2$sampleId)) %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxlab1[,1:2], by = c("ASV" = "OTU")) %>%
  dplyr::select(-ASV,) %>%
  pivot_longer(cols = -label, names_to = "sampleId", values_to = "abund")  %>%
  group_by(label) %>%
  summarise(median_value = median(abund, na.rm = TRUE))

m <- max(c(max(dat1$median_value),max(dat2$median_value)))

p2 <- p1 + 
  new_scale_fill() + 
  geom_fruit(
    data=dat1,
    geom=geom_tile,
    mapping=aes(y=label, fill=median_value),
    offset = 0.05,   
    pwidth = 0.4 
  ) + 
  scale_fill_carto_c(name = "Never smoker\n(control)\nmedian rel. abund.",
                     type = "sequential", palette = "DarkMint", direction = 1,
                     limits = c(0, m),
                     breaks = seq(0, m, by = 0.02)) +
  new_scale_fill() +
  geom_fruit(
    data=dat2,
    geom=geom_tile,
    mapping=aes(y=label, fill=median_value),
    offset = 0.06,
    pwidth = 0.4 
    ) + 
  scale_fill_carto_c(name = "Smoker baseline\nmedian rel. abund.",
                     type = "sequential", palette = "DarkMint", direction = 1,
                     limits = c(0, m),
                     breaks = seq(0, m, by = 0.02))

#### Smoking association #----

dat4 <- readRDS("./4-output/list_saliva.Rds")[[7]]%>%
  left_join(taxlab1[,1:2], by = c("feature" = "OTU")) %>%
  dplyr::select(-feature)

p3 <- p2 +
  geom_fruit(
    data = dat4,
    geom=geom_text,
    mapping=aes(y = label, label=association_smoking),
    offset = 0.06
  )

#### Relative abundance heatmaps baseline change #----

# M6 smoking cessation
metadat1 <- readRDS("./4-output/list_saliva.Rds")[[5]] %>%
  dplyr::filter(group == "M6 smoking cessation") 
dat1 <- readRDS("./4-output/list_saliva.Rds")[[4]] %>%
  dplyr::select(c(metadat1$sampleId)) %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxlab1[,1:2], by = c("ASV" = "OTU")) %>%
  dplyr::select(-ASV,) %>%
  pivot_longer(cols = -label, names_to = "sampleId", values_to = "abund")  %>%
  group_by(label) %>%
  summarise(mean_value = mean(abund, na.rm = TRUE))

# M6 no smoking cessation
metadat2 <- readRDS("./4-output/list_saliva.Rds")[[5]] %>%
  dplyr::filter(group == "M6 no smoking cessation") 
dat2 <- readRDS("./4-output/list_saliva.Rds")[[4]] %>%
  dplyr::select(c(metadat2$sampleId)) %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxlab1[,1:2], by = c("ASV" = "OTU")) %>%
  dplyr::select(-ASV,) %>%
  pivot_longer(cols = -label, names_to = "sampleId", values_to = "abund")  %>%
  group_by(label) %>%
  summarise(mean_value = mean(abund, na.rm = TRUE))

m <- max(c(max(dat1$mean_value),max(dat2$mean_value)))
  
p4 <- p3 + 
  new_scale_fill() + 
  geom_fruit(
    data=dat1,
    geom=geom_tile,
    mapping=aes(y=label, fill=mean_value),
    offset = 0.06,
    pwidth = 0.4 
  ) + 
  scale_fill_carto_c(name = "M6 smoking cessation\nmean baseline change",
                     type = "diverging", palette = "Earth", direction = -1,
                     limits = c(-m, m)) +
  new_scale_fill() +
  geom_fruit(
    data=dat2,
    geom=geom_tile,
    mapping=aes(y=label, fill=mean_value),
    offset = 0.061,
    pwidth = 0.4 
  ) + 
  scale_fill_carto_c(name = "M6 no smoking cessation\nmean baseline change",
                     type = "diverging", palette = "Earth", direction = -1,
                     limits = c(-m, m)) 

p4 <- p4 +
  theme(legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=10),
        legend.position = "bottom",
        legend.direction = "horizontal") 

ggsave("./5-output/saliva_ggtree.pdf", plot = p4, device = "pdf", width = 20, height = 15, units = "cm")

### Stool ###----

#### relabel tree and plot #----
taxlab1 <- readRDS("./4-output/list_stool.Rds")[[1]]
taxlab1 <- taxlab1 %>%
  mutate(label = case_when(
    is.na(Genus) ~ "unclassified",
    is.na(Species) & !is.na(Genus) ~ paste0(Genus," sp."),
    grepl("/",Species) ~ paste0(Genus," sp."),
    TRUE ~ paste0(Genus," ",Species)
  ))

# mark recovery yes/no
recovery <- readRDS("./4-output/list_stool.Rds")[[6]] %>%
  mutate(recovery = case_when(
    association_cessation %in% c(".","*","**","***") ~ "yes", TRUE ~ "no"
  ))
cols_recovery <- c("black","white")
names(cols_recovery) <- c("yes","no")
shapes_recovery <- c(21,25)
names(shapes_recovery) <- c("yes","no")
taxlab1 <- left_join(taxlab1,recovery, by = c("OTU" = "ASV")) 

phy_tree1 <- read.tree("./5-output/stool_ASV.nwk")
taxlab1 <- taxlab1[match(phy_tree1$tip.label, taxlab1$OTU),]
identical(phy_tree1$tip.label,c(taxlab1$OTU))
phy_tree1$tip.label <- c(taxlab1$label)
taxlab1 <- taxlab1 %>% relocate(label)

p1 <- ggtree(phy_tree1,
             branch.length="none",
             layout = "rectangular", 
             size = 0.2) %<+%
  taxlab1 +
  geom_tippoint(aes(shape = recovery,
                    fill = recovery)) +
  scale_fill_manual(values=cols_recovery) +
  scale_shape_manual(values=shapes_recovery) +
  theme(plot.margin=margin(0,0,0,0),
        legend.position="bottom") +
  guides(fill=guide_legend(ncol=3),
         shape=guide_legend(ncol=1))+
  geom_tiplab(size=5,offset = 2.1) + 
  coord_cartesian(clip = "off") + # allow plotting outside of plot panel (see: https://github.com/tidyverse/ggplot2/issues/2536)
  theme(plot.margin = margin(0,250,20,0, unit = "pt")) # increase margins

#### Relative abundance heatmaps #----
# Never smokers
metadat1 <- readRDS("./4-output/list_stool.Rds")[[3]] %>%
  dplyr::filter(group == "Never smoker (control)") 
dat1 <- readRDS("./4-output/list_stool.Rds")[[2]] %>%
  dplyr::select(c(metadat1$sampleId)) %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxlab1[,1:2], by = c("ASV" = "OTU")) %>%
  dplyr::select(-ASV,) %>%
  pivot_longer(cols = -label, names_to = "sampleId", values_to = "abund") %>%
  group_by(label) %>%
  summarise(median_value = median(abund, na.rm = TRUE))

# Smoker baseline
metadat2 <- readRDS("./4-output/list_stool.Rds")[[3]] %>%
  dplyr::filter(group == "Smoker baseline") 
dat2 <- readRDS("./4-output/list_stool.Rds")[[2]] %>%
  dplyr::select(c(metadat2$sampleId)) %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxlab1[,1:2], by = c("ASV" = "OTU")) %>%
  dplyr::select(-ASV,) %>%
  pivot_longer(cols = -label, names_to = "sampleId", values_to = "abund")  %>%
  group_by(label) %>%
  summarise(median_value = median(abund, na.rm = TRUE))

m <- max(c(max(dat1$median_value),max(dat2$median_value)))

p2 <- p1 + 
  new_scale_fill() + 
  geom_fruit(
    data=dat1,
    geom=geom_tile,
    mapping=aes(y=label, fill=median_value),
    offset = 0.065,   
    pwidth = 0.35 
  ) + 
  scale_fill_carto_c(name = "Never smoker\n(control)\nmedian rel. abund.",
                     type = "sequential", palette = "DarkMint", direction = 1,
                     limits = c(0, m),
                     breaks = seq(0, m, by = 0.005)) +
  new_scale_fill() +
  geom_fruit(
    data=dat2,
    geom=geom_tile,
    mapping=aes(y=label, fill=median_value),
    offset = 0.075,   
    pwidth = 0.35 
  ) + 
  scale_fill_carto_c(name = "Smoker baseline\nmedian rel. abund.",
                     type = "sequential", palette = "DarkMint", direction = 1,
                     limits = c(0, m),
                     breaks = seq(0, m, by = 0.005))


#### Smoking association #----

dat4 <- readRDS("./4-output/list_stool.Rds")[[7]]%>%
  left_join(taxlab1[,1:2], by = c("feature" = "OTU")) %>%
  dplyr::select(-feature)

p3 <- p2 +
  geom_fruit(
    data = dat4,
    geom=geom_text,
    mapping=aes(y = label, label=association_smoking),
    offset=0.075
  )

#### Relative abundance heatmaps baseline change #----

# M6 smoking cessation
metadat1 <- readRDS("./4-output/list_stool.Rds")[[5]] %>%
  dplyr::filter(group == "M6 smoking cessation") 
dat1 <- readRDS("./4-output/list_stool.Rds")[[4]] %>%
  dplyr::select(c(metadat1$sampleId)) %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxlab1[,1:2], by = c("ASV" = "OTU")) %>%
  dplyr::select(-ASV,) %>%
  pivot_longer(cols = -label, names_to = "sampleId", values_to = "abund")  %>%
  group_by(label) %>%
  summarise(mean_value = mean(abund, na.rm = TRUE))

# M6 no smoking cessation
metadat2 <- readRDS("./4-output/list_stool.Rds")[[5]] %>%
  dplyr::filter(group == "M6 no smoking cessation") 
dat2 <- readRDS("./4-output/list_stool.Rds")[[4]] %>%
  dplyr::select(c(metadat2$sampleId)) %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxlab1[,1:2], by = c("ASV" = "OTU")) %>%
  dplyr::select(-ASV,) %>%
  pivot_longer(cols = -label, names_to = "sampleId", values_to = "abund")  %>%
  group_by(label) %>%
  summarise(mean_value = mean(abund, na.rm = TRUE))

m <- max(c(max(dat1$mean_value),max(dat2$mean_value)))

p4 <- p3 + 
  new_scale_fill() + 
  geom_fruit(
    data=dat1,
    geom=geom_tile,
    mapping=aes(y=label, fill=mean_value),
    offset = 0.075,
    pwidth = 0.35 
  ) + 
  scale_fill_carto_c(name = "M6 smoking cessation\nmean baseline change",
                     type = "diverging", palette = "Earth", direction = -1,
                     limits = c(-m, m)) +
  new_scale_fill() +
  geom_fruit(
    data=dat2,
    geom=geom_tile,
    mapping=aes(y=label, fill=mean_value),
    offset = 0.076,
    pwidth = 0.35 
  ) + 
  scale_fill_carto_c(name = "M6 no smoking cessation\nmean baseline change",
                     type = "diverging", palette = "Earth", direction = -1,
                     limits = c(-m, m)) 

p4 <- p4 +
  theme(legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=10),
        legend.position = "bottom",
        legend.direction = "horizontal") 

ggsave("./5-output/stool_ggtree.pdf", plot = p4, device = "pdf", width = 20, height = 15, units = "cm")
