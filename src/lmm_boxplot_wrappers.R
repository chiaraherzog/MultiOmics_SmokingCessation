
# boxplot wrappers for plotting DA detected with lmm models
# wrappers depends on paired_longitudinal_compliance.R
# extra function for grouped boxplots for comparing smokers to control group 
#(boxplot_smokers_ctrl for taxa, boxplot_smokers_ctrl2 for functional data)
# Author: Charlotte Vavourakis
# version for smoking study arm

source("src/paired_longitudinal_compliance.R")

cols <- c("#1b69a1", "#48a0af", "#71b5a9", "#ec6669", "#f39668", "#bd647d", "#832c9b", "#5f70a8")

cols_comp <- cols[c(4, 1)]
names(cols_comp) <- c("lower compliance", "higher compliance")

cols_group <- cols[c(8, 6, 4, 1, 4, 1, 4, 1)]
names(cols_group) <- c("Never smoker (control)","Smoker baseline",
                       "M2 no smoking cessation","M2 smoking cessation",
                       "M4 no smoking cessation","M4 smoking cessation",
                       "M6 no smoking cessation","M6 smoking cessation")

boxplot_wrap_interaction <- function(dat1,dat2,exp,pval,interaction_only=F, taxa = NULL) {
  
  require(tidyverse)
  
  if (!is.null(taxa)){
    
    toPlot = taxa
    
  } else {
    
    # select significant taxa
    if (interaction_only == T) {
      
      toPlot = dat1 %>%
        dplyr::filter(`p.value_visitIdM2:compliancehigher compliance` < pval | 
                        `p.value_visitIdM4:compliancehigher compliance` < pval |
                        `p.value_visitIdM6:compliancehigher compliance` < pval) %>%
        pull(x)
      
    } else {
      
      toPlot = dat1 %>%
        dplyr::filter(`p.value_visitIdM2:compliancehigher compliance` < pval | 
                        `p.value_visitIdM4:compliancehigher compliance` < pval |
                        `p.value_visitIdM6:compliancehigher compliance` < pval |
                        `p.value_visitIdM2` < pval |
                        `p.value_visitIdM4` < pval |
                        `p.value_visitIdM6` < pval) %>%
        pull(x)
    }
  }
  
  
  # generate plots
  listPlots = list()
  for (i in 1:length(toPlot)){
    listPlots[[i]] = paired_longitudinal_compliance(dat2, 
                                                    variable = toPlot[i],
                                                    ylab = paste0(toPlot[i]), 
                                                    colours = cols_comp,
                                                    p = 'p.signif') +
      theme(aspect.ratio=1)
  }
  
  return(listPlots)
  
}

boxplot_wrap_time <- function(dat1,dat2,exp,pval){
  
  require(tidyverse)
  
  # select significant taxa
  toPlot = dat1 %>%
    dplyr::filter(`p.value_visitIdM2` < pval |
                    `p.value_visitIdM4` < pval |
                    `p.value_visitIdM6` < pval) %>%
    pull(x)
  
  # generate plots
  listPlots = list()
  for (i in 1:length(toPlot)){
    listPlots[[i]] = paired_longitudinal_compliance(dat2, 
                                                    variable = toPlot[i],
                                                    ylab = paste0(toPlot[i]),
                                                    colours = cols_comp,
                                                    p = 'p.signif')  +
      theme(aspect.ratio=1)
  }
  
  return(listPlots)
  
}

boxplot_smokers_ctrl <- function(dat1,dat2,exp,pval,interaction_only=F, taxa = NULL){
  
  require(tidyverse)
  require(gghalves)
  require(patchwork)
  
  if (!is.null(taxa)){
    
    toPlot = taxa
    
  } else {
  
    # select significant taxa
    if (any(grepl("p.value_visitIdM2:compliancehigher compliance", colnames(dat1)))){
  
      if (interaction_only == T) {
        
        toPlot = dat1 %>%
          dplyr::filter(`p.value_visitIdM2:compliancehigher compliance` < pval | 
                          `p.value_visitIdM4:compliancehigher compliance` < pval |
                          `p.value_visitIdM6:compliancehigher compliance` < pval) %>%
          pull(x)
        
      } else {
        
        toPlot = dat1 %>%
          dplyr::filter(`p.value_visitIdM2:compliancehigher compliance` < pval | 
                          `p.value_visitIdM4:compliancehigher compliance` < pval |
                          `p.value_visitIdM6:compliancehigher compliance` < pval |
                          `p.value_visitIdM2` < pval |
                          `p.value_visitIdM4` < pval |
                          `p.value_visitIdM6` < pval) %>%
          pull(x)
      }
      
    } else{
      toPlot = dat1 %>%
        dplyr::filter(`p.value_visitIdM2` < pval |
                        `p.value_visitIdM4` < pval |
                        `p.value_visitIdM6` < pval) %>%
        pull(x)
    }
    
  }
  
  # generate plot

  dat2 <- dat2 %>%
    dplyr::filter(!is.na(group) & Family %in% toPlot) %>%
    dplyr::filter(group %in% c("Never smoker (control)", "Smoker baseline", "M6 smoking cessation", "M6 no smoking cessation")) 
  
  dat2$group <- factor(dat2$group, levels = names(cols_group))
  
  plot <- dat2 %>%
    ggplot(aes(x = group,
               y = rel.abund.)) +
    # geom_boxplot(alpha = 0.3,
    #              aes(fill = group)) +
    gghalves::geom_half_boxplot(aes(fill = group),
                                alpha = 0.3,
                                width = 0.5,
                                outlier.shape = NA,errorbar.length = 0) +
    gghalves::geom_half_point_panel(alpha = 0.6,
                                    size = 1,
                                    aes(colour = group)) +
    ggpubr::stat_compare_means(ref.group = "Never smoker (control)",
                               label = 'p.signif',
                               size = 2.7,
                               method = "wilcox.test",
                               paired = F,
                               label.y.npc = 0.93) +
    scale_fill_manual(values = cols_group, name = "", aesthetics = c('colour', 'fill')) +
    theme_bw() +
    theme(legend.position = 'right',
          aspect.ratio = 0.6,
          #axis.text.x = element_text(angle = 60, hjust = 1),
          axis.text.x = element_blank(),
          axis.title.y = ggtext::element_markdown(),
          strip.text.x = ggtext::element_markdown(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank()) +
    labs(x = "") +
    facet_wrap(~Family, scales = "free", ncol=2)
  
  
  return(plot)
  
}

boxplot_smokers_ctrl2 <- function(dat1,dat2,exp,pval,interaction_only=F){
  
  require(tidyverse)
  require(patchwork)
  
  # select significant features
  if (any(grepl("p.value_visitIdM2:compliancehigher compliance", colnames(dat1)))){
    
    if (interaction_only == T) {
      
      toPlot = dat1 %>%
        dplyr::filter(`p.value_visitIdM2:compliancehigher compliance` < pval | 
                        `p.value_visitIdM4:compliancehigher compliance` < pval |
                        `p.value_visitIdM6:compliancehigher compliance` < pval) %>%
        pull(x)
      
    } else {
      
      toPlot = dat1 %>%
        dplyr::filter(`p.value_visitIdM2:compliancehigher compliance` < pval | 
                        `p.value_visitIdM4:compliancehigher compliance` < pval |
                        `p.value_visitIdM6:compliancehigher compliance` < pval |
                        `p.value_visitIdM2` < pval |
                        `p.value_visitIdM4` < pval |
                        `p.value_visitIdM6` < pval) %>%
        pull(x)
    }
    
  } else{
    toPlot = dat1 %>%
      dplyr::filter(`p.value_visitIdM2` < pval |
                      `p.value_visitIdM4` < pval |
                      `p.value_visitIdM6` < pval) %>%
      pull(x)
  }
  
  # generate plot
  
  dat2 <- dat2 %>%
    dplyr::filter(!is.na(group) & rowname %in% toPlot) %>%
    dplyr::filter(group %in% c("Never smoker (control)", "Smoker baseline", "M6 smoking cessation", "M6 no smoking cessation")) 
  
  dat2$group <- factor(dat2$group, levels = names(cols_group))
  
  plot <- dat2 %>%
    ggplot(aes(x = group,
               y = value)) +
    geom_boxplot(alpha = 0.3,
                 aes(fill = group)) +
    ggpubr::stat_compare_means(ref.group = "Never smoker (control)",
                               label = 'p.signif',
                               size = 2.7,
                               method = "wilcox.test",
                               paired = F,
                               label.y.npc = 0.95) +
    scale_fill_manual(values = cols_group) +
    theme_bw() +
    theme(legend.position = 'right',
          aspect.ratio = 0.5,
          #axis.text.x = element_text(angle = 60, hjust = 1),
          axis.text.x = element_blank(),
          axis.title.y = ggtext::element_markdown(),
          strip.text.x = ggtext::element_markdown(),
          axis.ticks.x = element_blank()) +
    labs(x = "", y = "abundance") +
    facet_wrap(~rowname, scales = "free")
  
  
  return(plot)
  
}
