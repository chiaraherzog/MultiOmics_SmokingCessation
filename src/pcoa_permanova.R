# Author: Charlotte Vavourakis

cols <- c("#1b69a1", "#48a0af", "#71b5a9", "#ec6669", "#f39668", "#bd647d", "#832c9b", "#5f70a8")

cols_group <- cols[c(8, 6, 4, 1, 4, 1, 4, 1)]
shapes_group <- c(21,25,24,23,24,23,24,23)
names(shapes_group) <- names(cols_group) <- c("Never smoker (control)","Smoker baseline",
                                              "M2 no smoking cessation","M2 smoking cessation",
                                              "M4 no smoking cessation","M4 smoking cessation",
                                              "M6 no smoking cessation","M6 smoking cessation")

### helper function to create PCOA plots for a given timepoint ###----
plotPcoa <- function(data, visitId, type, var, varname, root_tr, my_title = NULL) {
  
  library(vegan)
  library(ecodist)
  library(tidyverse)
  
  outList <- list()
  
  if (type=="Hell_eucl") {
    if(is.null(my_title)){
      tit = paste0("Euclidian-Hellinger PCoA, ", visitId)
    }
    dist = vegan::vegdist(vegan::decostand(t(data$abund), method = "hellinger"), method="euclidean")
  }
  
  else if (type=="Hell_bray") {
    if(is.null(my_title)){
      tit = paste0("Bray-Curtis-Hellinger PCoA, ", visitId)
    }
    dist = vegan::vegdist(vegan::decostand(t(data$abund), method = "hellinger"), method="bray")
  }
  
  else if (type=="Aitchison") {
    if(is.null(my_title)){
      tit = paste0("Aitchison PCoA, ", visitId)
    }
    dist = vegan::vegdist(t(data$abund), method="robust.aitchison")
  }
  
  else if (type=="GUniFrac") {
    
    require(ape)
    require(MiSPU)
    require(Matrix)
    
    if(is.null(my_title)){
      tit = paste0("Generalized UniFrac PCoA",visitId)
    }
    tree = ape::read.tree(file = root_tr)
    dist2 = GUniFrac(t(data$abund), tree)
    dist2 = dist2$d5
    colnames(dist2) = colnames(data$abund)
    rownames(dist2) = colnames(data$abund)
    # make a lower-triangular dissimilarity matrix
    dist = ecodist::lower(dist2)
    attr(dist, "Labels") = colnames(data$abund)
    attr(dist, "Diag") = FALSE
    attr(dist, "Upper") = FALSE
    class(d) = "dist"
  }
  
  pcoa = ecodist::pco(dist)
  
  df = data.frame(pcoa1 = pcoa$vectors[,1], 
                  pcoa2 = pcoa$vectors[,2],
                  pcoa3 = pcoa$vectors[,3])
  
  df = cbind(df,varname = data$metadata %>% pull(var))
  
  df$varname = factor(df$varname, levels = names(cols_group))
  
  outList[[1]] = ggplot(data = df, 
                        aes(x=pcoa1, 
                            y=pcoa2,
                            fill=varname,
                            color=varname,
                            shape=varname)) +
    geom_point() +
    labs(x = "PC1",
         y = "PC2", 
         title = my_title,
         color=varname) +
    stat_ellipse(aes(color = varname), type="norm", linetype=2) +
    scale_fill_manual(values = cols_group) + 
    scale_color_manual(values = cols_group) + 
    scale_shape_manual(values = shapes_group) +
    theme_bw() +
    theme(aspect.ratio = 1,
          legend.position = "bottom") +
    labs(fill  = "", color="", shape = "")+
    guides(color=guide_legend(ncol =1),
           fill=guide_legend(ncol=1),
           shape=guide_legend(ncol=1))
    
  
  outList[[2]] = ggplot(data = df, 
                        aes(x=pcoa1, 
                            y=pcoa3,
                            color=varname,
                            shape=varname)) +
    geom_point() +
    labs(x = "PC1",
         y = "PC3", 
         title = my_title,
         color=varname) +
    stat_ellipse(aes(color = varname), type="norm", linetype=2) +
    scale_color_manual(values = cols_group) +  
    scale_shape_manual(values = shapes_group) +
    theme_bw() +
    theme(aspect.ratio = 1,
          legend.position = "bottom") +
    guides(color=guide_legend(ncol =1))
    
  
  if (type=="GUniFrac"){
    outList[[3]] = dist2
  }
  
  return(outList)
  
}

### helper function for pairwaise significance testing  ###----
pairwise.adonis <- function(x,factors, sim.method, p.adjust.m, root_tr){
  
  library(vegan)
  
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  if (sim.method=="GUniFrac"){
    library(ape)
    library(MiSPU)
    library(Matrix)
    
    tree = ape::read.tree(file = root_tr)
    
    for(elem in 1:ncol(co)){
      dist = GUniFrac(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),], tree)
      dist = dist$d5
      colnames(dist) = colnames(t(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),]))
      rownames(dist) = colnames(t(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),]))
      ad = vegan::adonis2(dist ~
                            factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]);
      pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
      F.Model =c(F.Model,ad$F[1]);
      R2 = c(R2,ad$R2[1]);
      p.value = c(p.value,ad$`Pr(>F)`[1])
    }   
    
    
  } else{
    for(elem in 1:ncol(co)){
      ad = vegan::adonis2(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                            factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method);
      pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
      F.Model =c(F.Model,ad$F[1]);
      R2 = c(R2,ad$R2[1]);
      p.value = c(p.value,ad$`Pr(>F)`[1])
    }
  }
  
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  
  return(pairw.res)
}

### helper function for making summary full permanova analysis###----
# Hellinger transformation + Bray-Curtis
# Aitchison
# Generalized Unifrac
# additional a posteriori pair-wise testing M2,M4,M6
summarize_permanova_full <- function(d, tree){
  
  library(ampvis2)
  library(vegan)
  library(tidyverse)
  library(ape)
  library(MiSPU)
  library(Matrix)
  
  d0 <- amp_filter_samples(d,visitId %in% c("Never smoker (control)","M0"))
  d2 <- amp_filter_samples(d,visitId %in% c("Never smoker (control)","M2"))
  d4 <- amp_filter_samples(d,visitId %in% c("Never smoker (control)","M4"))
  d6 <- amp_filter_samples(d,visitId %in% c("Never smoker (control)","M6"))
  
  ## distance measures without phylogeny #----
  
  list_in <- list(labels = c("Bray-Curtis after Hellinger","Aitchison"),
            decostand_methods = c("hellinger","none"),
            vegdist_methods = c("bray","robust.aitchison"))
  list_out <- c()
  
  # helper function for data transformation
  transform_h <- function(dat,decostand_method){
    if (decostand_method == "none"){
      return(dat)
    } else {
      return(vegan::decostand(dat, method = decostand_method))
    }
  }
  
  for (i in 1:length(list_in$labels)){
    
    ## Make a table for saving PCOA + permanova testing results
    permTable <- as.data.frame(matrix(nrow=13, ncol=4))
    colnames(permTable) <- c("Distance","Comparison","P value","Adjusted p value")
    
    # M0
    res <- vegan::adonis2(vegan::vegdist(transform_h(t(d0$abund), list_in$decostand_methods[i]), method=list_in$vegdist_methods[i]) ~ group,
                          data=d0$metadata, 
                          permutations = 9999)
    permTable[1,] <- c(list_in$labels[i], "Smoker baseline vs Never smoker (control)", res$`Pr(>F)`[1], "NA")
    
    # M2
    res <- vegan::adonis2(vegan::vegdist(transform_h(t(d2$abund), list_in$decostand_methods[i]), method=list_in$vegdist_methods[i]) ~ group,
                          data=d2$metadata,
                          permutations = 9999)
    permTable[2,] <- c(list_in$labels[i], "M2 no smoking cessation vs M2 smoking cessation vs Never smoker (control)", res$`Pr(>F)`[1], "NA")
    
    # pair-wise pair-wise a posteriori tests M2
    PW.Adonis <- pairwise.adonis(transform_h(t(d2$abund),list_in$decostand_methods[i]),
                                 d2$metadata$group,
                                 sim.method=list_in$vegdist_methods[i],
                                 p.adjust.m = "bonferroni", 
                                 NULL)
    
    permTable[3,] <- c(list_in$labels[i], PW.Adonis[1,1], PW.Adonis[1,4], PW.Adonis[1,5])
    permTable[4,] <- c(list_in$labels[i], PW.Adonis[2,1], PW.Adonis[2,4], PW.Adonis[2,5])
    permTable[5,] <- c(list_in$labels[i], PW.Adonis[3,1], PW.Adonis[3,4], PW.Adonis[3,5])
    
    # M4
    res <- vegan::adonis2(vegan::vegdist(transform_h(t(d4$abund), list_in$decostand_methods[i]), method=list_in$vegdist_methods[i]) ~ group,
                          data=d4$metadata,
                          permutations = 9999)
    permTable[6,] <- c(list_in$labels[i], "M4 no smoking cessation vs M4 smoking cessation vs Never smoker (control)", res$`Pr(>F)`[1], "NA")
    
    # pair-wise pair-wise a posteriori tests M4
    PW.Adonis <- pairwise.adonis(transform_h(t(d4$abund), list_in$decostand_methods[i]),
                                 d4$metadata$group,
                                 sim.method=list_in$vegdist_methods[i],
                                 p.adjust.m = "bonferroni",
                                 NULL)
    
    permTable[7,] <- c(list_in$labels[i], PW.Adonis[1,1], PW.Adonis[1,4], PW.Adonis[1,5])
    permTable[8,] <- c(list_in$labels[i], PW.Adonis[2,1], PW.Adonis[2,4], PW.Adonis[2,5])
    permTable[9,] <- c(list_in$labels[i], PW.Adonis[3,1], PW.Adonis[3,4], PW.Adonis[3,5])
    
    # M6
    res <- vegan::adonis2(vegan::vegdist(transform_h(t(d6$abund), list_in$decostand_methods[i]), method=list_in$vegdist_methods[i]) ~ group,
                          data=d6$metadata,
                          permutations = 9999)
    permTable[10,] <- c(list_in$labels[i], "M6 no smoking cessation vs M6 smoking cessation vs control", res$`Pr(>F)`[1], "NA")
    
    # pair-wise pair-wise a posteriori tests M6
    PW.Adonis <- pairwise.adonis(transform_h(t(d6$abund), list_in$decostand_methods[i]),
                                 d6$metadata$group,
                                 sim.method=list_in$vegdist_methods[i],
                                 p.adjust.m = "bonferroni",
                                 NULL)
    
    permTable[11,] <- c(list_in$labels[i], PW.Adonis[1,1], PW.Adonis[1,4], PW.Adonis[1,5])
    permTable[12,] <- c(list_in$labels[i], PW.Adonis[2,1], PW.Adonis[2,4], PW.Adonis[2,5])
    permTable[13,] <- c(list_in$labels[i], PW.Adonis[3,1], PW.Adonis[3,4], PW.Adonis[3,5])
    
    list_out[[i]] <- permTable
  }
  
  table <- do.call(rbind,list_out)
  
  ## Generalized Unifrac #----
  
  ## Make a table for saving PCOA + permanova testing results
  permTable <- as.data.frame(matrix(nrow=13, ncol=4))
  colnames(permTable) <- c("Distance","Comparison","P value","Adjusted p value")
 
  # M0
  plotsM0 <- plotPcoa(d0, "M0", "GUniFrac", "group", "", tree)
  
  res <- vegan::adonis2(plotsM0[[3]] ~ group, 
                        data=d0$metadata, 
                        permutations = 9999)
  permTable[1,] <- c("GUniFrac", "Smoker baseline vs Never smoker (control)", res$`Pr(>F)`[1], "NA")
  
  # M2
  plotsM2 <- plotPcoa(d2, "M2", "GUniFrac", "group", "", tree)
  
  res <- vegan::adonis2(plotsM2[[3]] ~ group, 
                        data=d2$metadata, 
                        permutations = 9999)
  permTable[2,] <- c("GUniFrac", "M2 no smoking cessation vs M2 smoking cessation vs Never smoker (control)", res$`Pr(>F)`[1], "NA")
  
  # pair-wise pair-wise a posteriori tests M2
  PW.Adonis <- pairwise.adonis(t(d2$abund),
                               d2$metadata$group,
                               sim.method="GUniFrac",
                               p.adjust.m = "bonferroni",
                               tree)
  
  permTable[3,] <- c("GUniFrac", PW.Adonis[1,1], PW.Adonis[1,4], PW.Adonis[1,5])
  permTable[4,] <- c("GUniFrac", PW.Adonis[2,1], PW.Adonis[2,4], PW.Adonis[2,5])
  permTable[5,] <- c("GUniFrac", PW.Adonis[3,1], PW.Adonis[3,4], PW.Adonis[3,5])
  
  # M4
  plotsM4 <- plotPcoa(d4, "M4", "GUniFrac", "group", "", tree)
  
  res <- vegan::adonis2(plotsM4[[3]] ~ group, 
                        data=d4$metadata, 
                        permutations = 9999)
  permTable[6,] <- c("GUniFrac", "M4 no smoking cessation vs M4 smoking cessation vs Never smoker (control)", res$`Pr(>F)`[1], "NA")
  
  # pair-wise pair-wise a posteriori tests M4
  PW.Adonis <- pairwise.adonis(t(d4$abund),
                               d4$metadata$group,
                               sim.method="GUniFrac",
                               p.adjust.m = "bonferroni",
                               tree)
  
  permTable[7,] <- c("GUniFrac", PW.Adonis[1,1], PW.Adonis[1,4], PW.Adonis[1,5])
  permTable[8,] <- c("GUniFrac", PW.Adonis[2,1], PW.Adonis[2,4], PW.Adonis[2,5])
  permTable[9,] <- c("GUniFrac", PW.Adonis[3,1], PW.Adonis[3,4], PW.Adonis[3,5])
  
  # M6
  plotsM6 <- plotPcoa(d6, "M6", "GUniFrac", "group", "", tree)
  
  res <- vegan::adonis2(plotsM6[[3]] ~ group, 
                        data=d6$metadata, 
                        permutations = 9999)
  permTable[10,] <- c("GUniFrac", "M6 no smoking cessation vs M6 smoking cessation vs Never smoker (control)", res$`Pr(>F)`[1], "NA")
  
  ## pair-wise pair-wise a posteriori tests M6
  PW.Adonis <- pairwise.adonis(t(d6$abund),
                               d6$metadata$group,
                               sim.method="GUniFrac",
                               p.adjust.m = "bonferroni",
                               tree)
  
  permTable[11,] <- c("GUniFrac", PW.Adonis[1,1], PW.Adonis[1,4], PW.Adonis[1,5])
  permTable[12,] <- c("GUniFrac", PW.Adonis[2,1], PW.Adonis[2,4], PW.Adonis[2,5])
  permTable[13,] <- c("GUniFrac", PW.Adonis[3,1], PW.Adonis[3,4], PW.Adonis[3,5])
  
  table <- rbind(table,permTable)
  
  return(table)
}

summarize_permanova_aitchison <- function(d, collapse = NULL){
  
  library(ampvis2)
  library(vegan)
  library(tidyverse)
  
  d0 <- amp_filter_samples(d,visitId %in% c("Never smoker (control)","M0"))
  d2 <- amp_filter_samples(d,visitId %in% c("Never smoker (control)","M2"))
  d4 <- amp_filter_samples(d,visitId %in% c("Never smoker (control)","M4"))
  d6 <- amp_filter_samples(d,visitId %in% c("Never smoker (control)","M6"))
  
  if (!is.null(collapse)) {
    
    d0$abund <- aggregate_abund(d0$abund,d0$tax,tax_aggregate = collapse,tax_add = NULL,format = "abund",calcSums = FALSE)
    d2$abund <- aggregate_abund(d2$abund,d2$tax,tax_aggregate = collapse,tax_add = NULL,format = "abund",calcSums = FALSE)
    d4$abund <- aggregate_abund(d4$abund,d4$tax,tax_aggregate = collapse,tax_add = NULL,format = "abund",calcSums = FALSE)
    d6$abund <- aggregate_abund(d6$abund,d6$tax,tax_aggregate = collapse,tax_add = NULL,format = "abund",calcSums = FALSE)
    
  }
  
  ## distance measures without phylogeny #----
  
  list_in <- list(labels = c("Aitchison"),
                  decostand_methods = c("none"),
                  vegdist_methods = c("robust.aitchison"))
  list_out <- c()
  
  # helper function for data transformation
  transform_h <- function(dat,decostand_method){
    if (decostand_method == "none"){
      return(dat)
    } else {
      return(vegan::decostand(dat, method = decostand_method))
    }
  }
  
  for (i in 1:length(list_in$labels)){
    
    ## Make a table for saving PCOA + permanova testing results
    permTable <- as.data.frame(matrix(nrow=13, ncol=4))
    colnames(permTable) <- c("Distance","Comparison","P value","Adjusted p value")
    
    # M0
    res <- vegan::adonis2(vegan::vegdist(transform_h(t(d0$abund), list_in$decostand_methods[i]), method=list_in$vegdist_methods[i]) ~ group,
                          data=d0$metadata, 
                          permutations = 9999)
    permTable[1,] <- c(list_in$labels[i], "Smoker baseline vs Never smoker (control)", res$`Pr(>F)`[1], "NA")
    
    # M2
    res <- vegan::adonis2(vegan::vegdist(transform_h(t(d2$abund), list_in$decostand_methods[i]), method=list_in$vegdist_methods[i]) ~ group,
                          data=d2$metadata,
                          permutations = 9999)
    permTable[2,] <- c(list_in$labels[i], "M2 no smoking cessation vs M2 smoking cessation vs Never smoker (control)", res$`Pr(>F)`[1], "NA")
    
    # pair-wise pair-wise a posteriori tests M2
    PW.Adonis <- pairwise.adonis(transform_h(t(d2$abund),list_in$decostand_methods[i]),
                                 d2$metadata$group,
                                 sim.method=list_in$vegdist_methods[i],
                                 p.adjust.m = "bonferroni", 
                                 NULL)
    
    permTable[3,] <- c(list_in$labels[i], PW.Adonis[1,1], PW.Adonis[1,4], PW.Adonis[1,5])
    permTable[4,] <- c(list_in$labels[i], PW.Adonis[2,1], PW.Adonis[2,4], PW.Adonis[2,5])
    permTable[5,] <- c(list_in$labels[i], PW.Adonis[3,1], PW.Adonis[3,4], PW.Adonis[3,5])
    
    # M4
    res <- vegan::adonis2(vegan::vegdist(transform_h(t(d4$abund), list_in$decostand_methods[i]), method=list_in$vegdist_methods[i]) ~ group,
                          data=d4$metadata,
                          permutations = 9999)
    permTable[6,] <- c(list_in$labels[i], "M4 no smoking cessation vs M4 smoking cessation vs Never smoker (control)", res$`Pr(>F)`[1], "NA")
    
    # pair-wise pair-wise a posteriori tests M4
    PW.Adonis <- pairwise.adonis(transform_h(t(d4$abund), list_in$decostand_methods[i]),
                                 d4$metadata$group,
                                 sim.method=list_in$vegdist_methods[i],
                                 p.adjust.m = "bonferroni",
                                 NULL)
    
    permTable[7,] <- c(list_in$labels[i], PW.Adonis[1,1], PW.Adonis[1,4], PW.Adonis[1,5])
    permTable[8,] <- c(list_in$labels[i], PW.Adonis[2,1], PW.Adonis[2,4], PW.Adonis[2,5])
    permTable[9,] <- c(list_in$labels[i], PW.Adonis[3,1], PW.Adonis[3,4], PW.Adonis[3,5])
    
    # M6
    res <- vegan::adonis2(vegan::vegdist(transform_h(t(d6$abund), list_in$decostand_methods[i]), method=list_in$vegdist_methods[i]) ~ group,
                          data=d6$metadata,
                          permutations = 9999)
    permTable[10,] <- c(list_in$labels[i], "M6 no smoking cessation vs M6 smoking cessation vs Never smoker (control)", res$`Pr(>F)`[1], "NA")
    
    # pair-wise pair-wise a posteriori tests M6
    PW.Adonis <- pairwise.adonis(transform_h(t(d6$abund), list_in$decostand_methods[i]),
                                 d6$metadata$group,
                                 sim.method=list_in$vegdist_methods[i],
                                 p.adjust.m = "bonferroni",
                                 NULL)
    
    permTable[11,] <- c(list_in$labels[i], PW.Adonis[1,1], PW.Adonis[1,4], PW.Adonis[1,5])
    permTable[12,] <- c(list_in$labels[i], PW.Adonis[2,1], PW.Adonis[2,4], PW.Adonis[2,5])
    permTable[13,] <- c(list_in$labels[i], PW.Adonis[3,1], PW.Adonis[3,4], PW.Adonis[3,5])
    
    list_out[[i]] <- permTable
  }
  
  table <- do.call(rbind,list_out)
  
  return(table)
  
}
