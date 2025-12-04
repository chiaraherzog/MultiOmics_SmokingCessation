# plot relative abundances from a ASVtable (counts)
# Author: Charlotte Vavourakis
# version for S

plot_taxa_heatmap <- function(ASVtable, pheno, agg, add, abu, type, colors, ctrl=FALSE){
 
  require(ampvis2)
  require(tidyverse)
  require(ComplexHeatmap)
  require(rcartocolor)
  
  theme_set(theme_minimal(base_size=8))
  
  # load data for analysis with ampvis2, convert to rel abundance
  d = amp_load(otutable=ASVtable, metadata=pheno)
  d = normaliseTo100(d)
  
  # aggregate
  d = aggregate_abund(
    d$abund,
    d$tax,
    tax_aggregate = agg,
    tax_add = add,
    format = "abund",
    calcSums = FALSE)
  
  # abundance filter
  d = d %>% 
    mutate(threshold_bool = apply((dplyr::select(., ends_with(type))),1,function(x) +(any(x > abu)))) %>%
    filter(threshold_bool == 1) %>%
    dplyr::select(-threshold_bool) %>%
    droplevels()

  if (ctrl==FALSE){
    
    # annotation
    vars = c("compliance",
             "nic_replacement",
             "visitId")
    
    annot = pheno %>%
      dplyr::select(all_of(c("sampleId",vars))) %>%
      column_to_rownames(var="sampleId")
    
    annot$compliance <- factor(annot$compliance, levels = c("higher compliance","lower compliance"))
    annot$nic_replacement <- factor(annot$nic_replacement, 
                                    levels = c("e-cigarette","gum (up until at least last month before intervention)","heated tobacco product","plaster"))
    
    mycolors1 = cols[c(2, 1, 4, 5)]
    names(mycolors1) <- unique(annot$visitId)
    
    mycolors2 = cols[c(1,4)]
    names(mycolors2) <- c("higher compliance","lower compliance")
    
    mycolors3 = cols[c(1, 8,6,4)]
    names(mycolors3) <- c("e-cigarette","gum (up until at least last month before intervention)","heated tobacco product","plaster")
    
    mycolors <- list(nic_replacement = mycolors3,
                     compliance = mycolors2,
                     visitId = mycolors1)
    
  } else if (ctrl == TRUE) {
    
    # annotation
    vars = c("visitId","compliance")
    
    annot = pheno %>%
      dplyr::select(all_of(c("sampleId",vars))) %>%
      column_to_rownames(var="sampleId")
    
    annot$visitId <- factor(annot$visitId, levels = c("Never smoker (control)","M0","M2","M4","M6"))
    annot$compliance <- factor(annot$compliance, levels = c("Never smoker (control)","higher compliance","lower compliance"))
    
    mycolors1 = cols[c(6, 2, 1, 4, 5)]
    names(mycolors1) <- unique(annot$visitId)
    
    mycolors2 = cols[c(6, 1,4)]
    names(mycolors2) <- c("Never smoker (control)","higher compliance","lower compliance")
    
    mycolors <- list(compliance = mycolors2,
                     visitId = mycolors1)
  }
  
  # heatmap
  plot = pheatmap(as.matrix(d), 
                  annotation_col=annot,
                  annotation_colors = mycolors,
                  cluster_cols = FALSE,
                  show_colnames=FALSE,
                  border_color = NA)
  
  
  out <- list(plot = plot,
              data = d, 
              annot = annot)
  
  return(out)
}