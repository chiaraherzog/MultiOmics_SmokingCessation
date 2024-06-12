#' @param lmm_data_time Minimal model that provides time coefficient (no interaction) (out_lmm$`Minimal model`)
#' @param lmm_data_compliance Basic model that provides time*compliance interaction (out_lmm$`Basic model with BMI`)
#' @param cols Colour palette (leave blank for default colours [recommended], but can be altered if desired)
#' @param relabel (optional). Data frame supplying labels for the plot. x = variable, label = relabelled variable
#' @param filter_relabel (optional) boolean: if supplying new labels, should only variables provided in that dataframe be plotted? default to TRUE
#' @param relabel_assay (optional). Set to TRUE if you would like to rename the assay. Assay names should be provided in the relabel parameter, under column name 'assay2'
#' @param colour_assays (optional) if not NULL, list of default colours will be used for assays. This only works for clinical variables so far.
#' @param methyl_sub default to NULL; filter for methylation data of one tissue only (cervical, buccal, or blood)
#' @param m6_sep Should M6 values be plotted separately to avoid blank spots? default to TRUE
#' @param age_cor should correlation with age be plotted on the side? default to F
#' @param mark_age_cor description
#' @param cluster Should rows be clustered? If left blank, will set to 'default' and rows will be ordered by M6 time estimate. If 'cluster' (or other), rows will be clustered.
#' @return ComplexHeatmap object

# main function
plot_lmm_heatmap_frag_v2_metab <- function(lmm_data_time,
                                           lmm_data_compliance,
                                           sampletype='urine',
                                           cols = c("#1b69a1",
                                                    "#48a0af",
                                                    "#f39668",
                                                    "#ec6669"),
                                           relabel = NULL,
                                           filter_relabel = T,
                                           colour_assays = list(assay = c("Blood haemogram" = "#ec6669",
                                                                          "Body weight and composition" = "#832c9b",
                                                                          "Spirometry" = "#bd647d",
                                                                          "Functional exercise capacity" = "#f39668",
                                                                          "Blood test" = '#71b5a9',
                                                                          "Vascular features" = '#0d49a1')),
                                           relabel_assay = NULL,
                                           methyl_sub = NULL,
                                           m6_sep = T,
                                           age_cor = F,
                                           spy_cor = T,
                                           mark_age_cor = F,
                                           cluster = 'default'){
  
  if(!require("gtools")){
    install.packages("gtools")
  }
  
  tmp1 <- lmm_data_time |> 
    dplyr::select(!contains(c("std.", "age_at_consent", "smoking_py"))) 
    # # set pval
    # dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), " ", gtools::stars.pval(.)))) |> 
    # dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", " ", .)))
  
  tmp2 <- lmm_data_compliance |> 
    dplyr::select(!contains(c("std",
                              "smoking_py", 
                              "age_at_consent",
                              "estimate_compliancehigher compliance", "p.value_compliancehigher compliance"))) |> 
    dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6"))) 
    
    # # set pval
    # dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), " ", gtools::stars.pval(.)))) |> 
    # dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", " ", .)))
  
  

  
  tmp <- tmp1 |> 
    dplyr::left_join(tmp2) 
  
  tmp <-  tmp%>% 
    filter(rowSums(across(starts_with("p.value"), ~ . < 0.05), na.rm = TRUE) > 0) %>% 
    dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), "", gtools::stars.pval(.)))) |> 
    dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", "", .)))
  

    
    tmp <- tmp |> 
      tidyr::separate(x, into = c('assay', "x"),
                      sep = "_", 
                      extra = 'merge')
    

    # Load the correct reference map  
    #
    if (sampletype=='urine') {
      RefMet <-  readRDS("out/RefMet_mapped_urine.Rds")
      exp <- "Urine nuclear magnetic resonance: normalized_"
    } else if (sampletype=='saliva') {
      RefMet <-  readRDS("out/RefMet_mapped_saliva.Rds")
      exp <- "Saliva nuclear magnetic resonance: normalized_"
    }
  

    
    if(age_cor == T){
      load("out/corrAgeBaseline.R") 
      corr <- corr |> 
        dplyr::filter(grepl(exp, rowname)) |> 
        dplyr::mutate(rowname = gsub(exp, "", rowname))
      
      tmp <- tmp |> 
        dplyr::left_join(dplyr::select(corr, rowname, cor, p), by = c('assay' = 'rowname'))
    }
    
    if(spy_cor == T){
      load("out/corrCigBaseline.R")
      corr <- corr |> 
        dplyr::filter(grepl(exp, rowname)) |> 
        dplyr::mutate(rowname = gsub(exp, "", rowname)) |> 
        dplyr::rename(cor.c = cor,
                      p.c = p)
      
      tmp <- tmp |> 
        dplyr::left_join(dplyr::select(corr, rowname, cor.c, p.c), by = c('assay' = 'rowname'))
    }
    


  

  
  if(age_cor == T & spy_cor == T){
    hr <- rowAnnotation('Correlation with age\n(baseline)' = anno_points(tmp$cor,
                                                                         pch = ifelse(tmp$p < 0.05, 19, 1),
                                                                         ylim = c(-1, 1),
                                                                         gp = gpar(col = ifelse(tmp$cor<0, cols[1], cols[4]))),
                        
                        'Opposite change\nwith\nintervention' = anno_text(ifelse(sign(tmp$cor) != sign(tmp$estimate_visitIdM6) & tmp$p.value_visitIdM6!=" ",
                                                                                 "←", "")),
                        'Correlation with spy\n(baseline)' = anno_points(tmp$cor.c,
                                                                         pch = ifelse(tmp$p.c < 0.05, 19, 1),
                                                                         ylim = c(-1, 1),
                                                                         gp = gpar(col = ifelse(tmp$cor.c<0, cols[1], cols[4]))),
                        'Opposite change\nwith\ncessation' = anno_text(ifelse(sign(tmp$cor.c) != sign(tmp$`estimate_visitIdM6:compliancehigher compliance`) & tmp$`p.value_visitIdM6:compliancehigher compliance`!= " ",
                                                                              "←", "")))
  } else if(age_cor == T & spy_cor != T){
    hr <- rowAnnotation('Correlation with age\n(baseline)' = anno_points(tmp$cor,
                                                                         pch = ifelse(tmp$p < 0.05, 19, 1),
                                                                         ylim = c(-1, 1),
                                                                         gp = gpar(col = ifelse(tmp$cor<0, cols[1], cols[4]))),
                        
                        'Opposite change\nwith\nintervention' = anno_text(ifelse(sign(tmp$cor) != sign(tmp$estimate_visitIdM6) & tmp$p.value_visitIdM6!=" ",
                                                                                 "←", "")))
    
  } else if(age_cor != T & spy_cor == T){
    hr <- rowAnnotation('Correlation with spy\n(baseline)' = anno_points(tmp$cor.c,
                                                                         pch = ifelse(tmp$p.c < 0.05, 19, 1),
                                                                         ylim = c(-1, 1),
                                                                         gp = gpar(col = ifelse(tmp$cor.c<0, cols[1], cols[4]))),
                        'Opposite change\nwith\ncessation' = anno_text(ifelse(sign(tmp$cor.c) != sign(tmp$`estimate_visitIdM6:compliancehigher compliance`) & tmp$`p.value_visitIdM6:compliancehigher compliance`!= " ",
                                                                              "←", "")))
  } else {
    hr = NULL
  }
    
    ## Class annotation ##----
    
    # Harmonize input names
    tmp$assay <- gsub("^X(\\d+)", "\\1", tmp$assay)
    
    # Metabolites that are unknown to RefMet - try more standardized name
    tmp$assay[tmp$assay=="Acetate.mM."] <- "Acetate"
    tmp$assay[tmp$assay=="Methyl.2.oxovalerate"] <- "3-Methyl-2-oxovaleric acid"
    tmp$assay[tmp$assay=="Aminopentanoate"] <- "5-Aminopentanoate"
    tmp$assay[tmp$assay=="TMA..N.oxide"] <- "Trimethylamine N-oxide"

    
    
    
    # Replace all "-" entries with "unknown" in merged_data dataframe
    # RefMet <- RefMet %>%
    #   mutate_all(~replace(., . == NA, "unknown"))
    
    
    tmp <- merge(RefMet, tmp, by.x = "Input.name", by.y = "assay", all.x = F, all.y = TRUE)
    tmp$assay <- tmp$Input.name
    
    
    cols_class <- c("Nicotinic acid alkaloids" = "#1b69a1", "Amino acids and peptides" = "#48a0af", "-" = "#71b5a9", 
                    "Tryptophan alkaloids" = "#ec6669", "Phenolic acids" = "#f39668", "Fatty acids" = "#bd647d", 
                    "Carbonyl compounds" = "#395262", "Monosaccharides" = "#5f70a8", "Cholines" = "#d4a34b", "TCA acids" = "#4da457", 
                    "Azoles" = "#d16fa8", "Organonitrogen compounds" = "#78b646", "Carboxylic acids" = "#3b738f", 
                    "Short-chain acids" = "#a36f2d", "Sulfonic acids" = "#d4a34b", "Amines" = "#6e2b57", "Amine oxides" = "#c97f2e", 
                    "Organic carbonic acids" = "#8e4555", "Hydrocarbons" = "#4f6d32", "Phenylpropanoids" = "#ec6669", 
                    "Primary alcohols" = "#9c572c", "Alcohols and polyols" = "#5a9d79", "Phenols" = "#d09636", 
                    "Fatty amines" = "#7d3b65", "Pyrimidines" = "#1b69a1", "Purines" = "#48a0af", "Keto acids"="#832c9b")
    
    
    ha = rowAnnotation(Main.class = tmp$Main.class,
                           col = list(Main.class = cols_class),
                           show_annotation_name = c(Main.class = FALSE))
    
    

    
    
    # Indices
    ind_est <- grep("estimate", colnames(tmp))  
    ind_p <- grep("p.value", colnames(tmp))
    
    # Column splits and labels
    col_split <- rep(c("time", "high\ncompliance"), each = 3)
    col_split <- factor(col_split, levels = c("time", "high\ncompliance"))
    col_labs <- rep(c("M2", "M4", "M6"), 2)
    
    
    # Row splits and labels
    if(m6_sep == T){
      row_split = ifelse(is.na(tmp$`estimate_visitIdM2`), 2, 1)
    } else {
      row_split = rep(1, nrow(tmp))
    }
    
    row_labs <- ifelse(tmp$`p.value_visitIdM6:compliancehigher compliance` == "*" |
                         tmp$`p.value_visitIdM6:compliancehigher compliance` == "**" |
                         tmp$`p.value_visitIdM6:compliancehigher compliance` == "***" |
                         tmp$`p.value_visitIdM6:compliancehigher compliance` == ".",
                       paste0("**", tmp$assay, "**"),
                       tmp$assay)
    

  
  if(cluster == 'default'){
    p <- Heatmap(as.matrix(tmp[,ind_est]),
                 name = 'estimate\n(scaled)',
                 
                 # Row details
                 row_labels = gt_render(row_labs),
                 row_names_side = 'left',
                 row_split = row_split,
                 row_order = order(tmp$estimate_visitIdM6, decreasing = T),
                 show_row_dend = F, 
                 row_title = NULL,
                 
                 # Row annotation
                 left_annotation = ha,
                 right_annotation = hr,
                 
                 # Column details
                 column_labels = col_labs,
                 column_title_gp = grid::gpar(fontsize = 10),
                 column_split = col_split,
                 cluster_columns = F,
                 show_column_dend = F,
                 
                 # Row side colours
                 
                 
                 # Annotation
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(as.matrix(tmp[,ind_p])[i, j], x, y, gp = gpar(fontsize = 7), just = "centre")
                 },
                 
                 # Colours
                 na_col = 'white',
                 col = circlize::colorRamp2(breaks = seq(-max(abs(na.omit(tmp[,ind_est]))),
                                                         max(abs(na.omit(tmp[,ind_est]))),
                                                         length.out = 30),
                                            colors = colorRampPalette(c(cols[c(1, 2)],
                                                                        "grey95", cols[c(3, 4)]))(30)
                                            #                            # colors = rev(color("batlow")(5)),
                                            #                            # colors = rev(cols[c(2,5,4, 6, 7)])
                 ),
                 
                 # Titles
                 row_names_gp = grid::gpar(fontsize = 9),
                 column_names_gp = grid::gpar(fontsize = 9),
                 
                 # Borders
                 border_gp = gpar(lwd = 0.5),
                 border = T)  
  } else {
    p <- Heatmap(as.matrix(tmp[,ind_est]),
                 name = 'estimate\n(scaled)',
                 
                 # Row details
                 row_labels = gt_render(row_labs),
                 row_names_side = 'left',
                 row_split = row_split,
                 cluster_rows = T,
                 cluster_row_slices = F,
                 show_row_dend = F, 
                 row_title = NULL,
                 
                 # Row annotation
                 left_annotation = ha,
                 right_annotation = hr,
                 
                 # Column details
                 column_labels = col_labs,
                 column_title_gp = grid::gpar(fontsize = 10),
                 column_split = col_split,
                 cluster_columns = F,
                 show_column_dend = F,
                 
                 # Row side colours
                 
                 
                 # Annotation
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(as.matrix(tmp[,ind_p])[i, j], x, y, gp = gpar(fontsize = 7), just = "centre")
                 },
                 
                 # Colours
                 na_col = 'white',
                 col = circlize::colorRamp2(breaks = seq(-max(abs(na.omit(tmp[,ind_est]))),
                                                         max(abs(na.omit(tmp[,ind_est]))),
                                                         length.out = 30),
                                            colors = colorRampPalette(c(cols[c(1, 2)],
                                                                        "grey95", cols[c(3, 4)]))(30)
                                            #                            # colors = rev(color("batlow")(5)),
                                            #                            # colors = rev(cols[c(2,5,4, 6, 7)])
                 ),
                 
                 # Titles
                 row_names_gp = grid::gpar(fontsize = 9),
                 column_names_gp = grid::gpar(fontsize = 9),
                 
                 # Borders
                 border_gp = gpar(lwd = 0.5),
                 border = T)  
    
    if(age_cor == T){
      p <- print(p) + decorate_annotation("Correlation with age\n(baseline)", {
        grid.lines(c(0.5), gp = gpar(col = "grey80",lty = 'dotted'))
      })
    }
    
  }
  
  return(p)
  
}