#' @param exp which experiment should be plotted
#' @param lmm_data_time Model that provides time coefficient (out_lmm$`Higher compliance only` or out_lmm$`Higher compliance only with packyears`)
#' @param cols Colour palette (leave blank for default colours [recommended], but can be altered if desired)
#' @param relabel (optional). Data frame supplying labels for the plot. x = variable, label = relabelled variable
#' @param age_cor should correlation with age be plotted on the side? default to F
#' @param spy_cor should correlation with smoking pack years be plotted on the side? default to F
#' @param FDR should p-values be corrected? default to F
#' @param cluster Should rows be clustered? If left blank, will set to 'default' and rows will be ordered by M6 time estimate. If 'cluster' (or other), rows will be clustered.
#' @return ComplexHeatmap object
#' 
# version made for microbiomes, smoking study arm
# option to plot age and/or smoking pack years correlations at baseline
# option to FDR correct p-values, for each term (and category) separately 
# all terms with p.corrected < 0.05 are plotted.
# for uncorrected p-values, all rows with p < 0.01 are plotted

# main function
plot_lmm_heatmap_time <- function(exp,
                                  lmm_data_time,
                                  cols = c("#1b69a1",
                                           "#48a0af",
                                           "#f39668",
                                           "#ec6669"),
                                  relabel = NULL,
                                  age_cor = F,
                                  spy_cor = F,
                                  FDR = F,
                                  cluster = 'default'){
  
  require(gtools)
  require(stringr)
  
  source("src/FDR_correct.R")
  
  ## load output lmm models and correct p-values if requested, keep only significant results ##----
  
  if (FDR == F) {
    
    tmp <- lmm_data_time %>%
      filter(str_detect(x, exp)) %>%
      dplyr::select(!contains(c("std.", "age_at_consent", "smoking_py"))) %>% 
      dplyr::select(!ends_with(c("compliance")))
    
    tmp$x = gsub(paste0(exp,"_"),"",tmp$x)
    
    # filter out only significant rows and replace pval with * symbols
    tmp <- tmp %>%
      filter(rowSums(across(starts_with("p.value"), ~ . < 0.01), na.rm = TRUE) > 0) %>% 
      dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), "", gtools::stars.pval(.)))) |> 
      dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", "", .)))
    
  } else if (FDR == T){
    
    # time model FDR correct
    tmp <- lmm_data_time %>%
      filter(str_detect(x, exp)) %>%
      dplyr::select(!contains(c("std.", "age_at_consent", "smoking_py"))) %>% 
      dplyr::select(!ends_with(c("compliance")))
    
    tmp_p <- tmp %>%
      dplyr::select(x,contains("p.value")) %>%
      pivot_longer(cols = contains("p.value"),names_to = "estimate", values_to = "p_value") %>%
      FDR_correct() %>%
      pivot_wider(id_cols = x, names_from = estimate, values_from = p_value)
    
    tmp <- tmp %>%
      dplyr::select(!contains("p.value")) %>%
      left_join(tmp_p, by = "x") %>%
      dplyr::select(!contains(c("std.", "age_at_consent", "smoking_py")))
    
    tmp$x = gsub(paste0(exp,"_"),"",tmp$x)
    
    # filter out only significant rows and replace pval with * symbols
    tmp <- tmp %>%
      filter(rowSums(across(starts_with("p.value"), ~ . < 0.05)) > 0, na.rm = TRUE) %>% 
      dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), "", gtools::stars.pval(.)))) |> 
      dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", "", .)))
  }
  
  ## Column labels ##---- 
  col_labs <- c("M2", "M4", "M6")
  
  # Filter
  
  if(!is.null(relabel)){
    
    tmp <- tmp |> 
      # relabel
      dplyr::filter(x %in% relabel$x)
  }
    
  
  ## baseline correlations ##----
  
  if (age_cor == T & spy_cor == F) {
    
    load("out/corrAgeBaseline.R")
    corr <- corr |> 
      tidyr::separate(rowname, "_", into = c("assay", "name"), extra = 'merge') |> 
      dplyr::slice(match(tmp$x, name))
    
    tmp <- tmp |> 
      dplyr::left_join(dplyr::select(corr, name, cor, p), by = c('x' = 'name'))
    
    hr <- rowAnnotation('Correlation with age\n(baseline)' = anno_points(tmp$cor_age,
                                                                         pch = ifelse(tmp$p < 0.05, 19, 1),
                                                                         ylim = c(-1, 1),
                                                                         gp = gpar(col = ifelse(tmp$cor<0, cols[1], cols[4]))))
  } else if (age_cor == F & spy_cor == T) {
    load("out/corrCigBaseline.R")
    corr <- corr |> 
      tidyr::separate(rowname, "_", into = c("assay", "name"), extra = 'merge') |> 
      dplyr::slice(match(tmp$x, name))
    
    tmp <- tmp |> 
      dplyr::left_join(dplyr::select(corr, name, cor, p), by = c('x' = 'name'))
    
    hr <- rowAnnotation('Correlation with spy\n(baseline)' = anno_points(tmp$cor_spy,
                                                                         pch = ifelse(tmp$p < 0.05, 19, 1),
                                                                         ylim = c(-1, 1),
                                                                         gp = gpar(col = ifelse(tmp$cor<0, cols[1], cols[4]))))
  } else if (age_cor == T & spy_cor == T) {
    load("out/corrAgeBaseline.R")
    corr <- corr |> 
      tidyr::separate(rowname, "_", into = c("assay", "name"), extra = 'merge') |> 
      dplyr::slice(match(tmp$x, name)) %>%
      dplyr::rename(corr_age = cor) %>%
      dplyr::rename(p_age = p)
    
    tmp <- tmp |> 
      dplyr::left_join(dplyr::select(corr, name, corr_age, p_age), by = c('x' = 'name'))
    
    load("out/corrCigBaseline.R")
    corr <- corr |> 
      tidyr::separate(rowname, "_", into = c("assay", "name"), extra = 'merge') |> 
      dplyr::slice(match(tmp$x, name)) %>%
      dplyr::rename(corr_spy = cor) %>%
      dplyr::rename(p_spy = p)
    
    tmp <- tmp |> 
      dplyr::left_join(dplyr::select(corr, name, corr_spy, p_spy), by = c('x' = 'name'))
    
    hr <- rowAnnotation('Age' = anno_points(tmp$corr_age,
                                            pch = ifelse(tmp$p_age < 0.05, 19, 1),
                                            ylim = c(-1, 1),
                                            gp = gpar(col = ifelse(tmp$corr_age<0, cols[1], cols[4]))),
                        
                        'Spy' = anno_points(tmp$corr_spy,
                                            pch = ifelse(tmp$p_spy < 0.05, 19, 1),
                                            ylim = c(-1, 1),
                                            gp = gpar(col = ifelse(tmp$corr_spy<0, cols[1], cols[4]))))
  } else {
    hr = NULL
  }
  
  # relabel
  if(!is.null(relabel)){
    
    tmp <- tmp |> 
      # relabel
      dplyr::left_join(relabel) |> 
      dplyr::mutate(x = label)
  }
  
  
  ## Draw heatmap ##---- 
  
  # Indices
  ind_est <- grep("estimate", colnames(tmp))  
  ind_p <- grep("p.value", colnames(tmp))
  
  
  if(cluster == 'default'){
    p <- Heatmap(as.matrix(tmp[,ind_est]),
                 name = 'estimate\n(scaled)',
                 
                 # Row details
                 row_labels = tmp$x,
                 row_names_side = 'left',
                 row_order = order(tmp$estimate_visitIdM6, decreasing = T),
                 show_row_dend = F, 
                 row_title = NULL,
                 
                 # Row annotation
                 right_annotation = hr,
                 
                 # Column details
                 column_labels = col_labs,
                 column_title_gp = grid::gpar(fontsize = 10),
                 cluster_columns = F,
                 show_column_dend = F,
                 
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
                                                                        "grey95", cols[c(3, 4)]))(30)),
                 
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
                 row_labels = tmp$x,
                 row_names_side = 'left',
                 cluster_rows = T,
                 cluster_row_slices = F,
                 show_row_dend = F, 
                 row_title = NULL,
                 
                 # Row annotation
                 right_annotation = hr,
                 
                 # Column details
                 column_labels = col_labs,
                 column_title_gp = grid::gpar(fontsize = 10),
                 cluster_columns = F,
                 show_column_dend = F,
                 
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
  }
  
  
  out <- list(plot = p,
              full_df = tmp)
  return(out)
  
}
