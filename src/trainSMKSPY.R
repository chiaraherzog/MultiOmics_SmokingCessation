trainSMKSPY <- function(data, exp, populations,
                     mainOnly = F){
  
  # Extract data ------
  tr <- as.data.frame(wideFormat(data[,data@colData$visitId == 'M0',
                                      exp],
                                 colData = c('smoking_py'))) |> tidyr::drop_na() 
  
  colnames(tr) <- gsub("Flow.cytometry..T.cell.staining_|Flow.cytometry..white.blood.cell.staining_", "", colnames(tr))
  
  trdata <- tr |> dplyr::select(-c('smoking_py', 'primary'))
  
  if(mainOnly == T){
    trdata <- trdata |> 
      dplyr::select(any_of(dplyr::filter(populations, main_analysis == 'yes')$name))
  }
  
  trpheno <- tr |> dplyr::select(c('smoking_py'))
  
  # We train the immAge using the agePred function and select maximum correlation with age: ------
  models <- list("elnet" = 0.5,
                 "lasso" = 1)
  
  res <- lapply(models, function(x){
    smkPredictor(data_tr = t(trdata),
                 cig_tr = as.numeric(trpheno$smoking_py),
                 data_val = t(trdata),
                 cig_val = as.numeric(trpheno$smoking_py),
                 alpha = x)
  })
  
  cor <- lapply(names(res), function(m){
    cor(trpheno$smoking_py, res[[m]]$val_predictor)
  })
  
  keep <- which.max(cor)
  cat(names(res)[keep], ' model has the largest cor\n', sep = '')
  
  res <- res[[keep]]
  intercept <- as.numeric(coef(res$fit.cv, s = "lambda.min")[1]) 
  
  coefs <- as.matrix(coef(res$fit.cv, s = "lambda.min")) |> 
    as.data.frame() |> 
    tibble::rownames_to_column('name') |> 
    dplyr::filter(s1 != 0 & name != '(Intercept)') |> 
    dplyr::arrange(desc(abs(s1)))
  
  dt <- DT::datatable(coefs |> 
                        dplyr::left_join(dplyr::select(populations, name, `population name`)) |> 
                        dplyr::select(`population name`, s1))
  
  # Compute values in the full dataset -------
  full <- as.data.frame(wideFormat(data[,,exp],
                                   colData = c('interventionId', 'subjectId', 'visitId', 'compliance', 'age_at_consent', 'bmi_at_consent', 'smoking_ever', 'cig_curr', 'smoking_py'))) |> 
    dplyr::mutate(smoking_current = ifelse(interventionId == 'S', 'yes', 'no'),
                  smoking_current = factor(smoking_current, levels = c('no', 'yes')))
  
  colnames(full) <- gsub("Flow.cytometry..T.cell.staining_|Flow.cytometry..white.blood.cell.staining_", "", colnames(full))
  
  dat_fu <- full |> 
    dplyr::select(coefs$name) |> 
    as.matrix()
  
  if(!identical(colnames(dat_fu), coefs$name)){
    stop("colnames of full dataset not identical with coefs. please investigate.\n")
  }
  
  full$ImmuneSMKSPY <- apply(dat_fu, 1, function(x){
    intercept + sum(x * coefs$s1, na.rm = T)
  })
  
  # Initial plots: ---------
  a <- full |> 
    dplyr::filter(visitId == 'M0') |> 
    ggplot(aes(x = smoking_current,
               y = ImmuneSMKSPY)) +
    geom_boxplot()
  
  b <- full |> 
    dplyr::filter(visitId == 'M0') |> 
    ggplot(aes(x = as.numeric(cig_curr),
               y = ImmuneSMKSPY)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    ggpubr::stat_cor()
  
  c <- full |> 
    dplyr::filter(visitId == 'M0') |> 
    ggplot(aes(x = smoking_py,
               y = ImmuneSMKSPY)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    ggpubr::stat_cor()
  
  d <- full |> 
    dplyr::filter(visitId == 'M0') |> 
    ggplot(aes(x = smoking_ever,
               y = ImmuneSMKSPY,
               fill = interventionId)) +
    geom_boxplot()
  
  # Follow-up visits 
  follow_all <- full |> 
    dplyr::mutate(compliance = factor(compliance, levels = c("low", 'medium', 'high', 'lower compliance', 'higher compliance'))) |> 
    ggplot(aes(x = visitId,
               y = ImmuneSMKSPY)) +
    geom_boxplot() +
    geom_line(aes(group = subjectId),
              alpha = 0.1) +
    facet_wrap(interventionId~compliance) +
    ggpubr::stat_compare_means(ref.group = 'M0',
                               label = 'p.signif',
                               label.y.npc = 0.95,
                               paired = F) +
    labs(subtitle = 'Raw values over time')
  
  # Follow-up visits (complete)
  follow_paired <- full |> 
    dplyr::mutate(compliance = factor(compliance, levels = c("low", 'medium', 'high', 'lower compliance', 'higher compliance'))) |> 
    tidyr::pivot_wider(id_cols = c(subjectId, compliance, smoking_py, smoking_current, smoking_ever, cig_curr, interventionId),
                       names_from = visitId,
                       values_from = ImmuneSMKSPY) |> 
    dplyr::rowwise() |> 
    dplyr::mutate(M2 = M2-M0,
                  M4 = M4-M0,
                  M6 = M6-M0,
                  M0 = 0) |>
    dplyr::ungroup() |> 
    dplyr::filter(!is.na(M0) & !is.na(M2) & !is.na(M6) & !is.na(M4)) |> 
    tidyr::pivot_longer(cols = M0:M6,
                        names_to = 'visitId',
                        values_to = 'ImmuneSMKSPY') |> 
    ggplot(aes(x = visitId,
               y = ImmuneSMKSPY)) +
    geom_boxplot() +
    geom_line(aes(group = subjectId),
              alpha = 0.1) +
    facet_wrap(interventionId ~ compliance) +
    ggpubr::stat_compare_means(ref.group = 'M0',
                               label = 'p.signif', paired = T,
                               label.y.npc = 0.95) +
    labs(y = '∆ ImmuneSMK\nfrom baseline')
  
  out <- list(data = full,
              res = res,
              coef = dt,
              coef_raw = coefs,
              plots = list("smoking_current" = a,
                           "cig_curr_corr" = b,
                           "spy_corr" = c,
                           "boxplot_intervention" = d,
                           "follow_raw" = follow_all,
                           "follow_paired" = follow_paired))
  
  return(out)
}
