prep_data_hima <- function(count_dat, pheno, eps = 1e-12) {
  
  # Relative abundance
  dat_rel = sweep(count_dat, 1, rowSums(count_dat), FUN = "/")
  
  # Merge rel abundance data with phenotype
  df = cbind(pheno, dat_rel[pheno$sampleId, ])  # order rows by pheno
  
  # Baseline changes
  df_m0 = df[df$visitId == "M0", ]
  df_m2 = df[df$visitId == "M2", ]
  df_m4 = df[df$visitId == "M4", ]
  df_m6 = df[df$visitId == "M6", ]
  
  common_subj = Reduce(
    intersect,
    list(
      df_m0$subjectId,
      df_m2$subjectId,
      df_m4$subjectId,
      df_m6$subjectId
    )
  )
  
  df_m2 = df_m2[match(df_m0$subjectId, df_m2$subjectId), ]
  df_m4 = df_m4[match(df_m0$subjectId, df_m4$subjectId), ]
  df_m6 = df_m6[match(df_m0$subjectId, df_m6$subjectId), ]
  
  # fold change: M_month / M_baseline
  fc_m2 = (df_m2[, -(1:19)] + eps) / (df_m0[, -(1:19)] + eps) + eps
  fc_m4 = (df_m4[, -(1:19)] + eps) / (df_m0[, -(1:19)] + eps) + eps
  fc_m6 = (df_m6[, -(1:19)] + eps) / (df_m0[, -(1:19)] + eps) + eps
  
  rownames(fc_m2) <- df_m0$subjectId
  rownames(fc_m4) <- df_m0$subjectId
  rownames(fc_m6) <- df_m0$subjectId
  
  out <- list(fc_m2, fc_m4, fc_m6)
  names(out) <- c("fc_m2", "fc_m4", "fc_m6")
  
  return(out)
  
}

run_hima_microbiome <- function(dir_out,
                                time_label,
                                tax_level,
                                cessation_time,
                                microbiome_dat,
                                metabolites_dat,
                                mtb_interest,
                                pheno_dat,
                                form_covar,
                                sigcut = 0.05) {
  
  # Create output directory if needed
  if(!dir.exists(dir_out)){
    dir.create(dir_out, recursive = T)
  }
  
  # Filter pheno and outcome data by time
  pheno_dat = pheno_dat %>% 
    dplyr::filter(visitId == cessation_time) %>% # cessation at Mx
    full_join(metabolites_dat, by = "subjectId") # metabolites at Mx
  
  # Order everything
  microbiome_dat = microbiome_dat[match(pheno_dat$subjectId, rownames(microbiome_dat)), ]
  if (!identical(pheno_dat$subjectId, rownames(microbiome_dat))) {
    stop("Something went wrong: cannot match subjectId in pheno_dat to rownames in microbiome_dat.")
  }
  
  microbiome_dat = as.matrix(microbiome_dat)
  
  # Loop over outcome vars and run hima_microbiome
  for (i in 1:length(mtb_interest)) {
    
    cat(paste0("Running ", time_label, ", ",tax_level,", ", mtb_interest[i],"\n"))
    
    # Create regression formula
    form = as.formula(paste(mtb_interest[i], form_covar))
    
    # Filter pheno so only has covars of interest
    pheno_x = pheno_dat %>% dplyr::select(all.vars(form))
    
    # Y = metabolite outcome, M = matrix of 10 CLR ASVs, X = Quit
    out_hima = hima(
      form,
      data.pheno = pheno_x,
      data.M = microbiome_dat,
      mediator.type = "compositional",
      penalty = "DBlasso",
      sigcut = sigcut,
      verbose = T
    )
    saveRDS(out_hima, file = file.path(dir_out, paste0("out_hima_",time_label,"_", tax_level, "_", mtb_interest[i],"_FWER", sigcut,".Rds")))
    
  }
}
