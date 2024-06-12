#' @name checkASV
#' @description
#' Based on rmcorr table, append detailed information to ASVs
#' @param variable Feature to check (grepl), note: ensure this is unique
#' @param corr correlation dataframe
#' @return data frame with appended information

checkASV <- function(variable, corr){
  
  tmp1 <- corr[grepl(variable, corr$measure1) | grepl(variable, corr$measure2),] |> 
    dplyr::mutate(feature = ifelse(grepl(variable, measure1), measure2, measure1)) |> 
    tidyr::separate(feature, "_", into = c(NA, "feature"), extra = 'merge') |> 
    dplyr::filter(grepl("ASV|famil", measure1) | grepl("ASV|famil", measure2))
  
  ASVtable_saliva <- readRDS("~/Dropbox/tg-data-prep/prep-scripts/11-preprocessed/2-microbiome/6-study-arm-subsets/smk1/0-output/ASVtable_saliva_S.Rds")
  ASVtable_stool <- readRDS("~/Dropbox/tg-data-prep/prep-scripts/11-preprocessed/2-microbiome/6-study-arm-subsets/smk1/0-output/ASVtable_stool_S.Rds")
  ASV <- plyr::rbind.fill(ASVtable_saliva, ASVtable_stool)
  
  tmp2 <- ASV[match(x$feature, ASV$OTU),] |> 
    dplyr::select(OTU, Phylum:Species)
  
  tmp <- tmp1 |> 
    dplyr::left_join(tmp2, by = c('feature' = 'OTU'))
  
  return(tmp)
  
}
