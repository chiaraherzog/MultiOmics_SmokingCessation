#' @name plotSMKReplacement_T
#' @param dat data frame (eg change)
#' @param variable Variable to plot
#' @param metadata metadata object as data.frame (with t)
#' @param ylab Y label

plotSMKReplacement_T <- function(dat, variable, metadata,
                     ylab = ''){
  
  dat <- dat |> 
    dplyr::left_join(metadata, by = 'colname') |> 
    dplyr::group_by(subjectId) |> 
    dplyr::mutate(t_relative = ifelse(is.na(cessation_date), t-max(t), 
                                      t-cessation_date)) |> 
    dplyr::ungroup()
  
  p <- dat |> 
    dplyr::filter(rowname == variable) |> 
    ggplot(aes(x =  t_relative,
               y = value)) +
    geom_vline(xintercept = 0,
               linetype = 'solid',
               colour = cols[4]) +
    geom_line(aes(group = subjectId,
                  colour = nic_replacement),
              alpha = 0.7,
              linewidth = 0.7) +
    geom_hline(yintercept = 0,
               linetype = 'dashed',
               colour = 'grey60') +
    scale_colour_manual(values = c(cols[c(5, 8, 6, 2)], "NA"),
                        na.value = 'grey80',
                        name = 'Nicotine\nreplacement',
                        breaks = c("e-cigarette", "gum (up until at least last month before intervention)", "heated tobacco product", "plaster", NA),
                        labels = c("e-cigarette", "gum", "heated tobacco", "plaster", "none reported")) +
    theme_bw() +
    theme(axis.title.y = ggtext::element_markdown(),
          panel.grid.major = element_blank()) +
    labs(x = 'days from cessation date',
         y = ylab)
  
  return(p)
  
  }