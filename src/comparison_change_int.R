comparison_change_int <- function(dat, variable = 'variable',
                                  ylab = '',
                                  p = "p.signif",
                                  colours = cols[c(4, 1)]){
    
  complete <- dat |> 
    dplyr::filter(rowname == variable) |> 
    dplyr::filter(!is.na(value)) |> 
    dplyr::group_by(subjectId) |> 
    dplyr::filter(n() == 4)
    
    plot <- dat |> 
      dplyr::filter(visitId != "M0" & subjectId %in% complete$subjectId) |> 
      dplyr::filter(rowname == variable) |> 
      ggplot(aes(x = compliance,
                 y = value)) +
      geom_boxplot(alpha = 0.3,
                   aes(fill = compliance)) +
      ggbeeswarm::geom_beeswarm(aes(colour= compliance),
                                size = 0.9, alpha = 0.6) +
      ggpubr::stat_compare_means(
        comparisons = list(c("lower compliance", "higher compliance")),
                                 label = p,
                                 size = 2.7,
                                 method = 'wilcox.test',
                                 paired = F,
                                 label.y.npc = 0.95) +
      theme_bw() +
      theme(legend.position = 'none',
            axis.title.y = ggtext::element_markdown(),
            strip.text.x = ggtext::element_markdown(),
            aspect.ratio = 2) +
      labs(x = "", y = ylab) +
      scale_colour_manual(values = colours,
                          aesthetics = c("fill", 'colour')) +
      facet_wrap(~visitId)
    
    return(plot)
    
  }
