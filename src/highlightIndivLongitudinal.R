#' @param dat Dataset in long format(e.g., change dataset)
#' @param variable variable that should be plotted
#' @param highlight_var variable based on which to highlight
#' @param highlight_val value for highlight_var that should be highlighted
#' @param col colour for the highlighting
#' @param ylab y axis label. can include html (richtext)
#' @param aspect aspect ratio (default is 1.2)
#' @param by_compliance should plot be split by compliance groups?
#' @return ggplot object with highlight by certain variable

highlightIndivLongitudinal <- function(dat,
                           variable,
                           highlight_var,
                           highlight_val,
                           label = 'highlight',
                           ylab,
                           col = 'red',
                           aspect = 1.2,
                           by_compliance = T){
  
  
  label = ifelse(label == 'highlight', highlight_val,
                 label)
  
  if(by_compliance == F){
    
    n <- dat |> 
      dplyr::filter(rowname == variable & .data[[highlight_var]] == highlight_val) |> 
      dplyr::select(subjectId) |> dplyr::distinct() |> nrow()
      
    col = rep(col, n)
    
    p <- dat |> 
      dplyr::filter(rowname == variable) |> 
      dplyr::mutate(indicator = ifelse(.data[[highlight_var]] == highlight_val & visitId == 'M6', label, NA)) |> 
      ggplot(aes(x = as.numeric(gsub("M", "", visitId)),
                 y = value)) +
      geom_line(aes(group = subjectId,
                    colour = subjectId),
                alpha = 0.5) +
      facet_wrap(~compliance) +
      geom_hline(yintercept = 0,
                 linetype = 'dashed',
                 colour = 'grey60') +
      gghighlight::gghighlight(.data[[highlight_var]] == highlight_val,
                               use_direct_label = F) +
      ggtext::geom_richtext(aes(label = indicator),
                            hjust = 1,
                            vjust = 1,
                            fill = NA, label.color = NA,
                            colour = col) +
      theme_bw() +
      theme(axis.title.y = ggtext::element_markdown(),
            legend.position = 'none',
            aspect.ratio = aspect) +
      labs(x = 'visit',
           y = ylab) +
      scale_colour_manual(values = col)
    
    } else {
  
    dat <- dat |> 
      dplyr::filter(rowname == variable)
    
    ylims = c(min(dat$value)*1.1, max(dat$value)*1.1)
    
    dat1 <- dat |> 
      dplyr::filter(compliance == 'lower compliance')
    
    if(any(highlight_val %in% dat1[[highlight_var]])){
      
      n <- dat1 |> 
        dplyr::filter(rowname == variable & .data[[highlight_var]] == highlight_val) |> 
        dplyr::select(subjectId) |> dplyr::distinct() |> nrow()
      
      col1 = rep(col, n)
      
      p1 <- dat1 |> 
        dplyr::mutate(indicator = ifelse(.data[[highlight_var]] == highlight_val & visitId == 'M6', label, NA)) |> 
        ggplot(aes(x = as.numeric(gsub("M", "", visitId)),
                   y = value)) +
        geom_line(aes(group = subjectId,
                      colour = subjectId),
                  alpha = 0.5) +
        facet_wrap(~compliance) +
        geom_hline(yintercept = 0,
                   linetype = 'dashed',
                   colour = 'grey60') +
        gghighlight::gghighlight(.data[[highlight_var]] == highlight_val,
                                 use_direct_label = F) +
        ggtext::geom_richtext(aes(label = indicator),
                              hjust = 1,
                              vjust = 1,
                              fill = NA, label.color = NA,
                              colour = col) +
        theme_bw() +
        theme(axis.title.y = ggtext::element_markdown(),
              legend.position = 'none',
              aspect.ratio = aspect) +
        labs(x = 'visit',
             y = ylab) +
        ylim(ylims) +
        scale_colour_manual(values = col1)
      
    } else {
      p1 <- dat1 |> 
        ggplot(aes(x = as.numeric(gsub("M", "", visitId)),
                   y = value)) +
        geom_line(aes(group = subjectId),
                  colour = 'grey60',
                  alpha = 0.5) +
        facet_wrap(~compliance) +
        geom_hline(yintercept = 0,
                   linetype = 'dashed',
                   colour = 'grey60') +
        theme_bw() +
        theme(axis.title.y = ggtext::element_markdown(),
              legend.position = 'none',
              aspect.ratio = aspect) +
        labs(x = 'visit',
             y = ylab) +
        ylim(ylims) +
        scale_colour_manual(values = col)
      
    }
    
    
    dat2 <- dat |> 
      dplyr::filter(rowname == variable & compliance == 'higher compliance')
    
    if(any(highlight_val %in% dat2[[highlight_var]])){
      
      n <- dat2 |> 
        dplyr::filter(rowname == variable & .data[[highlight_var]] == highlight_val) |> 
        dplyr::select(subjectId) |> dplyr::distinct() |> nrow()
      
      col2 = rep(col, n)
      
      p2 <- dat2 |> 
        dplyr::mutate(indicator = ifelse(.data[[highlight_var]] == highlight_val & visitId == 'M6', label, NA)) |> 
        ggplot(aes(x = as.numeric(gsub("M", "", visitId)),
                   y = value)) +
        geom_line(aes(group = subjectId,
                      colour = subjectId),
                  alpha = 0.5) +
        facet_wrap(~compliance) +
        geom_hline(yintercept = 0,
                   linetype = 'dashed',
                   colour = 'grey60') +
        gghighlight::gghighlight(.data[[highlight_var]] == highlight_val,
                                 use_direct_label = F) +
        ggtext::geom_richtext(aes(label = indicator),
                              hjust = 1,
                              vjust = 1,
                              fill = NA, label.color = NA,
                              colour = col) +
        theme_bw() +
        theme(axis.title.y = ggtext::element_markdown(),
              legend.position = 'none',
              aspect.ratio = aspect,
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) +
        labs(x = 'visit',
             y = '') +
        ylim(ylims) +
        scale_colour_manual(values = col2)
      
    } else {
      
      p2 <- dat2 |> 
        ggplot(aes(x = as.numeric(gsub("M", "", visitId)),
                   y = value)) +
        geom_line(aes(group = subjectId),
                  colour = 'grey60',
                  alpha = 0.5) +
        facet_wrap(~compliance) +
        geom_hline(yintercept = 0,
                   linetype = 'dashed',
                   colour = 'grey60') +
        theme_bw() +
        theme(axis.title.y = ggtext::element_markdown(),
              legend.position = 'none',
              aspect.ratio = aspect,
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) +
        labs(x = 'visit',
             y = '') +
        ylim(ylims) +
        scale_colour_manual(values = col)
      
    }
    
    p <- p1|p2
  
    }
  
  return(p)  
}
  