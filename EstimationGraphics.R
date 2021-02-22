# Clear global environment
rm(list = ls())


##############################################################################################
# 1.Configuration
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set current dir as WD

raw_ext = "csv"  # raw_ext: Raw datafiles extension ("csv", "tsv", "txt") WITHOUT the dot (.)
sep = ","   # column separator for raw data files (";" or "," or "\t")
dec = "."   # decimal character for raw data files("." or ",")
required_packages = c("tidyverse","tools","RColorBrewer","devtools", "rlang", "besthr","cowplot")
nits = 1000 # A number of bootstrap iteration (default = 1000)

#---------------------------------------------------------------------------------------------

# Graphics parameters
graph_unit = "cm"   # unit for graphs sizes ("cm", "in" or "mm" )
graph_dpi = 600     # graph resolution (typically between 100 and 1200)
plot_device = "png" # graph file extension

#---------------------------------------------------------------------------------------------


##############################################################################################
# the "dot_plot" and "plot_hrest" functions from https://rdrr.io/github/TeamMacLean/besthr/
# ref : Besthr, Maclean 2019, Zenodo
# have been slightly modified to adjust graphics requirements and renamed  "mydot_plot" and "myplot_hrest" respectively thereafter.

mydot_plot <- function(hrest, group_col){
  hrest$ranked_data %>%
    dplyr::group_by(!!group_col, rank) %>%
    dplyr::summarise(Count = dplyr::n() ) %>%
    ggplot2::ggplot() +
    
    ggplot2::aes(!!group_col, rank) +
    ggplot2::geom_point(ggplot2::aes(size = Count, colour = !!group_col,
                                     fill = !!group_col)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = mean, colour = !!group_col),
                        data = hrest$group_means, linetype = 3, size = 1) +
    ggplot2::theme_minimal() +
    ggplot2::scale_size(range = c(5,10)) +
    ylab("Rank")+
    theme(axis.title = element_text(size=20, face="bold"),
          axis.text = element_text(size=14,color = "black",face="bold"),
          axis.title.x = element_text(margin = margin(t = 40, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
    
    theme(legend.title = element_text( vjust = 0.5,hjust=0.7, size=18, face="bold"),
          legend.text = element_text(vjust = 0.5, hjust=0.5, size=14),
          legend.position = "bottom", legend.box = "horizontal") +
    
    ggplot2::xlab(xtitle) +
    guides(colour="none", size="none",fill="none")
}


myplot.hrest <- function(hrest, which = "rank_simulation"){
  group_col <- names(hrest$group_n)[ names(hrest$group_n) != "n" ][[1]]
  group_col <- rlang::sym(group_col)
  
  a <- NULL
  
  quo_group_col <- hrest$column_info[[2]]
  a <- mydot_plot(hrest, quo_group_col)
  
  # legend <- cowplot::get_legend(a + ggplot2::theme(legend.position = "bottom"))
  legend_color <- cowplot::get_legend(a + guides(color = guide_legend(nrow = 1)))
  legend_size  <- cowplot::get_legend(a + guides(size = guide_legend(nrow = 1)))
  
  
  c <- hrest$bootstraps %>%
    ggplot2::ggplot() + ggplot2::aes(mean, UQ(group_col), fill = factor(..quantile..)) +
    ggplot2::xlim(min(hrest$ranked_data$rank), max(hrest$ranked_data$rank)) +
    ggridges::stat_density_ridges(geom = "density_ridges_gradient",
                                  calc_ecdf = TRUE,
                                  quantiles = c(hrest$low, hrest$high)) +
    ggplot2::scale_fill_manual(values = c("#0000FFA0",  "#A0A0A0A0", "#0000FFA0") ) +
    ggplot2::ylab(xtitle) +
    ggplot2::coord_flip() + 
    ggplot2::theme_minimal() +xlab("Mean")+
    ggplot2::theme(axis.title = element_text(size=20, face="bold"),
                   axis.text = element_text(size=14,color = "black",face="bold"),
                   axis.title.x = element_text(margin = margin(t = 40, r = 0, b = 40, l = 0)),
                   axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
    ggplot2::theme(legend.title = element_text( vjust = 0.5,hjust=0.7, size=18, face="bold"),
                   legend.text = element_text(vjust = 0.5, hjust=0.5, size=14),
                   legend.position = "bottom", legend.box = "horizontal") 
  
  d <- cowplot::plot_grid(a  + ggplot2::theme(legend.position = "none"),
                          c  + ggplot2::theme(legend.position = "none"),
                          align = "vh",
                          nrow = 1, axis = c("b"))
  
  return(cowplot::plot_grid(d, legend_color, legend_size, ncol = 1, rel_heights = c(1, .05, .05)))
}






##############################################################################################
# Run script

## chargement (ou installation si nécessaire) des packages nécessaires
for (pkg in required_packages) {
  if(pkg %in% rownames(installed.packages())) {
    print(paste0(c(pkg, " - installé")))
    library(pkg, character.only = TRUE)
  } else {
    print(paste0(c(pkg, " - pas installé")))
    if(pkg=="besthr"){
      devtools::install_github("TeamMacLean/besthr",
                               upgrade="never"      )
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

##############################################################################################


# path to rawdata file  
datafiles  = c("RGA4", "AVR-PikD", "RGA4 Strong") 

for (f in datafiles) {
  
  if (file.exists(file.path("rawdata", f, "visual.csv")) == FALSE) 
  {
    next
  } else 
  {
  resultsfolder = file.path("results","visual", f)
  
  # create folders if necessary
  dir.create("results", showWarnings = FALSE)
  dir.create(file.path("results","visual"), showWarnings = FALSE)
  dir.create(resultsfolder, showWarnings = FALSE)
  
  # load data
  scores <- read.csv(file.path("rawdata", f, "visual.csv")) %>% 
      gather() %>% 
      drop_na() %>% 
      separate(col=key, into=c("Treatment","rep"), sep=1)
    
    # translate treatment names
    if (f=="RGA4 strong") {
      scores$Treatment <- str_replace_all(scores$Treatment, c("A"="0", "B"="0.1", "C"="0.2", "D"="0.3", "E"="0.4", "F" = "0.5", "G" = "0.6") )
    }  else {
      scores$Treatment <- str_replace_all(scores$Treatment, c("A"="0", "B"="0.02", "C"="0.05", "D"="0.1", "E"="0.2", "F" = "0.5") )
    }
    
    
    scores$Treatment = factor(scores$Treatment)
    scores$rep = factor(scores$rep)
    
    # save modified csv
    write.csv(scores, file.path(resultsfolder, paste0("scores.csv")))
    save(scores, file=file.path(resultsfolder, "scores.RData"))
    
    
    # get filename
    xtitle = paste0(ifelse(f=="AVR-PikD", f, "RGA4"), " (OD600)")

    
    # dotplots
    graph_width = 25   # graph width
    graph_height = 34   # graph height
    
    
    ggplot(data=scores,aes(x=Treatment, y=value,  color=Treatment)) +
      geom_count(alpha=0.3, aes(color=Treatment)) +
      geom_jitter(size=1.3, alpha=0.9, width=0.3, height=0.15, aes(pch=rep, color=rep))+
      theme_bw() +
      scale_size(range = c(5,40), breaks=c(5,10,15,20,30,40,50))+
      ylab("Score")+
      xlab(xtitle)+
      guides(color="none", pch="none") +
      theme(axis.title = element_text(color = "black",size=32, face="bold"),
            axis.text = element_text(color = "black", face="bold",size=24),
            axis.title.x = element_text(margin = margin(t = 40, r = 0, b = 40, l = 0)),
            axis.title.y = element_text(margin = margin(t = 0, r = 40, b = 0, l = 0))) +
      theme(legend.title = element_text( vjust = 0.5,hjust=1, size=22, face="bold"),
            legend.text = element_text(vjust = 0.5, hjust=0.5, size=22),
            legend.position = "bottom", legend.box = "vertical") +
      guides(size = guide_legend(nrow = 1,
                            title.position = "left",
                            label.position="top",
                            order=1))

    
    ggsave(filename= paste0("dotplot",".",plot_device),
           path = resultsfolder,
           device=plot_device,
           width=graph_width,
           height=graph_height,
           unit=graph_unit,
           dpi=graph_dpi)
    
    
    
    ##############################################################################################
    # calculate bootstrap HR object creation
    hr <- estimate(scores, value, Treatment,  rep, control=0, nits=nits)
    hr
    
    
    ##############################################################################################
    # ranks estimates graph

    graph_width = 22    # graph width
    graph_height = 30   # graph height
    
    hr %>% myplot.hrest(which = "rank_simulation")
    
    
    ggsave(filename= paste0("est_ranks",".",plot_device),
           path = resultsfolder,
           device=plot_device,
           width=graph_width,
           height=graph_height,
           unit=graph_unit,
           dpi=graph_dpi)
    

  }
}
