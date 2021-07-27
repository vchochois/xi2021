################################################################################
# Raw data file format guidelines
################################################################################
# Raw data files format must be ".CSV"

# the separator is set by default on "," but can be set to other separators
# by changing the variable 'sep' in the 'Configuration' section below

# Raw data files must be stored in a subfolder in the "rawdata" folder

# the name of this subfolder corresponds to the construct name
# (for example "RGA4", "AVR-PikD" ...in our case)

# Filenames must correspond to the string specified in the
# "assay_method" variable below :
# 'visual' for the visual scoring method data (EstimationGraphics.R)
# 'chlorophyll' for the chlorophyll assay data (QChlorophyll.R)
# 'ion_leakage' for the ion leakage assay data (QIonLeakage.R)
# 'fluorescence' for the red fluorescence assay data (QRedFluo.R)

# Letters (A, B, C...) correspond to cell concentrations or treatments (OD600).
# Combining a number represents each independent experiments (A1, A2...).
# In each column, each row represents a biological replicate in each experiment

# see the provided sample file for each method in the rawdata_samples section

################################################################################
# Initialisation
################################################################################

# Clear global environment
rm(list = ls())

# set current dir as working directory (requires RStudio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# import required libraries
required_packages = c("tidyverse","tools","RColorBrewer","devtools", "rlang", "besthr","cowplot")

# Load (and install if necessary) all required packages
for (pkg in required_packages) {
  if(pkg %in% rownames(installed.packages())) {
    print(paste0(c(pkg, " - installed")))
    library(pkg, character.only = TRUE)
  } else {
    print(paste0(c(pkg, " - not installed")))
    if(pkg=="besthr"){
      devtools::install_github("TeamMacLean/besthr",
                               upgrade="never"      )
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

################################################################################
# Experiment specific options
################################################################################

datafiles  = c("RGA4", "AVR-PikD", "RGA4 Strong")
assay_method = "visual"
labels <- c(0, 1, 2, 3, 4) # y axis labels

################################################################################
# 1.Configuration

raw_ext = "csv"  # raw_ext: Raw datafiles extension ("csv", "tsv", "txt") WITHOUT the dot (.)
sep = ","   # column separator for raw data files (";" or "," or "\t")
dec = "."   # decimal character for raw data files("." or ",")
nits = 1000 # number of bootstrap iteration (default = 1000)


# Graphics parameters
graph_unit = "cm"   # unit for graphs sizes ("cm", "in" or "mm" )
graph_dpi = 600     # graph resolution (typically between 100 and 1200)
plot_device = "png" # graph file extension

graph_dotplot_width = 25   # graph width
graph_dotplot_height = 34   # graph height

graph_rank_width = 22    # ranks graph width
graph_rank_height = 30   # ranks graph height

alpha=0.5  # dotplots/boxplots alpha

################################################################################
# the "dot_plot" and "plot_hrest" functions from https://rdrr.io/github/TeamMacLean/besthr/
# ref : Besthr, Maclean 2019, Zenodo
# have been slightly modified to adjust graphics requirements and renamed "mydot_plot" and "myplot_hrest" respectively thereafter.


mydot_plot <- function(hrest, group_col){

  hrest$ranked_data %>%
    group_by(!!group_col, rank) %>%
    summarise(Count = n() ) %>%
    ggplot() +
    aes(!!group_col, rank) +
    geom_point(aes(size = Count,
                   colour = !!group_col,
                   fill = !!group_col),
               alpha=alpha) +
    geom_hline(aes(yintercept = mean,
                   colour = !!group_col),
               data = hrest$group_means,
               linetype = 3,
               size = 1) +
    theme_minimal() +
    scale_size(range = c(5,10)) +
    scale_colour_manual(values = dotplot_palette) +

    theme(axis.title = element_text(size=20, face="bold"),
          axis.text = element_text(size=14, color = "black", face="bold"),
          axis.title.x = element_text(margin = margin(t = 40, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +

    theme(legend.title = element_text( vjust = 0.5,hjust=0.7, size=18, face="bold"),
          legend.text = element_text(vjust = 0.5, hjust=0.5, size=14),
          legend.position = "bottom", legend.box = "horizontal") +

    ylab("Rank") + xlab(xtitle) +
    guides(colour="none", size="none",fill="none")
}


myplot.hrest <- function(hrest, which = "rank_simulation"){
  group_col <- names(hrest$group_n)[ names(hrest$group_n) != "n" ][[1]]
  group_col <- rlang::sym(group_col)

  a <- NULL

  quo_group_col <- hrest$column_info[[2]]
  a <- mydot_plot(hrest, quo_group_col)

  # legend <- get_legend(a + theme(legend.position = "bottom"))
  legend_color <- get_legend(a + guides(color = guide_legend(nrow = 1)))
  legend_size  <- get_legend(a + guides(size = guide_legend(nrow = 1)))


  c <- hrest$bootstraps %>%
    ggplot() + aes(mean, UQ(group_col), fill = factor(..quantile..)) +
    xlim(min(hrest$ranked_data$rank), max(hrest$ranked_data$rank)) +
    ggridges::stat_density_ridges(geom = "density_ridges_gradient",
                                  calc_ecdf = TRUE,
                                  quantiles = c(hrest$low, hrest$high)) +
    scale_fill_manual(values = c("#0000FFA0",  "#A0A0A0A0", "#0000FFA0") ) +
    ylab(xtitle) + xlab("Mean") +
    coord_flip() +
    theme_minimal() +
    theme(axis.title = element_text(size=20, face="bold"),
          axis.text = element_text(size=14,color = "black",face="bold"),
          axis.title.x = element_text(margin = margin(t = 40, r = 0, b = 40, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
    theme(legend.title = element_text( vjust = 0.5,hjust=0.7, size=18, face="bold"),
          legend.text = element_text(vjust = 0.5, hjust=0.5, size=14),
          legend.position = "bottom", legend.box = "horizontal")

  d <- plot_grid(a  + theme(legend.position = "none"),
                          c  + theme(legend.position = "none"),
                          align = "vh",
                          nrow = 1, axis = c("b"))

  return(plot_grid(d, legend_color, legend_size, ncol = 1, rel_heights = c(1, .05, .05)))
}

################################################################################
# path to rawdata file

for (f in datafiles) {
  if (file.exists(file.path("rawdata", f, paste0(assay_method, ".csv"))) == FALSE)
  {
    warning(paste0(f, "skipped. No data file found."))
    next
  }

  # custom palette
  if (f =="RGA4") {
    palette_values = c("0", "0.02", "0.05", "0.1", "0.2", "0.5", "None","")
    custom_palette = c("#F8766D", "#D89000", "#A3A500", "#00BF7D", "#00BFC4", "#00B0F6", "#aaaaaa", "#000000")

  } else if (f=="AVR-PikD"){
    palette_values = c("0", "0.02", "0.05", "0.1", "0.2", "0.4", "None","")
    custom_palette = c("#F8766D", "#D89000", "#A3A500", "#00BF7D", "#00BFC4", "#00B0F6", "#aaaaaa", "#000000")

  } else if (f=="RGA4 Strong"){
    palette_values = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6","None","")
    custom_palette = c("#F8766D", "#D89000", "#A3A500", "#00BF7D", "#00BFC4", "#00B0F6", "#9590FF", "#aaaaaa", "#000000")
  }
  names(custom_palette) = palette_values
  dotplot_palette = custom_palette

  df <- read.csv(file.path("rawdata", f, paste0(assay_method, ".csv")))
  resultsfolder = file.path("results", assay_method, f)


  # create folders if necessary
  dir.create("results", showWarnings = FALSE)
  dir.create(file.path("results",assay_method), showWarnings = FALSE)
  dir.create(resultsfolder, showWarnings = FALSE)


  # load data
  scores <- read.csv(file.path("rawdata", f, "visual.csv")) %>%
      gather() %>%
      drop_na() %>%
      separate(col=key, into=c("Treatment","rep"), sep=1)

    # translate treatment names
    if (f =="RGA4") {
      scores$Treatment <- str_replace_all(scores$Treatment, c("A"="0", "B"="0.02", "C"="0.05", "D"="0.1", "E"="0.2", "F" = "0.5", "G"="None") )
    } else if (f=="AVR-PikD"){
      scores$Treatment <- str_replace_all(scores$Treatment, c("A"="0", "B"="0.02", "C"="0.05", "D"="0.1", "E"="0.2", "F" = "0.4", "G"="None") )
    } else if (f=="RGA4 Strong"){
      scores$Treatment <- str_replace_all(scores$Treatment, c("A"="0", "B"="0.1", "C"="0.2", "D"="0.3", "E"="0.4", "F" = "0.5", "G"="0.6", "H"="None") )
    }

    scores$Treatment = factor(scores$Treatment)
    scores$rep = factor(scores$rep)

    # add a "colour" column to dataframe for colour consistency across samples and experiments
    # scores = merge(scores, custom_palette, by.x="Treatment", by.y="val")

    # save modified csv
    write.csv(scores, file.path(resultsfolder, paste0("scores.csv")))
    save(scores, file=file.path(resultsfolder, "scores.RData"))

    # get filename
    xtitle = paste0(ifelse(f=="AVR-PikD", f, "RGA4"), " (OD600)")

    # dotplots
    ggplot(data=scores,aes(x=Treatment, y=value, color=Treatment)) +
      # scale_colour_manual(values = custom_palette) +
      geom_count(alpha=alpha, aes(color=Treatment)) +
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


    ggsave(filename= paste0(f,"_visual_dotplot",".",plot_device),
           path = resultsfolder,
           device=plot_device,
           width=graph_dotplot_width,
           height=graph_dotplot_height,
           unit=graph_unit,
           dpi=graph_dpi)

    ################################################################################
    # calculate bootstrap. creation of hr object
    hr <- estimate(scores, value, Treatment,  rep, control=0, nits=nits)
    ################################################################################

    # ranks estimates graph
    hr %>% myplot.hrest(which = "rank_simulation")

    ggsave(filename= paste0(f, "_visual_est_ranks",".",plot_device),
           path = resultsfolder,
           device=plot_device,
           width=graph_rank_width,
           height=graph_rank_height,
           unit=graph_unit,
           dpi=graph_dpi)


}
