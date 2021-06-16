################################################################################
# Initialisation
################################################################################

# Clear global environment
rm(list = ls())

# set current dir as working directory (requires RStudio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# import required libraries
library("tidyverse")

################################################################################
# Experiment specific options
################################################################################

datafiles  = c("RGA4", "AVR-PikD", "RGA4 Strong")
assay_method = "chlorophyll"
labels <- c(6, 8, 10, 12, 14, 16, 18, 20) # y axis labels

################################################################################
# General options
################################################################################

# Graphics options
graph_unit = "cm"   # unit for graphs' width and height ("cm", "in" or "mm" )
graph_width = 30    # graph width
graph_height = 25   # graph height
graph_dpi = 600     # graph resolution (typically between 100 and 1200)
plot_device = "png"

alpha=0.5  # dotplots/boxplots alpha


################################################################################
# Data Import and Transform
################################################################################

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


  df <- read.csv(file.path("rawdata", f, paste0(assay_method, ".csv")))
  resultsfolder = file.path(getwd(), "results", assay_method, f)

  # create folders if necessary
  dir.create("results", showWarnings = FALSE)
  dir.create(file.path("results",assay_method), showWarnings = FALSE)
  dir.create(resultsfolder, showWarnings = FALSE)

  # data transform
  if (f =="RGA4") {
    trans <- function(x) x^0.25
  } else if (f=="AVR-PikD"){
    trans <- function(x) x
  } else if (f=="RGA4 Strong"){
    trans <- function(x) x
  }

  # Transform to long dataframe
  df2 <- reshape(df,
                 timevar="construct.rep",
                 varying=list(1:ncol(df)),
                 times=names(df)[1:ncol(df)],
                 v.names="surface",direction="long")

  # Reset row labels
  rownames(df2) <- NULL

  # remove holes (discard rows with NA)
  df2 <- df2[!is.na(df2$surface),]

  # Split construct (the different agro infiltrated mixes) and rep (biological replicates) that are separated by "-"
  cs <- strsplit(df2$construct.rep,'')
  cs2 <- data.frame(do.call(rbind,cs),stringsAsFactors=FALSE)
  names(cs2) <- c("old.construct","rep")

  # combine cs2 and df2 to a new dataframe
  df3 <- cbind(df2,cs2)

  # Modify construct names using "nomenclature.csv"
  if (f =="RGA4") {
    df3$construct <- str_replace_all(df3$old.construct, c("A"="0", "B"="0.02", "C"="0.05", "D"="0.1", "E"="0.2", "F" = "0.5", "G"="None") )
  } else if (f=="AVR-PikD"){
    df3$construct <- str_replace_all(df3$old.construct, c("A"="0", "B"="0.02", "C"="0.05", "D"="0.1", "E"="0.2", "F" = "0.4", "G"="None") )
  } else if (f=="RGA4 Strong"){
    df3$construct <- str_replace_all(df3$old.construct, c("A"="0", "B"="0.1", "C"="0.2", "D"="0.3", "E"="0.4", "F" = "0.5", "G"="0.6", "H"="None") )
  }

  # Convert columns format from character to factor
  df3$construct <- factor(df3$construct, levels=unique(df3$construct)); levels(df3$construct)
  df3$rep <- factor(df3$rep) ; levels(df3$rep)
  df3$construct.rep <- factor(df3$construct.rep) ; levels(df3$construct.rep)

  # add a "colour" column to dataframe for colour consistency across samples and experiments
  # df3 = merge(df3, custom_palette, by.x="construct", by.y="val")

  # Export final dataframe to .RData file
  save(df3, file=file.path(resultsfolder, "df3.RData"))

  ################################################################################
  # Create boxplots from transformed data frame
  ################################################################################

  # Calculate means and variances (for "construct.rep") and create a data frame
  ag1 <- aggregate(df3$surface,df3["construct.rep"],mean) ; names(ag1)[2] <- "mean"
  ag2 <- aggregate(df3$surface,df3["construct.rep"],var) ; names(ag2)[2] <- "variance"
  ag <- merge(ag1,ag2)

  # plot variance_vs_mean
  ggplot(ag, aes(x=mean, y=variance))+
    geom_point(shape=1) +
    geom_smooth(method="lm", formula="y ~ x", se=FALSE) +
    ggtitle("variance") + xlab("mean") + ylab("variance") +
    theme_bw()

  # Save  variance_vs_mean plot
  ggsave(filename= paste0("variance_vs_mean.", plot_device),
         path = resultsfolder,
         device=plot_device,
         width=graph_width,
         height=graph_height,
         unit=graph_unit,
         dpi=graph_dpi)

  # Create a table that associates columns "construct" and the corresponding calculated mean
  ag <- aggregate(df3$surface,df3[c("construct")],mean)

  # plot title
  xtitle = paste0(ifelse(f=="AVR-PikD", f, "RGA4"), " (OD600)")

  # Create a vertical boxplot (several options available: violin etc. See ggplot plots in R supporting documentation)
  ggplot(df3, aes(x=construct,y=surface)) +
    theme_classic() +
    theme(legend.position="none",axis.text=element_text(color="black", face="bold", size=24),
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          axis.title.y=element_text(size=30, face="bold", margin=margin(r =20)),
          axis.title.x=element_text(size=30, face="bold", margin=margin(t =20)))  +
    scale_x_discrete(limits=levels(df3$construct)) +
    scale_y_continuous(breaks=labels, labels=labels) +
    labs(x=xtitle, y="Total chlorophyll (mg/L)") +
    scale_fill_manual(values=custom_palette) +
    geom_boxplot(outlier.color="white", aes(fill=construct), alpha=alpha) +
    geom_point(position=position_jitterdodge(jitter.width=0.0, dodge.width=0.4), cex=3, alpha=alpha,
               aes(color=factor(rep)), show.legend=F) +
    stat_summary(fun=mean, geom="point", shape=20, size=5, color="blue", fill="blue")

  # Save boxplot
  ggsave(filename= paste0(f, "_", assay_method,"_", "boxplot.", plot_device),
         path = resultsfolder,
         device=plot_device,
         width=graph_height,
         height=graph_width,
         unit=graph_unit,
         dpi=graph_dpi)

}
# end of file
