# Clear global environment
rm(list = ls())
library("tidyverse")

# Config
options(width=160)

# Graphics options
graph_unit = "cm"   # unit for graphs sizes ("cm", "in" or "mm" )
graph_width = 30    # graph width
graph_height = 25   # graph height
graph_dpi = 600     # graph resolution (typically between 100 and 1200)
plot_device = "png" # graph file extension

datafiles  = c("RGA4", "AVR-PikD", "RGA4 Strong")
assay_method = "chlorophyll"

##############################################################################################
# Data Import and Transform
##############################################################################################

# set current dir as WD
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

for (f in datafiles) {
  if (file.exists(file.path("rawdata", f, paste0(assay_method, ".csv"))) == FALSE) 
  {
    next
  } else {
    
  df <- read.csv(file.path("rawdata", f, paste0(assay_method, ".csv")))
  resultsfolder = file.path("results", assay_method, f)
  
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
  
  # Split contruct (the different agro infiltrated mixes) and rep (biological replicates) that are separated by "-"
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
    df3$construct <- str_replace_all(df3$old.construct, c("A"="0", "B"="0.02", "C"="0.05", "D"="0.1", "E"="0.2", "F" = "0.5", "G"="None") )
  
  }

  # Convert columns format from character to factor
  df3$construct <- factor(df3$construct, levels=unique(df3$construct)); levels(df3$construct)
  df3$rep <- factor(df3$rep) ; levels(df3$rep)
  df3$construct.rep <- factor(df3$construct.rep) ; levels(df3$construct.rep)
  
  # Export final dataframe to .RData file
  save(df3, file=file.path(resultsfolder, "df3.RData"))
  
  ##############################################################################################
  # Create boxplots from transformed data frame
  ##############################################################################################
  
  # Re-import data if necessary (if scripts are saved separately)
  if (!exists("df3")) {load(file.path(outfolder,"df3.RData"))}
  
  # Calculate means and variances (for "construct.rep") and create a data frame
  ag1 <- aggregate(df3$surface,df3["construct.rep"],mean) ; names(ag1)[2] <- "mean"
  ag2 <- aggregate(df3$surface,df3["construct.rep"],var) ; names(ag2)[2] <- "variance"
  ag <- merge(ag1,ag2)

  # plot variance/mean
  # visualize standard deviation/mean
  windows()
  ggplot(ag, aes(x=mean, y=variance))+
    geom_point(shape=1) +
    geom_smooth(method="lm", se=FALSE) +
    ggtitle("variance") + xlab("mean") + ylab("variance") +
    theme_bw()
  
  # Save last plot using graphics settings
  ggsave(filename= paste0("variance_vs_mean.", plot_device),
         path = resultsfolder,
         device=plot_device,
         width=graph_width,
         height=graph_height,
         unit=graph_unit,
         dpi=graph_dpi)
  

  # Create a table that associates columns "construct" and the corresponding calculated mean  
  ag <- aggregate(df3$surface,df3[c("construct")],mean)
  
  # Modify construct names using "nomenclature.csv"
  if (f =="RGA4") {
    df3$construct <- str_replace_all(df3$old.construct, c("A"="0", "B"="0.02", "C"="0.05", "D"="0.1", "E"="0.2", "F" = "0.5", "G"="None") )
  } else if (f=="AVR-PikD"){
    df3$construct <- str_replace_all(df3$old.construct, c("A"="0", "B"="0.02", "C"="0.05", "D"="0.1", "E"="0.2", "F" = "0.4", "G"="None") )
  } else if (f=="RGA4 Strong"){
    df3$construct <- str_replace_all(df3$old.construct, c("A"="0", "B"="0.1", "C"="0.2", "D"="0.3", "E"="0.4", "F" = "0.5", "G"="0.6", "H"="None") )
  }
  
  # plot labels
  labels <- c(6,8,10,12,14,16,18,20)
  xtitle = paste0(ifelse(f=="AVR-PikD", f, "RGA4"), " (OD600)")
  
  
  # Create a vertical boxplot (several options available: violin etc. See ggplot plots in R supporting documentation)
  ggplot(df3, aes(x=construct,y=surface)) +
    theme_classic() + 
    theme(legend.position="none",axis.text = element_text(color = "black",face="bold",size=24), 
          axis.text.x=element_text(angle=90, hjust=1,vjust = 0.5),
          axis.title.y = element_text(size = 30,face="bold",margin = margin(t = , r =20 , b = , l = )),
          axis.title.x = element_text(size = 30,face="bold",margin = margin(t =20 , r = , b =, l = )))  +
    scale_x_discrete(limits=levels(df3$construct)) +
    scale_y_continuous(breaks=labels, labels=labels) +
    labs(x=xtitle, y="Total chlorophyll (mg/L)") +
    geom_boxplot(outlier.color = "white", aes(fill=construct), alpha=0.3) +
    geom_point(position=position_jitterdodge(jitter.width=0.0, dodge.width = 0.3), cex=3, alpha=0.3,
               aes(color=factor(rep)), show.legend = F) +
    stat_summary(fun=mean, geom="point", shape=20, size=5, color="blue", fill="blue")
    
  # Save boxplot
  ggsave(filename= paste0("boxplot.", plot_device),
         path = resultsfolder,
         device=plot_device,
         width=graph_height,
         height=graph_width,
         unit=graph_unit,
         dpi=graph_dpi)
  
  }
} 
# end of file