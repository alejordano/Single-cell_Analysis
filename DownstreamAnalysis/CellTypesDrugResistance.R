library(ggplot2)
library(dplyr)
library(ggpubr)

#function multiplot
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#get oxphos and cell death data
terms <- c("oxidative phosphorylation",
          "proton transmembrane transport",
          "mitochondrial ATP synthesis coupled electron transport")
terms2 <- c("negative regulation of cell death")
data <- read.csv("summaryrrvgoweigthedfreq", row.names = 1)
oxphos <- data %>% filter(term %in% terms) %>% distinct(sample)
ncd <- data %>% filter(term %in% terms2) %>% distinct(sample)

#get number of cell type proportions in each subpopulation with high expression of MT respiration and negative regulation of cell death
data <- readRDS("./AllMerged_NormSubsetCluster.Rds")

data$pat <- sapply(strsplit(data$batch, split = "_", fixed = TRUE), function(x) (x[1]))
data$tp <- sapply(strsplit(data$batch, split = "_", fixed = TRUE), function(x) (x[2]))
data$sample <- paste0(data$pat, "_", data$clusters)
df1 <- as.data.frame(data@meta.data)
df1$id <- rownames(df1)

df <- df1 %>% group_by(sample) %>% mutate(NcellsCluster = n_distinct(id))
df <- df %>% group_by(sample, celltypefinal) %>% mutate(NcellsClsCT = n_distinct(id))
df$percCT <- (df$NcellsClsCT / df$NcellsCluster) * 100

#### mitochondrial respiration ####
df$oxphosyn <- ifelse(df$sample %in% oxphos$sample, "OXPHOS", "Others")

names <- c("1", "2", "3", "6")
pats <- c("1216", "1886", "3904", "5750")
plot_list <- list()
#Plot figure 7A:
for (i in 1:4){
  df2 <- df %>% filter(pat == pats[i])
  pl5 <- ggplot(df2, aes(x = oxphosyn, y = percCT, fill = oxphosyn)) + geom_boxplot() + facet_wrap(~ celltypefinal) +
    theme_bw() + theme(legend.position = "none") + xlab("") + ylab("Percentage of cells") + ggtitle(paste0("Patient ", names[i])) +
    stat_compare_means(method = "t.test", label = "p.signif", label.x = 1.5, label.y = 85)  + ylim(0, 100)
  plot_list[[i]] <- pl5
}
multiplot(plotlist = plot_list, cols = 2)

## Is up-regulation of oxphos at relapse due to different cell types between timepoints or due to relapse?
#plot FIgure 7B:
plot_list2 <- list()
df3 <- df %>% filter(tp != "MRD", (celltypefinal == "HSC.MPP" | celltypefinal == "GMP"))
for (i in 1:3){
pl6 <- ggplot(df3, aes(x = oxphosyn, y = percCT, fill = tp)) + geom_boxplot() + facet_wrap(~ celltypefinal) +
  theme_bw() + theme(legend.position = "none") + xlab("") + ylab("Percentage of cells") + ggtitle(paste0("Patient ", names[i])) +
  ylim(0, 100)
plot_list2[[i]] <- pl6
}
pl6 <- ggplot(df3, aes(x = oxphosyn, y = percCT, fill = tp)) + geom_boxplot() + facet_wrap(~ celltypefinal) +
  theme_bw() + xlab("") + ylab("Percentage of cells") + ggtitle(paste0("Patient ", names[4])) +
  ylim(0, 100)
plot_list2[[4]] <- pl6
multiplot(plotlist = plot_list2, cols = 2)

#### cell death ####
df$cd <- ifelse(df$sample %in% ncd$sample, "CellDeath", "Others")
df2 <- df %>% filter(sample %in% ncd$sample)

#Plot figure 8:
names <- c("1", "3", "4", "6", "7", "9")
pats <- c("1216", "3904", "40389", "5750", "6323", "3853")
plot_list <- list()
for (i in 1:6){
  df2 <- df %>% filter(pat == pats[i])
  pl5 <- ggplot(df2, aes(x = cd, y = percCT, fill = cd)) + geom_boxplot() + facet_wrap(~ celltypefinal) +
    theme_bw() + theme(legend.position = "none") + xlab("") + ylab("Percentage of cells") + ggtitle(paste0("Patient ", names[i])) +
    ylim(0, 100)
  plot_list[[i]] <- pl5
}
multiplot(plotlist = plot_list, cols = 3)