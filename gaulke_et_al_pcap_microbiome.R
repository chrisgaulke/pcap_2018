# Title and author information --------------------------------------------
#!/usr/bin/R


##################################
#                                #
# gaulke_et_al_pcap_microbiome.R #
#                                #
##################################

# Copyright (C) 2017-2018  Christopher A. Gaulke
# author contact: gaulkec@oregonstate.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.


# Purpose -----------------------------------------------------------------

# This code document outlines the analyses described in the manuscript that
# links to this repository. The modeling analyses are described in a separate
# code document for clarity.

# Versions Options and Packages -------------------------------------------

#Packages

getRversion()
library(RColorBrewer)
library(vegan)
library(reshape2)
library(labdsv)
library(plyr)
library(dplyr)
library(ggplot2)

#options
options("stringsAsFactors"=F)

# FUNCTIONS ---------------------------------------------------------------

# There are a significant number of functions that are used in this script.
# Most of these functions are part of a function library that I maintain.
# Because these functions will likely be modified over time I am including a
# static copy here.


#------------------------------------------------#
#          Function adjust_colnames              #
#------------------------------------------------#

adjust_colnames <- function(df){
  # adjust colnames of dada2 sequence table to be colnames, not nt sequences
  x <- paste0(rep("seq", times = ncol(df)), c(1:ncol(df)))
  colnames(df) <- x
  return(df)
}

#------------------------------------------------#
#              Function normalize                #
#------------------------------------------------#

normalize <- function(df,
                      method="rel",
                      depth=depth){
  # normalize counts either by relative abundance of rarefying
  # requires vegan
  # default method = relative abundance
  if(method == "rare"){
    if( is.null(depth)){
      depth <- min(rowSums(df))
    }
    ndf <- df[which(rowSums(df) > depth),,drop=F]
    ndf <- rrarefy(ndf, depth)
    ndf <- ndf[,which(colSums(ndf) > 0), drop =F]
  }else{
    ndf <- sweep(df,
                 1,
                 rowSums(df),
                 `/`)
  }
  return(ndf)
}


#------------------------------------------------#
#             Function run_beta                  #
#------------------------------------------------#


run_beta <- function(obj){
  # calculates beta-div between all samples
  bdiv <- as.matrix(vegdist(obj$data))
  obj$bdiv <- bdiv
  return(obj)
}

#------------------------------------------------#
#             Function beta_intra                #
#------------------------------------------------#

beta_intra <- function(obj, group){
  #given a dist like matrix and grouping vector this returns within-group
  #beta-diversity
  #group: names of samples in group
  bdiv <- obj$bdiv
  bdiv_in <- bdiv[group, group]
  bdiv_in <- bdiv_in[upper.tri(bdiv_in)]
  return(bdiv_in)
}

#------------------------------------------------#
#             Function beta_inter                #
#------------------------------------------------#

beta_inter <- function(obj,x,y){
  # given a dist like matrix and two grouping vectors this returns between-group
  # beta-diversity for x and y
  # x: names of samples in group 1
  # y: names of samples in group 2
  bdiv <- obj$bdiv
  bdiv_in <- as.numeric(bdiv[x, y])
  return(bdiv_in)
}


#------------------------------------------------#
#             Function make_beta                 #
#------------------------------------------------#

make_beta <- function(obj, c){
  # obj: microbiome data object
  # c: column number in the metadata slot

  intra <- NULL
  # get list of all rownames across all groups
  ug <- unique(obj$meta[,c])
  groups <- list()

  for(j in 1:length(ug)){
    rn <- rownames(obj$meta[which(obj$meta[,c] == ug[j]),])
    # can't be 0 or you get an error so must be charater names
    groups[[as.character(ug[j])]] <- rn
  }

  intra <- lapply(groups, beta_intra, obj= obj)
  inter <- list()
  for(k in 1:length(ug)){
    for(l in 1:length(ug)){
      if(l > k ) {
        x <- beta_inter(obj, unlist(groups[ug[k]]), unlist(groups[ug[l]]))
        name <- paste0(ug[k], "_", ug[l])
        inter[[name]] <- x
      }
    }
  }
  bdiv_intra.lens  <- sapply(intra, length)
  bdiv_intra.names <- names(intra)
  bdiv_intra.names <- rep(bdiv_intra.names, times = bdiv_intra.lens)
  bdiv_intra.df <- cbind(names = bdiv_intra.names,
                         value = unlist(intra))
  bdiv_intra.df <- as.data.frame(bdiv_intra.df)
  bdiv_intra.df$value <- as.numeric(bdiv_intra.df$value)
  name <- paste0("bdiv_intra_",colnames(obj$meta)[c])
  obj[[name]] <- bdiv_intra.df


  bdiv_inter.lens  <- sapply(inter, length)
  bdiv_inter.names <- names(inter)
  bdiv_inter.names <- rep(bdiv_inter.names, times = bdiv_inter.lens)

  bdiv_inter.df <- cbind(names = bdiv_inter.names,
                         value = unlist(inter))
  bdiv_inter.df <- as.data.frame(bdiv_inter.df)
  bdiv_inter.df$value <- as.numeric(bdiv_inter.df$value)
  name <- paste0("bdiv_inter_",colnames(obj$meta)[c])

  obj[[name]] <- bdiv_inter.df

  return(obj)
}


#------------------------------------------------#
#        Function diversity_analysis             #
#------------------------------------------------#

diversity_analysis <- function(obj){
  # this function performs alpha and beta diversity analysis on the
  # user provided data frame
  # requires vegan
  # obj: microbiome data object
  obj$shannon <- diversity(obj$data, index = "shannon", MARGIN = 1)
  obj$simpson <- diversity(obj$data, index = "simpson", MARGIN = 1)
  obj$invsimpson <- diversity(obj$data, index = "invsimpson", MARGIN = 1)
  return(obj)
}

#------------------------------------------------#
#              Function ordinate                 #
#------------------------------------------------#

ordinate <- function(obj,
                     dims = 5,
                     trys = 20,
                     dm = "bray",
                     scale =T,
                     center = T,
                     maxit = 200 ){
  # This function creates ordination objects for plotting later
  # obj: a microbiome object to be passes to metaMDS and prcomp
  # dims: The number of dims to return
  # trys: THe max number of trys for try max
  # dm: distance metric for metaMDS
  # scale: whether or not to use scale in prcomp
  # center: Whether of not to use center in prcomp

  obj_mmds <- metaMDS(obj$data,
                      k =dims,
                      distance = dm,
                      trymax   = trys,
                      maxit = maxit

  )

  obj_prcomp <- prcomp(obj$data,
                       scale = scale,
                       center = center
  )
  obj$mds     <- NULL
  obj$prcomp  <- NULL

  obj_mmds.df <- as.data.frame(obj_mmds$points)
  obj_mmds.df <- obj_mmds.df[rownames(obj$meta),]
  obj$mds$df  <- obj_mmds.df
  obj$mds$obj <- obj_mmds


  obj_prcomp.df  <- as.data.frame(obj_prcomp$x[,1:dims])
  obj_prcomp.df  <- obj_prcomp.df[rownames(obj$meta),]
  obj$prcomp$df  <- obj_prcomp.df
  obj$prcomp$obj <- obj_prcomp
  return(obj)
}

#------------------------------------------------#
#        Function phylotype_analysis             #
#------------------------------------------------#

phylotype_analysis <- function(obj, tax){
  # Create an object of phylotype slots containing abundance data for all
  # phylotypes
  # obj: microbiome object with at least 1 slot (data)
  # tax: a tax object (named list taxa as names values in the list are seq ids)

  obj.out <- NULL
  for(h in 1:length(tax)){
    df <- NULL
    for( i in 1:length(tax[[h]])){
      v1       <- obj$data[,unlist(tax[[h]][[i]])]
      v2       <- names(tax[[h]])[i]
      if(is.null(dim(v1))){
        df[[v2]] <- v1
      }else{
        df[[v2]] <- rowSums(v1)
      }
    }
    obj.out[[names(tax)[h]]] <- as.data.frame(df)
  }
  return(obj.out)
}

#------------------------------------------------#
#            Function group_means                #
#------------------------------------------------#

group_means <- function(df, mapping, t = F){
  # This function will calculate the group means across all columns in a
  # dataframe. Note all columns must be numeric
  # df : data frame with rows as sample IDs and columns as objects (e.g. OTUs)
  # mapping : a mapping df with rownames df as rownames and group id in col2
  # t : boolean, defaults to FALSE. Use if df has sample ids as colnames

   if(t == T){
    df <- t(df)
  }
  groups <- base::unique(x=mapping[,2])
  my_df <- data.frame(matrix(nrow = length(groups), ncol = ncol(df)))
  for(i in 1:ncol(df)){
    tgvec <- NULL
    for(j in 1:length(groups)){
      s <- base::rownames(mapping[base::which(mapping[,2] %in% groups[j]),
                                  ,drop=F])
      m <- base::mean(df[s,i])
      tgvec <- c(tgvec, m)
    }
    my_df[,i] <- tgvec
  }
  rownames(my_df) <- groups
  colnames(my_df) <- colnames(df)
  return(my_df)
}


# IMPORT DATA -------------------------------------------------------------

pcap_indir  <- "/Users/gaulkec/Chris/dev/R_projects/pseudocap_long_2017/data/"

df <- read.table(paste0(pcap_indir,
                  "2018_pcap_longit/analysis/dada2_out/sequence_table.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

df <- adjust_colnames(df)

biome_metadf <- read.table(paste0(pcap_indir,
                          "pcap_long_metadata_and_worm_counts_rformat.txt"),
                           sep = "\t",
                           header = T)


rownames(biome_metadf) <- biome_metadf[,1]

taxadf <- read.table(paste0(pcap_indir,
                         "2018_pcap_longit/analysis/translated_tax_dict.txt"),
                          sep = "\t",
                          header = F)

rownames(taxadf) <- taxadf[,1]

pcap_worm_counts <- read.table(paste0(pcap_indir,
                         "pcap_long_metadata_and_worm_counts_rformat.txt"),
                          sep = "\t",
                          header =T)

rownames(pcap_worm_counts) <- pcap_worm_counts$Sample

#Analysis of kit blanks -----------------------------------------------------

# We have included several kit blanks in the this study and push them into
# their own data frame for subsequent analysis. This will allow us to analyze
# the fish data sepeartely as combining kit blanks and fish data can skew
# some analysis (i.e., ordination).

kit_blanks <- c("lane1-s212-index-AGCGATGCCTTA-PCAP-KB03",
                "lane1-s116-index-TTATCACGTGCA-PCAP-KB02",
                "lane1-s044-index-CGATCCGTATTA-PCAP-KB01")

# splice out just the innoc and kb for downstream analyses

set.seed(731)

df.rare <- normalize(df, "rare", depth = 3000)
dim(df.rare)

df.rare <- df.rare[,which(colSums(df.rare) > 0 ),drop=F]
dim(df.rare)

# The original metadata file has only fish in it, so I need to add the KB ids
# to this table.
full_metadf <- biome_metadf

full_metadf[c("lane1-s212-index-AGCGATGCCTTA-PCAP-KB03",
              "lane1-s116-index-TTATCACGTGCA-PCAP-KB02",
              "lane1-s044-index-CGATCCGTATTA-PCAP-KB01"),] <- NA

full_metadf[c("lane1-s212-index-AGCGATGCCTTA-PCAP-KB03",
              "lane1-s116-index-TTATCACGTGCA-PCAP-KB02",
              "lane1-s044-index-CGATCCGTATTA-PCAP-KB01"),"Exposure"] <- "KB"

# remove samples not in metadata (just the innoc samples)
df.rare <- df.rare[which(rownames(df.rare) %in% rownames(full_metadf)),]

df_rare.dist <- vegdist(df.rare)

df_rare.dist <- as.matrix(df_rare.dist)

fish.names <- rownames(df.rare[-which(rownames(df.rare)
                      %in% c("lane1-s212-index-AGCGATGCCTTA-PCAP-KB03",
                      "lane1-s116-index-TTATCACGTGCA-PCAP-KB02",
                      "lane1-s044-index-CGATCCGTATTA-PCAP-KB01")),])

within_fish <- df_rare.dist[fish.names,fish.names]
within_fish[upper.tri(within_fish)]

with_kb <- df_rare.dist[kit_blanks,kit_blanks]
with_kb[upper.tri(with_kb)]

btwn_fish_kb <- df_rare.dist[kit_blanks,fish.names]

boxplot(with_kb[upper.tri(with_kb)],
        within_fish[upper.tri(within_fish)],
        unlist(btwn_fish_kb)
        )

kruskal.test(list(with_kb[upper.tri(with_kb)],
             within_fish[upper.tri(within_fish)],
             unlist(btwn_fish_kb)
              )
             )

pairwise.wilcox.test(c(with_kb[upper.tri(with_kb)],
                          within_fish[upper.tri(within_fish)],
                          unlist(btwn_fish_kb)
                       )


, g = c(rep("within kb",
            times =    length(with_kb[upper.tri(with_kb)])),
         rep("within fish",
             times = length(within_fish[upper.tri(within_fish)])),
         rep("btwn fish:kb",
             times = length(unlist(btwn_fish_kb))))
,
p.adjust.method = "bonf"

)

# group means tax

df.rare <- df.rare[c(fish.names, kit_blanks),]
df.rare <- as.data.frame(df.rare)


df_rare.metadata <- cbind(id = rownames(df.rare), group =
                   rep(c("fish", "kb"), times = c(237,3)))

rownames(df_rare.metadata) <- df_rare.metadata[,1]
df_rare.metadata <- as.data.frame(df_rare.metadata)

for(i in 1:ncol(df.rare)){
  df.rare[,i] <- as.numeric(df.rare[,i])
}

df_rare.gm <- group_means(df.rare, df_rare.metadata)

df_rare.gm <- df_rare.gm[,order(df_rare.gm["fish",],
                                decreasing =T)]

df_rare.order <- colnames(df_rare.gm)
df_rare.gm$group <- c("fish", "kb")
df_rare.gm <- melt(df_rare.gm)


# look at top 100 fish sequence variants and their abundance
# in kit blanks

df_rare.gm <- df_rare.gm[which(df_rare.gm$variable %in%
                                 df_rare.order[1:100]),]

df_rare.gm.plot <- ggplot(df_rare.gm,
                          aes(x = variable,
                              y = value,
                              group =group,
                              colour = group)
                          )

df_rare.gm.plot +
  geom_line() +
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank()
  ) +
  xlab("Sequence Variant") +
  ylab("Abundance")

# Few abundant fish taxa are found at high abundance in kit blanks

# MAKE MICROBIOME OBJECT -----------------------------------------------------

pcap_microbiome <- NULL

drops <- c("lane1-s212-index-AGCGATGCCTTA-PCAP-KB03",
           "lane1-s116-index-TTATCACGTGCA-PCAP-KB02",
           "lane1-s044-index-CGATCCGTATTA-PCAP-KB01",
           "lane1-s002-index-AATCCGTACAGC-neg-innoc",
           "lane1-s001-index-AATCAGTCTCGT-pcap-innoc")

# remove kit blanks and innoc
pcap_microbiome$data <- df[-which(rownames(df) %in% drops),,drop=F]

# remove 0 sum columns resulting from removal of drops above
pcap_microbiome$data <- pcap_microbiome$data[,
                                 which(colSums(pcap_microbiome$data) > 0),
                                 drop=FALSE]

set.seed(731)
pcap_microbiome$data <- normalize(pcap_microbiome$data,
                                  "rare",
                                  depth = 5000)

# add metadata to the object
pcap_microbiome$meta <- biome_metadf[rownames(pcap_microbiome$data),]


# MICROBIOME Alpha Diversity --------------------------------------------

pcap_microbiome <- diversity_analysis(pcap_microbiome)

pcap_microbiome.shannon <- data.frame(shannon = pcap_microbiome$shannon,
                            time = pcap_microbiome$meta[,"DaysPE"],
                            exposure = pcap_microbiome$meta[,"Exposure"],
                            burden = pcap_microbiome$meta[,"Total"],
                            hist_total = pcap_microbiome$meta$histo_total,
                            cond_fac  = pcap_microbiome$meta$condition_factor)

# To determine if micorbiome alpha diversity is influenced by any of the
# covariates.
pcap_shannon.lm <- lm(log(shannon) ~ time *
                        factor(exposure, levels = c("Unexposed", "Exposed")) +
                        burden +
                        hist_total +
                        cond_fac,
                      data = na.omit(pcap_microbiome.shannon))

step.shannon <- stepAIC(pcap_shannon.lm, direction="both")

# build the best model
pcap_shannon.lm2 <- lm(log(shannon) ~ time +
                         factor(exposure, levels = c("Unexposed", "Exposed")) +
                         burden +
                         time:factor(exposure,
                                     levels = c("Unexposed", "Exposed")),
                       data = na.omit(pcap_microbiome.shannon))

summary(pcap_shannon.lm2)

pcap_microbiome_shannon.plot <- ggplot( pcap_microbiome.shannon,
                                        aes(x = factor(time,
                                           levels =c(0,7,10,21,30,43,59,86)),
                                        y = shannon,
                                           fill = factor(exposure,
                                           levels = c("Unexposed", "Exposed")
                                           )
                                          )
                                        )

pcap_microbiome_shannon.plot +
  geom_boxplot() +
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank()
  ) +
  xlab("Days Post Exposure")+
  scale_fill_manual(values = c("#217CA3", "#D24136"))


# MICROBIOME Beta Diversity -----------------------------------------------

# add a parameter to the metadata for a figure later on

pcap_microbiome$meta$daype_exp <- paste(pcap_microbiome$meta$DaysPE,
                                        pcap_microbiome$meta$Exposure,sep =
                                          "_")
# diversity
# should also add a dist object to the object
pcap_microbiome <- run_beta(pcap_microbiome)

pcap_microbiome <- make_beta(pcap_microbiome,4)
pcap_microbiome <- make_beta(pcap_microbiome,5)
pcap_microbiome <- make_beta(pcap_microbiome,3)
pcap_microbiome <- make_beta(pcap_microbiome,17)
pcap_microbiome <- make_beta(pcap_microbiome,18)

# no split out days post exposure and exposure status
x <- sapply(unlist(pcap_microbiome$bdiv_intra_daype_exp$names),
            function(x) strsplit(x,split = "_")[[1]])

# split into components
day <- x[1,]
names(day) <- NULL
exposure <- x[2,]
names(exposure) <- NULL
pcap_microbiome$bdiv_intra_daype_exp$day <- day
pcap_microbiome$bdiv_intra_daype_exp$exposure <- exposure


pcap_time_exp.plot <- ggplot(pcap_microbiome$bdiv_intra_daype_exp,
                             aes(x = factor(day,
                                            levels =c(0,7,10,21,30,43,59,86)),
                                 y = value,
                                 fill = factor(exposure,
                                            levels = c("Unexposed", "Exposed")
                                 )
                             )
)


pcap_time_exp.plot +
  geom_boxplot()+
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black", angle =45, hjust =1),
        legend.key = element_blank(),
        legend.title = element_blank()
  ) +
  xlab("Days Post Exposure")+
  ylab("Beta Diversity Within (Bray Curtis)")+
  scale_fill_manual(values = c("#217CA3", "#D24136"))



# Test if beta dispersion is different between exposed and unexposed fish

pcap_bdiv.dist <- vegdist(pcap_microbiome$data)
pcap_bdiv.labels <- labels(pcap_bdiv.dist)
pcap_bdiv.betadisp <-
  betadisper( pcap_bdiv.dist,
              pcap_microbiome$meta[pcap_bdiv.labels, "Exposure"])

anova(pcap_bdiv.betadisp)

# MICROBIOME ORDINATION: NMDS ----------------------------------------------

pcap_microbiome <- ordinate(pcap_microbiome, dims = 5, trys = 20, maxit =400,
                            scale = F, center = F)# no mds convergence


pcap_microbiome$mds$df$group <- pcap_microbiome$meta[,4]
pcap_microbiome$mds$df$day   <- pcap_microbiome$meta[,3]
pcap_microbiome$mds$df$worm <- pcap_microbiome$meta[,8]
pcap_microbiome$mds$df$hist <- pcap_microbiome$meta[,17]


pcap_microbiome$mds$df$rescale_worm <- cut(pcap_microbiome$mds$df$worm,
                                           breaks=c(-1,0,3,5,8,10,100),
                                           labels=c("0",
                                                    "1-3",
                                                    "3-5",
                                                    "5-8",
                                                    "8-10",
                                                    ">10"
                                           )
)


pcap_microbiome$mds$df$rescale_hist <- cut(pcap_microbiome$mds$df$hist,
                                           breaks=c(-1,0,1,2,3,4,5,6),
                                           labels=c("0",
                                                    "1",
                                                    "2",
                                                    "3",
                                                    "4",
                                                    "5",
                                                    "6"
                                           )
)


# Define color set for nmds
darkcols <- brewer.pal(8, "Dark2")

# Oridnation by day

pcap_mds <- ggplot(pcap_microbiome$mds$df, aes(x = MDS1,
                                               y = MDS2,
                                               color = as.factor(day),
                                               shape = group
)
)

pcap_mds +
  geom_point(size =4)+
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )+
  scale_colour_manual(values=c(darkcols))+
  coord_fixed()


# Ordination by hist

my_blues <- colorRampPalette(colors = c("white", "blue", "black"))
pcap_mds <- ggplot(pcap_microbiome$mds$df, aes(x = MDS1,
                                               y = MDS2,
                                               color = rescale_hist,
                                               shape = group
)
)

pcap_mds +
  geom_point(size =4)+
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )+
  coord_fixed()+
  scale_colour_manual(values=c(my_blues(100)[c(100,10,20,30,35,40,55)]),
                      na.value="grey")

pcap_mds <- ggplot(pcap_microbiome$mds$df, aes(x = MDS1,
                                               y = MDS2,
                                               color = rescale_worm,
                                               shape = group
)
)

pcap_mds +
  geom_point(size =4)+
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )+
  coord_fixed()+
  scale_colour_manual(values=c(my_blues(100)[c(100,10,20,30,40,55)]),
                      na.value="grey")


# Quantify the associations between components of microbiome diversity and worm
# burden and histopath.

cor.test(pcap_microbiome$mds$df$MDS1 ,
         pcap_microbiome$mds$df$worm,
         method = "kendall")
cor.test(pcap_microbiome$mds$df$MDS1 ,
         pcap_microbiome$mds$df$hist,
         method = "kendall")

# Plot worm burden by MDS1
mds1_worm.plot <- ggplot(pcap_microbiome$mds$df, aes(x = MDS1,
                                                     y = worm)
)

mds1_worm.plot +
  geom_point(size = 3)+
  geom_smooth(method = "loess")+
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank()
  )+
  ylab("Total Parasite Burden")

# Plot total histopathology score by MDS1

mds1_hist.plot <- ggplot(pcap_microbiome$mds$df, aes(x = MDS1,
                                                     y = hist)
)

mds1_hist.plot +
  geom_point(size = 3)+
  geom_smooth(method = "loess")+
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank()
  )+
  ylab("Total Histopathology")

# MICROBIOME DIVERSITY: Adonis ----------------------------------------------

# Adonis does not handle NAs well. Because there are NAs at the 0 dpe time
# point we need to remove this time point and other samples that contain NAs
# in any of the parameters we are testing.

drops <- rownames(pcap_microbiome$meta[
  which(pcap_microbiome$meta$DaysPE == 0),])

# Eliminate NAs
drops <- c(drops,
           "lane1-s204-index-TGCACAATGGCG-Z0574",
           "lane1-s210-index-AGTTAGTGCGTC-Z0580",
           "lane1-s215-index-TCTTAACGACTG-Z0584")

tdata <- pcap_microbiome$data[-which(rownames(pcap_microbiome$data)
                                     %in% drops),
                              ,drop=F]
tmetadata <- pcap_microbiome$meta[-which(rownames(pcap_microbiome$meta)
                                         %in% drops),
                                  ,drop=F]


# For clarity, a list of the column IDs of parameters can be found below.

# 3 Days
# 4 Exposure
# 5 Sex
# 8 Total
# 15 Condition Factor
# 16 Tank
# 17 Total histopath


# Adonis is sensitive to the order of parameters. However, it may be
# inappropriate to omit some parameters when considering certain parameters.
# For example, tank is confounded by exposure status (i.e.,
# tank 1,2,3 = unexposed, while tank 4,5,6 = exposed) thus exposure should
# be included in any model in which tank appears. Since it is difficult, or
# impossible, to know which parameters are may impact the others influence on
# diversity ahead of time I adopt an approach where I generate individual \
# models and also a full model which includes all parameters.

# single parameter models

days.adonis <- adonis(tdata ~ tmetadata[,3]  ,
                      permutations = 5000)

exposure.adonis <- adonis(tdata ~ tmetadata[,4]  ,
                          permutations = 5000)

sex.adonis      <- adonis(tdata ~ tmetadata[,5]  ,
                          permutations = 5000)

burden.adonis   <- adonis(tdata ~ tmetadata[,8]  ,
                          permutations = 5000)

cf.adonis       <- adonis(tdata ~ tmetadata[,15]  ,
                          permutations = 5000)

hist.adonis     <- adonis(tdata ~ tmetadata[,17]  ,
                          permutations = 5000)

# necessarily need exposure in addition to tank in this model as explained above

tank.adonis     <- adonis(tdata ~ tmetadata[,4] + tmetadata[,16] ,
                          permutations = 5000)

# A full model

pcap_adonis_f8 <- adonis(tdata ~ tmetadata[,4] *
                           tmetadata[,3] +
                           tmetadata[,5] +
                           tmetadata[,15] +
                           tmetadata[,8]+
                           tmetadata[,16]+
                           tmetadata[,17],
                         permutations = 5000)

# Full model

pcap_adonis_f8

#holm correction

p.adjust(pcap_adonis_f8$aov.tab$`Pr(>F)`, method = "holm")

pcap_adonis_f8.adj <- cbind(pcap_adonis_f8$aov.tab,
      p.adj = p.adjust(pcap_adonis_f8$aov.tab$`Pr(>F)`, method = "holm")
    )

rownames(pcap_adonis_f8.adj) <- c(
  "Exposure",
  "DPE",
  "Sex",
  "Condition Factor",
  "Total Burden",
  "Tank",
  "Total histopath",
  "Exposure:DPE",
  "Residuals",
  "Total"
)

pcap_adonis_f8.adj

# Individual models

days.adonis
exposure.adonis
sex.adonis
burden.adonis
cf.adonis
hist.adonis
tank.adonis


# MICROBIOME: Dispersion index --------------------------------------------

# Next I quantify the dispersion of burden data. Dispersion index D is defined
# as variance / mean

dispersion.index <-
  var(tmetadata[which(tmetadata$Exposure=="Exposed"), "Total"]) /
  mean(tmetadata[which(tmetadata$Exposure=="Exposed"), "Total"])

# MICROBIOME: Parasite prevalence ------------------------------------------

# Here we are calculating parasite prevalence by day

tmetadata.exp <- tmetadata[which(tmetadata$Exposure == "Exposed"),]
parasite_prev <- cbind(days= tmetadata.exp$DaysPE,
                       total = as.numeric(as.logical(tmetadata.exp$Total))
)

parasite_prev <- as.data.frame(parasite_prev)


x <- aggregate(parasite_prev$total ~ parasite_prev$days,
               FUN = length)$`parasite_prev$total`

y <- aggregate(parasite_prev$total ~ parasite_prev$days,
               FUN = sum)$`parasite_prev$total`

z <- aggregate(parasite_prev$total ~ parasite_prev$days,
               FUN = sum)$`parasite_prev$days`

x.y <- y/x

parasite_prev <- as.data.frame(cbind(days = z, prev = x.y) )


parasite_prev.plot <- ggplot(parasite_prev, aes(x = factor(days),
                                                y = prev,
                                                group = 1))

parasite_prev.plot +
  geom_line(size = 2, linetype = 1)  +
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black")
  ) +
  xlab("Days PE")+
  ylab("Prevalence (Exposed)") +
  scale_y_continuous(labels = scales::percent)

# Microbiome sex and burden -----------------------------------------------

# Quantify the effect of sex on burden.
wilcox.test(tmetadata.exp$Total ~ tmetadata.exp$Sex)


# Next we determine if sex is evenly distributed amongst treatment groups and
# days.
aggregate(Sex ~ Exposure + DaysPE, data = tmetadata, table)

aggregate(Sex ~ Exposure, data = tmetadata, table)

#Same question for length of fish

aggregate(Length_.mm. ~ Exposure, data = tmetadata, mean)

aggregate(Length_.mm. ~ Exposure, data = tmetadata, sd)



# FISH PARAM ANALYSIS -----------------------------------------------------

# Fish Weight by day
fish_wt <- ggplot(subset.data.frame(pcap_worm_counts, DaysPE != 0),
                  aes(x = factor(DaysPE,
                                 # levels =c(0,7,10,21,30,43,59,86)),
                                 levels =c(7,10,21,30,43,59,86)),
                      y = Weight_.mg.,
                      fill = factor(Exposure,
                                    levels = c("Unexposed", "Exposed")
                      )
                  )
)

fish_wt +
  geom_boxplot()+
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank()
  ) +
  xlab("Days PE")+
  ylab("Weight(mg)")+
  scale_fill_manual(values = c("#217CA3", "#D24136"))+
  scale_y_continuous(expand = c(0.05,0))

# Fish Weight by day and sex

fish_wt <- ggplot(na.omit(subset.data.frame(pcap_worm_counts, DaysPE != 0)),
                  aes(x = factor(DaysPE,
                                 levels =c(7,10,21,30,43,59,86)),
                      y = Weight_.mg.,
                      fill = factor(Exposure,
                                    levels = c("Unexposed", "Exposed"))
                  )

)

fish_wt +
  geom_boxplot()+
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank()
  ) +
  xlab("Days PE")+
  ylab("Weight(mg)")+
  scale_fill_manual(values = c("#217CA3", "#D24136")) +
  facet_grid(.~Sex)+
  scale_y_continuous(expand = c(0,2))


# Total burden by day as quantified by wet mount

fish_total <- ggplot(na.omit(subset.data.frame(pcap_worm_counts, DaysPE != 0 &
                                                 Exposure == "Exposed")),
                     aes(x = factor(DaysPE,
                                    levels =c(7,10,21,30,43,59,86)),
                         y = Total
                      )
)

fish_total +
  geom_boxplot(fill = "#4D648D")+
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank()
  ) +
  xlab("Days PE")+
  ylab("Total Worms (wet mount)")+
  scale_y_continuous(expand = c(0,0.25))


mature_worms_by_time.plot <- ggplot(
  tmetadata[which(tmetadata$Exposure == "Exposed"),],
  aes(x = factor(DaysPE),
      y = Mature))

# Fish total worm burden by sex as quantified by wet mount

fish_total <- ggplot(na.omit(subset.data.frame(pcap_worm_counts, DaysPE != 0 &
                                                 Exposure == "Exposed")),
                     aes(x = factor(DaysPE,
                                    levels =c(7,10,21,30,43,59,86)),
                         y = Total,
                         fill = factor(Sex)
                     )
)


fish_total +
  geom_boxplot()+
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank()
  ) +

  xlab("Days PE")+
  ylab("Total Worms (wet mount)")+
  scale_fill_manual(values = c("#75B1A9", "#F2C057"))+
  scale_y_continuous(expand = c(0,0.2))

# Fish mature worm burden

mature_worms_by_time.plot +
  geom_boxplot(fill = "#4D648D") +
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black")
  ) +
  xlab("Days PE")+
  ylab("Mature Worms")

# Fish immature worm burden

immature_worms_by_time.plot <- ggplot(
  tmetadata[which(tmetadata$Exposure == "Exposed"),],
  aes(x = factor(DaysPE),
      y = Immature))

immature_worms_by_time.plot +
  geom_boxplot(fill = "#4D648D") +
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black")
  ) +
  xlab("Days PE")+
  ylab("Immature Worms")


# Burden distribution

fish_total <- ggplot(na.omit(subset.data.frame(pcap_worm_counts, DaysPE != 0 &
                                                 Exposure == "Exposed")),
                     aes(x = Total#,
                      )
)

fish_total +
  geom_histogram(binwidth = 1, fill = "#75B1A9", colour = "black")+
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  xlab("Total Worms")+
  ylab("Count")+
  scale_y_continuous(expand = c(0,0), limits = c(0,25))+
  scale_x_continuous(expand = c(0,0),
                     breaks = round(seq(0,35, by = 5), 1 ))


# Total worm burden by sex as quantified by histology

fish_hw <- ggplot(na.omit(subset.data.frame(pcap_worm_counts, DaysPE != 0 &
                                              Exposure == "Exposed")),
                  aes(x = factor(DaysPE,
                                 levels =c(7,10,21,30,43,59,86)),
                      y = Worms.histo,
                      fill = factor(Sex)
                  )
)


fish_hw +
  geom_boxplot()+
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank()
  ) +
  xlab("Days PE")+
  ylab("Total Worms (Histology)")+
  scale_fill_manual(values = c("#75B1A9", "#F2C057"))+
  scale_y_continuous(expand = c(0,.25))


# Gut hyperplasia by day
fish_hyper <- ggplot(na.omit(subset.data.frame(pcap_worm_counts, DaysPE != 0 &
                                                 Exposure == "Exposed")),
                     aes(x = factor(DaysPE,
                                    levels =c(7,10,21,30,43,59,86)),
                         y = Hyperplasia
                     )
)

fish_hyper +
  geom_boxplot(fill = "#4D648D")+
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank()
  ) +
  xlab("Days PE")+
  ylab("Hyperplasia")+
  scale_y_continuous(expand = c(0,0.05))

#Gut inflammation by day

fish_inf <-  ggplot(na.omit(subset.data.frame(pcap_worm_counts, DaysPE != 0 &
                                                Exposure == "Exposed")),
                    aes(x = factor(DaysPE,
                                   levels =c(7,10,21,30,43,59,86)),
                        y = Inflammation
                    )
)

fish_inf +
  geom_boxplot(fill = "#4D648D")+
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank()
  ) +
  xlab("Days PE")+
  ylab("Inflammation")+
  scale_y_continuous(expand = c(0,0.05))

# Total histopathology score by day

total_histo.df <- tmetadata[
  which(tmetadata$Exposure == "Exposed"),]

total_histo.plot <- ggplot(total_histo.df,
                           aes(x = factor(DaysPE),
                               y = histo_total))

total_histo.plot +
  geom_boxplot(fill = "#4D648D")+
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank()
  ) +
  xlab("Days PE")+
  ylab("Total Histopathology Score")+
  scale_y_continuous(expand = c(0,0.25))


# Fish condition factor by day

fish_cf <- ggplot(subset.data.frame(pcap_worm_counts, DaysPE != 0),
                  aes(x = factor(DaysPE,
                                 levels =c(7,10,21,30,43,59,86)),
                      y = condition_factor,
                      fill = factor(Exposure,
                                    levels = c("Unexposed", "Exposed"))
                  )
)

fish_cf +
  geom_boxplot()+
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank()
  ) +
  xlab("Days PE")+
  ylab("Condition Factor") +
  scale_fill_manual(values = c("#217CA3", "#D24136"))+
  scale_y_continuous(expand = c(0,0.05))



# Phylotyping -------------------------------------------------------------

# Next we aggregate counts of sequence variants that share the same label.

taxadf <- taxadf[which(rownames(taxadf) %in% colnames(pcap_microbiome$data)),
                 ,drop=F]

kingdom.df <- replicate(length(unique(taxadf[,2])), c())
names(kingdom.df) <- unique(taxadf[,2])
phylum.df  <- replicate(length(unique(taxadf[,3])), c())
names(phylum.df) <- unique(taxadf[,3])
class.df   <- replicate(length(unique(taxadf[,4])), c())
names(class.df) <- unique(taxadf[,4])
order.df   <- replicate(length(unique(taxadf[,5])), c())
names(order.df) <- unique(taxadf[,5])
family.df  <- replicate(length(unique(taxadf[,6])), c())
names(family.df) <- unique(taxadf[,6])
genus.df   <- replicate(length(unique(taxadf[,7])), c())
names(genus.df) <- unique(taxadf[,7])

for(i in 1:nrow(taxadf)){

  kingdom.df[[taxadf[i, 2]]] <- c(kingdom.df[[taxadf[i, 2]]], taxadf[i,1])
  phylum.df[[taxadf[i, 3]]]  <- c(phylum.df[[taxadf[i, 3]]], taxadf[i,1])
  class.df[[taxadf[i, 4]]]   <- c(class.df[[taxadf[i, 4]]], taxadf[i,1])
  order.df[[taxadf[i, 5]]]   <- c(order.df[[taxadf[i, 5]]], taxadf[i,1])
  family.df[[taxadf[i, 6]]]  <- c(family.df[[taxadf[i, 6]]], taxadf[i,1])
  genus.df[[taxadf[i, 7]]]   <- c(genus.df[[taxadf[i, 7]]], taxadf[i,1])

}

tax.obj <- NULL
tax.obj$kingdom <- kingdom.df
tax.obj$phylum  <- phylum.df
tax.obj$class   <- class.df
tax.obj$order   <- order.df
tax.obj$family  <- family.df
tax.obj$genus   <- genus.df



# aggregate phylotype counts (not really a df, actually an obj)
phylotype.df <- phylotype_analysis(pcap_microbiome, tax.obj)


