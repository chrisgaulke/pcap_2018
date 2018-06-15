# Title and author information --------------------------------------------
#!/usr/bin/R


#########################################
#                                       #
# gaulke_et_al_pcap_microbiome_models.R #
#                                       #
#########################################


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
# links to this repository.

# Versions Options and Packages -------------------------------------------

#Packages

getRversion()
library(ggplot2)
library(MASS)
library(R.utils)
library(glmmTMB)
library(bbmle)
library(reshape2)
library(randomForest)

# options
options("stringsAsFactors"=F)

# FUNCTIONS ---------------------------------------------------------------

#------------------------------------------------#
#           Function adjust_colnames            #
#------------------------------------------------#

adjust_colnames <- function(df){
  # adjust colnames of dada2 sequence table to be colnames, not nt sequences
  x <- paste0(rep("seq", times = ncol(df)), c(1:ncol(df)))
  colnames(df) <- x
  return(df)
}

#------------------------------------------------#
#             Function get_tank                  #
#------------------------------------------------#


get_tank <- function(x){
  # this function will parse the sample id of the model data frame
  # and return a numeric value of tank. NC tanks will be 1-3, PC tanks
  # 4-6 to avoid confusion in the model.
  if(length(grep("PC", x)) == 0){
    y <- as.numeric(strsplit(strsplit(x, "C")[[1]][2], "-")[[1]][1])
  }else{
    y <- as.numeric(strsplit(strsplit(x, "C")[[1]][2], "-")[[1]][1]) + 3
  }
  return(y)
}

#------------------------------------------------#
#             Function get_coeffs                #
#------------------------------------------------#

get_coeffs <- function(obj,
                       exclude = NULL,
                       k1= 1,
                       k2 = 4,
                       length=FALSE,
                       intercept =FALSE){
  # this function will take a micro model object and find extract the
  # coefficients and p values for those coeffs. exclude is a list of fit names
  # that you wish to be exclude (e.g. "fit0")
  # k1: column of coeff (estimate)
  # k2: column of pval
  # length: the number of parameters in model (exlcuding intercept). Use when
  # using multiple model fits with different numbers of parameters
  # intercept: should coeff and pvals be collected for intercept?

  df <- NULL

  if(length){
    if(intercept){
      for(i in 1:length(obj$fit)){
        if(obj$bfs[[i]] %in% exclude){
          next
        }
        sum.fit   <- summary(obj$fit[[i]][[obj$bfs[[i]]]])
        name      <- names(obj$fit)[i]
        fit.coeff <- sum.fit$coefficients[1:nrow(sum.fit$coefficients),k1]
        coef.pval <- sum.fit$coefficients[1:nrow(sum.fit$coefficients),k2]
        if(length(fit.coeff) < length){
          add.x <- length - length(fit.coeff)
          fit.coeff <- c(fit.coeff, rep(NA, times = add.x))
          coef.pval <- c(coef.pval, rep(NA, times = add.x))
        }
        df[[name]]  <- c(fit.coeff, coef.pval)
      }
    }else{
      for(i in 1:length(obj$fit)){
        if(obj$bfs[[i]] %in% exclude){
          next
        }
        sum.fit   <- summary(obj$fit[[i]][[obj$bfs[[i]]]])
        name      <- names(obj$fit)[i]
        fit.coeff <- sum.fit$coefficients[2:nrow(sum.fit$coefficients),k1]
        coef.pval <- sum.fit$coefficients[2:nrow(sum.fit$coefficients),k2]
        if(length(fit.coeff) < length){
          add.x <- length - length(fit.coeff)
          fit.coeff <- c(fit.coeff, rep(NA, times = add.x))
          coef.pval <- c(coef.pval, rep(NA, times = add.x))
        }
        df[[name]]   <- c(fit.coeff, coef.pval)
      }
    }
  }else{
    if(intercept){
      for(i in 1:length(obj$fit)){
        if(obj$bfs[[i]] %in% exclude){
          next
        }
        sum.fit   <- summary(obj$fit[[i]][[obj$bfs[[i]]]])
        name      <- names(obj$fit)[i]
        fit.coeff <- sum.fit$coefficients[1:nrow(sum.fit$coefficients),k1]
        coef.pval <- sum.fit$coefficients[1:nrow(sum.fit$coefficients),k2]
        df[[name]]   <- c(fit.coeff, coef.pval)
      }
    }else{
      for(i in 1:length(obj$fit)){
        if(obj$bfs[[i]] %in% exclude){
          next
        }
        sum.fit   <- summary(obj$fit[[i]][[obj$bfs[[i]]]])
        name      <- names(obj$fit)[i]
        fit.coeff <- sum.fit$coefficients[2:nrow(sum.fit$coefficients),k1]
        coef.pval <- sum.fit$coefficients[2:nrow(sum.fit$coefficients),k2]
        df[[name]]   <- c(fit.coeff, coef.pval)
      }
    }
  }
  return(as.data.frame(df))
}

#------------------------------------------------#
#        Function phylotype_analysis             #
#------------------------------------------------#


phylotype_analysis <- function(obj, tax){
  #obj: microbiome object with at least 1 slot (data)
  #tax: a tax object (named list taxa as names values in the list are seq ids)
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
#             Function make_dict                 #
#------------------------------------------------#

#This function makes a lookup table which maps IDs
#generated as in adjust_colnames with the sequence variant names


make_dict <- function(df){
  x <- cbind(id = paste0(rep("seq", times = ncol(df)), c(1:ncol(df))),
             sequence = colnames(df)
  )
  return(x)
}


#------------------------------------------------#
#           Function get_best_model              #
#------------------------------------------------#

get_best_model <- function(l) {
  #this function will select the best model for list "l" using AICctab
  #this function is intended to be applied to an object with many models
  #(see test data for examples).
  #note: the list of objects should be named, or else AICctab will rename them
  #it is doubtful that this behavoir would be desired.
  #Value: the name of the fit with the lowest AIC
  l <- l
  aic.obj <- AICctab(l, mnames = names(l))
  fit <- attr(aic.obj, "row.names")[1]
  return(fit)
}

#------------------------------------------------#
#          Function get_prevalence               #
#------------------------------------------------#

get_prevalence <- function(df, id.col, value.col){
  #this function will take an long form data frame and
  #calculate the prevalence for each taxa in id.col
  #id.col: the name of the column that contains taxa
  #value.col: the column that contains the taxa abd information

  prev.vec <- NULL
  seq_id <- unique(df[[id.col]])
  for(i in 1:length(seq_id)){
    sub.vec <- unlist(df[which(df[,id.col] == seq_id[i]),value.col])
    #turn into a logical (i.e., T/F). Since T=1 and F=0 summing gives us an
    #count of the individuals with non-zero abd for this taxa
    seq.prev <- sum(as.logical(sub.vec)) / length(sub.vec)
    prev.vec <- c(prev.vec, seq.prev)
  }
  names(prev.vec) <- seq_id
  return(prev.vec)
}

#------------------------------------------------#
#          Function check_converge               #
#------------------------------------------------#

check_converge <- function(obj, model.name, j, k) {
  #this function is to be use in conjuction with get_glmm_bf to ensure that
  #the model with the best AIC also converges.

  #check to make sure that the model is not NA, that sdr obj exists
  # and that the model is not of type null
  if(!is.na(obj$fit[[j]]$AICctab$dAICc[k]) &
     !is.null(obj$fit[[j]][[model.name]]$sdr$pdHess) &
     !is.null(obj$fit[[j]][[model.name]])){
    #then check that the model converges (pdHess and checking that there aren't
    #NAs in the pvals)
    if(obj$fit[[j]][[model.name]]$sdr$pdHess &
       !is.na(summary(obj$fit[[j]][[model.name]])$coefficients$cond[1,4])){
      name.select <- TRUE

    }else{
      name.select <- FALSE
    }

  }else{
    name.select <- FALSE
  }
  #return a true or false
  return(name.select)
}

#------------------------------------------------#
#            Function get_glmm_bf                #
#------------------------------------------------#

get_glmm_bf <- function(obj){
  #this function will find the best fitting convergent model from glmm mixed
  #models and will return named vector with the name of taxa and the best_fit
  #model or an NA
  #obj: should be a list of lists of glmmtmb objects

  obj <- obj
  fit.names <- NULL
  fit.vec   <- NULL

  for(i in 1:length(obj$fit)){
    print(i)
    fit.aic <- obj$fit[[i]]$AICctab
    name <- names(obj$fit)[i]
    fit.names <- c(fit.names, name)

    if(length(fit.aic) == 0){
      #i.e., no models were able to be fit
      fit.vec   <- c(fit.vec, NA)
    }else if(length(obj$fit[[i]]$AICctab) == 1){
      #i.e., only one model was able to be fit so AIC could not be computed
      #because the length is one, we know that only one name will be associated
      #with the AICctab
      model.name <- obj$fit[[i]]$AICctab
      if(obj$fit[[i]][[model.name]]$sdr$pdHess){
        #model converges
        fit.vec  <- c(fit.vec, model.name)
      }else{
        #model does not converge
        fit.vec  <- c(fit.vec, NA)
      }
    }else{
      if(all(is.na(obj$fit[[i]]$AICctab$dAICc))){
        #no models converge
        fit.vec <- c(fit.vec, NA)
      }else{
        #check for convergence. Note that in some instances even though the
        #AIC is the best for a model it fails to converge. So we go through
        #each entry in AICctab to look for the best convergent model starting
        #with the model with the best AIC and working towards that with the
        #worst AIC. Check convergence is kinda a hacky solution, but looking
        #at the glmmtmb docs it seems like there aren't obvious methods for
        #checking convergence of models, likely because it is doubtful that the
        #creators of the package foresaw some applying this procedure to
        #hundreds of taxa. Regardless, it works.

        m1 <- attr(obj$fit[[i]]$AICctab, "row.names")[1]
        m2 <- attr(obj$fit[[i]]$AICctab, "row.names")[2]
        m3 <- attr(obj$fit[[i]]$AICctab, "row.names")[3]
        m4 <- attr(obj$fit[[i]]$AICctab, "row.names")[4]

        if(check_converge(obj = obj, model.name = m1, j = i, k = 1)) {
          fit.vec <- c(fit.vec, m1)
        }else if(check_converge(obj = obj, model.name = m2, j = i, k = 2)){
          fit.vec <- c(fit.vec, m2)
        }else if(check_converge(obj = obj, model.name = m3, j = i, k = 3)){
          fit.vec <- c(fit.vec, m3)
        }else if(check_converge(obj = obj, model.name = m4, j = i, k = 4)){
          fit.vec <- c(fit.vec, m4)
        }else {
          fit.vec <- c(fit.vec, NA)
        }
      }
    }
  }

  df <- cbind(fit.names,fit.vec)
  return(df)
}

#------------------------------------------------#
#            Function make_tax_obj               #
#------------------------------------------------#

make_tax_obj <- function(df, obj){
  #This function will generate lists of seqs at each taxonomic level for
  #use with phylotype analysis
  #df is taxa df
  #obj is microbiome object should have data as a slot

  taxadf <- df[which(rownames(df) %in% colnames(obj$data)),
               ,drop=F]

  kingdom.df <- replicate(length(unique(taxadf[,1])), c())
  names(kingdom.df) <- unique(taxadf[,1])
  phylum.df  <- replicate(length(unique(taxadf[,2])), c())
  names(phylum.df) <- unique(taxadf[,2])
  class.df   <- replicate(length(unique(taxadf[,3])), c())
  names(class.df) <- unique(taxadf[,3])
  order.df   <- replicate(length(unique(taxadf[,4])), c())
  names(order.df) <- unique(taxadf[,4])
  family.df  <- replicate(length(unique(taxadf[,5])), c())
  names(family.df) <- unique(taxadf[,5])
  genus.df   <- replicate(length(unique(taxadf[,6])), c())
  names(genus.df) <- unique(taxadf[,6])

  for(i in 1:nrow(taxadf)){

    kingdom.df[[taxadf[i, 1]]] <- c(kingdom.df[[taxadf[i, 1]]],
                                    rownames(taxadf)[i])
    phylum.df[[taxadf[i, 2]]]  <- c(phylum.df[[taxadf[i, 2]]],
                                    rownames(taxadf)[i])
    class.df[[taxadf[i, 3]]]   <- c(class.df[[taxadf[i, 3]]],
                                    rownames(taxadf)[i])
    order.df[[taxadf[i, 4]]]   <- c(order.df[[taxadf[i, 4]]],
                                    rownames(taxadf)[i])
    family.df[[taxadf[i, 5]]]  <- c(family.df[[taxadf[i, 5]]],
                                    rownames(taxadf)[i])
    genus.df[[taxadf[i, 6]]]   <- c(genus.df[[taxadf[i, 6]]],
                                    rownames(taxadf)[i])

  }

  tax.obj <- NULL
  tax.obj$kingdom <- kingdom.df
  tax.obj$phylum  <- phylum.df
  tax.obj$class   <- class.df
  tax.obj$order   <- order.df
  tax.obj$family  <- family.df
  tax.obj$genus   <- genus.df

  return(tax.obj)
}


# IMPORT DATA -------------------------------------------------------------
wd <- "/Users/gaulkec/Chris/dev/R_projects/pseudocap_long_2017/"

all_seqs_model.df <- read.table(paste0(wd,
              "analysis/flat_files/model_dfs/all_model_2018_06_06.txt"),
                                header = T, row.names = 1 )

all_genus_model.df <- read.table(paste0(wd,
            "analysis/flat_files/model_dfs/all_genus_model_2018_06_06.df"),
                                 header = T, row.names = 1 )


# BURDEN by GENUS ----------------------------------------------------

# Our prior work has demonstrated that the distribution of parasite burden in
# zebrafish follows a negative binomial distribution. Thus we used negative
# binomial generalized linear models to quantify the relationship between
# parasite burden and microbial abundance for each genus. I will build two
# models and then use AIC to determine the model that is likely the best fit.
# The first model will model burden as a function of parasite exposure and time
# post exposure. The second model will incorporate the abundance of an
# individual genus as an additional parameter. We reasons that If the "best"
# model is the fit it suggest that this taxa may influence  parasite burden.

#initialize a multi-model object
burden_by_genus <- NULL
burden_by_genus$fit <- NULL

#collect a list of genera names
seq_id <- unique(all_genus_model.df$variable)

for(i in 1:length(seq_id)){
  #subset the data and remove all rows with NAs
  df <- na.omit(subset(all_genus_model.df, variable == seq_id[i]))

  if(sum(df$value) == 0){
    next
  }
  fit0 <- glm.nb(Total ~
                   factor(Exposure, levels = c("Unexposed", "Exposed"))+
                   as.numeric(DaysPE),
                 data = df )
  fit1 <- glm.nb(Total ~ scale(value)+
                   factor(Exposure, levels = c("Unexposed", "Exposed"))+
                   as.numeric(DaysPE),
                 data = df )
  burden_by_genus$fit[[seq_id[i]]] <- list(fit0=fit0,fit1=fit1)
}

#Get name of best models (i.e., the model with the lowest AIC)

burden_by_genus$bfs <- lapply( burden_by_genus$fit, get_best_model)

#grab the coefficients of the best model
burden_by_genus$coef <- get_coeffs(burden_by_genus,
                                   exclude = c("fit0"),
                                   length=3,#changed from 4
                                   intercept = FALSE )


rownames(burden_by_genus$coef) <- c("abd_est",
                                    "exposure_ex_est",
                                    "days_pe_est",
                                    "abd_pval",
                                    "exposure_ex_pval",
                                    "days_pe_pval"
)

#how many slopes are positive (1), and negative (0)
table(as.numeric(ifelse(burden_by_genus$coef[1,] > 0 , TRUE, FALSE)))

coef.df <- cbind(melt(t(burden_by_genus$coef[1:3,])),
                 melt(t(burden_by_genus$coef[4:6,]))[,2:3])

colnames(coef.df) <- c("taxa",
                       "estimate_type",
                       "estimate_value",
                       "pval_type",
                       "pval_value" )

coef.df$qval <- qvalue::qvalue(coef.df$pval_value)$qvalue


#Next we will remove low prevalence taxa.
genus.prev <- get_prevalence(all_genus_model.df,
                             id.col = "variable",
                             value.col = "value")

coef.names <- colnames(burden_by_genus$coef[,
                       order(burden_by_genus$coef["abd_est",])])
coef.df <- coef.df[which(coef.df$taxa
                         %in% names(genus.prev[which(genus.prev > .1)])),]

coef.names <-
  coef.names[which(coef.names %in% names(genus.prev[which(genus.prev > .1)]))]


coef.df$rescale_estimate <- cut(coef.df$estimate_value,
                                breaks=c(-10,-5,-1,-.1, -0.01,
                                         0, 0.01, 0.1, 1, 5,10),
                                labels=c("<-5",
                                         "-5 - -1",
                                         "-1 - -0.1",
                                         "-0.1 - -0.01",
                                         "-0.01 - 0",
                                         "0 - 0.01",
                                         "0.01 - 0.1",
                                         "0.1 - 1",
                                         "1 - 5",
                                         ">5")
)



coef.df$rescale_pval <- cut(coef.df$pval_value,
                            breaks=c(1e-100, 1e-3, 1e-2, 5e-2, 1),
                            labels=c("***",
                                     "**",
                                     "*",
                                     "")
)

coef.df$rescale_qval <- cut(coef.df$qval,
                            breaks=c(1e-100, 0.15,1),
                            labels=c("*",
                                      "")
)


#reverse levels so that this will show up right in scale
coef.df$rescale_estimate <- factor(as.character(coef.df$rescale_estimate),
                             levels = rev(levels(coef.df$rescale_estimate)))


coef.df$taxa <- factor(coef.df$taxa,
                       levels = coef.names)

#Remove the "NONE" genus which simply represents all taxa for which we did not
#have enough confidence to call a genus.

coef.df <- coef.df[-which(coef.df$taxa == "NONE"),]

my_reds  <- colorRampPalette(colors = c("white", "red"))
my_blues <- colorRampPalette(colors = c("blue", "white"))

#generate a heatmap

burden_by_genus.tile <- ggplot(coef.df,
                               aes(x= estimate_type,
                                   y = taxa,
                                   fill = rescale_estimate))

burden_by_genus.tile +
  geom_tile()+
  scale_fill_manual(values=c(my_blues(100)[c(1,50,70,85,95)],
                             my_reds(100)[c(5,15,30,50,99)]),
                    na.value="grey")+
  scale_x_discrete(expand = c(0,0),
                   labels = c("Abd",
                              "Exp(exposed)",
                              "Days PE",
                              "Abd:Exp(exposed)"))+
  scale_y_discrete(expand = c(0,0))+
  geom_text(aes(label=rescale_qval),
            nudge_y = -.2)+
  theme(text = element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = .9, hjust = 1, size =16),
        axis.text = element_text(color = "black"),
        legend.title=element_blank()
  )+
  ylab("")+
  xlab("")

# GENUS by BURDEN and HISTO -----------------------------------------------

# To determine if taxonomic abundance is assocaited with parasite burden or
# histopathological changes that accompany parasite infection. Here I use
# zero inflatted negative binomial generalized linear mixed effects models to
# quantify these assocaitions. The glmmtmb package has two unique negative
# binomial distributions that differ in their expectation of mean-variance
# relationship. So I fit a single model using each of these distributions and
# use AIC as above to determine which is the "best" model. As individuals in
# a tank may share more similar gut microbiome structure I incorporate tank
# as a random effect.

genus.nbglmm <- NULL
genus.nbglmm$fit <- NULL

seq.vec <- unique(all_genus_model.df$variable)

for(i in 1:length(seq.vec)){

  df <- all_genus_model.df[ which(all_genus_model.df$variable %in%
                                    c(seq.vec[i])), ]

  if(sum(df$value) == 0){
    next
  }

  fit0 <- NULL
  fit2 <- NULL

  try(fit0 <- glmmTMB(value ~ factor(Exposure,
                                     levels = c("Unexposed", "Exposed")) +
                        DaysPE +
                        Total +
                        histo_total +
                        (1|tank),
                      data = na.omit(df),
                      ziformula = ~1,
                      family = nbinom2), silent = T
  )

  try(fit2 <- glmmTMB(value ~ factor(Exposure,
                                     levels = c("Unexposed", "Exposed")) +
                        DaysPE +
                        Total +
                        histo_total +
                        (1|tank),
                      data = na.omit(df),
                      ziformula = ~1,
                      family = nbinom1), silent = T
  )

  genus.nbglmm$fit[[seq.vec[i]]] <- list(nb2.0=fit0,
                                         nb1.0=fit2
  )

}

# Now the goal is to 1) calculate AICctab for each set of models, 2) Id the
# best model for each taxa, 3) collect the coefficients for each of the best
# models, 4) collect the mixed effects and zero inflation parameters for each
# best fit model.

#Calculate AIC

for(i in 1:length(genus.nbglmm$fit)){
  print(i)
  nb2.0 <- genus.nbglmm$fit[[i]][[1]]
  nb1.0 <- genus.nbglmm$fit[[i]][[2]]

  #check if any of the models are null (no fit)
  if(any(is.null(nb2.0),is.null(nb1.0))){
    my_list <- list(nb2.0 = nb2.0,
                    nb1.0 = nb1.0
    )
    #this will return a logical vector. TRUE !is.null, FALSE is.null
    my_ex   <- !sapply(list(nb2.0,nb1.0),is.null)
    my_list <- my_list[my_ex]

    if(length(my_list) == 0){
      genus.nbglmm$fit[[names(genus.nbglmm$fit)[i]]]$AICctab <- NULL
    }else if(length(my_list) == 1){
      genus.nbglmm$fit[[names(genus.nbglmm$fit)[i]]]$AICctab <-
        names(my_list)[[1]]
    }
  }else{
    genus.nbglmm$fit[[names(genus.nbglmm$fit)[i]]]$AICctab <-
      AICctab(nb2.0,nb1.0)
  }
}


#now i get the best fits for each model
genus_fits.df <- get_glmm_bf(genus.nbglmm)
genus_fits.df <- as.data.frame(genus_fits.df)

genus.nbglmm$conv_fixed_coef <- NULL
genus.nbglmm$raneff          <- NULL
genus.nbglmm$conv_zero_coef  <- NULL

rownames(genus_fits.df) <- genus_fits.df$fit.names
fit.vec <- genus_fits.df$fit.vec
converged.vec <- rownames(na.omit(genus_fits.df))

for(i in 1:length(fit.vec)){
  print(i)
  #check to make sure that the model converged, if not move on
  if(names(genus.nbglmm$fit[i]) %in% converged.vec){
    fit   <- genus.nbglmm$fit[[i]][[fit.vec[i]]]
    fit.s <- summary(fit)
    #check nrows of coeff to determine if an NA needs to be added for the
    #interaction term
    fit.cond_est <- fit.s$coefficients$cond[,1]
    fit.cond_p   <- fit.s$coefficients$cond[,4]
    fit.zero_est <- fit.s$coefficients$zi[,1]
    fit.zero_p   <- fit.s$coefficients$zi[,4]
    fit.ranef    <- ranef(fit)
    nterms       <- nrow(fit.s$coefficients$cond)
    name         <- names(genus.nbglmm$fit[i])
    if(nterms == 5){
      fit.cond_est <- c(fit.cond_est, NA)
      fit.cond_p   <- c(fit.cond_p, NA)
    }
    genus.nbglmm$conv_fixed_coef[[name]] <- c(fit.cond_est,fit.cond_p)
    genus.nbglmm$raneff[[name]]          <- fit.ranef
    genus.nbglmm$conv_zero_coef[[name]]  <- c(fit.zero_est,fit.zero_p)
  }else{
    next
  }
}

#clean up and add rownames
genus.nbglmm$conv_fixed_coef <- as.data.frame(genus.nbglmm$conv_fixed_coef)
rownames(genus.nbglmm$conv_fixed_coef) <- c("intercept_est",
                                            "exposure_est",
                                            "days_pe_est",
                                            "burden_est",
                                            "histopath_est",
                                            "histo_burden_int_est",
                                            "intercept_pval",
                                            "exposure_pval",
                                            "days_pe_pval",
                                            "burden_pval",
                                            "histopath_pval",
                                            "histo_burden_int_pval")

#melt into a df for ggplot
genus_coef.df <- cbind(melt(t(genus.nbglmm$conv_fixed_coef[1:6,])),
                       melt(t(genus.nbglmm$conv_fixed_coef[7:12,]))[,2:3])
colnames(genus_coef.df) <- c("taxa",
                             "estimate_type",
                             "estimate_value",
                             "pval_type",
                             "pval_value" )

#remove the intercept pvalues
genus_coef.df <- genus_coef.df[-which(genus_coef.df$estimate_type %in%
                                        c( "intercept_est",
                                           "histo_burden_int_est")),]

genus_coef.df$qval <- qvalue::qvalue(genus_coef.df$pval_value)$qvalue

# Some of these taxa are likely to be those that are very rare (i.e., low
# prevalence) in the fish. As above we can filter these taxa in the interest of
# reducing noise and improving visualization

genus_coef.df <-  genus_coef.df[which(genus_coef.df$taxa
                %in% names(genus.prev[which(genus.prev > .1)])),]

#remove the "NONE" genus, because it really doesn't mean anything

genus_coef.df <- genus_coef.df[-which(genus_coef.df$taxa == "NONE"),]

#As above, in order to plot something that looks nice I need to rescale the
#estimates into ranges

genus_coef.df$rescale_estimate <- cut(genus_coef.df$estimate_value,
                                      breaks=c(-1000,-5,-3,-1,-0.1,
                                               0,
                                               0.1, 1, 3, 5, 1000),
                                      labels=c("<-5",
                                               "-5 - -3",
                                               "-3 - -1",
                                               "-1 - -0.1",
                                               "-0.1 - 0",
                                               "0 - 0.1",
                                               "0.1 - 1",
                                               "1 - 3",
                                               "3 - 5",
                                               ">5")
)

genus_coef.df$rescale_pval <- cut(genus_coef.df$pval_value,
                                  breaks=c(1e-100, 1e-3, 1e-2, 5e-2, 1),
                                  labels=c("***",
                                           "**",
                                           "*",
                                           "")
)
genus_coef.df$rescale_qval <- cut(genus_coef.df$qval,
                                  breaks=c(1e-100,.15, 1),
                                  labels=c("*",
                                           "")
)


# reverse levels so that this will show up right in scale
genus_coef.df$rescale_estimate <-
  factor(as.character(genus_coef.df$rescale_estimate),
         levels = rev(levels(genus_coef.df$rescale_estimate)))

#now I will re level the taxa so they will print nicely

genus_coef.names <- genus_coef.df[which(
  genus_coef.df$estimate_type == "exposure_est"),]
genus_coef.names <- genus_coef.names[order(
  genus_coef.names$estimate_value, decreasing = T),]
genus_coef.names <- as.character(unlist(genus_coef.names$taxa))

genus_coef.df$taxa <- factor(genus_coef.df$taxa,
                             levels = genus_coef.names)


my_reds  <- colorRampPalette(colors = c("white", "red"))
my_blues <- colorRampPalette(colors = c("blue", "white"))

#Plot a heat map of coefficients

genus_by_exposure.tile <- ggplot(genus_coef.df,
                                 aes(x= estimate_type,
                                     y = taxa,
                                     fill = rescale_estimate))

genus_by_exposure.tile +
  geom_tile()+
  scale_fill_manual(values=c(my_blues(100)[c(1,50,70,85,95)],
                             my_reds(100)[c(5,15,30,50,99)]),
                    na.value="grey")+
  scale_x_discrete(expand = c(0,0),
                   labels = c("Exposure",
                              "DaysPE",
                              "Burden",
                              "Histopath",
                              "Burden:Histopath"))+
  scale_y_discrete(expand = c(0,0))+
  geom_text(aes(label=rescale_qval),
            nudge_y = -.2)+
  theme(text = element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 30,
                                   vjust = .9,
                                   hjust = 1,
                                   size =16),
        axis.text = element_text(color = "black"),
        legend.title=element_blank()
  )+
  ylab("")+
  xlab("")

# RANDOM FOREST MICROBIOME: Exposure -------------------------------------

#Filter data for prevalent taxa

genus_rf.df <-  all_genus_model.df[which(all_genus_model.df$variable
                         %in% names(genus.prev[which(genus.prev > .1)])),]


genus_rf_exp.df <- dcast(genus_rf.df, Seq_ID ~ variable)
rownames(genus_rf_exp.df) <- genus_rf_exp.df$Seq_ID

genus_rf_exp.df$exp <- genus_rf.df$Exposure[1:237]
genus_rf_exp.df$Seq_ID <- NULL
genus_rf_exp.df$NONE   <- NULL


set.seed(731)
genus_rf_exp.rf <- randomForest(factor(exp) ~.,data=genus_rf_exp.df,
                                na.action = na.omit,
                                importance = T,
                                ntree =10000)

#What are the most important variables
varImpPlot(genus_rf_exp.rf)

exp_rf_imp.df <- genus_rf_exp.rf$importance
exp_rf_imp.df <- as.data.frame(exp_rf_imp.df)

exp_rf_imp.df <- exp_rf_imp.df[order(exp_rf_imp.df$MeanDecreaseAccuracy,
                                     decreasing = T),]

exp_rf_imp.df <- exp_rf_imp.df[1:15, 3, drop =F]

exp_rf_imp.df$names <- rownames(exp_rf_imp.df)
exp_rf_imp.df$sd <- genus_rf_exp.rf$importanceSD[exp_rf_imp.df$names, 3]



exp_rf_imp.df$mda_scale <- exp_rf_imp.df$MeanDecreaseAccuracy /
  exp_rf_imp.df$sd

exp_rf_imp.df <- exp_rf_imp.df[order(exp_rf_imp.df$mda_scale,decreasing = T),]

exp_rf_imp.df$names <- factor(exp_rf_imp.df$names,
                              levels = rev(exp_rf_imp.df$names))

#Plot the importance of variables

exp_rf_imp.plot <- ggplot(exp_rf_imp.df, aes(x = names,
                                             y= mda_scale ))

exp_rf_imp.plot +
  geom_point(size = 5) +
  coord_flip() +
  xlab("") +
  ylab("Mean Decrease in Accuracy")+
  ylim(c(0,150)) +
  theme(text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        axis.text = element_text(color = "black"),
        legend.title=element_blank()
  )


# RANDOM FOREST MICROBIOME: Exposure test -----------------------------------
pcap_2015_seqtable <-
  paste0("/Users/gaulkec/Chris/dev/R_projects/pseudocap_long_2017/",
         "data/pcap_2015_dada2_rerun/out/sequence_table.txt")

pcap_2015_taxtable <-
  paste0("/Users/gaulkec/Chris/dev/R_projects/pseudocap_long_2017/",
         "data/pcap_2015_dada2_rerun/out/tax_final.txt")


pcap_2015_seqtable.df <- read.table(pcap_2015_seqtable,
                                    sep = "\t",
                                    row.names = 1,
                                    header = T)

#create a mapping of short seqvar id to sequence

seq_id_to_seq_var.df <- make_dict(pcap_2015_seqtable.df)
rownames(seq_id_to_seq_var.df) <- seq_id_to_seq_var.df[,2]
seq_id_to_seq_var.df <- as.data.frame(seq_id_to_seq_var.df)

pcap_2015_seqtable.df <- adjust_colnames(pcap_2015_seqtable.df)



pcap_2015_taxtable.taxadf <- read.table(pcap_2015_taxtable,
                                        sep = "\t",
                                        header = T)

#align colnames for lookup and tax table
seq_id_to_seq_var.df <-
  seq_id_to_seq_var.df[rownames(pcap_2015_taxtable.taxadf),]

#swap out the seq var names so we can do phylotyping
rownames(pcap_2015_taxtable.taxadf) <- seq_id_to_seq_var.df$id

#make new microbiome obj

pilot_parasite <- NULL
pilot_parasite$data <- pcap_2015_seqtable.df
pilot_parasite$data <- normalize(pilot_parasite$data, "rare", depth = 5000)



pcap_2015.taxobj <- make_tax_obj(pcap_2015_taxtable.taxadf,pilot_parasite)

pcap_2015.phylo  <- phylotype_analysis(pilot_parasite, pcap_2015.taxobj)


#since Random forest is going to expect that all of the classifiers are present
#so I have to grab the colnames from the rf data frame and then filter the
#above colnames to conform to that df

#grab the genus
pcap_2015_genus.rf <-  pcap_2015.phylo$genus

#grab colnames to keep
keeps.rf <- colnames(genus_rf_exp.df)
keeps.rf <- keeps.rf[-which(keeps.rf == "exp")]


missing.rf <- keeps.rf[-which(keeps.rf %in% colnames(pcap_2015_genus.rf))]

#filter genera not present in the rf model df
pcap_2015_genus.rf <- pcap_2015_genus.rf[,which(colnames(pcap_2015_genus.rf) %in% keeps.rf)]

#add missing columns that are in model df but not in test df
for(i in 1:length(missing.rf)){
  pcap_2015_genus.rf[[missing.rf[i]]] <- 0
}

genus_rf_exp.predict <- predict(genus_rf_exp.rf, pcap_2015_genus.rf)

# all of these animals were exposed to P. tomentosa so perfect accuracy
# would be 65 exposed 0 unexposed
table(genus_rf_exp.predict)


# Weight and Condition factor ~ exposure ----------------------------------

#grab histopath data for one genus so we can do some stats on it
#note this would be the same no matter what taxa is select as this is a long
#form data table.

wt_mod.df <- all_genus_model.df[which(
  all_genus_model.df$variable == "Megamonas"),]

wt.lm <- lm(
  Weight_.mg. ~ Sex*factor(Exposure, levels = c("Unexposed", "Exposed")) +
    Total + histo_total + Length_.mm.,
  data =na.omit( wt_mod.df))

step.wt <- stepAIC(wt.lm, direction="both")

#The optimal stepAIC model
nwt.lm <- lm(
  Weight_.mg. ~ Sex + histo_total + Length_.mm.,
  data =na.omit( wt_mod.df))

summary(nwt.lm)


cf.lm <- lm(
  condition_factor ~ Sex*factor(Exposure, levels = c("Unexposed", "Exposed")) +
    Total + histo_total, data = na.omit(wt_mod.df))

step.cf <- stepAIC(cf.lm, direction="both")

#The optimal stepAIC model
ncf.lm <- lm(
  condition_factor ~ Sex + histo_total, data = na.omit(wt_mod.df))

summary(ncf.lm)

# Association between Mycoplasma abundance, exposure and hyperplasia ------

mycoplasma_seqs.df <- all_seqs_model.df[which(all_seqs_model.df$variable %in%
                                                c("seq3100",
                                                  "seq1607",
                                                  "seq3845",
                                                  "seq47",
                                                  "seq56",
                                                  "seq3321",
                                                  "seq3106",
                                                  "seq2782",
                                                  "seq21")),]

mycoplasma_seqs.sum <- 0

for( i in 1:length(unique(mycoplasma_seqs.df$variable))) {
  seq <- unique(mycoplasma_seqs.df$variable)[i]
  mycoplasma_seqs.sum <- mycoplasma_seqs.sum +
    mycoplasma_seqs.df[which(mycoplasma_seqs.df$variable == seq), "value"]

}


mycoplasma_clade_abd.df <- cbind(sample_id =
                                   unique(mycoplasma_seqs.df$Sample.ID),
                                 abundance = mycoplasma_seqs.sum )
mycoplasma_clade_abd.df <- as.data.frame(mycoplasma_clade_abd.df)

highlight.clade <- ifelse( mycoplasma_clade_abd.df$sample_id
                           %in% c("PC1-35","PC2-33"),
                           "highlight", "normal")

mycoplasma_clade_abd.df$highlight <- highlight.clade

mycoplasma_clade_abd.df$hyperplasia <- wt_mod.df$Hyperplasia
mycoplasma_clade_abd.df$exposure <- wt_mod.df$Exposure

# Wilcox test between exposed and unexposed
wilcox.test(as.numeric(mycoplasma_clade_abd.df$abundance) ~
              mycoplasma_clade_abd.df$exposure)

# Examine correaltion between mycoplasma and hyperplasia

mycoplasma_clade_abd.kendall <-
  cor.test(as.numeric(mycoplasma_clade_abd.df$abundance),
           as.numeric(mycoplasma_clade_abd.df$hyperplasia),
           method = "kendall" )

# Set highlight colors
mycolours <- c("highlight" = "red", "normal" = "black")

#plot abundance
mycoplasma_clade_abd.plot <- ggplot(mycoplasma_clade_abd.df, aes(x =
                          factor(exposure, levels = c("Unexposed", "Exposed")),
                                           y = as.numeric(abundance)
))

mycoplasma_clade_abd.plot +
  geom_boxplot(outlier.shape = NA,fill = "#4D648D")+
  geom_point(
    size = 2,
    aes(colour = highlight), data =mycoplasma_clade_abd.df
  )+
  scale_colour_manual(guide=FALSE, values = mycolours) +

  theme(text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = "black",size =16),
        legend.title=element_blank()
  ) +
  ylab("Mycoplasma Combined Sequence Abundance") +
  xlab("")

# Now get the summed abundance of the seqs in the burns clade
burns_seqs.sum <- 0

burns.seqs <- c("seq21", "seq2782")

for( i in 1:length(burns.seqs)) {
  seq <- unique(mycoplasma_seqs.df$variable)[i]
  burns_seqs.sum <- burns_seqs.sum +
    mycoplasma_seqs.df[which(mycoplasma_seqs.df$variable == seq), "value"]

}

burns_abd.df <- cbind(sample_id = unique(mycoplasma_seqs.df$Sample.ID),
                                 abundance = burns_seqs.sum )
burns_abd.df <- as.data.frame(burns_abd.df)

highlight.clade <- ifelse( burns_abd.df$sample_id %in% c("PC1-35","PC2-33"),
                           "highlight", "normal")

burns_abd.df$highlight <- highlight.clade

burns_abd.df$hyperplasia <- wt_mod.df$Hyperplasia
burns_abd.df$exposure <- wt_mod.df$Exposure

# Wilcox test between exposed and unexposed
wilcox.test(as.numeric(burns_abd.df$abundance) ~
              burns_abd.df$exposure)

# Examine correaltion between burns clade mycoplasma and hyperplasia

burns_abd.kendall <-
  cor.test(as.numeric(burns_abd.df$abundance),
           as.numeric(burns_abd.df$hyperplasia),
           method = "kendall" )


# Plot abundance

burns_abd.plot <- ggplot(burns_abd.df, aes(x =
                        factor(exposure, levels = c("Unexposed", "Exposed")),
                        y = as.numeric(abundance)
))

burns_abd.plot +
  geom_boxplot(outlier.shape = NA,fill = "#4D648D")+
  geom_point(
  size = 2,
  aes(colour = highlight), data =burns_abd.df
    )+
  scale_colour_manual(guide=FALSE, values = mycolours) +

  theme(text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = "black",size =16),
        legend.title=element_blank()
  ) +
  ylab("Burns Clade Abundance") +
  xlab("")

