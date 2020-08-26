# load packages and data
library(tidyverse)
library(knitr)
library(quickpsy)
library(ez)
library(ggplot2)
library(ggsignif)
library(PairedData)
library(cowplot)
library(dplyr)
library(gtools)
library(ggpubr)

theme_set(bayesplot::theme_default())

## define theme for whole report
mytheme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme_classic() + 
  theme(strip.background = element_blank()) 

## definition of fuction to calulate geometric mean
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


## load data Exp1, Exp2, and Exp3

# Exp 1
dat_exp1 = read.csv('../data/AllData_Exp1.csv') %>% filter(valid  == 1)
dat_exp1$cond = factor(dat_exp1$cond, labels = c("PS",  "NS"))
dat_exp1$curDur = dat_exp1$curDur * 1000;
dat_exp1$NB <- NULL

# Exp 2
dat_exp2 = read.csv('../data/AllData_Exp2.csv') %>% filter(valid  == 1)
dat_exp2$cond = factor(dat_exp2$cond, labels = c("DF", "AF"))
dat_exp2$curDur = dat_exp2$curDur * 1000;
dat_exp2$moda <- NULL
dat_exp2$NB <- NULL

# Exp 3
dat_exp3 = read.csv('../data/AllData_Exp3.csv') %>% filter(valid  == 1)
dat_exp3$cond = factor(dat_exp3$cond, labels = c("U-shaped",  "I T-shaped"))
dat_exp3$curDur = dat_exp3$curDur * 1000;
dat_exp3$moda <- NULL


#merge 3 experiment data together 
dat_all <- rbind(data.frame(rbind(dat_exp1,dat_exp2)), dat_exp3)


# Linear regression of data
lm_model <- function(df) {
    lm(x ~ y, data = df)
}

