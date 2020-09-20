# load R packages
library(parallel)
library(doParallel)
library(rstan)
library(tidyverse)
library(rlist)
library(ggpubr)
library(gtools)
library(reshape)

# Setting
options(mc.cores = parallel::detectCores())
rstan_options (auto_write=TRUE)
# flag for merge model results 
needmerge= 1
# flag for running rstan model and saving the results
runRstanModels = FALSE
# flag for running rstan model locally
runRstanModelslocally = FALSE

threshold = 0.75
#threshold = 0.84
midthreshold = 0.5

### model number 1-3
model123 <- "
  data {
    int nsubjs;
    int nstim;
    int n[2*nsubjs,nstim];
    int r[2*nsubjs,nstim];
    int x[2*nsubjs,nstim];
    real xmean[2,nsubjs];
    real xsd[2,nsubjs];
  }
  parameters {
    real mua;
    real mub;
    real<lower=0,upper=3000> sigmaa;
    real<lower=0,upper=3000> sigmab;
    vector[nsubjs] alpha;
    vector[nsubjs] beta;
  }

  model {
    // Priors
    mua ~ normal(0, inv_sqrt(.001));
    mub ~ normal(0, inv_sqrt(.001));
    alpha ~ normal(mua, sigmaa);
    beta ~ normal(mub, sigmab);
    
    for (i in 1:nsubjs) {
      for (j in 1:nstim) {
        real theta[2];
        theta[1] = inv_logit(alpha[i] + beta[i] * (x[i,j]/xmean[1, i] - 1));  // * wf[i]
        r[i,j] ~ binomial(n[i,j], theta[1]);
        
        theta[2] = inv_logit(alpha[i] + beta[i] * (x[i+nsubjs,j]/xmean[2, i] - 1)); //* wf[nsubjs+i]
        r[i+nsubjs,j] ~ binomial(n[i+nsubjs,j], theta[2]);
      }
    }
  }
"

# compile models
stanmodel123 <- stan_model(model_code = model123, model_name="stanmodel123")





### model number 4
model4 <- "
// Logistic Psychophysical Function
data {
  int nsubjs;
  int nstim;
  int n[2*nsubjs,nstim];
  int r[2*nsubjs,nstim];
  int x[2*nsubjs,nstim];
  real xmean[2,nsubjs];
  real xsd[2,nsubjs];
}
parameters {
  real mua;
  real mub;
  real<lower=0,upper=3000> sigmaa;
  real<lower=0,upper=3000> sigmab;
  real<lower=0,upper=200> sigmac;
  vector[nsubjs] alpha;
  vector[nsubjs] beta;
  real<lower=500,upper=1300>  xmeanhat[2,nsubjs];
}

model {
  mua ~ normal(0, inv_sqrt(.001));
  mub ~ normal(0, inv_sqrt(.001));
  alpha ~ normal(mua, sigmaa);
  beta ~ normal(mub, sigmab);
  for (k in 1:2){
    for (i in 1:nsubjs) {
      xmeanhat[k,i] ~ normal(xmean[k, 1], sigmac);
      for (j in 1:nstim) {
          real theta;
          theta = inv_logit(alpha[i] + beta[i] * (x[i+nsubjs*(k-1),j]/xmeanhat[k, i] - 1)); 
          r[i+nsubjs*(k-1),j] ~ binomial(n[i+nsubjs*(k-1),j], theta);
      }
    }
  }
}"
stanmodel4 <- stan_model(model_code = model4, model_name="stanmodel4")


### model number 5 EDA
model5 <- "
// Logistic Psychophysical Function
data {
  int nsubjs;
  int nstim;
  int n[2*nsubjs,nstim];
  int r[2*nsubjs,nstim];
  int x[2*nsubjs,nstim];
  real xmean[2,nsubjs];
  real xsd[2,nsubjs];
}
parameters {
  real mua[2];
  real mub[2];
  real<lower=0,upper=3000> sigmaa[2];
  real<lower=0,upper=3000> sigmab[2];
  real<lower=0,upper=200> sigmac[2];
  real alpha[2, nsubjs];
  real beta[2, nsubjs];
  real<lower=500,upper=1300>  xmeanhat[2,nsubjs];
}

model {
  for (k in 1:2){
    mua[k] ~ normal(0, inv_sqrt(.001));
    mub[k] ~ normal(0, inv_sqrt(.001));
    alpha[k,] ~ normal(mua[k], sigmaa[k]);
    xmeanhat[k,]~ normal(xmean[k, 1], sigmac[k]);
    beta[k,] ~ normal(mub[k]*xmean[k, 1]/xsd[k,1], sigmab[k]);
    for (i in 1:nsubjs) {
      for (j in 1:nstim) {
          real theta;
          theta = inv_logit(alpha[k,i] + beta[k, i] * (x[i+nsubjs*(k-1),j]/xmeanhat[k, i] - 1));  
          r[i+nsubjs*(k-1),j] ~ binomial(n[i+nsubjs*(k-1),j], theta);
      }
    }
  }
}"
# compile models
stanmodel5 <- stan_model(model_code = model5, model_name="stanmodel5")


### model number 6 EDA (less variables)
model6 <- "
// Logistic Psychophysical Function
data {
  int nsubjs;
  int nstim;
  int n[2*nsubjs,nstim];
  int r[2*nsubjs,nstim];
  int x[2*nsubjs,nstim];
  real xmean[2,nsubjs];
  real xsd[2,nsubjs];
}
parameters {
  real mua;
  real mub;
  real<lower=0,upper=3000> sigmaa;
  real<lower=0,upper=3000> sigmab;
  real<lower=0,upper=200> sigmac;
  vector[nsubjs] alpha;
  real beta[2, nsubjs];
  real<lower=500,upper=1300>  xmeanhat[2,nsubjs];
}

model {
  mua ~ normal(0, inv_sqrt(.001));
  mub ~ normal(0, inv_sqrt(.001));
  alpha ~ normal(mua, sigmaa);
  
  for (k in 1:2){
    for (i in 1:nsubjs) {
      xmeanhat[k,i] ~ normal(xmean[k, 1], sigmac);
      beta[k,i] ~ normal(mub*xmean[k, 1]/xsd[k,1], sigmab);   
      for (j in 1:nstim) {
          real theta;
          theta = inv_logit(alpha[i] + beta[k,i]* (x[i+nsubjs*(k-1),j]/xmeanhat[k, i] - 1)); 
          r[i+nsubjs*(k-1),j] ~ binomial(n[i+nsubjs*(k-1),j], theta);
      }
    }
  }
}"
stanmodel6 <- stan_model(model_code = model6, model_name="stanmodel6")


gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

############################## PSYCHOMETRIC FUNCTIONS ##########################
# only the MAP estimate; use this to plot psychometric functions
F1 <- function(pX, s, x_mean, aMAP, bMAP)
{
  exp(aMAP[s] + bMAP[s]*(pX - x_mean[s]))/
    (1+exp(aMAP[s] + bMAP[s]*(pX - x_mean[s])))
}

F1inv <- function(Y,s, aMAP, bMAP)
{
  (log(-Y/(Y-1))-aMAP[s])/bMAP[s]
}

# function for all the posterior alpha/beta values; use this to calculate JND
# posterior
F2 <- function(X, s, x_mean, a, b)
{
  exp(a[,s] + b[,s]*(X - x_mean[s]))/
    (1+exp(a[,s] + b[,s]*(X - x_mean[s])))
}

F2inv <- function(Y, s, a, b)
{
  (log(-Y/(Y-1))-a[,s])/(b[,s])
}

F3invMAP <- function(Y,s, aMAP, bMAP, x_meanhat)
{
  (log(-Y/(Y-1))-aMAP[s])* x_meanhat/bMAP[s] + x_meanhat
}

F3inv <- function(Y, s, a, b, x_meanhat)
{
  (log(-Y/(Y-1))-a[,s])*x_meanhat/(b[,s]) + x_meanhat
}

# function for 20 grabbed posterior alpha/beta values; use this to plot
# overlapping sigmoids to visualize variance
F3 <- function(X, s, g, x_mean, a_sel,  b_sel)
{
  exp(a_sel[g,s] + b_sel[g,s]*(X - x_mean[s]))/
    (1+exp(a_sel[g,s] + b_sel[g,s]*(X - x_mean[s])))
}

pltaicbicjnd<- function(dat, condlist, modellist){
  newdat <- dat %>% filter(cond %in% condlist)%>%filter(model %in% modellist)
  plt_aicjnd<- newdat %>%
    ggplot(aes(model, aicJnd, color = model, fill = model)) +
    geom_bar(stat = 'identity', width = 0.5) +
    theme_classic()+
    facet_wrap(~cond)
  
  plt_bicjnd<- newdat %>% 
    ggplot(aes(model, bicJnd, color = model, fill = model)) +
    geom_bar(stat = 'identity', width = 0.5) +
    theme_classic()+
    facet_wrap(~cond)
  
  pltaicbic_jnd<- ggarrange(plt_aicjnd, plt_bicjnd, common.legend = TRUE, ncol=2, nrow=1)
  return(pltaicbic_jnd)
}


pltaicbicpse<- function(dat, condlist, modellist){
  newdat <- dat %>% filter(cond %in% condlist)%>%filter(model %in% modellist)
  plt_aicpse<- newdat %>%
    ggplot(aes(model, aicPse, color = model, fill = model)) +
    geom_bar(stat = 'identity', width = 0.5) +
    theme_classic()+
    facet_wrap(~cond)
  
  plt_bicjnd<- newdat %>% 
    ggplot(aes(model, bicPse, color = model, fill = model)) +
    geom_bar(stat = 'identity', width = 0.5) +
    theme_classic()+
    facet_wrap(~cond)
  
  pltaicbic_pse<- ggarrange(plt_aicpse, plt_bicjnd, common.legend = TRUE, ncol=2, nrow=1)
  return(pltaicbic_pse)
}

### function runRstanModel
runRstanModelonlocal<- function(Expdata){
  modelname <-  unique(Expdata$model)
  filename <- unique(Expdata$Exp)
  
  print(paste0('current time: ', Sys.time()))
  print(paste0('Start run rstan model ', modelname,' on ', filename))
  RStanModel <- stanmodel123
  rprop <- {}
  x <- {}
  n <- {}
  r <- {}
  condlist <- unique(Expdata$cond)
  nsubjs <- length(unique(Expdata$NSub))
  nstim <- length(unique(Expdata$curDur))
  xmeanfit <- {}
  xsdfit <- {}
  
  #construct the data to be passed to Stan
  m_expdata<- Expdata %>% 
    dplyr::group_by(cond,curDur,NSub) %>%
    dplyr::summarise(N= n(), y = sum(RP))
  m_expdata$rprop <- m_expdata$y / m_expdata$N
  
  xmean <-  matrix(c(rep(0, length(condlist))), length(condlist), nsubjs)  
  xsd <- matrix(c(rep(0, length(condlist))), length(condlist), nsubjs)  
  nstim <- 0
  wf <- c(rep(1, 2*nsubjs))
  
  myinits <- list(
    list(alpha=rep(0, nsubjs), beta=rep(0, nsubjs),
         mua=0, mub=0, sigmaa=1, sigmab=1, sigmac=1),  
    list(alpha=rep(0, nsubjs), beta=rep(0, nsubjs),
         mua=0, mub=0, sigmaa=1, sigmab=1, sigmac=1))  
  
  for(j in 1: length(condlist)){
    for(i in unique(m_expdata$NSub)){
      subdat <- Expdata %>% filter(NSub == i & cond == condlist[j])
      xsd[j, i] <- sd(subdat$curDur)
      nstim  <- length(unique(subdat$curDur))
      dat <- m_expdata %>% filter(NSub == i & cond == condlist[j])
      x <- rbind(x, t(dat$curDur))
      r <- rbind(r, t(dat$y))
      n <- rbind(n, t(dat$N))
      rprop <- rbind(rprop, t(dat$rprop))
      
      if(modelname == 'Bisection Model'){ #simple bisection AM
        xmean[j, i] <- (max(dat$curDur)+min(dat$curDur))/2
      }
      else if (modelname == 'Spacing Model'){ #spacing model AM
        xmean[j, i] <- sum(unique(dat$curDur))/length(unique(dat$curDur))
      }else if(modelname == 'Ensemble Mean'){ #Ensemble Mean AM
        xmean[j, i] <- (dat$curDur  %*% dat$N)/sum(dat$N)
      }else if(modelname == 'Two-stage Ensemble Mean'){#Two stage ensemble AM
        xmean[j, i] <- (dat$curDur  %*% dat$N)/sum(dat$N)
        RStanModel <- stanmodel4
      }else if(modelname == 'Bisection Model GM'){ #simple bisection GM
        xmean[j, i] <- exp(log(max(dat$curDur))/2+ log(min(dat$curDur))/2)
      }else if(modelname == 'Spacing Model GM'){ #spacing model GM
        xmean[j, i] <- exp(sum(log(unique(dat$curDur)))/length(unique(dat$curDur)))
      }else if(modelname == 'Ensemble Mean GM'){ #ensembel mean model GM
        xmean[j, i] <- exp(sum(log(dat$curDur)*dat$N) /sum(dat$N))
      }else if(modelname == 'Two-stage Ensemble Mean GM'){ #Two stage ensemble GM
        xmean[j, i] <- exp(sum(log(dat$curDur)*dat$N) /sum(dat$N))
        RStanModel <- stanmodel4
      }else if(modelname == 'EDA'){ #EDA AM
        xmean[j, i] <- (dat$curDur  %*% dat$N)/sum(dat$N)
        RStanModel <- stanmodel6
        myinits <- list(
          list(alpha=rep(0, nsubjs), beta=matrix(0, 2, nsubjs),
               mua=0, mub=0, sigmaa=1, sigmab=1, sigmac=1),
          list(alpha=rep(0, nsubjs), beta=matrix(0, 2, nsubjs),
               mua=0, mub=0, sigmaa=1, sigmab=1, sigmac=1))
      }else if(modelname == 'EDA GM'){ #EDA GM
        xmean[j, i] <- exp(sum(log(dat$curDur)*dat$N) /sum(dat$N))
        RStanModel <- stanmodel6
        myinits <- list(
          list(alpha=rep(0, nsubjs), beta=matrix(0, 2, nsubjs),
               mua=0, mub=0, sigmaa=1, sigmab=1, sigmac=1),
          list(alpha=rep(0, nsubjs), beta=matrix(0, 2, nsubjs),
               mua=0, mub=0, sigmaa=1, sigmab=1, sigmac=1))
      } 
    }
  }
  
  data<- list("x"=x, "xmean"=xmean, "xsd"=xsd, "n"=n, "r"=r, "nsubjs"= nsubjs,  "nstim"=nstim, "rprop" =rprop)
  
  
  
  parameters <- c("alpha", "beta")  # Parameters to be monitored
  
  # calls Stan with specific options
  samples <- sampling(RStanModel,
                      data=data,
                      init=myinits,
                      #pars=parameters,
                      iter=8000,
                      chains=2,
                      thin=1,
                      seed = 12345,  # Setting seed; Default is random seed
                      control = list(adapt_delta = 0.99, 
                                     max_treedepth = 15)
  )
  
  
  # Constructing MAP-estimates and alpha/beta range
  beta1 <- matrix(0, ncol = nsubjs, nrow =8000)
  beta2 <- matrix(0, ncol = nsubjs, nrow =8000)
  alpha1 <- matrix(0, ncol = nsubjs, nrow =8000)
  alpha2 <- matrix(0, ncol = nsubjs, nrow =8000)
  
  alphaMAP  = matrix(0, 2, nsubjs)
  betaMAP  = matrix(0, 2, nsubjs)
  
  # Extracting the parameters
  #wf = data.frame(rstan::extract(samples)$wf)
  
  
  alpha      = data.frame(rstan::extract(samples)$alpha)
  beta       = data.frame(rstan::extract(samples)$beta)
  
  mua = data.frame(rstan::extract(samples)$mua)
  mub = data.frame(rstan::extract(samples)$mub)
  mua1 = 0
  mua2 = 0
  mub1 = 0 
  mub2 = 0
  
  
  if(modelname %in% c('Bisection Model','Spacing Model','Ensemble Mean', 'Two-stage Ensemble Mean','Bisection Model GM','Spacing Model GM', 'Ensemble Mean GM','Two-stage Ensemble Mean GM')){
    for (i in 1:nsubjs)
    {
      beta1[,i] <- beta[, i]
      beta2[,i] <- beta[, i]
      alpha1[,i] <- alpha[, i]
      alpha2[,i] <- alpha[, i]
    }
    mua1 = mean(mua[,1])
    mua2 = mua1
    mub1 = mean(mub[,1])
    mub2 = mub1
  }else if(modelname %in% c('EDA', 'EDA GM')){
    for (i in 1:nsubjs)
    {
      beta1[,i] <- beta[, 2*i-1]
      beta2[,i] <- beta[, 2*i]
      alpha1[,i] <- alpha[, i]
      alpha2[,i] <- alpha[, i]
    }
    mua1 = mean(mua[,1])
    mua2 = mua1
    mub1 = mean(mub[,1])
    mub2 = mub1
  }
  
  
  Beta <- list('beta1' = beta1,'beta2' = beta2)
  Alpha <- list('alpha1' = alpha1,'alpha2' = alpha2)
  mu <- list('mua1' = mua1,'mua2' = mua2, 'mub1' =mub1, 'mub2'=mub2)
  for (i in 1:nsubjs)
  {
    for(j in 1: length(condlist)){
      alphaMAP[1,i]   <- density(alpha1[,i])$x[which(density(alpha1[,i])$y ==
                                                       max(density(alpha1[,i])$y))]
      betaMAP[1,i]    <- density(beta1[,i])$x[which(density(beta1[,i])$y ==
                                                      max(density(beta1[,i])$y))]
      alphaMAP[2,i]   <- density(alpha2[,i])$x[which(density(alpha2[,i])$y ==
                                                       max(density(alpha2[,i])$y))]
      betaMAP[2,i]    <- density(beta2[,i])$x[which(density(beta2[,i])$y ==
                                                      max(density(beta2[,i])$y))]
    }
  }
  
  
  xmeanhatMAP = matrix(c(rep(0, length(condlist))), length(condlist), nsubjs)  
  if(modelname %in% c('Two-stage Ensemble Mean', 'Two-stage Ensemble Mean GM', 'EDA', 'EDA GM', 'EDA1', 'EDA1 GM')){
    xmeanhat = data.frame(rstan::extract(samples)$xmeanhat)
    for (i in 1:nsubjs){
      for(j in 1: length(condlist)){
        k =(i-1)*2+j
        xmeanhatMAP[j, i]  <- density(xmeanhat[,k])$x[which(density(xmeanhat[,k])$y ==
                                                              max(density(xmeanhat[,k])$y))]
      }
    }
    data$xmeanhat = xmeanhatMAP
  }else{
    data$xmeanhat = data$xmean
  }
  
  # Constructing theta-estimates
  theta = matrix(NA,nsubjs*2, data$nstim)
  for (k in 1:2){
    for (i in 1:nsubjs)
    {
      for (j in 1 : data$nstim){
        # theta[i+(k-1)*nsubjs, j] = inv.logit(alphaMAP[k,i] + betaMAP[k, i] * (data$x[i+(k-1)*nsubjs,j] - data$xmeanhat[k,i]));
        theta[i+(k-1)*nsubjs, j] = inv.logit(alphaMAP[k,i] + betaMAP[k, i] * (data$x[i+(k-1)*nsubjs,j]/data$xmeanhat[k,i] - 1));
        
      }
    }
  }
  
  newtheta = data.frame(theta)
  durs <- length(newtheta[1,])
  colnames(newtheta) <-c(paste0('Prob', 1:durs))
  newtheta$dur1 = c(rep(x[1,1], nsubjs), rep(x[nsubjs+1,1], nsubjs))
  newtheta$dur2 = c(rep(x[1,2], nsubjs), rep(x[nsubjs+1,2], nsubjs))
  newtheta$dur3 = c(rep(x[1,3], nsubjs), rep(x[nsubjs+1,3], nsubjs))
  newtheta$dur4 = c(rep(x[1,4], nsubjs), rep(x[nsubjs+1,4], nsubjs))
  newtheta$dur5 = c(rep(x[1,5], nsubjs), rep(x[nsubjs+1,5], nsubjs))
  newtheta$dur6 = c(rep(x[1,6], nsubjs), rep(x[nsubjs+1,6], nsubjs))
  newtheta$dur7 = c(rep(x[1,7], nsubjs), rep(x[nsubjs+1,7], nsubjs))
  newtheta$model =  modelname
  newtheta$exp =  filename
  newtheta$cond =  c(rep(condlist[1], nsubjs), rep(condlist[2], nsubjs))
  newtheta$cond = factor(newtheta$cond, labels = condlist)
  newtheta$NSub =  c(rep(1:nsubjs, 2))
  gen_theta = NULL
  if(durs ==7 ){
    gen_theta<- newtheta %>% 
      tidyr::unite(DurProb1,  Prob1, dur1)%>%
      tidyr::unite(DurProb2, Prob2, dur2)%>%
      tidyr::unite(DurProb3, Prob3, dur3)%>%
      tidyr::unite(DurProb4, Prob4, dur4)%>%
      tidyr::unite(DurProb5, Prob5, dur5)%>%
      tidyr::unite(DurProb6, Prob6, dur6)%>%
      tidyr::unite(DurProb7, Prob7, dur7)%>%
      melt(id.vars = c('NSub','model','exp','cond'))%>%
      dplyr::rename(
        DurID = variable,
        DurProb = value
      )%>%
      separate(DurProb, c("Prob", "Dur"), sep = "_")
  }else if(length(x[1,]) == 8){
    newtheta$dur8 = c(rep(x[1,8], nsubjs), rep(x[nsubjs+1,8], nsubjs))
    gen_theta<- newtheta %>% 
      tidyr::unite(DurProb1,  Prob1,dur1)%>%
      tidyr::unite(DurProb2, Prob2,dur2)%>%
      tidyr::unite(DurProb3, Prob3,dur3)%>%
      tidyr::unite(DurProb4, Prob4,dur4)%>%
      tidyr::unite(DurProb5, Prob5,dur5)%>%
      tidyr::unite(DurProb6, Prob6,dur6)%>%
      tidyr::unite(DurProb7, Prob7,dur7)%>%
      tidyr::unite(DurProb8, Prob8,dur8)%>%
      melt(id.vars = c('NSub','model','exp','cond'))%>%
      dplyr::rename(
        DurID = variable,
        DurProb = value
      )%>%
      separate(DurProb, c("Prob", "Dur"), sep = "_")
  }
  
  
  #JND/PSE calculation 
  PSE1 <- F3inv(midthreshold,c(1:nsubjs), Alpha$alpha1, Beta$beta1, data$xmeanhat[1,])
  PSE2 <- F3inv(midthreshold,c(1:nsubjs), Alpha$alpha2, Beta$beta2, data$xmeanhat[2,])
  JND1 <- F3inv(threshold,c(1:nsubjs), Alpha$alpha1, Beta$beta1, data$xmeanhat[1,])- F3inv(midthreshold,c(1:nsubjs), Alpha$alpha1, Beta$beta1, data$xmeanhat[1,] )
  JND2 <- F3inv(threshold,c(1:nsubjs), Alpha$alpha2, Beta$beta2, data$xmeanhat[2,])- F3inv(midthreshold,c(1:nsubjs), Alpha$alpha2, Beta$beta2, data$xmeanhat[2,])
  
  PSEmap = matrix(0, 2, nsubjs)
  JNDmap = matrix(0, 2, nsubjs)
  for (k in 1:2){
    PSEmap[k,] <- F3invMAP(midthreshold,c(1:nsubjs), alphaMAP[k,], betaMAP[k,], data$xmeanhat[k,])
    JNDmap[k,] <- F3invMAP(threshold,c(1:nsubjs), alphaMAP[k,], betaMAP[k,], data$xmeanhat[k,])-F3invMAP(midthreshold,c(1:nsubjs), alphaMAP[k,], betaMAP[k,], data$xmeanhat[k,])
  }
  
  resultList <- {}
  for (k in 1:2){
    result <- {}
    result <- data.frame(cbind(as.numeric(data$xmean[k,]), as.numeric(data$xmeanhat[k, ])))
    result <- data.frame(cbind(result, as.numeric(alphaMAP[k,])))
    result <- data.frame(cbind(result, as.numeric(betaMAP[k,])))
    result <- data.frame(cbind(result, as.numeric(PSEmap[k,])))
    result <- data.frame(cbind(result, as.numeric(JNDmap[k,])))
    colnames(result) <- c("xmean","xmeanhat", "alphaMAP", "betaMAP", "PSEmap", "JNDmap")
    result$NSub <- 1:nsubjs
    result$model <- modelname
    result$cond <- condlist[k]
    result$Exp <- filename
    resultList <- rbind(resultList, result)
  }
  
  rlt <- list("mu" =mu,  "theta" = gen_theta, "data" = data, "modelname" = modelname, "exp" = filename, "result" = resultList)
  #returns the model result
  return(rlt)
}



#Parallel compute parameters for each subject each model
runModelcluster <- function(sub_exp_dat){
  cl <- makeCluster(detectCores() - 1)
  clusterEvalQ(cl, {library(dplyr, rstan) })
  clusterExport(cl, c('runRstanModelParl', 'F1', 'F1inv', 'F2', 'F2inv', 'F3','F3inv','F3invMAP'))
  t0 = proc.time()
  resultlist <- clusterMap(cl, runRstanModelParl, sub_exp_dat)
  stopCluster(cl)
  print(proc.time()-t0)
  return(resultlist)
}


