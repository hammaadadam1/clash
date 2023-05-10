usePackage <- function(p){
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
usePackage('dplyr')
usePackage('tidyr')
usePackage('gsDesign')
usePackage('grf')
usePackage('ranger')
usePackage('glmnet')
usePackage('survival')

# Simulate data with specified arguments for Gaussian outcomes
simDataXb <- function(args){
  
  set.seed(args$seed)
  min.size <- args$min.size
  N <- floor(args$N/2)*2
  W <- rep(c(0,1), N / 2)
  ncolX <- args$ncolX
  
  X <- matrix(rbinom(N/2*ncolX, 1, 0.5), nrow=N/2)
  X <- X[rep(seq(1,nrow(X)),each=2),]
  
  if(args$eff.cols > 1){
    tau <- apply(X[,1:args$eff.cols],1,prod)
  } else if(args$eff.cols==1){
    tau <- X[,1]
  } else if(args$eff.cols==0){
    tau <- rep(1, nrow(X))
  }
  G <- (tau > 0)*1
  
  if(args$distr == 'normal'){
    y <- rnorm(N, args$base.mean, sqrt(args$base.var)) + 
      args$effect.maj*W*(1-G) + args$effect.min*W*G
  }
  
  if(args$distr == 'binary'){
    p <- args$base.mean + args$effect.maj*W*(1-G) + args$effect.min*W*G
    y <- rbinom(N, 1, p)
  }
  
  df <- as.data.frame(cbind(rep(seq(1,N/2,1), each=2),seq(1,N,1), W, G, X, y))
  names(df) <- c("time","id", "treatment", "group", paste0("X", seq(1,ncolX)),"outcome")
  
  return(df)
}

# Simulate data with specified arguments for TTE outcomes
simTteXb <- function(args){
  
  set.seed(args$seed)
  
  N <- floor(args$N/2)*2
  W <- rep(c(0,1), N/2)
  lambda0 <- args$lambda0
  lambda1_min <- args$lambda1_min
  lambda1_maj <- args$lambda1_maj
  Maxtime <- args$Maxtime
  ncolX <- args$ncolX
  
  X <- matrix(rbinom(N*ncolX, 1, 0.5), nrow=N)
  
  if(args$eff.cols > 1){
    tau <- apply(X[,1:args$eff.cols],1,prod)
  } else if(args$eff.cols==1){
    tau <- X[,1]
  } else if(args$eff.cols==0){
    tau <- rep(1, nrow(X))
  }
  G <- (tau > 0)*1
  
  Etime <- runif(N, 0, 12)
  
  trt <- W*G
  
  ttime <- rexp(N)
  
  ttime[G==0] <- ttime[G==0]*(W[G==0]/lambda1_maj[1]+(1-W[G==0])/lambda0)
  ttime[G==1] <- ttime[G==1]*(W[G==1]/lambda1_min[1]+(1-W[G==1])/lambda0)
  
  ctime <- rexp(N)/0.014*12
  y <- pmin(ttime, ctime, Maxtime-Etime)
  delta <- 1*(y==ttime)
  
  df <- as.data.frame(cbind(seq(1,N,1), Etime, W, G, X, y, delta))
  names(df) <- c("id", "accrual", "treatment", "group", paste0("X", seq(1,ncolX)), "outcome", "noncensored")
  
  return(df)
}

# Prepare data as differences in treatment and control outcomes
prepData <- function(df){
  df <- df %>% select(time, treatment, group, outcome) %>% 
    pivot_wider(id_cols = c(time,group), names_from = treatment, 
                names_prefix = 'y',values_from = outcome) %>% 
    mutate(y = y1 - y0, 
           y2 = y^2,
           ymean = cumsum(y) / time,
           y0mean = cumsum(y0) / time,
           y1mean = cumsum(y1) / time,
           pmean = (y0mean+ y1mean)/2)
  return(df)
}

# Estimate CLASH weights with cross validation
estGWeightedCV <- function(df, cv, cutoff, est_time, method, sim_args, nthreads){
  ncolX <- sim_args$ncolX
  
  # Only consider data observed thus far
  df <- df %>% filter(time <= est_time)
  
  # Split into folds
  folds <- rep(sample(cut(seq(1,nrow(df)/2),breaks=cv,labels=FALSE)),each=2)
  w <- rep(NA, nrow(df))
  
  df$pred <- 1
  df$w <- 1
  
  for(f in 1:cv){
    train_ids <- which(folds!=f)
    pred_ids <- which(folds==f)
    X <- df[train_ids,] %>% select(paste0("X", seq(1,ncolX))) %>% as.matrix()
    # Train causal forest on train fold
    cf <- causal_forest(X=X, 
                         Y=df$outcome[train_ids],
                         W=df$treatment[train_ids],
                         W.hat = rep(0.5, nrow(X)),
                         num.trees = 2000, 
                         num.threads = nthreads)
    
    Xn <- df[pred_ids,] %>% select(paste0("X", seq(1,ncolX))) %>% as.matrix()
    
    # Get predicted treatment effects from CF
    pred <- predict(cf, Xn, estimate.variance = TRUE, num.threads = nthreads)
    
    # Compute CLASH weights
    w[pred_ids] <- 1 - pnorm((cutoff-pred$predictions) / sqrt(pred$variance.estimates))
  }
  return(w)
}

# Wrapper function to run different stopping tests with Gaussian data
stopping_time <- function(data, weights, method, method_params){
  
  h0 <- method_params$h0
  h1 <- method_params$h1
  sigmasq <- method_params$base.var
  sim_args <- method_params$args
  interims <- method_params$interims
  
  if(method == 'sprt'){
    st <- wsprt(data, weights, interims, h0, h1, sigmasq, sim_args)
  } 
  else if(method=='maxsprt'){
    st <- wmaxsprt(data, weights, interims, h0, sigmasq, method_params$maxsprt_bound, sim_args)
  }
  else if(method %in% c('Pocock', 'OF')){
    st <- wfrequentistStop(data, weights, interims, sigmasq, method, sim_args)
  }
  else if(method == 'mixsprt'){
    st <- wmixsprt(data, weights, interims, h0, method_params$tausq, sigmasq, sim_args)
  } 
  else if(method =='bayesian'){
    st <- bayesian(data, weights, interims, h0, method_params$tausq, sigmasq, method_params$bayes_bound, sim_args)
  }
  return(st)
}

# Wald's SPRT with Gaussian outcomes
wsprt <- function(df, weights, interims, h0, h1, sigmasq, sim_args){
  if(sim_args$distr=='normal'){
    sprt_df <- df %>% 
      mutate(w = weights, 
             llh0 = dnorm(y, mean = h0, sd = sqrt(2*sigmasq),log=TRUE), 
             llh1 = dnorm(y, mean = h1, sd = sqrt(2*sigmasq),log=TRUE), 
             llr = w*(llh1 - llh0), 
             llr = cumsum(llr))
  }
  llr_bound <- 2.995732 #waldBoundary(type1 = sim_args$alpha, type2 = 1-sim_args$beta)
  
  max_time <- max(sprt_df$time)
  sprt_calc <- sprt_df %>% filter(time %in% interims) %>%
    mutate(decision = 1*((llr <= llr_bound))) %>%
    filter(decision == 0)
  
  stop_time <- ifelse(nrow(sprt_calc)==0, max(sprt_df$time), sprt_calc$time[1])
  stop_type <- ifelse(nrow(sprt_calc)==0, "None", 
                      ifelse(h1 > h0, "Increase", "Decrease"))
  return(cbind(stop_type, stop_time))
}

# maxSPRT (with weights)
wmaxsprt <- function(df, weights, interims, h0, sigmasq, llr_bound, sim_args){
  
  if(sim_args$distr=='normal'){
    sprt_df <- df %>%
      mutate(w = weights, 
             wy = w*y, 
             wymean = cumsum(wy) / cumsum(w), 
             llr = (wymean/(4*sigmasq))*(2*cumsum(wy) - 
                                           wymean*cumsum(w)))
  }
  
  max_time <- max(sprt_df$time)
  sprt_calc <- sprt_df %>% filter(time %in% interims) %>%
    mutate(decision = 1*(llr <= llr_bound)) %>%
    filter(decision == 0, ymean > 0)
  
  stop_time <- ifelse(nrow(sprt_calc)==0, max(sprt_df$time), sprt_calc$time[1])
  stop_ymean <- ifelse(nrow(sprt_calc)==0, h0, 
                       sprt_calc$ymean[1])
  stop_type <- ifelse(stop_ymean==h0, "None", 
                      ifelse(stop_ymean > h0, "Increase", "Decrease"))
  return(cbind(stop_type, stop_time))
}

# Frequentist tests (with weights): supports Pocock and OBF
wfrequentistStop <- function(df, weights, interims, sigmasq, method, sim_args){
  if(sim_args$distr=='normal'){
    sprt_df <- df %>%  
      mutate(w = weights,# / sum(weights), 
             w2 = w^2,
             wy = w*y, 
             wysum = cumsum(wy),
             zscore = cumsum(wy)/ sqrt(2*sigmasq*cumsum(w2)))
    
    max_time <- max(sprt_df$time)
    
    sprt_df <- sprt_df %>% filter(time %in% c(interims, max_time))
    k = length(interims) + 1
    
    sprt_df$zbound <- gsDesign(k = k, test.type = 1, 
                               alpha = sim_args$alpha,
                               sfu = method)$upper$bound
    
    sprt_calc <- sprt_df %>%
                  mutate(decision = 1*(zscore <= zbound)) %>%
                    filter(decision == 0, time %in% interims)
    
    stop_time <- ifelse(nrow(sprt_calc)==0, max_time, sprt_calc$time[1])
    stop_type <- ifelse(nrow(sprt_calc)==0, "None", "Increase")
    # return(sprt_calc)
    return(cbind(stop_type, stop_time))
  }
}

# mSPRT (with weights), based on Johari et al (2017)
wmixsprt <- function(df, weights, interims, h0, tausq, sigmasq, sim_args){
  if(sim_args$distr == 'normal'){
    sprt_df <- df %>%
                mutate(w = weights, 
                       wy = w*y, 
                       wsum = cumsum(w),
                       wymean = cumsum(wy) / cumsum(w))
    
    sprt_df$mult1 <- 2*sigmasq / (2*sigmasq + sprt_df$wsum*tausq)
    sprt_df$exp1 <- sprt_df$wsum^2 * tausq * (sprt_df$wymean - h0)^2 / (4*sigmasq*(2*sigmasq + sprt_df$wsum*tausq))
    
    sprt_df$llr <- 0.5*log(sprt_df$mult1) + sprt_df$exp1
  }
  llr_bound <- -log(sim_args$alpha)
  max_time <- max(sprt_df$time)
  
  sprt_calc <- sprt_df %>% filter(time %in% interims) %>%
    mutate(decision = 1*((llr <= llr_bound))) %>%
    filter(decision == 0)
  
  stop_time <- ifelse(nrow(sprt_calc)==0, max_time, sprt_calc$time[1])
  stop_type <- ifelse(nrow(sprt_calc)==0, "None", 
                      ifelse(sprt_calc$wymean[1] > 0, "Increase", "Decrease"))
  return(cbind(stop_type, stop_time))
}

# Bayesian estimation-based test for early stopping
bayesian <- function(df, weights, interims, h0, tausq, sigmasq, loss_bound, sim_args){
  if(sim_args$distr == 'normal'){
    sprt_df <- df %>% 
                mutate(w = weights, 
                       wy = w*y, 
                       wsum = cumsum(w),
                       wymean = cumsum(wy) / cumsum(w))
    
    sprt_df$post_mean <- (h0*2*sigmasq + sprt_df$wsum*tausq*sprt_df$wymean) / (sprt_df$wsum*tausq + 2*sigmasq)
    sprt_df$post_var <- (2*sigmasq*tausq)/(sprt_df$wsum*tausq + 2*sigmasq)
    
    sprt_df$loss1 <- bayesloss(sprt_df$post_mean, sqrt(sprt_df$post_var))
    sprt_df$loss2 <- bayesloss(-sprt_df$post_mean, sqrt(sprt_df$post_var))
    
    sprt_df$loss <- sprt_df$loss1
  }
  
  sprt_calc <- sprt_df %>% filter(time %in% interims) %>%
    mutate(decision = 1*(loss <= loss_bound)) %>%
      filter(decision == 0)
  
  stop_time <- ifelse(nrow(sprt_calc)==0, max(sprt_df$time), sprt_calc$time[1])
  stop_postmean <- ifelse(nrow(sprt_calc)==0, h0, 
                          sprt_calc$post_mean[1])
  stop_type <- ifelse(stop_postmean==h0, "None", 
                      ifelse(stop_postmean > h0, "Increase", "Decrease"))
  
  return(cbind(stop_type, stop_time))
}

bayesloss <- function(mu, sigma){
  -(sigma/sqrt(2*pi))*exp(-mu^2/(2*sigma^2)) + mu*pnorm(-mu/sigma)
}

# Wrapper for standard group sequential testing 
gsTTE <- function(df, interims, method, metric, sim_args){
  
  if(interims[length(interims)]==sim_args$Maxtime){gsk = length(interims)} else{gsk = length(interims)+1}
  
  design <- gsDesign(k=gsk, test.type = 1, 
                     alpha=sim_args$alpha, sfu=method, 
                     timing = interims / sim_args$Maxtime)
  
  p_bounds <- 1-pnorm(design$upper$bound)
  
  for(check_time in interims){
    df_check <- df[check_time > df$accrual,]
    df_check$outcome <- pmin(df_check$outcome, check_time-df_check$accrual)
    df_check$noncensored <- df_check$noncensored * (df_check$outcome < check_time-df_check$accrual)
    
    if(metric == 'RMST'){
      mod <- rmst2(df_check$outcome, df_check$noncensored, df_check$treatment)$unadjusted[1,]
      interim_est <- mod[1]
      interim_p <- mod[4]
    } 
    else if(metric=='PH'){
      mod <- coef(summary(coxph(Surv(outcome, noncensored)~treatment, df_check)))[1,]
      interim_est <- mod[1]
      interim_p <- mod[5]
    }
    if(interim_p < p_bounds[which(interims==check_time)]){
      if(interim_est < 0){stop_type <- "Decrease"} else if(interim_est > 0) {stop_type <- "Increase"}
      return(c(stop_type, check_time))
    }
  }
  return(c("None", sim_args$Maxtime))
}

# Wrapper for group sequential testing with weights
wgsTTE <- function(df, interims, wmethod, cutoff, nfolds, method, metric, sim_args){
  
  if(interims[length(interims)]==sim_args$Maxtime){gsk = length(interims)} else{gsk = length(interims)+1}
  
  design <- gsDesign(k=gsk, test.type = 1, 
                     alpha=sim_args$alpha, sfu=method, 
                     timing = interims / sim_args$Maxtime)
  
  p_bounds <- 1-pnorm(design$upper$bound)
  
  if(wmethod=='known'){
    return(gsTTE(df[df$group==1,], interims, method, metric, sim_args))
  }
  
  for(check_time in interims){
    df_check <- df[check_time > df$accrual,]
    df_check$outcome <- pmin(df_check$outcome, check_time-df_check$accrual)
    df_check$noncensored <- df_check$noncensored * (df_check$outcome < check_time-df_check$accrual)
    
    if(wmethod=='est'){
      # determine weights
      df_check$w <- pmax(estWTTE(df, cutoff, check_time, nfolds, sim_args, NULL),0.01)
    } 
    else if(wmethod == 'unknown'){
      df_check$w <- 1
    }
    
    
    if(metric=='PH'){
      mod <- coef(summary(coxph(Surv(outcome, noncensored)~treatment, df_check, 
                                weights = df_check$w, robust = FALSE)))[1,]
      interim_est <- mod[1]
      interim_p <- mod[length(mod)]
    }
    if(interim_p < p_bounds[which(interims==check_time)]){
      if(interim_est < 0){stop_type <- "Decrease"} else if(interim_est > 0) {stop_type <- "Increase"}
      return(c(stop_type, check_time))
    }
  }
  return(c("None", sim_args$Maxtime))
}

# Calculate weights for TTE data
estWTTE <- function(df, cutoff, est_time, nfolds, sim_args, nthreads){
  
  ncolX <- sim_args$ncolX
  
  df_est <- df[est_time > df$accrual,]
  df_est$outcome <- pmin(df_est$outcome, est_time-df_est$accrual)
  df_est$noncensored <- df_est$noncensored * (df_est$outcome < est_time-df_est$accrual)
  
  folds <- sample(cut(seq(1,nrow(df_est)),breaks=nfolds,labels=FALSE))
  
  w <- rep(NA, nrow(df_est))
  
  for(f in 1:nfolds){
    train_ids <- which(folds!=f)
    pred_ids <- which(folds==f)
    X <- df_est[train_ids,] %>% select(paste0("X", seq(1,ncolX))) %>% as.matrix()
    csf <- causal_survival_forest(X=X, 
                                  Y=df_est$outcome[train_ids],
                                  W=df_est$treatment[train_ids],
                                  D=df_est$noncensored[train_ids],
                                  W.hat = rep(0.5, nrow(X)), 
                                  horizon = est_time-6, # heuristic to stop instability
                                  target = "RMST", 
                                  num.threads = nthreads)
    
    Xn <- df_est[pred_ids,] %>% select(paste0("X", seq(1,ncolX))) %>% as.matrix()
    pred <- predict(csf, Xn, estimate.variance = TRUE, num.threads = nthreads)
    w[pred_ids] <- 1 - pnorm((cutoff-pred$predictions) / sqrt(pred$variance.estimates))
  }
  return(w)
}

# Implementation of SUBTLE (Yu et al 2021)
subtle <- function(df, interims, ibatch, tbatch, tausq, sim_args, nthreads){
  
  ncolX <- sim_args$ncolX
  rf_form <- as.formula(paste0('outcome ~ ',
                               paste(paste0("X", seq(1,ncolX)), collapse=' + ')))
  tau <- sqrt(tausq)
  dmeans <- c()
  dsigmas <- c()
  
  check_times <- seq(ibatch, max(interims) - tbatch, tbatch)
  
  for(check_time in check_times){
    stop_time <- check_time + tbatch
    data_prev <- df %>% filter(time <= check_time)
    data_batch <- df %>% filter(time > check_time, 
                                time <= check_time + tbatch)
    
    rf_c <- ranger(rf_form, data_prev %>% filter(treatment == 0), classification=FALSE, num.threads = nthreads)
    rf_t <- ranger(rf_form, data_prev %>% filter(treatment == 1), classification=FALSE, num.threads = nthreads)
    
    # previous batch prediction
    data_prev$pred_c  <- predict(rf_c, data_prev, num.threads = nthreads)$prediction
    data_prev$pred_t  <- predict(rf_t, data_prev, num.threads = nthreads)$prediction
    
    # new batch prediction
    data_batch$pred_c  <- predict(rf_c, data_batch, num.threads = nthreads)$prediction
    data_batch$pred_t  <- predict(rf_t, data_batch, num.threads = nthreads)$prediction
    
    # propensity score
    p = sum(data_prev$treatment) / length(data_prev$treatment)
    
    # d in new batch
    data_batch <- data_batch %>% 
      mutate(mu = pred_c, 
             theta =  pred_t - pred_c, 
             opt = 1*(theta>0),
             d1 = outcome*(treatment==opt) / p - ((treatment==opt)/p - 1)*(mu+theta*opt), 
             d0 = outcome*(1-treatment)/(1-p) - ((1-treatment)/(1-p) - 1)*mu, 
             d = d1-d0)
    # d in old batch
    data_prev <- data_prev %>% 
      mutate(mu = pred_c, 
             theta =  pred_t - pred_c, 
             opt = 1*(theta>0),
             d1 = outcome*(treatment==opt) / p - ((treatment==opt)/p - 1)*(mu+theta*opt), 
             d0 = outcome*(1-treatment)/(1-p) - ((1-treatment)/(1-p) - 1)*mu, 
             d = d1-d0)
    
    dmean <- mean(data_batch$d)
    dsigma <- sqrt(var(data_prev$d) / nrow(data_batch))
    dsigma <- max(dsigma, 0.01)
    
    dmeans <- c(dmeans, dmean)
    dsigmas <- c(dsigmas, dsigma)
    
    k <- length(dmeans)
    R_k <- sum(dmeans / dsigmas) / sqrt(k)
    
    Fmean = sqrt(k)*sum(R_k/dsigmas) / (k + sum(1/dsigmas)^2)
    Fsd = k*tausq / (k + tausq*sum(1/dsigmas)^2)
    
    lambda_k <- 2*sqrt(k/(k + tausq*sum(1/dsigmas)^2))*
                exp((0.5*tausq*sum(R_k/dsigmas)^2 ) / (k + tausq*sum(1/dsigmas)^2))*
                (1-pnorm(0, mean = Fmean, sd = Fsd))
    
    if(lambda_k > 1/sim_args$alpha){
      if(stop_time %in% interims){
        return(cbind("Increase", stop_time))
      }
    }
  }
  return(cbind("None", max(df$time)))
}

