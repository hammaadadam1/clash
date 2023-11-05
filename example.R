# This notebook provides a simple example on how to use CLASH for early stopping

## We first import a few packages we'll need and set a random seed

library(grf)
library(gsDesign)
set.seed(100)

#############################################################################
############################ Generate Data ##################################
#############################################################################

## The experiment recruits 4,000 participants (two at each time point)
N <- 4000
participant_id <- 1:N
time_id <- rep(seq(1,N/2,1), each=2)

## Assign treatment and control groups (50-50 split)
D <- rep(c(0,1), N / 2) # assign treatment and control groups

## Generate covariates and group membership
ncolX <- 5 # Number of covariates
X <- matrix(rbinom(ncolX*N, 1, 0.5), nrow=N) # Simulate binary covariates
G <- X[,1]*X[,2]*X[,3] # assign minority group membership based on three covariates

## Generate outcomes
effect_G <- 1            # Minority group is harmed with effect size 1
effect_notG <- -0.1      # Majority group is benefited with effect size 0.1
y <- rnorm(N) + effect_notG*D*(1-G) + effect_G*D*G      # Simulate outcomes

# ## Gather data into data frame
# df <- as.data.frame(cbind(time_id, participant_id, D, G, X, y))
# names(df) <- c("time","id", "treatment", "group", paste0("X", 1:ncolX),"outcome")

#############################################################################
############## Conduct interim checkpoint  ##################################
#############################################################################

## Here, we'll conduct one interim checkpoint halfway through the trial

X_half <- X[1:(N/2),]
D_half <- D[1:(N/2)]
y_half <- y[1:(N/2)]

################## CLASH Stage 1: Estimate weights ###########################

# At this checkpoint, we first use the interim data to compute the CLASH weights 

estWeightsCLASH <- function(y_interim, D_interim, X_interim, cv, delta_cutoff){
  ## This function computes the CLASH weights at an interim checkpoint
  ## It has five arguments: 
  ##  1. y_interim: outcomes observed at interim checkpoint
  ##  2. D_interim: treatment assignments at interim checkpoint
  ##  3. X_interim: covariates collected at interim checkpoint
  ##  3. cv: the number of cross-validation folds to use during estimation
  ##. 4. delta_cutoff: hyperparameter which should be set to minimum effect size of interest
  
  ### Split data into cross-validation folds and initialize weights vector
  folds <- sample(cut(1:nrow(X_interim), breaks=cv, labels=FALSE))
  w <- rep(NA, nrow(X_interim))
  
  ### Conduct estimation in each fold separately
  for(f in 1:cv){
    
    # Get train data for this fold
    train_ids <- which(folds!=f)
    X_fold <- X_interim[train_ids,]
    y_fold <- y_interim[train_ids]
    D_fold <- D_interim[train_ids]
    
    # Fit causal forest on training data
    cf <- causal_forest(X=X_fold, 
                        Y=y_fold,
                        W=D_fold)
    
    # Predict CATE on held out data
    pred_ids <- which(folds==f)
    X_test <- X_interim[pred_ids, ]
    pred <- predict(cf, X_test, estimate.variance = TRUE)
    
    # Compute and store CLASH weights
    w_fold <- 1 - pnorm((delta_cutoff-pred$predictions) / sqrt(pred$variance.estimates))
    w[pred_ids] <- w_fold
  }
  return(w)
}

## We estimate the CLASH weights using 5-fold CV and delta=0.1 (as in the paper)
w_half <- estWeightsCLASH(y_half, D_half, X_half, cv=5, delta_cutoff = 0.1)

################## CLASH Stage 2: Early Stopping #############################

## We now use the estimate weights for early stopping
## In this example, we'll use a two-sample z-test with the O'Brien-Fleming adjustment
## We assume variance of y is known and equal to 1, and conduct a one-sided test for harm

sigmasq <- 1

## Compute unweighted z-statistic (i.e., the homogeneous approach)
treatment_mean <- mean(y_half[D_half==1])
control_mean <- mean(y_half[D_half==0])
treatment_n <- sum(D_half)
control_n <- sum(1-D_half)

zscore <- (treatment_mean - control_mean) / sqrt(sigmasq/treatment_n + sigmasq / control_n)

## Compute weighted z-statistic (i.e., CLASH)
wy_half <- w_half*y_half
treatment_wmean <- sum(wy_half[D_half==1]) / sum(w_half[D_half==1])
control_wmean <- sum(wy_half[D_half==0]) / sum(w_half[D_half==0])
treatment_wn <- sum(w_half*D_half)
control_wn <- sum(w_half*(1-D_half))

wzscore <- (treatment_wmean - control_wmean) / sqrt(sigmasq/treatment_wn + sigmasq / control_wn)

## Get O'Brien Fleming z-bound with one interim analysis

OF_design <- gsDesign(k = 2,         # one interim analysis + one final analysis
                      test.type = 1, # one-sided test for harm
                      alpha = 0.05,  # Type I error control: 5%
                      sfu = 'OF')    # O'Brien Fleming

OF_bound <- OF_design$upper$bound[1] # z-bound at first analysis

## Compare CLASH z-score to the bound to make stopping decision

Homogeneous_stop <- (zscore > OF_bound)
CLASH_stop <- (wzscore > OF_bound)

if(CLASH_stop){
  print("CLASH suggests stopping the experiment")
} else{
  print("CLASH suggests continuing the experiment")
} 

if(Homogeneous_stop){
  print("The homogeneous approach suggests stopping the experiment")
} else{
  print("The homogeneous approach suggests continuing the experiment")
}
