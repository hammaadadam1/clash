start <- Sys.time()
source('../utils.R')
args <- commandArgs(trailingOnly=TRUE)

seed <- as.integer(args[1])
resfolder <- args[2]
resname <- paste0(resfolder, '/sim', seed, '.csv')

##############################################################################################
####################################### Settings #############################################
##############################################################################################

N <- 2000   # Number of participants
maxT <- 30  # Follow up time
interims <- c(0.4, 0.6, 0.8) * maxT  # Timing of interim checkpoints

lambda1s <- c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1) # Survival function lambda1 for minority group
maj_lambda1s <- c(0.1, 0.12)                     # Survival function lambda1 for majority group
knowntypes <- c("known", "unknown", "est")       # known: oracle, unknown: homogeneous baseline, est: CLASH
ncolXs <- c(3,5,10)        # Number of binary covariates
eff.columns <- c(0,1,2)    # Number of covariates that drive heterogeneity (0: 100% group, 1: 50% minority group, 2: 25% minority group)

delta.eff <- 0.05          # Choice of delta (CLASH hyperparameter)
cv <- 5                    # k for k-fold cross validation

methods <- c('OF', 'Pocock')

# Simulation arguments
tte_args <- list(N=N,
                 Maxtime = maxT,
                 lambda0 = 0.1,
                 lambda1_min = 0.1,
                 lambda1_maj = 0.12,
                 alpha=0.05,
                 beta=1,
                 seed = seed,
                 ncolX = 5, 
                 eff.cols = 2)

##############################################################################################
####################################### Estimation ###########################################
##############################################################################################

i <- 1

result_names <- c('type', 'ncolX','lambda1', 'lambda1_maj', 'eff.cols', 'HR',
                  'method', 'sim', 'stop_type', 'stop_time')

results <- data.frame(matrix(NA,nrow=20000, ncol = length(result_names)))
names(results) <- result_names

for(lambda1 in lambda1s){
  print(paste0('Lambda1: ', lambda1))
  for(eff.cols in eff.columns){
    for(ncolX in ncolXs){
      for(maj_lambda1 in maj_lambda1s){
        
        # Simulate data
        tte_args$lambda1_min <- lambda1
        tte_args$lambda1_maj <- maj_lambda1
        tte_args$eff.cols <- eff.cols
        tte_args$ncolX <- ncolX
        data <- simTteXb(tte_args)
        
        # Determine true hazard ratio (HR)
        hr <- NA
        if(eff.cols==0 & tte_args$lambda1_min != tte_args$lambda0){
          hr <- exp(coxph(Surv(outcome, noncensored)~treatment, data)$coefficients[1])
        } 
        else if (eff.cols==0 & tte_args$lambda1_min == tte_args$lambda0){
          hr <- 1
        } 
        
        # Run CLASH, Homogeneous, and Oracle to determine stopping time
        for(knowntype in knowntypes){
          for(method in methods){
            setting <- cbind(knowntype, ncolX, lambda1, maj_lambda1, eff.cols, hr, method, seed)
            
            if(knowntype %in% c('unknown', 'known')){
              stop_time <- wgsTTE(data, interims, knowntype, delta.eff, 5, method, 'PH', tte_args)
              results[i,] <- c(setting, stop_time)
              i <- i+1
            }
            else if (knowntype == 'est'){
              stop_time <- wgsTTE(data, interims, knowntype, delta.eff, 5, method, 'PH', tte_args)
              results[i,] <- c(setting, stop_time)
              i <- i+1
            }
          }
        }
      }
    }
  }
}
results <- results %>%
  mutate(ncolX = as.numeric(ncolX),
         lambda1=as.numeric(lambda1),
         lambda1_maj = as.numeric(lambda1_maj),
         HR=as.numeric(HR),
         eff.cols = as.numeric(eff.cols),
         sim=as.numeric(sim),
         stop_time = as.numeric(stop_time)) %>% 
          filter(!is.na(stop_time))
if(!dir.exists(resfolder)){
  dir.create(resfolder)
}
write.csv(results, resname, row.names = FALSE)
end <- Sys.time()
print(start-end)


