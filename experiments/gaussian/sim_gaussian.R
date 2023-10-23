start <- Sys.time()
source('../utils.R')
args <- commandArgs(trailingOnly=TRUE)

seed <- as.integer(args[1])
resfolder <- args[2]
resname <- paste0(resfolder, '/sim', seed, '.csv')
nthreads <- 2

##############################################################################################
####################################### Settings #############################################
##############################################################################################

N <- 4000                            # Number of participants
interims <- c(0.25, 0.5, 0.75) * N/2 # Timing of interim checkpoints

effect_sizes <- seq(0,1,0.2)               # Effect sizes for the minority group
maj_effect_sizes <- c(-0.1, 0)             # Effect sizes for the majority group
knowntypes <- c("known", "unknown", "est") # known: oracle, unknown: homogeneous baseline, est: CLASH
ncolXs <- c(3,5,10)                        # Number of binary covariates
eff.columns <- c(1,2,3)                    # Number of covariates that drive heterogeneity
                                           # (1: 50% minority group, 2: 25% minority group, 3: 11% minority group)

delta.eff <- 0.1                           # Choice of delta (CLASH hyperparameter)
cv <- 5                                    # k for k-fold cross validation

methods <- c('maxsprt', 'mixsprt',  'bayesian', 'OF')

# Simulation arguments
args_normal <- list(N=N,
                    distr='normal',
                    min.size=0.1,
                    base.mean= 0,
                    base.var = 1, 
                    effect.maj = -0.1,
                    effect.min = 1,
                    alpha=0.05,
                    beta=1,
                    min.stop.time = 0,
                    seed = seed,
                    ncolX = 3,
                    eff.cols = 3)

# Required parameters for stopping tests
method_params <- list(h0=0, # Null hypothesis
                      h1=1, # Alternate hypothesis (only used for SPRT)
                      interims = interims,
                      base.var = args_normal$base.var, # Known variance for OF/Pocock
                      tausq = 1, # for mSPRT
                      maxsprt_bound = 3.701775,    # calibrated from prior simulations
                      bayes_bound = -0.0001858695, # calibrated from prior simulations
                      args=args_normal)

# Hyperparameters for SUBTLE
subtle_init_batch <- 250
subtle_batch_size <- 50
subtle_tausq <- 1

##############################################################################################
####################################### Estimation ###########################################
##############################################################################################
  
i <- 1
result_names <- c('type','ncolX','effect', 'effect_maj', 'eff.cols',
                  'method', 'sim','stop_type', 'stop_time')
results <- data.frame(matrix(NA,nrow=2000, ncol = length(result_names)))
names(results) <- result_names

for(effect in effect_sizes){
  print(paste0('Effect: ', effect))
  for(eff.cols in eff.columns){
    for(ncolX in ncolXs){
      for(maj_effect in maj_effect_sizes){
        
        # Simulate and prepare data
        args_normal$effect.min <- effect
        args_normal$effect.maj <- maj_effect
        args_normal$eff.cols <- eff.cols
        args_normal$ncolX <- ncolX
        data_raw <- simDataXb(args_normal)
        data <- prepData(data_raw)
        
        # Run CLASH, Homogeneous, and Oracle to determine stopping time
        for(knowntype in knowntypes){
          # First compute weights. 
          # Weights are a matrix: each column contains the weights estimated at one interim check
          
          # For Oracle, weights are true group indicator
          if(knowntype=="known"){
            w_final <- matrix(rep(data$group, length(interims)), ncol=length(interims))
            h1 <- max(effect, delta.eff)
          } 
          # For Homogeneous, weights are all 1
          else if(knowntype=="unknown"){ 
            w_final <- matrix(1, nrow=nrow(data), ncol=length(interims))
            h1 <- sum(data$group) / nrow(data) * effect
            h1 <- max(h1, delta.eff)
          }
          # For CLASH, weights are estimated by causal forest with CV
          else if(knowntype=="est"){
            w_final <- matrix(0, nrow=nrow(data), ncol=length(interims))
            for(check in 1:length(interims)){
              w <- estGWeightedCV(data_raw, cv, delta.eff, interims[check],
                                  'cf', args_normal, nthreads = nthreads)
              w <- w[seq(1, length(w), 2)]
              w_final[1:interims[check], check] <- w
            }
            
            pg <- sum(data$group) / nrow(data)
            h1 <- 2*pg*effect / (pg + 1)
            h1 <- (h1==0)*0.1 + h1
          }
          method_params$h1 <- h1
          method_params$args <- args_normal
          
          for(method in methods){
            setting <- cbind(knowntype, ncolX, effect, maj_effect, eff.cols, method, seed)
            stop_time <- stopping_time(data, w_final, method, method_params)
            results[i,] <- c(setting, stop_time)
            i <- i+1
            
            # Run SUBTLE to compare with CLASH
            if(knowntype == 'est' & method=='mixsprt'){
              setting <- cbind('subtle', ncolX, effect, maj_effect, eff.cols, method, seed)
              stop_time <- subtle(data_raw, interims, subtle_init_batch, 
                                  subtle_batch_size, subtle_tausq, args_normal,nthreads = nthreads) 
              results[i,] <- c(setting, stop_time)
              i <- i+1
            }
          }
        }
      }
    }
  }
}

results <- na.omit(results)
results <- results %>%
  mutate(ncolX = as.numeric(ncolX),
         effect=as.numeric(effect),
         effect_maj=as.numeric(effect_maj),
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
