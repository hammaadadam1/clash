## Causal Latent Analysis for Stopping Heterogeneously (CLASH)

This repository contains replication code for the experiments described in *Should I Stop or Should I Go: Early Stopping with Heterogeneous Populations*. All experiments were run using R version 4.2.2.

### Gaussian Outcomes

The files contained in the `experiments/gaussian` folder run the simulation experiments with Gaussian outcomes (as in Figures 2-3). To run one replication of the experiment (i.e., all considered simulation settings with a single random seed), navigate to the `experiments/gaussian` folder and run the following line of code from the command line (replacing <seed> with the desired random seed):
  
``` Rscript --vanilla sim_gaussian.R <seed> ./results```
  
The results will be stored in `experiments/gaussian/results`. Running this line of code with seeds 1-1000 will generate the full set of simulation results. Note that we used a computing cluster with the Slurm workload manager to run all 1,000 replications in parallel. To do the same, navigate to the `experiments/gaussian` folder and run the following lines of code from the command line:
  
```
  mkdir output
  sbatch tte.slurm
  ```
 
Finally, `experiments/gaussian/plot_gaussian.R` uses the generated results to create the plots from the Paper. Running this file interactively (through RStudio) will produce Figures 2, 3, S7, S8, and S9. 
  
### TTE Outcomes

The files contained in the `experiments/tte` folder run the simulation experiments with time-to-event (TTE) outcomes. To run one replication of the experiment (i.e., all considered simulation settings with a single random seed), navigate to the `experiments/tte` folder and run the following line of code from the command line (replacing <seed> with the desired random seed):
  
``` Rscript --vanilla sim_tte.R <seed> ./results```
  
The results will be stored in `experiments/tte/results`. Running this line of code with seeds 1-1000 will generate the full set of simulation results. Note that we used a computing cluster with the Slurm workload manager to run all 1,000 replications in parallel. To do the same, navigate to the `experiments/tte` folder and run the following lines of code from the command line:
  
```
  mkdir output
  sbatch tte.slurm
  ```
 
Finally, `experiments/gaussian/plot_tte.R` uses the generated results to create the plots from the Paper. Running this file interactively (through RStudio) will produce Figures S10, S11, and S12.
  
 ### Real-world Application
  
 The data from our real-world application is proprietary, and thus cannot be provided as part of this replication package.
