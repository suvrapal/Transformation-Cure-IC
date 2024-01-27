# Transformation-Cure-IC

Code description:

Sim-bc-gen-code-SM.R: Simulation code containing EM-SM algorithm 
 contains functions 
1. LR_int - generates interval
2. data_gen_BC - generates data from Box-Cox model
3. BC_gen_EM_Wei - EM-SM algorithm under assumption of Weibull distribution for baseline survival
4. std - Observed log-likelihood function
5. MC.Sim.log.like - Monte Carlo simulation for parameter estimation 
   using direct max of observed log-likelihood function
6. MC.Sim.BC.EM.IC Monte Carlo simulation for parameter estimation using EM Algorithm

Sim-bc-gen-code-grid.R: Simulation code containing EM-PL algorithm
contains functions
1. LR_int - generates interval
2. data_0_BC - generates data from Box-Cox model when alpha=0
3. data_gen_BC - generates data from Box-Cox model when alpha is in (0, 1]
4. std0 - Observed log-likelihood function when alpha=0
5. stdgen - Observed log-likelihood function when alpha is in (0, 1]
6. BC_0_EM_Wei - EM algorithm for alpha=0
7. BC_gen_EM_Wei - EM-PL algorithm for alpha in (0, 1]
   Note: alpha is treated as fixed in this function because a loop later in the code
   (beginning line 549) initiates a grid containing each permissible value of alpha
