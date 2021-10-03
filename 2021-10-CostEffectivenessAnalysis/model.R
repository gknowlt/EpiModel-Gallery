
##
## Simple Cost-effectiveness Model (Simple SI model with cost and utility tracking)
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Gregory Knowlton (University of Minnesota)
## Date: October 2021
##

# Load EpiModel
suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

if (interactive()) {
  nsims <- 5
  ncores <- 5
  nsteps <- 256
} else {
  nsims <- 1
  ncores <- 1
  nsteps <- 104
}


# Vital Dynamics Setup ----------------------------------------------------

# Yearly mortality rates by age: (1-84: 0.02, 85-99: 1.00, 100+: 2.00)
ages <- 0:100
dr <- 0.02
death_rate_pp_pw <- c(dr / 52, 50*dr / 52, 100*dr / 52)

# Build out a mortality rate vector
age_spans <- c(85, 15, 1)
dr_vec <- rep(death_rate_pp_pw, times = age_spans)
data.frame(ages, dr_vec)

# Network Model Estimation ------------------------------------------------

# Initialize the network
n <- 500
init_ages <- 16:85
ageVec <- sample(init_ages, n, replace = TRUE)

# Throughout simulation, individuals over age 65 are sexually inactive
activeVec <- ifelse(ageVec <= 65, 1, 0)

n_active_s <- length(which(ageVec <= 65))
n_inactive_s <- length(which(ageVec > 65))

nw <- network_initialize(n)

# Set up ages
nw <- set_vertex_attribute(nw, "age", ageVec)
nw <- set_vertex_attribute(nw, "active_s", activeVec)

# Define the formation model: edges
formation <- ~edges + absdiff("age") + nodefactor("active_s", levels = 1)

# Input the appropriate target statistics for each term
mean_degree <- 0.8
edges <- mean_degree * (n_active_s/2)
avg.abs.age.diff <- 1.5
absdiff <- edges * avg.abs.age.diff
inactive_s.edges <- 0

target.stats <- c(edges, absdiff, inactive_s.edges)

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(~offset(edges), 60, mean(dr_vec))
coef.diss


# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps,
            nwstats.formula = ~edges + absdiff("age") + isolates + degree(0:5)
            + nodefactor("active_s", levels = 1))
print(dx)
plot(dx)


# Epidemic model simulation -----------------------------------------------

# Epidemic model parameters
param <- param.net(inf.prob = 0.15,
                   death.rates = dr_vec,
                   end.horizon = 52,
                   arrival.rate = dr / 52,
                   # Weekly costs by health status
                   sus_cost = 150,
                   inf_cost = 300,
                   # QALYs by health status
                   sus_qaly = 1.00,
                   inf_qaly = 0.75,
                   # 
                   age_decrement = -0.003,
                   disc_rate = 0.03)

# Initial conditions
init <- init.net(i.num = 50)

# Read in the module functions

if (interactive()) {
  source("2021-10-CostEffectivenessAnalysis/module-fx.R", echo = TRUE)
} else {
  source("module-fx.R")
}

# Control settings
control <- control.net(type = NULL,
                       nsims = nsims,
                       ncores = ncores,
                       nsteps = nsteps,
                       aging.FUN = aging,
                       departures.FUN = dfunc,
                       arrivals.FUN = afunc,
                       infection.FUN = ifunc,
                       cea.FUN = costeffect,
                       resim_nets.FUN = resimfunc,
                       resimulate.network = TRUE,
                       verbose = TRUE)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)

# Plot outcomes
par(mfrow = c(1, 3))
plot(sim, y = "d.flow", mean.smooth = TRUE, qnts = 1, main = "Departures")
plot(sim, y = "a.flow", mean.smooth = TRUE, qnts = 1, main = "Arrivals")
plot(sim, y = "si.flow", mean.smooth = TRUE, qnts = 1, main = "Infections")

par(mfrow = c(2, 2))
plot(sim, y = "cost", mean.smooth = TRUE, qnts = 1, main = "Costs (undiscounted)")
plot(sim, y = "qaly", mean.smooth = TRUE, qnts = 1, main = "QALYs (undiscounted)")
plot(sim, y = "cost.disc", mean.smooth = TRUE, qnts = 1, main = "Costs (discounted)")
plot(sim, y = "qaly.disc", mean.smooth = TRUE, qnts = 1, main = "QALYs (discounted)")

