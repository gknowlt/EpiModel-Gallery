
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

nsims <- 1
ncores <- 1
nsteps <- 200


# Vital Dynamics Setup ----------------------------------------------------

ages <- 0:100

# death rate per person, per week
dr <- 0.01 / 52
death_rate_pp_pw <- c(dr / 52, Inf)

# Build out a mortality rate vector
age_spans <- c(100, 1)
dr_vec <- rep(death_rate_pp_pw, times = age_spans)
data.frame(ages, dr_vec)


# Network Model Estimation ------------------------------------------------

# Initialize the network
n_active <- 500
n_inactive <- 100
n <- n_active + n_inactive
nw <- network_initialize(n)

# Set up ages
ageVec <- sample(ages, n, replace = TRUE)
nw <- set_vertex_attribute(nw, "age", ageVec)
nw <- set_vertex_attribute(nw, "active", c(rep(1, n_active), rep(0, n_inactive)))
nw <- set_vertex_attribute(nw, "alive", rep(1, n))

# Define the formation model: edges
formation <- ~edges + absdiff("age") + nodefactor("active", levels = 1)

# Input the appropriate target statistics for each term
mean_degree <- 0.8
edges <- mean_degree * (n_active/2)
avg.abs.age.diff <- 1.5
absdiff <- edges * avg.abs.age.diff
inactive.edges <- 0

target.stats <- c(edges, absdiff, inactive.edges)

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(~offset(edges), 60, dr)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps,
            nwstats.formula = ~edges + absdiff("age") + isolates + degree(0:5) 
            + nodefactor("active", levels = 1))
print(dx)
plot(dx)


# Epidemic model simulation -----------------------------------------------

# Epidemic model parameters
param <- param.net(inf.prob = 0.15,
                   death.rates = dr_vec,
                   arrival.rate = dr)

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
                       infection.FUN = infection.net,
                       resim_nets.FUN = resim_nets,
                       resimulate.network = TRUE,
                       verbose = TRUE)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)

# Plot outcomes
par(mfrow = c(1,2))
plot(sim, main = "State Prevalences", popfrac = TRUE)
plot(sim, main = "State Sizes", sim.lines = TRUE,
     qnts = FALSE, mean.smooth = FALSE)

par(mfrow = c(1, 2))
plot(sim, y = "num", main = "Population Size", qnts = 1, ylim = c(450, 550))
plot(sim, y = "meanAge", main = "Mean Age", qnts = 1, ylim = c(35, 50))

par(mfrow = c(1, 2))
plot(sim, y = "d.flow", mean.smooth = TRUE, qnts = 1, main = "Departures")
plot(sim, y = "a.flow", mean.smooth = TRUE, qnts = 1, main = "Arrivals")

# Examine the data
df <- as.data.frame(sim, out = "mean")
head(df, 25)
