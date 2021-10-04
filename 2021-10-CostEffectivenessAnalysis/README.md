# Simple Cost-effectiveness Model (Simple SI model with cost and utility tracking)

## Description
This example shows how to model a relatively simple SI epidemic over a dynamic network, but with vital dynamic processes for aging, births, and deaths. The modules implement a different approach to mortality, in which the stochastic death process is a function of both an age-specific mortality rate and disease-induced mortality. The network model uses an `absdiff` term to specify a parametric form of age mixing, and that age attribute that is initialized on the network is then pulled into the epidemic modules.

Within the module that simulates the sexual network(s), you must manually remove the partnerships that involve egos who are still alive but no longer sexually active. Ordinarily, egos are automatically deleted by the departure module when they reach age 65, but because those individuals are being retained for the purposes of tracking quality of life and health care spending, these inactive individuals must be removed from the edge list. In this module and several others, an additional chunk of code is added to instruct EpiModel to skip over this module when the time-step has exceeded the end of the main analytic time window/the beginning of the end horizon. 

Within the departure module, we add more code to differentiate between individuals that become sexually inactive and individuals who die. We add in a new indicator for egos who have reached age 65 at this time-step (named `ids65`) and set the `active` attribute to be equal to `0`. Unlike model code in other versions of EpiModelHIV, those in `ids65` are not included in the vector of egos to delete from the simulation, `idsDepAll`. 

Deeper in the `departure_msm` module, another change is made to facilitate the simulation of the end horizon. When the time-step of the simulation reaches the beginning of the end-horizon, all ego's `active` attribute is set to `0`.

During the simulation of the end horizon, certain modules are disabled so that transmission of HIV/STIs cease but other processes related to testing, treatment, disease progression, and outcome tracking continue.

If QALYs are age-adjusted additively, then we can easily calculate QALYs externally. The age decrement would be equal to the following at each time-step: (age.mean) * (num.alive) * age.adj
If QALYs are age-adjusted multiplicatively, then the picture becomes more complicated. We would need to know the mean (or sum) of ages within each separate group for each health outcome.
The number of individuals in each health group will also need to be tracked whether or not age depreciation is accounted for. Within the prevalence module, new trackers should be added akin to `dat$epi$num.nohiv.prep[at] <- length(ids.nohiv.prep)` seen at the end of `track_util_msm()` within the vignette `cea_module_development.Rmd`.

### Modules
The **cost-effectiveness module** (function = `costeffect`) calculates the costs and QALYs accrued by individuals along with relevant intervention costs at each time-step. Individual costs and QALYs are functions of attributes within the model.

The **aging module** (function = `aging`) sets up the age attribute at the initial time step by pulling from the network object. Then it will subsequently update the age attribute as updating the age in increments of a week at each time step.

The **death module** (function = `dfunc`)  simulates mortality as a function of an age-specific mortality rate and disease-induced mortality. 

The **birth module** (function = `bfunc`) implements a more simplified birth process compared to the built-in birth module, mainly for greater legibility. Note that the `age` attribute must be set on both the `attr` list and the `nw` object for new incoming nodes. 

### Parameters
The epidemic model parameters are basic here because we're not changing any of the core epidemiology from a simple SI model.

* `cea.start`: sets the time-step at which the accounting for costs and effects begins to take place.

* `end.horizon`: When we reach the end of the analytic time window, we enter the "end horizon". During the end horizon, further HIV transmission, PrEP usage, and new entry to the population are disabled, and the health status and treatment behaviors of individuals remaining in the population are simulated until death. The purpose of extending the time horizon in this way is to capture the full costs and effects accrued by individuals remaining at the end of the analytic time horizon. If the time horizon is abruptly truncated, then the QALYs gained from averting HIV infections towards the end of the window will not be fully captured. This is especially important when the interventions have a significant impact on mortality.

* `d_r`: The discount rate (`d_r`) is used to deprecate costs accrued in the simulation as a function of how much time has passed since the beginning of the analytic time horizon. The general rationale for discounting in CEA is that costs and benefits that are deferred have lower value than those that are realized immediately. The discount rate (`d_r`) is used to deprecate costs accrued in the simulation as a function of how much time has passed since the beginning of the analytic time horizon. The general rationale for discounting in CEA is that costs and benefits that are deferred have lower value than those that are realized immediately.

* `inf.prob`: the probability that an infection will occur given an act between a susceptible and infected node in the absence of the intervention. 
* `mortality.rates`: a vector of mortality rates, where the position on the vector correspond to the mortality rate of persons of `age+1` (for example, the first rate is for persons less than one year of age). 
* `birth.rate`: a scalar for the rate of births per person per week.


end.horizon <- 52 + 14
  sus.cost <- 150
  inf.cost <- 300
  sus.qaly <- 1.00
  inf.qaly <- 0.75
  age.decrement <- -0.003
  disc.rate <- 0.03
  cea.start <- 14
  nsteps <- 118
  inter.eff <- 0.50
  inter.start <- 14
  inter.cost <- 500000

## Author
Greg Knowlton, University of Minnesota
