
##
## Simple Cost-effectiveness Model (Simple SI model with cost and utility tracking)
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Gregory Knowlton (University of Minnesota)
## Date: October 2021
##

# Cost/Utility Tracking Module ------------------------------------------------

costeffect <- function(dat, at) {
  
  # Import attributes
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  status <- get_attr(dat, "status")
  
  # Import parameters
  sus_cost <- get_param(dat, "sus_cost")
  inf_cost <- get_param(dat, "inf_cost")
  sus_qaly <- get_param(dat, "sus_qaly")
  inf_qaly <- get_param(dat, "inf_qaly")
  age_decrement <- get_param(dat, "age_decrement")
  disc_rate <- get_param(dat, "disc_rate")
  
  # Identify relevant sub-populations
  idsSus <- which(active == 1 & status == "s")
  idsInf <- which(active == 1 & status == "i")
  
  # Account for costs of each sub-population
  pop_sus_cost <- length(idsSus) * sus_cost
  pop_inf_cost <- length(idsInf) * inf_cost
  
  # Account for effects (QALYs) of each sub-population
  # QALY parameters by health state and age decrements must be converted to weekly time-step
  pop_sus_qaly <- sum((age[idsSus] * age_decrement + sus_qaly) / 52, na.rm = TRUE)
  pop_inf_qaly <- sum((age[idsInf] * age_decrement + inf_qaly) / 52, na.rm = TRUE)
  
  # Aggregate costs and effects
  pop_cost <- pop_sus_cost + pop_inf_cost
  pop_qaly <- pop_sus_qaly + pop_inf_qaly
  
  # Discount aggregated costs and effects
  pop_cost_disc <- pop_cost * (1 - disc_rate) ^ (at / 52)
  pop_qaly_disc <- pop_qaly * (1 - disc_rate) ^ (at / 52)
  
  ## Summary statistics ##
  dat <- set_epi(dat, "cost", at, pop_cost)
  dat <- set_epi(dat, "qaly", at, pop_qaly)
  dat <- set_epi(dat, "cost.disc", at, pop_cost_disc)
  dat <- set_epi(dat, "qaly.disc", at, pop_qaly_disc)
  
  return(dat)
}

# Updated Aging Module --------------------------------------------------------

aging <- function(dat, at) {

  # Update age on attr and also the network
  active <- get_attr(dat, "active")
  idsActive <- which(active == 1)
  age <- get_attr(dat, "age")
  age[idsActive] <- age[idsActive] + 1/52
  dat <- set_attr(dat, "age", age)

  ## Summary statistics ##
  dat <- set_epi(dat, "meanAge", at, mean(age[idsActive], na.rm = TRUE))

  return(dat)
}


# Updated Departure Module -----------------------------------------------------

dfunc <- function(dat, at) {
  
  ## Attributes
  active <- get_attr(dat, "active")
  active_s <- get_attr(dat, "active_s")
  exitTime <- get_attr(dat, "exitTime")
  age <- get_attr(dat, "age")
  status <- get_attr(dat, "status")

  ## Parameters
  death.rates <- get_param(dat, "death.rates")
  end.horizon <- get_param(dat, "end.horizon")
  
  ## Query active
  idsElig <- which(active == 1)
  nElig <- length(idsElig)
  nDeaths <- 0

  if (nElig > 0) {

    ## Calculate age-specific departure rates for each eligible node
    ## Everyone older than 85 gets the final mortality
    whole_ages_of_elig <- pmin(ceiling(age[idsElig]), 101)
    death_rates_of_elig <- death.rates[whole_ages_of_elig]
    
    ## Simulate departure process
    vecDeaths <- which(rbinom(nElig, 1, 1 - exp(-death_rates_of_elig)) == 1)
    idsDeaths <- idsElig[vecDeaths]
    nDeaths <- length(idsDeaths)

    ## Update nodal attributes
    if (nDeaths > 0) {
      active_s[idsDeaths] <- 0
      active[idsDeaths] <- 0
      exitTime[idsDeaths] <- at
    }
    
    ## 65+ who did not die this time-step
    idsRetire <- setdiff(which(age >= 65 & active_s == 1), idsDeaths)
    active_s[idsRetire] <- 0
    
  }
  
  if (at == end.horizon) {
    active_s <- rep(0, length(active_s))
  }

  ## Reset attr
  dat <- set_attr(dat, "active_s", active_s)
  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)

  ## Summary statistics
  dat <- set_epi(dat, "d.flow", at, nDeaths)
  
  return(dat)
}


# Updated Arrivals Module ----------------------------------------------------

afunc <- function(dat, at) {
  # browser()
  
  end.horizon <- get_param(dat, "end.horizon")
  
  if (at >= end.horizon) {
    return(dat)
  }

  ## Parameters
  n <- sum(get_attr(dat, "active") == 1)
  a.rate <- get_param(dat, "arrival.rate")

  ## Process
  nArrivalsExp <- n * a.rate
  nArrivals <- rpois(1, nArrivalsExp)

  ## Update attributes
  if (nArrivals > 0) {
    dat <- append_core_attr(dat, at, nArrivals)
    dat <- append_attr(dat, "status", "s", nArrivals)
    dat <- append_attr(dat, "infTime", NA, nArrivals)
    dat <- append_attr(dat, "age", 16, nArrivals)
    dat <- append_attr(dat, "active_s", 1, nArrivals)
  }

  ## Summary statistics
  dat <- set_epi(dat, "a.flow", at, nArrivals)

  return(dat)
}

# Updated Infection Module ---------------------------------------------------

ifunc <- function (dat, at) {
  
  end.horizon <- get_param(dat, "end.horizon")
  
  # 
  if (at >= end.horizon) {
    return(dat)
  }
  
  #
  active_s <- get_attr(dat, "active_s")
  #
  idsActive_s <- which(active_s == 1)
  
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat, "act.rate")
  inter.eff <- get_param(dat, "inter.eff", override.null.error = TRUE)
  inter.start <- get_param(dat, "inter.start", override.null.error = TRUE)
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)
  nInf <- 0
  if (nElig > 0 && nElig < nActive) {
    del <- discord_edgelist(dat, at, network = 1)
    if (!(is.null(del))) {
      #
      del <- del[which(del$sus %in% idsActive_s & del$inf %in% idsActive_s),]
      del$infDur <- at - infTime[del$inf]
      del$infDur[del$infDur == 0] <- 1
      linf.prob <- length(inf.prob)
      del$transProb <- ifelse(del$infDur <= linf.prob, 
                              inf.prob[del$infDur], inf.prob[linf.prob])
      if (!is.null(inter.eff) && at >= inter.start) {
        del$transProb <- del$transProb * (1 - inter.eff)
      }
      lact.rate <- length(act.rate)
      del$actRate <- ifelse(del$infDur <= lact.rate, act.rate[del$infDur], 
                            act.rate[lact.rate])
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]
      idsNewInf <- unique(del$sus)
      status <- get_attr(dat, "status")
      status[idsNewInf] <- "i"
      dat <- set_attr(dat, "status", status)
      infTime[idsNewInf] <- at
      dat <- set_attr(dat, "infTime", infTime)
      nInf <- length(idsNewInf)
    }
  }
  if (nInf > 0) {
    dat <- set_transmat(dat, del, at)
  }
  dat <- set_epi(dat, "si.flow", at, nInf)
  return(dat)
}

# Updated Network Resimulation Module ----------------------------------------

resimfunc <- function (dat, at) {
  
  # This function is identical to the default of EpiModel::resim_nets()
  # except for the following line.
  # Skipping over this function after the start of the end horizon greatly
  # improves computational speed.
  end.horizon <- get_param(dat, "end.horizon")
  if (at >= end.horizon) {
    return(dat)
  }
  
  tergmLite <- get_control(dat, "tergmLite")
  isTERGM <- get_control(dat, "isTERGM")
  save.nwstats <- get_control(dat, "save.nwstats")
  resimulate.network <- get_control(dat, "resimulate.network")
  nwstats.formula <- get_control(dat, "nwstats.formula")
  set.control.stergm <- get_control(dat, "set.control.stergm")
  tergmLite.track.duration <- get_control(dat, "tergmLite.track.duration")
  dat <- edges_correct(dat, at)
  active <- get_attr(dat, "active")
  idsActive <- which(active == 1)
  anyActive <- ifelse(length(idsActive) > 0, TRUE, FALSE)
  if (dat$param$groups == 2) {
    group <- get_attr(dat, "group")
    groupids.1 <- which(group == 1)
    groupids.2 <- which(group == 2)
    nActiveG1 <- length(intersect(groupids.1, idsActive))
    nActiveG2 <- length(intersect(groupids.2, idsActive))
    anyActive <- ifelse(nActiveG1 > 0 & nActiveG2 > 0, TRUE, 
                        FALSE)
  }
  nwparam <- get_nwparam(dat)
  if (anyActive == TRUE & resimulate.network == TRUE) {
    if (tergmLite == FALSE) {
      if (isTERGM == TRUE) {
        suppressWarnings(dat$nw[[1]] <- simulate(dat$nw[[1]], 
                                                 formation = nwparam$formation, dissolution = nwparam$coef.diss$dissolution, 
                                                 coef.form = nwparam$coef.form, coef.diss = nwparam$coef.diss$coef.adj, 
                                                 constraints = nwparam$constraints, time.start = at, 
                                                 time.slices = 1, time.offset = 0, monitor = nwstats.formula, 
                                                 control = set.control.stergm))
      }
      else {
        dat$nw[[1]] <- simulate(object = nwparam$formation, 
                                basis = dat$nw[[1]], coef = nwparam$coef.form, 
                                constraints = nwparam$constraints, dynamic = FALSE, 
                                monitor = nwstats.formula, nsim = 1)
      }
      if (save.nwstats == TRUE) {
        new.nwstats <- tail(attributes(dat$nw[[1]])$stats, 
                            1)
        keep.cols <- which(!duplicated(colnames(new.nwstats)))
        new.nwstats <- new.nwstats[, keep.cols, drop = FALSE]
        dat$stats$nwstats[[1]] <- rbind(dat$stats$nwstats[[1]], 
                                        new.nwstats)
      }
    }
    if (tergmLite == TRUE) {
      dat <- tergmLite::updateModelTermInputs(dat)
      if (isTERGM == TRUE) {
        rv <- tergmLite::simulate_network(state = dat$p[[1]]$state, 
                                          coef = c(nwparam$coef.form, nwparam$coef.diss$coef.adj), 
                                          control = dat$control$mcmc.control[[1]], save.changes = TRUE)
        dat$el[[1]] <- rv$el
        if (tergmLite.track.duration == TRUE) {
          dat$p[[1]]$state$nw0 %n% "time" <- rv$state$nw0 %n% 
            "time"
          dat$p[[1]]$state$nw0 %n% "lasttoggle" <- rv$state$nw0 %n% 
            "lasttoggle"
        }
      }
      else {
        rv <- tergmLite::simulate_ergm(state = dat$p[[1]]$state, 
                                       coef = nwparam$coef.form, control = dat$control$mcmc.control[[1]])
        dat$el[[1]] <- rv$el
      }
      if (save.nwstats == TRUE) {
        nwL <- tergmLite::networkLite(dat$el[[1]], dat$attr)
        if (tergmLite.track.duration == TRUE) {
          nwL %n% "time" <- dat$p[[1]]$state$nw0 %n% 
            "time"
          nwL %n% "lasttoggle" <- dat$p[[1]]$state$nw0 %n% 
            "lasttoggle"
        }
        nwstats <- summary(dat$control$nwstats.formulas[[1]], 
                           basis = nwL, term.options = dat$control$mcmc.control[[1]]$term.options, 
                           dynamic = isTERGM)
        keep.cols <- which(!duplicated(names(nwstats)))
        dat$stats$nwstats[[1]] <- rbind(dat$stats$nwstats[[1]], 
                                        nwstats[keep.cols])
      }
    }
  }
  return(dat)
}
