
##
## Simple Cost-effectiveness Model (Simple SI model with cost and utility tracking)
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Gregory Knowlton (University of Minnesota)
## Date: October 2021
##

# Cost/Utility Tracking Module ------------------------------------------------

ce <- function(dat, at) {
  
  # #
  # age <- get_attr()
  # 
  # ## Summary statistics ##
  # dat <- set_epi()
  
  return(dat)
}

# Updated Aging Module --------------------------------------------------------

aging <- function(dat, at) {

  # Update age on attr and also the network
  age <- get_attr(dat, "age")
  age <- age + 1/52
  dat <- set_attr(dat, "age", age)

  ## Summary statistics ##
  dat <- set_epi(dat, "meanAge", at, mean(age, na.rm = TRUE))

  return(dat)
}


# Updated Departure Module -----------------------------------------------------

dfunc <- function(dat, at) {
  # browser()
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
  nDepts <- 0

  if (nElig > 0) {

    ## Calculate age-specific departure rates for each eligible node
    whole_ages_of_elig <- ceiling(age[idsElig])
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
    idsRetire <- setdiff(which(age >= 65), idsDeaths)
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
  dat <- set_epi(dat, "d.flow", at, nDepts)
  
  return(dat)
}


# Updated Arrivals Module ----------------------------------------------------

afunc <- function(dat, at) {

  ## Parameters
  n <- sum(get_attr(dat, "active_s") == 1)
  a.rate <- get_param(dat, "arrival.rate")

  ## Process
  nArrivalsExp <- n * a.rate
  nArrivals <- rpois(1, nArrivalsExp)

  ## Update attributes
  if (nArrivals > 0) {
    dat <- append_core_attr(dat, at, nArrivals)
    dat <- append_attr(dat, "status", "s", nArrivals)
    dat <- append_attr(dat, "infTime", NA, nArrivals)
    dat <- append_attr(dat, "age", 0, nArrivals)
    dat <- append_attr(dat, "active_s", 1, nArrivals)
  }

  ## Summary statistics
  dat <- set_epi(dat, "a.flow", at, nArrivals)

  return(dat)
}

simnet_msm <- function(dat, at) {
  # NEW CODE
  ##########
  if (at > dat$param$end.horizon) {
    return(dat)
  }
  
  inactive_ids <- which(dat$attr$active == 0)
  for (i in 1:3) {
    rm <- union(which(dat$el[[i]][, 1] %in% inactive_ids), which(dat$el[[i]][, 2] %in% inactive_ids))
    if (length(rm > 0)) {
      dat$el[[i]] <- as.edgelist(x = dat$el[[i]][-rm, ], n = length(dat$attr$age))
    }
  }
  
  # ## Edges correction
  # dat <- edges_correct_msm(dat, at)
  # 
  # ## Main network
  # nwparam.m <- EpiModel::get_nwparam(dat, network = 1)
  # 
  # dat$attr$deg.casl <- get_degree(dat$el[[2]])
  # dat <- tergmLite::updateModelTermInputs(dat, network = 1)
  # 
  # dat$el[[1]] <- tergmLite::simulate_network(p = dat$p[[1]],
  #                                            el = dat$el[[1]],
  #                                            coef.form = nwparam.m$coef.form,
  #                                            coef.diss = nwparam.m$coef.diss$coef.adj,
  #                                            save.changes = TRUE)
  # 
  # plist1 <- update_plist(dat, at, ptype = 1)
  # 
  # 
  # ## Casual network
  # nwparam.p <- EpiModel::get_nwparam(dat, network = 2)
  # 
  # dat$attr$deg.main <- get_degree(dat$el[[1]])
  # dat <- tergmLite::updateModelTermInputs(dat, network = 2)
  # 
  # dat$el[[2]] <- tergmLite::simulate_network(p = dat$p[[2]],
  #                                            el = dat$el[[2]],
  #                                            coef.form = nwparam.p$coef.form,
  #                                            coef.diss = nwparam.p$coef.diss$coef.adj,
  #                                            save.changes = TRUE)
  # 
  # plist2 <- update_plist(dat, at, ptype = 2)
  # 
  # dat$temp$plist <- rbind(plist1, plist2)
  # if (dat$control$truncate.plist == TRUE) {
  #   to.keep <- which(is.na(dat$temp$plist[, "stop"]))
  #   dat$temp$plist <- dat$temp$plist[to.keep, ]
  # }
  # 
  # ## One-off network
  # nwparam.i <- EpiModel::get_nwparam(dat, network = 3)
  # 
  # dat$attr$deg.tot <- pmin(dat$attr$deg.main + get_degree(dat$el[[2]]), 3)
  # dat <- tergmLite::updateModelTermInputs(dat, network = 3)
  # 
  # dat$el[[3]] <- tergmLite::simulate_ergm(p = dat$p[[3]],
  #                                         el = dat$el[[3]],
  #                                         coef = nwparam.i$coef.form)
  # 
  # if (dat$control$save.nwstats == TRUE) {
  #   dat <- calc_nwstats(dat, at)
  # }
  
  return(dat)
}
