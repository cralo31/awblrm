library(rstan)
rstan_options(auto_write = TRUE)
mod = stan_model('model.stan', verbose = FALSE)

#'Next Recommended Dose
#'
#'Function to determine the next recommended dose according to the input data.
#'
#'@export
#'@param dose_levels A vector of length d, where d is the number of dose levels considered.
#'@param max_cycle Integer; Number of total treament cycles.
#'@param cycle_time Double; The length of a single treatment cycle.
#'@param mode Integer; Choose the decision rule to implement.
#'@param target The defined target probability of toxicity. For point-estimate-based decision rules, please input a single numeric value.
#'For interval-estimate-based decision rules, please input a vector indicating the target toxicity range, e.g. c(0.2, 0.33).
#'@param od_control Double(0-1); Escalation with overdose control.
#'@param cohort Vector; The cohort label of all patients currently observed.
#'@param dlt Vector; If toxic response was observed during the patient's evaluation window.
#'@param dose Vector;The dose level administered to the patients.
#'@param dlt_cycle Vector;Cycle in which DLT is observed for each patient. If no DLT please input -1.
#'@param obs_time Vector; How long has each patient been in the trial.
#'@param mcmc_chains Integer; Number of chains to run for MCMC sampling.
#'@param mcmc_samples Integer; Number of posterior samples to obtain
#'
#'@return
#'
#'A list with two elements. The first element is the next recommended dose. The second element is a summary information
#'matrix of the implemented decision rule.
#'
#'@examples
#'#################################################################################
#'# Consider a 3-cycle treatment scenario with observed information of 9 patients #
#'
#'all_doses = c(1, 2, 3, 4, 5)
#'max_cycle = 3
#'cycle_time = 30
#'mode = 1
#'target = 0.25
#'od_control = 1
#'cohort = c(1, 1, 1, 2, 2, 2, 3, 3, 3)
#'dlt = c(0, 0, 0, 1, 0, 0, 1, 0, 1)
#'dose = c(1, 1, 1, 2, 2, 2, 2, 2, 2)
#'dlt_cycle = c(-1, -1, -1, 2, -1, -1, 3, -1, 2)
#'obs_time = c(90, 45, 90, 35, 90, 75, 80, 90, 55)
#'
#'test = decision.dose(all_doses, max_cycle, cycle_time, mode, target, od_control, cohort, dlt,
#'dose, dlt_cycle, obs_time)
#'


decision.dose = function(dose_levels, max_cycle, cycle_time, mode, target, od_control, cohort, dlt, dose, dlt_cycle, obs_time,
                         mcmc_chains = 4, mcmc_samples = 3000) {

  # Calculate weight of each observation according to observed information
  weight = find_weight2(length(cohort), dlt, obs_time, cycle_time, max_cycle, dlt_cycle)

  # Combine information into a matrix
  mat = data.frame(cbind(cohort, dlt, dose, dlt_cycle, obs_time, weight))

  # Run stan to find the betas
  fit = run_stan2(mat, 500, mcmc_samples, mcmc_chains, 1, 0.95)

  # Determine next recommended dose depending on implemented decision rule
  decision = decision_rule(dose_levels, mat$dose, tail(mat$dose, 1), target, fit[,1], fit[,2], mode, od_control)

  if (decision$next_dose == -1) {
    return("All doses are identified to be too toxic. No MTD found.")
  }

  # Final information table to return depending on
  if (mode == 1) {  # mode 1
    vital = rbind(decision$info[c(3,5),], decision$safety)
    rownames(vital) = c("Prob DLT", "Prob Overdose", "Patient Protection")
    return(list(decision = decision$next_dose, information = round(vital, 4)))

  } else if (mode == 2) {  # mode 2
    vital = rbind(decision$info[c(3,4),], decision$safety)
    rownames(vital) = c("Dist to d_hat", "Prob Overdose", "Patient Protection")
    return(list(decision = decision$next_dose, information = round(vital, 4), d_hat = decision$d_hat))

  } else if (mode == 3) {  # mode 3
    vital = rbind(decision$info[c(3,4), decision$safety])
    rownames(vital) = c("Target Probability", "Prob Overdose", "Patient Protection")
    return(list(decision = decision$next_dose, information = round(vital, 4)))

  } else {  # mode 4
    vital = rbind(decision$info[3,], decision$ratio[,2], decision$info[4,], decision$safety)
    rownames(vital) = c("Prob DLT", "Proportion of Overlap", "Prob Overdose", "Patient Protection")
    return(list(decision = decision$next_dose, information = round(vital, 4)))
  }
}

############################################################################
# Weight function

find_weight2 = function(size, obs_DLT, obs_time, cycle_time, K, cycle) {

  # Linear Weight (1 if DLT==1 OR actual observation time reaches the full assessment time)
  Gt = numeric(size)
  for (i in 1:size) {
    if (obs_DLT[i] || (obs_time[i] %% cycle_time == 0)) {
      Gt[i] = 1
    } else {
      Gt[i] = (obs_time[i] %% cycle_time) / cycle_time
    }
  }

  # Update cycle information
  alpha = rep(1, K) # Prior for dirichlet distribution (assume all a's are 1)

  # Update posterior distribution of P given observed cycle information
  for (i in 1:size) {
    if (obs_DLT[i] & cycle[i] >= 1) {
      alpha[cycle[i]] = alpha[cycle[i]] + 1
    }
  }
  p_hat = alpha / (sum(alpha))  # Posterior mean of P|Z (Updated probaility of cycle information)

  cyc_assess = ceiling(obs_time / cycle_time)
  # Combine P_hat and G(t) to calculate the final adaptive cycle weight of the observed data
  weight = numeric(size)
  for (i in 1:size) {

    # Weight is equal to 1 if DLT is observed during assessment window, else the weight is calculated linearly
    if (obs_DLT[i] == 1) {
      weight[i] = 1
      next
    }
    if (cyc_assess[i] > 1) {
      for (j in 1:(cyc_assess[i]-1)) {
        weight[i] = weight[i] + p_hat[j]
      }
    }
    weight[i] = weight[i] + p_hat[cyc_assess[i]] * Gt[i]
  }
  info = cbind(Gt, weight)

  return(info)
}

###############################################################################################
### Part 2: Run MCMC via rStan
### model.stan saved in a seperate file

#' Bayesian Logistic Regression with Stan
#'
#' @export
#' @param mat Current information matrix.
#' @param warmup Number of burn-in samples.
#' @param iter Number of posterior samples to obtain.
#' @param chains Number of chains to run.
#' @param samp_cores Number of cores to perform sampling. Recommend to 1 if iter < 10000.
#' @param cont Step size of the path search algorithm.
#' @return The posterior samples of the slope and intercept.

run_stan2 = function(mat, warmup, iter, chains, cores, cont) {

  # Save into rstan required format
  stan_dat = list(N = length(mat$dlt), y = mat$dlt, dose = mat$dose, W = mat$weight)

  sink(tempfile())
  # Run the MCMC
  fit = sampling(mod, data = stan_dat, warmup = warmup, iter = iter, chains = chains,
                 cores = cores, verbose = FALSE, refresh = -1, control = list(adapt_delta = cont))
  sink(NULL)

  # Extract posterior samples
  extract_fit = extract(fit)
  alpha = extract_fit$alpha; beta = extract_fit$beta

  return(data.matrix(cbind(alpha, beta)))

}



################################################################################################
### Rule 1: Dose Level that gives the closest probability of toxicity from target toxicity.
#
# all_dose: Predefined doses considered for the experiment.
# dose: Actual dose information of the patients.
# current: Current dose level.
# target: Target MTD probability.
# beta: Posterior mean of the estimated coeffcient of dose-toxicity curve.

rule_1 = function(all_dose, dose, current, target, alpha, beta, od_control, od_prob) {

  # Add levels and indiator function
  levels = 1:length(all_dose)
  tested = ifelse(all_dose %in% dose, 1, 0)

  # Calculate estimated toxicity probability given dose and beta
  prob_hat = exp(beta * all_dose + alpha) / (1 + exp(beta * all_dose + alpha))

  # Calculate distance between estimated and targeted probabilty
  dist = abs(prob_hat - rep(target, length(all_dose)))

  # Store info within a matrix
  info = rbind(all_dose, levels, tested, prob_hat, dist, od_prob)
  colnames(info) = info[1,]
  info = info[-1,]

  # Escalation / Expansion / De-escalation (no dose skipping if dose level has not been tested before)
  safe_dose = as.matrix(info[,which(info[5,] <= od_control)]) # Apply overdose control
  if (length(safe_dose) == 0) {
    return(list(info = info, next_dose = -1))
  }

  cand = as.numeric(safe_dose[1, which.min(safe_dose[4,])])
  level_diff = cand - current
  if (level_diff > 1 & info[2, which(colnames(info) == cand)] == 0) {
    cand = as.numeric(colnames(info)[which(colnames(info) == current) + 1])
  }

  return(list(info = info, next_dose = cand))
}

############################################################################
### Rule 2: CRM-based. Dose level computed from target closest to dose level

rule_2 = function(all_dose, dose, current, target, alpha, beta, od_control, od_prob) {

  # Set up info matrix, levels and tested indicator
  levels = 1:length(unique(all_dose))
  tested = ifelse(all_dose %in% dose, 1, 0)

  # Estimate the "target" dose given the target toxicity rate and beta
  target_dose = (log(target / (1 - target)) - alpha) / beta

  # Estimate the maximum "safe" dose before crossing overdose threshold (33%)
  safe_thres = (log(0.33 / (1 - 0.33)) - alpha) / beta

  # Difference between target dose and pre-defined all_dose
  diff = abs(target_dose - all_dose)

  # Bind info together
  info = rbind(all_dose, levels, tested, diff, od_prob)
  colnames(info) = all_dose
  info = info[-1,]

  # Escalation / Expansion / De-escalation (no dose skipping if dose level has not been tested before)
  safe_dose = as.matrix(info[,which(info[4,] <= od_control)]) # Apply overdose control
  if (length(safe_dose) == 0) {
    return(list(info = info, next_dose = -1))
  }

  cand = as.numeric(safe_dose[1, which.min(safe_dose[3,])])
  level_diff = cand - current
  if (level_diff > 1 & info[2, which(colnames(info) == cand)] == 0) {
    cand = as.numeric(colnames(info)[which(colnames(info) == current) + 1])
  }

  return(list(info = info, next_dose = cand, d_hat = target_dose))
}

###########################
### Rule 3: Interval-based

#target = c(0.2, 0.3)

rule_3 = function(all_dose, dose, current, target, alpha, beta, od_control, od_prob) {

  # Stop algorithm and throw error message if target is not a length-2 vector (interval)
  if (length(target) != 2) {
    stop("Target toxicity is not a interval. A 2-tuple with the lower and upper bound of the target interval is required.")
  }

  # Set up info matrix, levels and tested indicator
  levels = 1:length(unique(all_dose))
  tested = ifelse(all_dose %in% dose, 1, 0)
  info = rbind(all_dose, levels, tested)
  colnames(info) = all_dose
  info = info[-1,]

  # Calculate posterior interval
  res = c()
  for (i in 1:length(all_dose)) {
    post_y = exp(beta * all_dose[i] + alpha) / (1 + exp(beta * all_dose[i] + alpha))
    below = mean(post_y < min(target))
    bingo = mean(post_y > min(target) & post_y < max(target))
    above = mean(post_y > 0.33)
    output = c(all_dose[i], below, bingo, above)
    res = rbind(res, output)
  }
  colnames(res) = c("Dose", "Below Target", "On Target", "Above Target")

  info = rbind(info, res[,3], res[,4]);
  rownames(info) = c("levels", "tested", "target", "overdose")

  # Escalation / Expansion / De-escalation (no dose skipping if dose level has not been tested before)
  safe_dose = as.matrix(info[,which(info[4,] <= od_control)]) # Apply overdose control
  if (length(safe_dose) == 0) {
    return(list(info = info, next_dose = -1))
  }
  cand = as.numeric(safe_dose[1, which.max(safe_dose[3,])])
  level_diff = cand - current
  if (level_diff > 1 & info[2, which(colnames(info) == cand)] == 0) {
    cand = as.numeric(colnames(info)[which(colnames(info) == current) + 1])
  }

  return(list(info = info, post_interval = res, next_dose = cand))

}

##########################
### Rule 4: Interval Ratio

# Function to calculate midpoint of two points
midpoint = function(vec) {
  res = c()
  for (i in 1:(length(vec)-1)) {
    mid = (vec[i] + vec[i+1]) / 2
    res = c(res, mid)
  }
  res = c(vec[1], res, vec[length(vec)])
  return(res)
}

# Ratio Calculation
ratio = function(midpoints, target) {
  res = c()
  for (i in 1:(length(midpoints) - 1)) {
    range = Overlap(c(midpoints[i], midpoints[i+1]), target)
    ratio = range / (midpoints[i + 1] - midpoints[i])
    #ratio = range / (max(target) - min(target))
    res = rbind(res, c(range, ratio))
  }
  colnames(res) = c("range", "ratio")
  return (data.frame(res))
}


# The actual rule function
rule_4 = function(all_dose, dose, current, target, alpha, beta, od_control, od_prob) {

  # Stop algorithm and throw error message if target is not a length-2 vector (interval)
  if (length(target) != 2) {
    stop("Target toxicity is not a interval. A 2-tuple with the lower and upper bound of the target interval is required.")
  }

  temp_dose = c(min(all_dose) - 0.5, all_dose, max(all_dose) + 0.5)

  # Set up levels and tested indicator
  levels = 1:length(unique(all_dose))
  tested = ifelse(all_dose %in% dose, 1, 0)

  # Create info matrix
  info = rbind(all_dose, levels, tested)
  colnames(info) = all_dose
  info = info[-1,]

  # Calculate estimated toxicity probability given dose and beta
  prob_hat = exp(beta * temp_dose + alpha) / (1 + exp(beta * temp_dose + alpha))
  temp_prob = prob_hat[-c(1, length(prob_hat))]
  info = rbind(info, temp_prob, od_prob)

  # Calculate midpoint of estimated toxicity probabilities
  intervals = midpoint(temp_prob)
  intervals[c(1,length(intervals))] = prob_hat[c(1, length(prob_hat))]

  # Decision rule for candidate dose
  if (max(target) < min(prob_hat)) {
    ratio_tab  = 0; safe_dose = info
    cand = min(all_dose)
  } else if (min(target) > max(prob_hat)) {
    ratio_tab = 0; safe_dose = info
    cand = max(all_dose)
  } else {
    ratio_tab = ratio(intervals, target)
    index = which.max(ratio_tab$ratio)
    safe_dose = as.matrix(info[,which(info[4,] <= od_control)]) # Apply overdose control
    if (length(safe_dose) == 0) {
      return(list(info = info, next_dose = -1))
    }
    cand = ifelse(index <= ncol(safe_dose), safe_dose[1,index], current)
  }

  level_diff = cand - current
  if (level_diff > 1 & info[2, which(colnames(info) == cand)] == 0) {
    cand = as.numeric(colnames(info)[which(colnames(info) == current) + 1])
  }
  else if (level_diff < 0){
    cand = as.numeric(colnames(info)[which(colnames(info) == current) - 1])
  }

  return(list(info = info, ratio = ratio_tab, next_dose = cand))
}

#############################################################################
# Wrapper function for decision rule

decision_rule = function(all_dose, dose, current, target, alpha, beta, mode, od_control) {

  # od_prob: calculates the posterior probability of the probability of toxicity falling above 33%
  # safety_prob: calculates the posterior probability of the probability of toxicity falling above target MTD
  safety_prob = c(); od_prob = c()
  for (i in 1:length(all_dose)) {
    post_y = exp(beta * all_dose[i] + alpha) / (1 + exp(beta * all_dose[i] + alpha))
    od = mean(post_y > 0.33); od_prob = c(od_prob, od)
    safe = mean(post_y > target); safety_prob = c(safety_prob, safe)
  }

  # If the current dose is extremely toxic, de-escalate immediately and remove any dose levels equal to or higher.
  # In particular, if the first dose is extremely toxic, then the trial is terminated.

  if (safety_prob[current] > 0.95) {
    if (all_dose[1] == current) {
      decision = list(dose_tab = 0, next_dose = -1)
    } else {
      decision = list(dose_tab = all_dose[-c(current:length(all_dose))], next_dose = current - 1)
    }

    # Special case when the only available dose level is the starting dose. This happens after applying the
    # safety precautions in certain trial scenarios.
  } else if (length(all_dose) == 1){
    decision = list(dose_tab = all_dose[current], next_dose = all_dose[current])
  } else {

    # Rule 3 (Interval based BLRM)
    if (mode == 3) {
      decision = rule_3(all_dose, dose, current, target, alpha, beta, od_control = od_control, od_prob)
    } else {

      alpha = mean(alpha); beta = mean(beta)
      # Rule 1
      if (mode == 1) {
        decision = rule_1(all_dose, dose, current, target, alpha, beta, od_control = od_control, od_prob)

        # Rule 2
      } else if (mode == 2) {
        decision = rule_2(all_dose, dose, current, target, alpha, beta, od_control = od_control, od_prob)

        # Rule 4
      } else {
        decision = rule_4(all_dose, dose, current, target, alpha, beta, od_control = od_control, od_prob)
      }
    }
    decision$dose_tab = all_dose
  }

  decision$safety = safety_prob
  return(decision)
}
