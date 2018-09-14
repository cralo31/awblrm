library(rstan)
library(LaplacesDemon)
library(DescTools)
library(doParallel)
library(foreach)

#setwd("C:/Users/yuton/Dropbox/Takeda/Code")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
mod = stan_model('model.stan', verbose = FALSE)

##################################################################################################
# Parallel computing for 1 trial of one model

#' Simulated Dose-finding Trial
#'
#' This function simulates a user defined number of dose-finding trials with the pre-specified design in parallel.
#'
#' @param true Matrix; Pre-defined true dose-toxicity relationship. The input should be a 2-by-p matrix, with the 1st row being the
#' p dose levels and the 2nd row being the corresponding probability of DLT.
#' @param n Integer; Number of trials to be simulated. Minimum of 1.
#' @param cores Integer; Number of cores to utilize for parallel computing. Defaults to 1.
#' @param K Integer; Number of treatment cycles considered.
#' @param size Integer; Number of patients to enrolled in each cohort.
#' @param cycle_prob K-simplex; A vector that contains the true probability of DLT happening in each of the cycles. Sums up to one.
#' @param cycle_time Integer; Duration of a single treatment cycle (NOT entire treatment window!). Units are arbitrary.
#' @param target The defined target probability of toxicity. For point-estimate-based decision rules, please input a single numeric value.
#' For interval-estimate-based decision rules, please input a vector indicating the target toxicity range, e.g. c(0.2, 0.33).
#' @param samples Integer; The maximum number of patients to be enrolled into the trial. Also serves as the default stopping rule.
#' @param assessment Integer; The cycle in which the model is assessed to determine the dose level of the next cohort (1 means the model is assessed
#' after the completion of the 1st cycle of each cohort, 2 means 2nd cycle...).
#' @param fixed Boolean; TRUE if patients in the same cohort are enrolled at the same time. FALSE if the enrollment time of patients in the same
#' cohort follows a Uniform (0, cycle_time).
#' @param mode Integer; Determine which decision rule to implement to perform dose escalation / de-escalation. 1 = Original CRM, 2 = Inverse Point-Estimate,
#' 3 = BLRM, 4 = Ratio. Defaults to 1. For details on the decision rule please refer to the report.
#' @param stopping Integer; Defined the stopping criteria. 1 = stop the trial when maximum enrollment is achieved, 2 = stop the trial if at least 12
#' patients are enrolled and 9 patients have already been treated at the next recommended dose. Defaults to 1.
#' @param least Integer; Least amount of patients to enrolled if using the second stopping rule. Defaults to 12, directly impacts "stopping".
#' @param dropout Boolean; TRUE if patient dropout time is assumened to follow an exponential distribution with a cumulative probability of 15% at the end
#' of the treatment window.
#' @param od_control Double (0 - 1); Escalation with overdose control. Escalate to the next dose if the probability of overdosing does not exceed the
#' input probability. No overdose control if od_control is set to 1.
#' @param mcmc_chains Integer; Number of chains to consider during MCMC sampling. Defaults to 3.
#' @param mcmc_samples Integer; Numbers of samples from the posterior distribution of the betas. Default to 3000 with a burn-in period of 500
#' samples.
#' @return Returns a list with two sublist.
#'
#' ###########################################################################
#'
#' The first sublist contains the two information tables of the simulated trials.
#'
#' The first table provides the allocation of patients and number of DLTs averaged over all simulations.
#' The second table gives the summary information of the simulations. Below are the meaning of the summary information.
#'
#' # The first four probabilities sums up to 1.
#'
#' P(Choose True MTD): Probability of choosing the true MTD as the final recommended dose.
#'
#' P(Choose Over MTD): Probability of choosing a dose level above the true MTD as the final recommended dose.
#'
#' P(Choose Under MTD): Probability of choosing a dose level below the true MTD as the final recommended dose.
#'
#' P(No Selection): Probability of terminating the trial due to all dose levels deemed toxic.
#'
#' Mean Final Dose: Numeric average of the final recommended dose.
#'
#' Allo Under / True / Over MTD: Average allocation of patients under / at / over the true MTD during trials.
#'
#' Total / Prop DLT: Average total number and proportion of patients who had DLT.
#'
#' Patients Enrolled / Trial Duration: Average number of patients enrolled before meeting stopping criteria
#' and duration of the trials.
#'
#' Computation Time: Average computation time of a trial.
#'
#' ###########################################################################
#'
#' The second sublist gives the detail generated information table of each trial. The last element within the sublist
#' gives the total computation time of simulated trials.
#'
#'
#' @examples
#' ################################# Example 1 ################################
#' # Simple uniform 2-cycle treatment design simulation with 4 trials and 4 cores. #
#'
#' set.seed(1031)
#' dose1 = rbind(1:5, c(0.05, 0.1, 0.25, 0.4, 0.6))
#' n1 = 4
#' cores1 = 4
#' K1 = 2
#' size1 = 3
#' cycle_prob1 = c(0.5, 0.5)
#' cycle_time1 = 30
#' samples1 = 36
#' assessment1 = 1
#' fixed1 = TRUE
#' mode1 = 1
#' stopping1 = 1
#' least1 = 12
#' dropout1 = FALSE
#' od_control1 = 0.8
#' target1 = 0.25
#' mcmc_chains = 4
#' mcmc_samples = 1000
#'
#' # Simulate 4 trials in parallel with 4 cores according to the defined parameters. #
#' test = parallel.trial(dose1, n1, cores1, K1, size1, cycle_prob1, cycle_time1, target1, samples1, assessment1,
#' fixed1, mode1, stopping1, least1, dropout1, od_control1, mcmc_chains, mcmc_samples)
#'
#' ################################# Example 2 ################################
#' # 3-cycle treatment design simulation with 4 trials and 4 cores. #
#' # BLRM Decision Rule
#' # Early Dropout
#'
#' set.seed(1031)
#' dose2 = rbind(1:5, c(0.05, 0.1, 0.25, 0.4, 0.6))
#' n2 = 4
#' cores2 = 4
#' K2 = 3
#' size2 = 3
#' cycle_prob2 = c(0.25, 0.25, 0.5)
#' cycle_time2 = 30
#' samples2 = 36
#' assessment2 = 2
#' fixed2 = TRUE
#' mode2 = 3
#' stopping2 = 2
#' least2 = 12
#' dropout2 = TRUE
#' od_control2 = 0.75
#' target2 = c(0.2, 0.33)
#' mcmc_chains = 3
#' mcmc_samples = 2000
#'
#' # Simulate 4 trials in parallel with 4 cores according to the defined parameters. #
#' test = parallel.trial(dose2, n2, cores2, K2, size2, cycle_prob2, cycle_time2, target2, samples2, assessment2,
#' fixed2, mode2, stopping2, least2, dropout2, od_control2, mcmc_chains, mcmc_samples)
#'
#' @export

parallel.trial = function(true, n, cores = 1, K, size, cycle_prob, cycle_time, target, samples, assessment, fixed, mode,
                    stopping, least, dropout, od_control, mcmc_chains = 4, mcmc_samples = 3000) {
  if (cores > detectCores()) {
    cores = detectCores()
    warning("Input number of clusters is greater than the number of clusters on the current machine.")
  } else if (cores < 1) {
    cores = 1
    warning("Input number of clusters is less than 1. Uses 1 core instead.")
  }

  cl = makeCluster(cores)
  registerDoParallel(cl)
  ptm = proc.time()
  run = foreach (i = 1:n, .packages = c("rstan", "LaplacesDemon", "DescTools", "doParallel", "foreach"),
                 .errorhandling = 'remove') %dopar% {
                   trial(true, K, size, cycle_prob, cycle_time, target, samples, assessment, fixed, mode,
                         stopping, least, dropout, od_control, mcmc_chains, mcmc_samples) }
  time = proc.time() - ptm
  stopCluster(cl)
  run[[n + 1]] = time[3]

  result = mean.res(run, true, target)
  result$allocation = round(result$allocation, 4); result$summary = round(result$summary, 4)

  return(list(result = result, all_trials = run))
}

##################################################################################################

#' Simulated Dose-finding Trial
#'
#' This function simulates one dose-finding trial with the pre-specified design.
#' @param true Matrix; Pre-defined true dose-toxicity relationship. The input should be a 2-by-p matrix, with the 1st row being the
#' p dose levels and the 2nd row being the corresponding probability of DLT.
#' @param K Integer; Number of treatment cycles considered.
#' @param size Integer; Number of patients to enrolled in each cohort.
#' @param cycle_prob K-simplex; A vector that contains the true probability of DLT happening in each of the cycles. Sums up to one.
#' @param cycle_time Integer; Duration of a single treatment cycle (NOT entire treatment window!). Units are arbitrary.
#' @param target The defined target probability of toxicity. For point-estimate-based decision rules, please input a single numeric value.
#' For interval-estimate-based decision rules, please input a vector indicating the target toxicity range, e.g. c(0.2, 0.33).
#' @param samples Integer; The maximum number of patients to be enrolled into the trial. Also serves as the default stopping rule.
#' @param assessment Integer; The cycle in which the model is assessed to determine the dose level of the next cohort (1 means the model is assessed
#' after the completion of the 1st cycle of each cohort, 2 means 2nd cycle...).
#' @param fixed Boolean; TRUE if patients in the same cohort are enrolled at the same time. FALSE if the enrollment time of patients in the same
#' cohort follows a Uniform (0, cycle_time).
#' @param mode Integer; Determine which decision rule to implement to perform dose escalation / de-escalation. 1 = Original CRM, 2 = Inverse Point-Estimate,
#' 3 = BLRM, 4 = Ratio. Defaults to 1. For details on the decision rule please refer to the report.
#' @param stopping Integer; Defined the stopping criteria. 1 = stop the trial when maximum enrollment is achieved, 2 = stop the trial if at least 12
#' patients are enrolled and 9 patients have already been treated at the next recommended dose. Defaults to 1.
#' @param least Integer; Least amount of patients to enrolled if using the second stopping rule. Defaults to 12, directly impacts "stopping".
#' @param dropout Boolean; TRUE if patient dropout time is assumened to follow an exponential distribution with a cumulative probability of 15% at the end
#' of the treatment window.
#' @param od_control Numeric (0 - 1); Escalation with overdose control. Escalate to the next dose if the probability of overdosing does not exceed the
#' input probability. No overdose control if od_control is set to 1.
#' @param mcmc_chains Integer; Number of chains to consider during MCMC sampling. Defaults to 3.
#' @param mcmc_samples Integer; Numbers of samples from the posterior distribution of the betas. Default to 3000 with a burn-in period of 500
#' samples.
#' @return Returns a list with three information tables.
#'
#' The first table gives the dose allocation of the patients during the trial(number of patients treated at each dose), including the corresponding
#' rate of DLT at each dose level.
#'
#' The second table provides a summary table of the results after the trial is finished. This table includes the
#' number of patients treated, the duration of the trial, the number / proportion of DLTs, the proportion of patients allocated in dose levels
#' that are under / on / over the true MTD, and the final recommended dose.
#'
#' The third table entails a full detailed information matrix of the trial.
#'
#'
#' @examples
#' # Set up the parameters for a simple simulation #
#' dose1 = rbind(1:5, c(0.05, 0.1, 0.25, 0.4, 0.6))
#' K1 = 2
#' size1 = 3
#' cycle_prob1 = c(0.5, 0.5)
#' cycle_time1 = 30
#' samples1 = 36
#' assessment1 = 1
#' fixed1 = TRUE
#' mode1 = 1
#' stopping1 = 1
#' least1 = 12
#' dropout1 = FALSE
#' od_control1 = 0.8
#' target1 = 0.25
#'
#' # Simulate one trial according to the defined parameters. #
#' trial(dose1, K1, size1, cycle_prob1, cycle_time1, target1, samples1, assessment1,
#' fixed1, mode1, stopping1, least1, dropout1, od_control1)
#' @export

trial = function(true, K, size, cycle_prob, cycle_time, target, samples,
                 assessment, fixed = TRUE, mode = 1, stopping = 1, least = 12, dropout, od_control,
                 mcmc_chains = 4, mcmc_samples = 3000) {

  # Basic initialization of parameters
  cohort = 1; index = 1 # dose index
  dose = true[1,]; true_prob = true[2,]

  # Administer the drug to the first cohort
  new_cohort = generate(cohort, 3, dose[index], true_prob[index], K, cycle_prob, cycle_time, assessment, cumu_days = 0, dropout = dropout)
  full_tab = new_cohort$info; cumulate = new_cohort$cumu_days;

  # Deterministic dose-escalation until first DLT is witnessed.
  while (sum(full_tab$obs_y) < 1 & nrow(full_tab) < samples) {

    # Special case to stop
    if (stopping == 2 & length(which(full_tab$dose == max(dose))) > 6) {
      break
    }

    # New cohort of patients comes into the study
    cohort = cohort + 1; index = ifelse(index == max(dose), max(dose), index + 1)
    new_cohort = generate(cohort, 3, dose[index], true_prob[index], K, cycle_prob,
                          cycle_time, assessment, cumu_days = cumulate, dropout = dropout)
    full_tab = rbind(full_tab, new_cohort$info)
    cumulate = new_cohort$cumu_days

    # Continously update previous cohorts' information as new patients are enrolled
    if (nrow(full_tab) > size) {
      full_tab = update(full_tab, cohort, cumulate, K, cycle_time, final = FALSE)
    }
  }

  # Enter modeling loop when first DLT is observed
  if (sum(full_tab$obs_y) > 0) {
    # Primary loop depending on the total sample size. Stop generating new patient info after sample size has been reached.
    while(nrow(full_tab) < samples) {

      # Run MCMC to obtain posterior samples of slope and intercept of dose toxicity curve
      fit = run_stan(full_tab, 500, mcmc_samples, mcmc_chains, 1, 0.975)

      # Current allocation of patients at each different dose level
      allocation = table(full_tab$dose)

      # Implement dose decision rule to find the next recommended dose
      curr_dose = full_tab[nrow(full_tab),]$dose;
      recommend = decision_rule(dose, full_tab$dose, curr_dose, target, fit[,1], fit[,2], mode, od_control)
      new_dose = as.numeric(recommend$next_dose); dose = recommend$dose_tab

      # Stop current trial if starting dose level is too toxic.
      if (new_dose == -1) {
        return(list(data = full_tab, new_dose = -1))
      }

      # Additional stopping rule.
      # Default: stopping = 1, the if statements below will not be visited. Stopping rule purely based on exhausting the sample size.
      # stopping = 2: Stop the trial if "max" patients are at the current dose and the next recommended dose is the same as current dose.
      if (stopping == 2 & length(full_tab$y) >= least) {
        if (as.numeric(allocation[new_dose]) >= 9) {
          break
        }
      }

      # Administer recommended dose to next cohort
      cohort = cohort + 1
      new_cohort = generate(cohort, size, new_dose, true_prob[new_dose], K,
                            cycle_prob, cycle_time, assessment, cumu_days = cumulate, dropout = dropout)
      full_tab = rbind(full_tab, new_cohort$info)
      cumulate = new_cohort$cumu_days
      full_tab = update(full_tab, cohort, cumulate, K, cycle_time, final = FALSE)

    }
  }

  # Update final cohort if assessment is not at the end of cycle
  if (K != assessment) {
    full_tab = update(full_tab, cohort, cumulate, K, cycle_time, final = TRUE)
  }

  if (sum(full_tab$y) > 0) {
    fit = run_stan(full_tab, 500, mcmc_samples, mcmc_chains, 1, 0.975)
    allocation = table(full_tab$dose)
    final_MTD = decision_rule(dose, full_tab$dose, full_tab[nrow(full_tab),]$dose, target, fit[,1], fit[,2], mode, od_control)
    if (allocation[as.numeric(final_MTD$next_dose)] < 9) {
      final_MTD$next_dose = full_tab[nrow(full_tab), 3]
    }
  } else {
    final_MTD = max(full_tab$dose)
  }

  ### Report Results
  res = results(full_tab, true, final_MTD, target)

  # Merge all results into one list
  final_res = res
  final_res$data = full_tab

  return(final_res)
}

##########################################################################
### Report the results of the simulation
# Parameters:
# tab: Generated data table.
# true: True dose probability relationship

results = function(tab, true, final_MTD, target) {

  low_dose = true[1, which(true[2,] < 0.16)]
  MTD_dose = true[1, which(true[2,] >= 0.16 & true[2,] < 0.33)]
  high_dose = true[1, which(true[2,] > 0.33)]

  # Basic summary statistics of results
  num_patients = nrow(tab)
  num_dose = table(factor(tab$dose, levels = 1:ncol(true)))
  prop_dose = prop.table(num_dose)
  num_dlt_dose = table(factor(tab$dose, levels = 1:ncol(true)), tab$y)[,2]
  prop_dlt_dose = ifelse(num_dose > 0, num_dlt_dose / num_dose, 0)
  dose_tab = rbind(num_dose, prop_dose, num_dlt_dose, prop_dlt_dose)

  # Duration of trial
  duration = tail(tab$cumu_arrival, n=1) + tail(tab$obs_time, n=1)

  # Number of total DLT's
  total_dlt = sum(tab$y); prop_dlt = total_dlt / num_patients

  # Proportion of patients treated at different dose category(low, MTD, toxic)
  prop_low = sum(prop.table(num_dose)[as.numeric(names(num_dose)) %in% low_dose])
  prop_MTD = sum(prop.table(num_dose)[as.numeric(names(num_dose)) %in% MTD_dose])
  prop_high = sum(prop.table(num_dose)[as.numeric(names(num_dose)) %in% high_dose])

  final_MTD = ifelse(sum(tab$y >0), as.numeric(final_MTD$next_dose), final_MTD)

  info = cbind(num_patients, total_dlt, prop_dlt, duration, prop_low, prop_MTD, prop_high, final_MTD)

  return(list(table = dose_tab, info = info))
}

###############################################################################################
### Generate cohort information for simulation
# Parameters for the function
# cohort: Current cohort number
# size: Size of cohort at a dose level
# x: dose level
# prob: true dose skeleton
# K: number of considered cycles
# cycle_prob: true probability of DLT happening in any of the cycles
# cycle_time: time of observation for a full cycle
# assessment: assessment period by clinicians, set as the first cycle after DLT is observed, then similarly for future cycles
# fixed: control the time of arrival between patients. If fixed is true then the time between arrival is set as constant
# for all patients. If fixed=FALSE then time between arrival is randomly generated.

generate = function(cohort, size, x, prob, K, cycle_prob, cycle_time, assessment, fixed = TRUE, cumu_days, dropout) {

  # Generate True Info #
  # 1. If patient has DLT.
  # 2. For patients that have DLTs, when did the DLT happen.
  # 3. Given patient's cycle information of DLT, the time to toxicity within the cycle
  # 4. Time to dropout using an exponential hazard function
  # 5. Time between arrival of patients

  cohort = rep(cohort, size)
  y = rbinom(size, 1, prob) # True y, DLT or no DLT after administering dose
  m = sum(y) # Number of DLTs given this dose
  cycle = ifelse(y == 1, sample(1:K, m, replace = TRUE, prob = cycle_prob), -1) # True cycle when DLT happens
  ttt = ifelse(y == 1, sample(1:cycle_time, m, replace = TRUE), -1) # True time-to-toxicity
  dose = rep(x, size) # Dose level for current cohort
  total = cycle_time * K
  rate = 0.15

  # Time of dropout of patients using a hazard function
  # Adds a corresponding dropout indicator if the dropout time is less than the total time (dropout is true)
  if (dropout) {
    lambda = log(1 - rate) / -total
    drop_time = round(rexp(size, lambda), 3)
    drop_time = ifelse(drop_time > total, total, drop_time)
  } else {
    drop_time = total
  }
  drop = ifelse(drop_time == total, 0, 1)


  # Time between arrivals of patients
  if (fixed) {
    tba = rep(0, size)
  } else {
    tba = c(0, sample((cycle_time/2):cycle_time, size-1, replace = TRUE))
  }

  # Time of patient entering study
  time_arrival = numeric(size); time_arrival[1] = 0
  for (i in 2:size) {
    time_arrival[i] = time_arrival[i-1] + tba[i]
  }

  # Calculate observed time during clinical assessement
  time_as = time_arrival[size] + assessment * cycle_time
  obs_time = time_as - time_arrival
  obs_time[which(obs_time > total)] = total

  # Set observed time as dropout out time
  obs_time = ifelse(drop_time < obs_time, drop_time, obs_time)

  # Calculate time of toxicity
  toxic_time = ifelse(y == 1, time_arrival + (cycle - 1) * cycle_time + ttt, -1)

  # cumu Days of arrival
  if (cumu_days == 0) {
    cumu_arrival = time_arrival
  } else {
    cumu_arrival = time_arrival + cumu_days
  }

  # Cumulative counter to keep track of total days
  cumu_days = cumu_days + time_arrival[size] + cycle_time * assessment

  # Indicator if DLT has been observed at time of assessment
  obs_DLT = (y == 1) & ((toxic_time - time_arrival) < obs_time)

  # Cycle of patient at assessment time
  cyc_assess = ceiling(obs_time / cycle_time)

  # Calculate weight using the weight function
  weight_info = find_weight(size = size, obs_DLT = obs_DLT, obs_time = obs_time, cycle_time = cycle_time,
                       K = K, cycle = cycle, y = y, toxic_time = toxic_time, time_arrival = time_arrival,
                       cyc_assess = cyc_assess)

  obs_y = ifelse(weight_info[,2] == 1 & y == 1, 1, 0)

  # Complete information matrix
  mat = data.frame(cbind(cohort, y, dose, cycle, time_arrival, cumu_arrival, ttt, toxic_time,
                         obs_time, cyc_assess, weight_info, drop_time, drop, obs_y))

  # Output: A listed format of the data to be pass into Stan for MCMC computation.
  return(list(info = mat, cumu_days = cumu_days))
}

###########################

# Function to update observed time and weight after inclusion of new cohorts

update = function(mat, cohort, cumulate, K, cycle_time, final) {

  # Update all cohorts excluding the lastest one
  if (!final) {
    # Update observed time
    prev = mat[which(mat$cohort != cohort),]  # Extract table of previous cohort
    prev$obs_time = ifelse(prev$drop == 1, prev$drop_time ,cumulate - prev$cumu_arrival)
    prev[which(prev$obs_time > (K * cycle_time)),9] = cycle_time * K

    # Update cycle at assessment and corresponding weight
    prev$cyc_assess = ceiling(prev$obs_time / cycle_time)
    temp = rbind(prev, mat[which(mat$cohort == cohort),])
    weight_info = find_weight(size = nrow(temp), temp$obs_y == 1, obs_time = temp$obs_time, cycle_time = cycle_time, K = K,
                              cycle = temp$cycle, y = temp$y, toxic_time = temp$toxic_time, time_arrival = temp$time_arrival,
                              cyc_assess = temp$cyc_assess)
    temp$obs_y = ifelse(temp$toxic_time > temp$obs_time + temp$time_arrival, temp$obs_y, temp$y)
    temp[,11:12] = weight_info

    # Merge back together
    new_mat = temp

  # Update only the last cohort
  } else {
    last = mat[which(mat$cyc_assess != K),]
    last$obs_time = last$drop_time
    last$cyc_assess = ceiling(last$obs_time / cycle_time)
    weight_info = find_weight(size = nrow(last), last$obs_y == 1, obs_time = last$obs_time, cycle_time = cycle_time, K = K,
                              cycle = last$cycle, y = last$y, toxic_time = last$toxic_time, time_arrival = last$time_arrival,
                              cyc_assess = last$cyc_assess)
    last$obs_y = ifelse(last$toxic_time > last$obs_time + last$time_arrival, last$obs_y, last$y)
    last[,11:12] = weight_info

    # Merge back together
    new_mat = rbind(mat[which(mat$cyc_assess == K),], last)
  }

  return(new_mat)
}

#############################

# The weight function

find_weight = function(size, obs_DLT, obs_time, cycle_time, K, cycle, y, toxic_time, time_arrival, cyc_assess) {

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

  # Combine P_hat and G(t) to calculate the final adaptive cycle weight of the observed data
  weight = numeric(size)
  for (i in 1:size) {

    # Weight is equal to 1 if DLT is observed during assessment window, else the weight is calculated linearly
    if (y[i] == 1 & toxic_time[i] < (time_arrival[i] + obs_time[i])) {
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
#'

run_stan = function(mat, warmup, iter, chains, samp_cores = 1, cont) {

  # Save into rstan required format
  stan_dat = list(N = length(mat$y), y = mat$obs_y, dose = mat$dose, W = mat$weight)

  sink(tempfile())
  # Run the MCMC
  fit = sampling(mod, data = stan_dat, warmup = warmup, iter = iter, chains = chains,
                 cores = samp_cores, verbose = FALSE, refresh = -1, control = list(adapt_delta = cont))
  sink(NULL)

  # Extract posterior samples
  extract_fit = extract(fit)
  alpha = extract_fit$alpha; beta = extract_fit$beta

  return(data.matrix(cbind(alpha, beta)))

}

###########################################################################################
### Part 3: Decision Rule for assigning dose level for the next patient

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
    return(list(info = info, next_dose = current))
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
    return(list(info = info, next_dose = current))
  }

  cand = as.numeric(safe_dose[1, which.min(safe_dose[3,])])
  level_diff = cand - current
  if (level_diff > 1 & info[2, which(colnames(info) == cand)] == 0) {
    cand = as.numeric(colnames(info)[which(colnames(info) == current) + 1])
  }

  return(list(info = info, next_dose = cand))
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
    return(list(info = info, next_dose = current))
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
      return(list(info = info, next_dose = current))
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

###############################################################################
# Function to calculate the results of the simulated trials.

mean.res = function(list, true, target_prob) {

  if (length(target_prob) == 2) {
    target = true[1,which(true[2,] < max(target_prob) & true[2,] > min(target_prob))]
  } else {
    target = true[1,which(true[2,] < 0.33 & true[2,] > 0.2)]
  }

  # Dummy variable to store list
  temp = list

  # Computation time of trial
  time = temp[[length(temp)]] / (length(temp) - 1)

  # Number of no selection
  no.select = length(Filter(function(x) {length(x) == 2}, temp)) / (length(temp) - 1)

  # All simulation trials that did not end pre-maturely
  full = Filter(function(x) {length(x) == 3 | length(x) == 1}, temp)

  allo = full[[1]]$table
  for (i in 2:(length(full)-1)) {
    allo = allo + full[[i]]$table
  }

  sum = full[[1]]$info
  for (i in 2:(length(full)-1)) {
    sum = sum + full[[i]]$info
  }
  prob_right_MTD = 0
  prob_over_MTD = 0
  prob_under_MTD = 0
  for (i in 1:(length(full)-1)) {
    if (full[[i]]$info[8] < target) {
      prob_under_MTD = prob_under_MTD + 1
    } else if (full[[i]]$info[8] > target) {
      prob_over_MTD = prob_over_MTD + 1
    } else {
      prob_right_MTD = prob_right_MTD + 1
    }
  }
  sum = cbind(sum, prob_right_MTD, prob_over_MTD, prob_under_MTD)
  mean_allo = allo / (length(full) - 1)
  mean_info = sum / (length(full)-1)
  prob_true_mtd = mean_info[,9] - no.select
  mean_info = cbind(mean_info, no.select, prob_true_mtd, time)
  rownames(mean_info) = NULL
  info = mean_info[,c(13, 10, 11, 12, 8, 5, 6, 7, 2, 3, 1, 4, 14)]
  names(info) = c("P(Choose True MTD)", "P(Choose Over MTD)", "P(Choose Under MTD)", "P(No Selection)",
                  "Mean Final Dose", "Allo Under MTD", "Allo True MTD","Allo Over MTD", "Total DLT",
                  "Prop DLT", "Patients Enrolled" ,"Trial Duration", "Computation Time")
  rownames(mean_allo) = c("Num Patients Treated", "Prop Patients Treated", "Number of DLTs", "Prop of DLTs")

  return(list(allocation = mean_allo, summary = info))
}
