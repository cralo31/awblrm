# Adaptive Time-to-event Weight Bayesian Logistic Regression Model Incorporating Cycle Information for Phase 1 Dose Escalation Studies

This package is for the talk on Bayesian Dose-Finding Model with Adaptive Time-to-event Weight
Incorporating Cycle Information for Immuno-oncology Compounds. The package provides a flexible framework for phase 1 dose-escalation trial designs, including patient enrollment, treatment cycles, toxicity assessment window, decision rule, stopping criteria, and escalation with overdose control. 

To install this on R/Rstudio, simply install directly from github with the _devtools_ package.
```{r}
library(devtools)
install_github("cralo31/awblrm")
```
There are 2 main functions. One for recommending the next dose according to the current trial information, and the other for repeated random trial simulations for hypothetical scenarios.

1. Recommend the next dose according to the current patient's information. 
```{r}
?decision.dose()
```

2. Set up a simulated trial with prespecifed number of cores to compute in parallel. There are many parameters available to tune depending on the nature of the study and patient information. 
```{r}
# For running trials in parallel
?parallel.trial() 

# For running a single simulated trial
?trial()
```

Run the following example. Set up a 3-cycle treatment design with 5 different dose levels. The are currently 9 patients (3 cohorts of 3) observed at the assessment time window. The function outputs the list with two elements. The first element is the next recommended dose. The second element is a summary information
#'matrix of the implemented decision rule.
```{r}
# Consider a 3-cycle treatment scenario with observed information of 9 patients #

# Current trial information
all_doses = c(1, 2, 3, 4, 5)
max_cycle = 3
cycle_time = 30
mode = 1
target = 0.25
od_control = 1
cohort = c(1, 1, 1, 2, 2, 2, 3, 3, 3)
dlt = c(0, 0, 0, 1, 0, 0, 1, 0, 1)
dose = c(1, 1, 1, 2, 2, 2, 2, 2, 2)
dlt_cycle = c(-1, -1, -1, 2, -1, -1, 3, -1, 2)
obs_time = c(90, 45, 90, 35, 90, 75, 80, 90, 55)

# Output the next recommended dose
test = decision.dose(all_doses, max_cycle, cycle_time, mode, target, od_control, cohort, dlt,
dose, dlt_cycle, obs_time)
test
```
