# Adaptive Time-to-event Weight Bayesian Logistic Regression Model Incorporating Cycle Information for Dose Finding 

This package is for the talk on Bayesian Dose-Finding Model with Adaptive Time-to-event Weight
Incorporating Cycle Information for Immuno-oncology Compounds. The package provides a flexible framework for phase 1 dose-escalation trial designs, including patient enrollment, treatment cycles, toxicity assessment window, decision rule, stopping criteria, and escalation with overdose control. 

To install this on R/Rstudio, simply install directly from github with the _devtools_ package.
```{r}
library(devtools)
install_github("cralo31/awblrm")
```
There are 2 main functions. 

1. Recommend the next dose according to the current patient's information. 
```{r}
?next.dose()
```

2. Set up a simulated trial with prespecifed number of cores to compute in parallel. There are many parameters available to tune depending on the nature of the study and patient information. 
```{r}
?simulation()
```

Run the following example. 
```{r}
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
```
