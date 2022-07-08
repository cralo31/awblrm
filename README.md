# awblrm
Dose-finding package

This package is for the Bayesian Dose-Finding Model with Adaptive Time-to-event Weight
Incorporating Cycle Information for Immuno-oncology Compounds. 

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

*true Matrix; Pre-defined true dose-toxicity relationship. The input should be a 2-by-p matrix, with the 1st row being the
*levels and the 2nd row being the corresponding probability of DLT.
*n Integer; Number of trials to be simulated. Minimum of 1.
*cores Integer; Number of cores to utilize for parallel computing. Defaults to 1.
*K Integer; Number of treatment cycles considered.
*size Integer; Number of patients to enrolled in each cohort.
*cycle_prob K-simplex; A vector that contains the true probability of DLT happening in each of the cycles. Sums up to one.
*cycle_time Integer; Duration of a single treatment cycle (NOT entire treatment window!). Units are arbitrary.
*target The defined target probability of toxicity. For point-estimate-based decision rules, please input a single numeric value. For interval-estimate-based decision *rules, please input a vector indicating the target toxicity range, e.g. c(0.2, 0.33).
*samples Integer; The maximum number of patients to be enrolled into the trial. Also serves as the default stopping rule.
*assessment Integer; The cycle in which the model is assessed to determine the dose level of the next cohort (1 means the model is assessed
*after the completion of the 1st cycle of each cohort, 2 means 2nd cycle...).
*fixed Boolean; TRUE if patients in the same cohort are enrolled at the same time. FALSE if the enrollment time of patients in the same
*cohort follows a Uniform (0, cycle_time).
*mode Integer; Determine which decision rule to implement to perform dose escalation / de-escalation. 1 = Original CRM, 2 = Inverse Point-Estimate, 3 = BLRM, 4 = Ratio. Defaults to 1. For details on the decision rule please refer to the report.
stopping Integer; Defined the stopping criteria. 1 = stop the trial when maximum enrollment is achieved, 2 = stop the trial if at least 12
patients are enrolled and 9 patients have already been treated at the next recommended dose. Defaults to 1.
least Integer; Least amount of patients to enrolled if using the second stopping rule. Defaults to 12, directly impacts "stopping".
dropout Boolean; TRUE if patient dropout time is assumened to follow an exponential distribution with a cumulative probability of 15% at the end of the treatment window.
od_control Double (0 - 1); Escalation with overdose control. Escalate to the next dose if the probability of overdosing does not exceed the input probability. No overdose control if od_control is set to 1.
mcmc_chains Integer; Number of chains to consider during MCMC sampling. Defaults to 3.
mcmc_samples Integer; Numbers of samples from the posterior distribution of the betas. Default to 3000 with a burn-in period of 500 samples.


To run an example. 
```
