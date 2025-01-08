# Case study on scheduling MRI facilities 

This case study aims to evaluate the efficiency of MRI scheduling policies at a hospital with scarce MRI resources. 
Two patient types with distinct scan durations and arrival patterns currently use separate MRI facilities, but management is considering combining these facilities to optimize utilization and reduce waiting times. 

## Part 1 - Statistical Analysis of MRI scan data 

To aid this decision, the study will first estimate the mean and variance metrics related to scan durations and inter-arrival rates of both patient types. Subsequently, parametric and non-parametric bootstrap methods are applied to measure the uncertainty of the estimates. Finally, Monte Carlo simulations are performed to validate the model specifications. 

## Part 2 - Optimizing Scheduling Policies 

A discrete event simulation is used to model scheduling under the current and proposed new system. The study evaluates the performance indicators  overtime, patient waiting time, system utilization and total number of patients served and considers the uncertainty around the estimates from the statistical analysis to offer an advice on optimal time slot lengths that aim to minimize overtime while meeting patient needs.
