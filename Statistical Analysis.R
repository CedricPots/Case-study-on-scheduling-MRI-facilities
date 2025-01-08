
#===============================================================================
#Statistical Analysis Part - Computational Research Skills 
#===============================================================================

#Course code: EBS4043
#Authors: Jorn Meessen, Pauline Hanssen, Thijn de Vries, Cedric Pots. 
#Topic: A case study on scheduling MRI facilities 

#===============================================================================
#DATA EXPLORATION 
#===============================================================================

library(tseries)
library(goft)
library(stats)
library(ggplot2)
library(MASS)

#Two types of patients (Type 1 and Type 2)
type1.data <- scan.records.data[scan.records.data$PatientType == "Type 1", ]
type2.data <- scan.records.data[scan.records.data$PatientType == "Type 2", ]

#-------------------------------------------------------------------------------
#TYPE 1 PATIENTS 
#-------------------------------------------------------------------------------

mean(type1.data$Duration)*60 #=>26 min
sd(type1.data$Duration)*60 #=>6 min

#-------------------------------------------------------------------------------
#Informal test: histogram against normal distribution 
#-------------------------------------------------------------------------------

hist(type1.data$Duration, , col="#104E8B", prob=TRUE, 
     main="MRI Scan Duration of Type 1 Patients", xlab="MRI Scan Duration", 
     breaks=20)
x.axis.1 <- seq(min(type1.data$Duration), max(type1.data$Duration), 0.01)
lines(x.axis.1, dnorm(x.axis.1, mean=mean(type1.data$Duration), sd=sd(type1.data$Duration)), col="#CC0000", 
      lwd=2.4)
legend("topright", legend="Fitted Normal Distribution", col="", lwd=2.4, cex=0.8)

#-------------------------------------------------------------------------------
#Formal test: test for normality (Jarque-Bera test)
#-------------------------------------------------------------------------------

jarque.bera.test(type1.data$Duration) #p-value: 0.388

#Conclusion, fail to reject H0 that data follows a normal distribution. 
#Hence, we accept that data type 1 may follow a normal distribution. 

#-------------------------------------------------------------------------------
#Informal test: #patients per day Poisson distributed. 
#-------------------------------------------------------------------------------
#If #patients per day is Poisson distributed, then the inter-arrival times a are exponentially distributed.

# Opening hours: 8:00 - 17:00
n1 <- nrow(type1.data)
X.inter.arrival <- rep(NA, times = n1 - 1)  # Initialize vector for inter-arrival times

# Compute inter-arrival times
for (i in 1:(n1 - 1)) {
  if (type1.data$Date[i] == type1.data$Date[i + 1]) {
    # Same day arrival
    X.inter.arrival[i] <- type1.data$Time[i + 1] - type1.data$Time[i]
  } else {
    # Different days (crossing the closing time)
    X.inter.arrival[i] <- (17 - type1.data$Time[i]) + (type1.data$Time[i + 1] - 8)
  }
}

# Check the histogram with a fitted exponential distribution
hist.arrival.t1 <- ggplot(data.frame(X.inter.arrival), aes(x = X.inter.arrival)) +
  geom_histogram(bins = 12, fill = "#104E8B", color = "black", alpha = 0.7, aes(y = ..density..)) +
  stat_function(fun = dexp, args = list(rate = 1 / mean(X.inter.arrival)), 
                aes(color = "Fitted Exponential Distribution"), size = 1.5) +  # Add legend label
  labs(title = "Inter-Arrival Time - Type 1 Patients",
       x = "Inter-Arrival Time (hours)",
       y = "Density") +
  scale_color_manual(name = "", values = c("Fitted Exponential Distribution" = "#FF3030")) +  # Define legend color
  theme_bw() +  # White background
  theme(legend.position = "top")  # Position the legend at the top
print(hist.arrival.t1)

#-------------------------------------------------------------------------------
#Formal test: #patients per day Poisson distributed (Kolmogorov-Smirnov K-S test) 
#-------------------------------------------------------------------------------

avg.inter.arrival <- mean(X.inter.arrival)
ks.test.result <- ks.test(X.inter.arrival, "pexp", rate = 1 / avg.inter.arrival)
print(ks.test.result) #p-value: 0.55 
#Conclusion: fail to reject H0 that inter-arrival times follow an exponential distribution

#-------------------------------------------------------------------------------
#TYPE 2 PATIENTS 
#-------------------------------------------------------------------------------

mean(type2.data$Duration)*60 #=>40 min
sd(type2.data$Duration)*60 #=>11 min

#-------------------------------------------------------------------------------
#Informal test: histogram against normal distribution 
#-------------------------------------------------------------------------------

hist(type2.data$Duration, col="#104E8B", prob=TRUE, 
     main="MRI Scan Duration of Type 2 Patients", xlab="MRI Scan Duration", 
     breaks=20)
x.axis.2 <- seq(min(type2.data$Duration), max(type2.data$Duration), 0.01)
lines(x.axis.2, dnorm(x.axis.2, mean=mean(type2.data$Duration), sd=sd(type2.data$Duration)), col="#FF3030", 
      lwd=2.4)
#Data seems to look more like a gamma distribution. 

#-------------------------------------------------------------------------------
#Formal test: test for normality (Jarque-Bera test)
#-------------------------------------------------------------------------------

jarque.bera.test(type2.data$Duration) #p-value: 0.027

#Conclusion: reject H0 that data is normally distributed.


#-------------------------------------------------------------------------------
#Informal test 1: histogram against gamma distribution 
#-------------------------------------------------------------------------------

beta.gamma <- sd(type2.data$Duration)^2/mean(type2.data$Duration)
alpha.gamma <- mean(type2.data$Duration)/beta.gamma
lines(x.axis.2, dgamma(x.axis.2, shape=alpha.gamma, rate=1/beta.gamma), type="l",
      col="#EEC900", lwd=2)
legend("topright", legend=c("Fitted Normal Distribution", "Fitted Gamma Distribution"), 
       col=c("#FF3030", "#EEC900"), lwd=2, cex=0.8)

#-------------------------------------------------------------------------------
#Informal test 2: empirical distribution function (EDF) gamma vs. normal  
#-------------------------------------------------------------------------------
plot(ecdf(type2.data$Duration), col="deepskyblue4", lwd=2, main="EDF of Type 2 Patients", ylab="(Empirical) CDF")
lines(x.axis.2, pnorm(x.axis.2, mean=mean(type2.data$Duration), sd=sd(type2.data$Duration)), col="#CC0000", lwd=2)
lines(x.axis.2, pgamma(x.axis.2, shape=mean(type2.data$Duration) / (sd(type2.data$Duration)^2 / mean(type2.data$Duration)), rate=1/(sd(type2.data$Duration)^2 / mean(type2.data$Duration))), col="#FF9900", lwd=2)

#-------------------------------------------------------------------------------
#Formal test: test for gamma distribution (Villasenor & Gonzales-Estrada, 2015)
#-------------------------------------------------------------------------------

gamma_test(type2.data$Duration) #p-value 0.55 

#Conclusion: fail to reject H0 of gamma distribution. 

#-------------------------------------------------------------------------------
#Informal test: #patients per day normally distributed. 
#-------------------------------------------------------------------------------

n2 <- nrow(type2.data)  # Number of Type 2 records
X.inter.arrival.t2 <- rep(NA, times = n2 - 1)

# Compute inter-arrival times
for (i in 1:(n2 - 1)) {
  if (type2.data$Date[i] == type2.data$Date[i + 1]) {
    # Same day arrival
    X.inter.arrival.t2[i] <- type2.data$Time[i + 1] - type2.data$Time[i]
  } else {
    # Different days (crossing the closing time)
    X.inter.arrival.t2[i] <- (17 - type2.data$Time[i]) + (type2.data$Time[i + 1] - 8)
  }
}

# Create the histogram with the fitted normal distribution
hist.arrival.t2.normal <- ggplot(data.frame(X.inter.arrival.t2), aes(x = X.inter.arrival.t2)) +
  geom_histogram(bins = 12, fill = "#104E8B", color = "black", alpha = 0.7, aes(y = ..density..)) +
  stat_function(fun = dnorm, 
                args = list(mean = mean(X.inter.arrival.t2, na.rm = TRUE), 
                            sd = sd(X.inter.arrival.t2, na.rm = TRUE)), 
                aes(color = "Fitted Normal Distribution"), linewidth = 1) +  # Add legend label for normal distribution
  labs(title = "Inter-Arrival Time - Type 2 Patients",
       x = "Inter-Arrival Time (hours)",
       y = "Density") +
  scale_color_manual(name = "", values = c("Fitted Normal Distribution" = "#FF3030")) +  # Define legend color
  theme_minimal() +  # Minimal theme
  theme(legend.position = "top")  # Position the legend at the top

print(hist.arrival.t2.normal)

#-------------------------------------------------------------------------------
#Formal test: #patients per day normal distributed (Jarque Bera test) 
#-------------------------------------------------------------------------------

print(jarque.bera.test(X.inter.arrival.t2))

#Conclusion: fail to reject H0 of normally distributed inter-arrival times.

#===============================================================================
#TYPE 1 PATIENTS 
#===============================================================================

#Note: only a parametric approach is conducted. 

#-------------------------------------------------------------------------------
#BOOTSTRAP TYPE 1 PATIENTS (PARAMETRIC APPROACH)
#-------------------------------------------------------------------------------

'''
Quantities of interest: 
1. Mean MRI (with confidence intervals calculated by bootstrap)
Motivation: to determine average utilization 
2. Standard deviation MRI (with confidence intervals calculated by bootstrap)
Motivation: buffer times. High standard deviation mwould require more buffer time
between appointments to ensure flexibility and avoid delays. 
3. Quantiles of crossing threshold 
Motivation: handle extreme cases to ensure that schedules can accommodate 
unusually long scans without causing excessive waiting times. 
'''
set.seed(515)

B1 <- 10000 #number of bootstrap replications 
alpha.t1 <- 0.05 #significance level 

n1 <- length(type1.data$Duration)
X.bar.t1 <- mean(type1.data$Duration)      # Mean of Type 1 durations
St.Dev.t1 <- sd(type1.data$Duration)        # Standard deviation of Type 1 durations

# Initialize vectors to store bootstrap statistics
X.star.bar.t1 <- rep(NA, B1)                # Bootstrap means
X.star.sd.t1 <- rep(NA, B1)                 # Bootstrap standard deviations
Q.star.t1 <- rep(NA, B1)                    # T-statistics

#-------------------------------------------------------------------------------
# Q1: Bootstrap CI for Mean
#-------------------------------------------------------------------------------
# Motivation: to determine average utilization

for(b in 1:B1) {
  J1 <- sample.int(n1, n1, replace = TRUE)  # Resample with replacement
  X.star.t1 <- type1.data$Duration[J1]      # Bootstrap sample
  X.star.bar.t1[b] <- mean(X.star.t1)       # Mean of bootstrap sample
  X.star.sd.t1[b] <- sd(X.star.t1)          # SD of bootstrap sample
  Q.star.t1[b] <- sqrt(n1) * (X.star.bar.t1[b] - X.bar.t1) / X.star.sd.t1[b]  # T-statistic
}

# Main bootstrap results for mean duration
t1.avg.duration <- mean(X.star.bar.t1)

# Critical values from the bootstrap distribution
cv.lower1 <- quantile(Q.star.t1, probs = alpha.t1 / 2)
cv.upper1 <- quantile(Q.star.t1, probs = 1 - alpha.t1 / 2)

# 95% CI for the mean duration
CI.dur.t1.lower <- X.bar.t1 - cv.upper1 * St.Dev.t1 / sqrt(n1)
CI.dur.t1.upper <- X.bar.t1 - cv.lower1 * St.Dev.t1 / sqrt(n1)

cat("Bootstrap Mean MRI Duration (Type 1):", t1.avg.duration, "which is", t1.avg.duration * 60, "minutes\n")
cat("Bootstrap 95% CI for Mean Duration:", CI.dur.t1.lower, "-", CI.dur.t1.upper, "which is between", CI.dur.t1.lower * 60, "-", CI.dur.t1.upper * 60, "minutes\n")

# Plot histogram of bootstrap means
hist(X.star.bar.t1, main="Bootstrap Distribution of MRI Mean Duration (Type 1)",
     xlab="Mean MRI Duration", col="skyblue", breaks=30)
abline(v=c(CI.dur.t1.lower, CI.dur.t1.upper), col="red", lwd=2, lty=2)
legend("topright", legend=c("95% CI"), col="red", lty=2, lwd=2)

#-------------------------------------------------------------------------------
# Q2: Bootstrap CI for Standard Deviation 
#-------------------------------------------------------------------------------
# Motivation: To plan buffer times between appointments
t1.sd.duration <- mean(X.star.sd.t1)  # Calculate the mean of bootstrap standard deviations

# Bootstrap confidence intervals for the standard deviation
CI.sd.dur.t1.lower <- quantile(X.star.sd.t1, probs = alpha.t1 / 2)
CI.sd.dur.t1.upper <- quantile(X.star.sd.t1, probs = 1 - alpha.t1 / 2)

cat("Bootstrap SD MRI Duration (Type 1):", t1.sd.duration, "which is", t1.sd.duration * 60, "minutes\n")
cat("Bootstrap 95% CI for SD Duration:", CI.sd.dur.t1.lower, "-", CI.sd.dur.t1.upper, "which is between", CI.sd.dur.t1.lower * 60, "-", CI.sd.dur.t1.upper * 60, "minutes\n")

# Plot histogram of bootstrap SDs
hist(X.star.sd.t1, main="Bootstrap Distribution of MRI Standard Deviation (Type 1)",
     xlab="Standard Deviation of MRI Duration", col="skyblue", breaks=30)
abline(v=c(CI.sd.dur.t1.lower, CI.sd.dur.t1.upper), col="red", lwd=2, lty=2)
legend("topright", legend=c("95% CI"), col="red", lty=2, lwd=2)

#-------------------------------------------------------------------------------
#Q3: Quantiles with Confidence Intervals
#-------------------------------------------------------------------------------
# Motivation: To handle extreme cases and avoid excessive waiting times. 
# Quantiles indicate the range within which specific percentiles of the data are likely to fall. 

# Vector of quantiles of interest for Type 1
quantiles.interest.t1 <- c(0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975)
observed.quantiles.t1 <- quantile(type1.data$Duration, quantiles.interest.t1)

# Matrix of differences between bootstrap quantiles and observed quantiles
bootstrap.differences.t1 <- replicate(B1, 
                                      quantile(type1.data$Duration[sample(n1, replace = TRUE)], probs = quantiles.interest.t1) - observed.quantiles.t1)

# Confidence intervals for quantiles
ci.upper.quantiles.t1 <- apply(bootstrap.differences.t1, 1, quantile, 0.975)
ci.lower.quantiles.t1 <- apply(bootstrap.differences.t1, 1, quantile, 0.025)

# Confidence intervals for each quantile
ci.quantile.intervals.t1 <- cbind(observed.quantiles.t1 - ci.upper.quantiles.t1, 
                                  observed.quantiles.t1 - ci.lower.quantiles.t1)

# Print quantile confidence intervals (in minutes)
print(ci.quantile.intervals.t1 * 60)
#Interpretation: 95% of cases are likely to be below 35.8 minutes. 

#-------------------------------------------------------------------------------
#TYPE 1: SIMULATION STUDY 
#-------------------------------------------------------------------------------

# Parameters for normal distribution (duration)
mean.duration.t1 <- mean(type1.data$Duration)
sd.duration.t1 <- sd(type1.data$Duration)

# Parameters for exponential distribution (inter-arrival times)
lambda.t1 <- 1 / mean(X.inter.arrival, na.rm = TRUE) # lambda is the reciprocal of the mean inter-arrival time

bootstrap.t1 <- function(duration.t1, inter.arrivals.t1, B.t1, alpha.t1) { 
  n1.t1 <- length(duration.t1)
  n2.t1 <- length(inter.arrivals.t1)
  
  mean.duration.t1 <- mean(duration.t1)
  sd.duration.t1 <- sd(duration.t1)
  
  mean.inter.arrivals.t1 <- mean(inter.arrivals.t1)
  sd.inter.arrivals.t1 <- sd(inter.arrivals.t1)
  
  # Storing the statistics in vectors:
  # T-statistics for CI surrounding the mean:
  T.stat.dur.t1 <- rep(NA, B.t1)
  T.stat.time.t1 <- rep(NA, B.t1)
  
  # Statistics for the CI surrounding the variances:
  Q.stat.dur.t1 <- rep(NA, B.t1)
  Q.stat.time.t1 <- rep(NA, B.t1)
  
  # Bootstrap procedure:
  for (b in 1:B.t1) {
    # Resampling the durations with replacement
    J1.t1 <- sample(1:n1.t1, replace = TRUE)  
    duration.star.t1 <- duration.t1[J1.t1]  # Bootstrap sample for durations
    mean.duration.star.t1 <- mean(duration.star.t1)
    sd.duration.star.t1 <- sd(duration.star.t1)
    
    # Calculate T-statistic for duration
    T.stat.dur.t1[b] <- sqrt(n1.t1) * (mean.duration.star.t1 - mean.duration.t1) / sd.duration.star.t1
    # Calculate Q-statistic for duration variance
    Q.stat.dur.t1[b] <- (n1.t1 - 1) * sd.duration.star.t1^2 / (sd.duration.t1^2)
    
    # Resampling the inter-arrival times with replacement
    J2.t1 <- sample(1:n2.t1, replace = TRUE)  
    inter.arrivals.star.t1 <- inter.arrivals.t1[J2.t1]  # Bootstrap sample for inter-arrivals
    mean.inter.arrivals.star.t1 <- mean(inter.arrivals.star.t1)
    sd.inter.arrivals.star.t1 <- sd(inter.arrivals.star.t1)
    
    # Calculate T-statistic for inter-arrival times
    T.stat.time.t1[b] <- sqrt(n2.t1) * (mean.inter.arrivals.star.t1 - mean.inter.arrivals.t1) / sd.inter.arrivals.star.t1
    # Calculate Q-statistic for inter-arrival variance
    Q.stat.time.t1[b] <- (n2.t1 - 1) * sd.inter.arrivals.star.t1^2 / (sd.inter.arrivals.t1^2)
  }
  
  # Critical values
  cv.T.dur.t1 <- quantile(T.stat.dur.t1, probs = c(alpha.t1 / 2, 1 - alpha.t1 / 2))
  cv.T.time.t1 <- quantile(T.stat.time.t1, probs = c(alpha.t1 / 2, 1 - alpha.t1 / 2))
  
  cv.Q.dur.t1 <- quantile(Q.stat.dur.t1, probs = c(alpha.t1 / 2, 1 - alpha.t1 / 2))
  cv.Q.time.t1 <- quantile(Q.stat.time.t1, probs = c(alpha.t1 / 2, 1 - alpha.t1 / 2))
  
  # Two-sided 1-alpha-confidence-intervals
  CIs <- matrix(NA, nrow = 4, ncol = 2, 
                dimnames = list(c("CI Mean Duration", "CI Variance Duration",
                                  "CI Mean Inter-Arrival Time", "CI Variance Inter-Arrival Time"), 
                                c("Lowerbound", "Upperbound")))
  
  # Confidence intervals for duration mean and variance
  CIs[1, 1] <- mean.duration.t1 - cv.T.dur.t1[2] * sd.duration.t1 / sqrt(n1.t1)
  CIs[1, 2] <- mean.duration.t1 - cv.T.dur.t1[1] * sd.duration.t1 / sqrt(n1.t1)  # Corrected to add
  CIs[2, 1] <- (n1.t1 - 1) * sd.duration.t1^2 / cv.Q.dur.t1[2]
  CIs[2, 2] <- (n1.t1 - 1) * sd.duration.t1^2 / cv.Q.dur.t1[1]
  
  # Confidence intervals for inter-arrival mean and variance
  CIs[3, 1] <- mean.inter.arrivals.t1 - cv.T.time.t1[2] * sd.inter.arrivals.t1 / sqrt(n2.t1)
  CIs[3, 2] <- mean.inter.arrivals.t1 - cv.T.time.t1[1] * sd.inter.arrivals.t1 / sqrt(n2.t1)  # Corrected to add
  CIs[4, 1] <- (n2.t1 - 1) * sd.inter.arrivals.t1^2 / cv.Q.time.t1[2]
  CIs[4, 2] <- (n2.t1 - 1) * sd.inter.arrivals.t1^2 / cv.Q.time.t1[1]
  
  return(CIs)
}

#-------------------------------------------------------------------------------
# TYPE 1: Monte Carlo Simulation (parametric approach)
#-------------------------------------------------------------------------------
no.sims.t1 <- 1000                    # Number of simulations
B.t1 <- 1000                          # Number of bootstrap draws
alpha.t1 <- 0.05                      # Rejection level
n.t1 <- nrow(type1.data)              # Sample size

# Initialize rejection probability matrix
rej.prob.matrix.t1 <- matrix(NA, nrow=0, ncol=5,
                             dimnames = list(NULL, c("Sample Size", 
                                                     "Rej.Prob.Mean.Duration", 
                                                     "Rej.Prob.Var.Duration",
                                                     "Rej.Prob.Mean.Inter-Arrival.Time", 
                                                     "Rej.Prob.Var.Inter-Arrival.Time")))

for (i.t1 in 1:length(n.t1)) {
  n.i.t1 <- n.t1[i.t1]  # Current sample size
  
  # Store rejection results
  Reject.mat.t1 <- matrix(NA, nrow=no.sims.t1, ncol=4, 
                          dimnames = list(NULL, c("Reject Mean Duration", 
                                                  "Reject Variance Duration",
                                                  "Reject Mean Inter-Arrival Time", 
                                                  "Reject Variance Inter-Arrival Time")))
  
  for (j.t1 in 1:no.sims.t1) {
    # Generate data
    X.dur.t1 <- rnorm(n.i.t1, mean.duration.t1, sd.duration.t1)  # Normal distribution for duration
    X.time.t1 <- rexp(n.i.t1 - 1, rate = lambda.t1)              # Exponential distribution for inter-arrival times
    
    # Get confidence intervals from bootstrap
    Conf.Ints.t1 <- bootstrap.t1(X.dur.t1, X.time.t1, B.t1, alpha.t1)
    
    # Record rejection results
    Reject.mat.t1[j.t1, 1] <- as.integer(mean(type1.data$Duration) < Conf.Ints.t1[1, 1]) + 
      as.integer(mean(type1.data$Duration) > Conf.Ints.t1[1, 2])
    Reject.mat.t1[j.t1, 2] <- as.integer(sd(type1.data$Duration)^2 < Conf.Ints.t1[2, 1]) + 
      as.integer(sd(type1.data$Duration)^2 > Conf.Ints.t1[2, 2])
    Reject.mat.t1[j.t1, 3] <- as.integer(mean(X.inter.arrival) < Conf.Ints.t1[3, 1]) + 
      as.integer(mean(X.inter.arrival) > Conf.Ints.t1[3, 2])
    Reject.mat.t1[j.t1, 4] <- as.integer(sd(X.inter.arrival)^2 < Conf.Ints.t1[4, 1]) + 
      as.integer(sd(X.inter.arrival)^2 > Conf.Ints.t1[4, 2])
  }
  
  # Calculate rejection probabilities
  rej.prob.t1 <- c(n.i.t1, colMeans(Reject.mat.t1))
  rej.prob.matrix.t1 <- rbind(rej.prob.matrix.t1, rej.prob.t1)
}

print(rej.prob.matrix.t1)

'''
Interpretation: 
All rejection probabilities are close to rejection level of 0.05. Hence, this builds 
confidence that the parametric model with duration being normally distributed
and inter-arrival times being exponentially distributed is well-specified. 
'''

print(Conf.Ints.t1*60)

#===============================================================================
#TYPE 2 PATIENTS 
#===============================================================================

# Two approaches will be studied: 
# 1. Parametric approach (from known gamma and normal distributions)
# 2. Non-parametric approach (estimation of empirical distribution function)

#------------------------------------------------------------------------------
#Quantiles 
#------------------------------------------------------------------------------
B2 <- 10000
# Vector of quantiles of interest for Type 2
quantiles.interest.t2 <- c(0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975)
observed.quantiles.t2 <- quantile(type2.data$Duration, quantiles.interest.t2)

# Matrix of differences between bootstrap quantiles and observed quantiles
bootstrap.differences.t2 <- replicate(B2, 
                                      quantile(type2.data$Duration[sample(n2, replace = TRUE)], probs = quantiles.interest.t2) - observed.quantiles.t2)

# Confidence intervals for quantiles
ci.upper.quantiles.t2 <- apply(bootstrap.differences.t2, 1, quantile, 0.975)
ci.lower.quantiles.t2 <- apply(bootstrap.differences.t2, 1, quantile, 0.025)

# Confidence intervals for each quantile
ci.quantile.intervals.t2 <- cbind(observed.quantiles.t2 - ci.upper.quantiles.t2, 
                                  observed.quantiles.t2 - ci.lower.quantiles.t2)

# Print quantile confidence intervals (in minutes)
print(ci.quantile.intervals.t2 * 60)

#===============================================================================
#PARAMETRIC APPROACH: 
#===============================================================================
#Key assumptions: 
#1. Duration is gamma distributed
#2. Inter-arrival times are normally distributed. 

# Parameters for gamma distribution (duration)
beta.gamma.t2 <- sd(type2.data$Duration)^2 / mean(type2.data$Duration)
alpha.gamma.t2 <- mean(type2.data$Duration) / beta.gamma.t2

# Parameters for normal distribution (inter-arrival times)
mean.normal.t2 <- mean(X.inter.arrival.t2, na.rm = TRUE) #inter-arrival times were already calculated
sd.normal.t2 <- sd(X.inter.arrival.t2, na.rm = TRUE)

bootstrap.t2 <- function(duration.t2, inter.arrivals.t2, B.t2, alpha.t2) { 
  n1.t2 <- length(duration.t2)
  n2.t2 <- length(inter.arrivals.t2)
  
  mean.duration.t2 <- mean(duration.t2)
  sd.duration.t2 <- sd(duration.t2)
  
  mean.inter.arrivals.t2 <- mean(inter.arrivals.t2)
  sd.inter.arrivals.t2 <- sd(inter.arrivals.t2)
  
  # Storing the statistics in vectors:
  # T-statistics for CI surrounding the mean:
  T.stat.dur.t2 <- rep(NA, B.t2)
  T.stat.time.t2 <- rep(NA, B.t2)
  
  # Statistics for the CI surrounding the variances:
  Q.stat.dur.t2 <- rep(NA, B.t2)
  Q.stat.time.t2 <- rep(NA, B.t2)
  
  # Bootstrap procedure:
  for (b in 1:B.t2) {
    # Resampling the durations with replacement
    J1.t2 <- sample(1:n1.t2, replace = TRUE)  
    duration.star.t2 <- duration.t2[J1.t2]  # Bootstrap sample for durations
    mean.duration.star.t2 <- mean(duration.star.t2)
    sd.duration.star.t2 <- sd(duration.star.t2)
    
    # Calculate T-statistic for duration
    T.stat.dur.t2[b] <- sqrt(n1.t2) * (mean.duration.star.t2 - mean.duration.t2) / sd.duration.star.t2
    # Calculate Q-statistic for duration variance
    Q.stat.dur.t2[b] <- (n1.t2 - 1) * sd.duration.star.t2^2 / (sd.duration.t2^2)
    
    # Resampling the inter-arrival times with replacement
    J2.t2 <- sample(1:n2.t2, replace = TRUE)  
    inter.arrivals.star.t2 <- inter.arrivals.t2[J2.t2]  # Bootstrap sample for inter-arrivals
    mean.inter.arrivals.star.t2 <- mean(inter.arrivals.star.t2)
    sd.inter.arrivals.star.t2 <- sd(inter.arrivals.star.t2)
    
    # Calculate T-statistic for inter-arrival times
    T.stat.time.t2[b] <- sqrt(n2.t2) * (mean.inter.arrivals.star.t2 - mean.inter.arrivals.t2) / sd.inter.arrivals.star.t2
    # Calculate Q-statistic for inter-arrival variance
    Q.stat.time.t2[b] <- (n2.t2 - 1) * sd.inter.arrivals.star.t2^2 / (sd.inter.arrivals.t2^2)
  }
  
  # Critical values
  cv.T.dur.t2 <- quantile(T.stat.dur.t2, probs = c(alpha.t2 / 2, 1 - alpha.t2 / 2))
  cv.T.time.t2 <- quantile(T.stat.time.t2, probs = c(alpha.t2 / 2, 1 - alpha.t2 / 2))
  
  cv.Q.dur.t2 <- quantile(Q.stat.dur.t2, probs = c(alpha.t2 / 2, 1 - alpha.t2 / 2))
  cv.Q.time.t2 <- quantile(Q.stat.time.t2, probs = c(alpha.t2 / 2, 1 - alpha.t2 / 2))
  
  # Two-sided 1-alpha-confidence-intervals
  CIs <- matrix(NA, nrow = 4, ncol = 2, 
                dimnames = list(c("CI Mean Duration", "CI Variance Duration",
                                  "CI Mean Inter-Arrival Time", "CI Variance Inter-Arrival Time"), 
                                c("Lowerbound", "Upperbound")))
  
  # Confidence intervals for duration mean and variance
  CIs[1, 1] <- mean.duration.t2 - cv.T.dur.t2[2] * sd.duration.t2 / sqrt(n1.t2)
  CIs[1, 2] <- mean.duration.t2 - cv.T.dur.t2[1] * sd.duration.t2 / sqrt(n1.t2)  # Corrected to add
  CIs[2, 1] <- (n1.t2 - 1) * sd.duration.t2^2 / cv.Q.dur.t2[2]
  CIs[2, 2] <- (n1.t2 - 1) * sd.duration.t2^2 / cv.Q.dur.t2[1]
  
  # Confidence intervals for inter-arrival mean and variance
  CIs[3, 1] <- mean.inter.arrivals.t2 - cv.T.time.t2[2] * sd.inter.arrivals.t2 / sqrt(n2.t2)
  CIs[3, 2] <- mean.inter.arrivals.t2 - cv.T.time.t2[1] * sd.inter.arrivals.t2 / sqrt(n2.t2)  # Corrected to add
  CIs[4, 1] <- (n2.t2 - 1) * sd.inter.arrivals.t2^2 / cv.Q.time.t2[2]
  CIs[4, 2] <- (n2.t2 - 1) * sd.inter.arrivals.t2^2 / cv.Q.time.t2[1]
  
  return(CIs)
}

#-------------------------------------------------------------------------------
# TYPE 2: Monte Carlo Simulation (parametric approach)
#-------------------------------------------------------------------------------
no.sims.t2 <- 1000                    # Number of simulations
B.t2 <- 1000                          # Number of bootstrap draws
alpha.t2 <- 0.05                      # Rejection level
n.t2 <- nrow(type2.data)              # Sample size

# Initialize rejection probability matrix
rej.prob.matrix.t2.1 <- matrix(NA, nrow=0, ncol=5,
                             dimnames = list(NULL, c("Sample Size", 
                                                     "Rej.Prob.Mean.Duration", 
                                                     "Rej.Prob.Var.Duration",
                                                     "Rej.Prob.Mean.Inter-Arrival.Time", 
                                                     "Rej.Prob.Var.Inter-Arrival.Time")))

for (i.t2 in 1:length(n.t2)) {
  n.i.t2 <- n.t2[i.t2]  # Current sample size
  
  # Store rejection results
  Reject.mat.t2 <- matrix(NA, nrow=no.sims.t2, ncol=4, 
                          dimnames = list(NULL, c("Reject Mean Duration", 
                                                  "Reject Variance Duration",
                                                  "Reject Mean Inter-Arrival Time", 
                                                  "Reject Variance Inter-Arrival Time")))
  
  for (j.t2 in 1:no.sims.t2) {
    # Generate data
    X.dur.t2 <- rgamma(n.i.t2, alpha.gamma.t2, 1/beta.gamma.t2)
    X.time.t2 <- rnorm(n.i.t2 - 1, mean.normal.t2, sd.normal.t2)
    
    # Get confidence intervals from bootstrap
    Conf.Ints.t2 <- bootstrap.t2(X.dur.t2, X.time.t2, B.t2, alpha.t2)
    
    # Record rejection results
    Reject.mat.t2[j.t2, 1] <- as.integer(mean(type2.data$Duration) < Conf.Ints.t2[1, 1]) + 
      as.integer(mean(type2.data$Duration) > Conf.Ints.t2[1, 2])
    Reject.mat.t2[j.t2, 2] <- as.integer(sd(type2.data$Duration)^2 < Conf.Ints.t2[2, 1]) + 
      as.integer(sd(type2.data$Duration)^2 > Conf.Ints.t2[2, 2])
    Reject.mat.t2[j.t2, 3] <- as.integer(mean.normal.t2 < Conf.Ints.t2[3, 1]) + 
      as.integer(mean.normal.t2 > Conf.Ints.t2[3, 2])
    Reject.mat.t2[j.t2, 4] <- as.integer(sd.normal.t2^2 < Conf.Ints.t2[4, 1]) + 
      as.integer(sd.normal.t2^2 > Conf.Ints.t2[4, 2])
  }
  
  # Calculate rejection probabilities
  rej.prob.t2 <- c(n.i.t2, colMeans(Reject.mat.t2))
  rej.prob.matrix.t2.1 <- rbind(rej.prob.matrix.t2.1, rej.prob.t2)
}

print(rej.prob.matrix.t2.1)

'''
Interpretation: 
All rejection probabilities are close to rejection level of 0.05. Hence, this build 
confidence that the parametric model with duration being gamma distributed
and inter-arrival times being normally distributed is well-specified. 
'''

print(Conf.Ints.t2*60)

#===============================================================================
#NON-PARAMETRIC APPROACH: 
#===============================================================================
#Key approach: draw bootstrap samples with replacement from EDF (empirical distribution function)

#-------------------------------------------------------------------------------
# TYPE 2: Monte Carlo Simulation (parametric approach)
#-------------------------------------------------------------------------------

no.sims.t2 <- 1000                    # Number of simulations
B.t2 <- 1000                          # Number of bootstrap draws
alpha.t2 <- 0.05                      # Rejection level
n.t2 <- nrow(type2.data)              # Sample size

# Initialize rejection probability matrix
rej.prob.matrix.t2.1 <- matrix(NA, nrow=0, ncol=5,
                               dimnames = list(NULL, c("Sample Size", 
                                                       "Rej.Prob.Mean.Duration", 
                                                       "Rej.Prob.Var.Duration",
                                                       "Rej.Prob.Mean.Inter-Arrival.Time", 
                                                       "Rej.Prob.Var.Inter-Arrival.Time")))

for (i.t2 in 1:length(n.t2)) {
  n.i.t2 <- n.t2[i.t2]  # Current sample size
  
  # Store rejection results in the following matrix:
  Reject.mat.t2.1 <- matrix(NA, nrow=no.sims.t2, ncol=4, 
                            dimnames = list(NULL, c("Reject Mean Duration", 
                                                    "Reject Variance Duration",
                                                    "Reject Mean Inter-Arrival Time", 
                                                    "Reject Variance Inter-Arrival Time")))
  
  for (j.t2 in 1:no.sims.t2) {
    # Sampling from EDF instead of gamma distribution
    EDF.duration <- sample(1: length(type2.data$Duration), n.i.t2, replace = TRUE) 
    X.duration <- type2.data$Duration[EDF.duration]
    #sampling from EDF instead of normal distribution 
    EDF.time <- sample(1:length(X.inter.arrival.t2), n.i.t2 - 1, replace = TRUE) #
    X.time <- X.inter.arrival.t2[EDF.time]
    # Get confidence intervals from bootstrap
    Conf.Ints.t2.np <- bootstrap.t2(X.duration, X.time, B.t2, alpha.t2)
    
    # Record rejection results
    Reject.mat.t2.1[j.t2, 1] <- as.integer(mean(type2.data$Duration) < Conf.Ints.t2.np[1, 1]) + 
      as.integer(mean(type2.data$Duration) > Conf.Ints.t2.np[1, 2])
    Reject.mat.t2.1[j.t2, 2] <- as.integer(sd(type2.data$Duration)^2 < Conf.Ints.t2.np[2, 1]) + 
      as.integer(sd(type2.data$Duration)^2 > Conf.Ints.t2.np[2, 2])
    Reject.mat.t2.1[j.t2, 3] <- as.integer(mean.normal.t2 < Conf.Ints.t2.np[3, 1]) + 
      as.integer(mean.normal.t2 > Conf.Ints.t2.np[3, 2])
    Reject.mat.t2.1[j.t2, 4] <- as.integer(sd.normal.t2^2 < Conf.Ints.t2.np[4, 1]) + 
      as.integer(sd.normal.t2^2 > Conf.Ints.t2.np[4, 2])
  }
  
  # Calculate rejection probabilities
  rej.prob.t2.1 <- c(n.i.t2, colMeans(Reject.mat.t2.1))
  rej.prob.matrix.t2.1 <- rbind(rej.prob.matrix.t2.1, rej.prob.t2.1)
}

print(rej.prob.matrix.t2.1)

'''
Interpretation: 
All rejection probabilities are close to rejection level of 0.05. Hence, this build 
confidence that the parametric model with duration being gamma distributed
and inter-arrival times being normally distributed is well-specified. 
'''
print(Conf.Ints.t2.np*60)


