###########################################
# Coxphf Benchmarking
# 2024-10-11

# load libraries
library(crrf)
library(coxphf)
library(survival)

# Load an example simulation result
res.dir <- "Y:/Anna/CRR/Results/"
load(paste0(res.dir, "N50umax25p30_test.RData"))
# extract an example dataset
dat <- res$data[[300]]
head(dat)

# combine events 1 and 2 so we can fit a Cox model
dat$status <- as.factor(ifelse(dat$stat==0, 0, 1))

# fit the coxphf model
fit1 <- coxphf(Surv(time, status)~x1+x2+x3, data=dat)
fit2 <- crrf(Surv(time, status)~x1 + x2 + x3, dset=dat, etype=1, CI=TRUE)

# Difference in coefficients
dif <- coef(fit1) - coef(fit2)
dif
# difference in exponentiated lower CI
dif.lower <- fit1$ci.lower - exp(fit2$CI.tbl[,2])
dif.upper <- fit1$ci.upper - exp(fit2$CI.tbl[,3])
dif.lower
dif.upper

# Fit for all of the data and save the difference
coef.dif <- data.frame()
lower.dif <- data.frame()
upper.dif <- data.frame()
for(i in 1:length(res$data)){
  print(i)
  # extract an example dataset
  dat <- res$data[[i]]
  #head(dat)
  
  # combine events 1 and 2 so we can fit a Cox model
  dat$status <- as.factor(ifelse(dat$stat==0, 0, 1))
  
  # fit the coxphf model
  fit1 <- coxphf(Surv(time, status)~x1+x2+x3, data=dat)
  fit2 <- crrf(Surv(time, status)~x1 + x2 + x3, dset=dat, etype=1, CI=TRUE)
  
  # Difference in coefficients
  dif <- coef(fit1) - coef(fit2)
  #dif
  coef.dif <- rbind.data.frame(coef.dif, dif)
  # difference in exponentiated lower CI
  dif.lower <- fit1$ci.lower - exp(fit2$CI.tbl[,2])
  dif.upper <- fit1$ci.upper - exp(fit2$CI.tbl[,3])
  #dif.lower
  #dif.upper
  lower.dif <- rbind.data.frame(lower.dif, dif.lower)
  upper.dif <- rbind.data.frame(upper.dif, dif.upper)
}


