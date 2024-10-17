##############################################
# Firth Fine-Gray Simulation
# AES
# 2024-04-19
#############################################

###########################################
# specify the result directory

message("Checking that result directory exists.")

res.dir="/home/aseffern/Firth-CRR/Results/2024-04-19/"

if (!dir.exists(res.dir))
  stop("result directory not found.")
message("Result directory exists.")

################################################
# load code libraries
library(devtools)
#devtools::install_github("annaSeffernick/crrf")
library(crrf)
library(fastcmprsk)
library(survival)
library(cmprsk)


#####################################
# Simulate 2-cause Fine-Gray Data
# Using fastcmprsk package

# parameters to change
seed1 <- 6466
N <- 200 # number of observations (50, 100, 200)
umax=0.25 # control censoring rate (increase censorign by decreasing umax) (umax=1 ~ 50% censoring; umax=0.55 ~65% censoring; umax=0.25 ~ 80% censoring)
p=0.3 # control event rate for event of interest (0.5, 0.3)

# constant params
# Create coefficient vector for event of interest and competing event
beta1 <- c(0.40, 0.60, -0.80)
beta2 <- -beta1
umin=0
B <- 1000 # number of simulation replicates

data.list <- list()
event.prop <- data.frame()
firth.bias <- data.frame()
crr.bias <- data.frame()
#beta.dif <- data.frame()
firth.coverage <- data.frame()
crr.coverage <- data.frame()
firth.width <- data.frame()
crr.width <- data.frame()
for(k in 1:B){
  #set.seed(209+10*i)
  set.seed(seed1+11*k)
  # Simulate design matrix
  Z <- matrix(rnorm(N*length(beta1)), nrow=N)
  # Generate data
  dat <- simulateTwoCauseFineGrayModel(N, beta1, beta2, Z,
                                       u.min=umin, u.max=umax, p=p)
  x1 <- Z[,1]
  x2 <- Z[,2]
  x3 <- Z[,3]
  # Event counts (0=censored; 1=event of interest; 2=competing event)
  #table(dat$fstatus)
  dat$status <- as.factor(dat$fstatus)
  # First 6 observed survival times
  #head(dat$ftime)
  dset <- cbind.data.frame(dat$ftime, dat$status, x1, x2, x3)
  colnames(dset)[1:2] <- c("time", "stat")
  data.list[[k]] <- dset
  event.tab <- as.vector(table(dset$stat)/nrow(dset))
  event.prop <- rbind.data.frame(event.prop, event.tab)
  
  # Fit models
  form <- Surv(time, stat)~x1+x2+x3
  failcode="1"
  # Firth
  #firth.fit <- firth.fg(dat,form, failcode, use.pen = T)
  #fgm(Surv(ftime,status)~x1+x2+x3,dat,T,T)
  firth.fit <- try(crrf(form,failcode,dset,T,T, eps=1e-10), silent=TRUE)
  if(inherits(firth.fit, "try-error")){
    print(table(dset$stat))
    f.bias <- rep(NA, times=length(beta1))
    firth.bias <- rbind.data.frame(firth.bias, f.bias)
    f.covg <- rep(NA, times=length(beta1))
    firth.coverage <- rbind.data.frame(firth.coverage, f.covg)
    #temp.dif <- rep(NA, times=length(beta1))
    #beta.dif <- rbind.data.frame(beta.dif, temp.dif)
    f.width <- rep(NA, times=length(beta1))
    firth.width <- rbind.data.frame(firth.width, f.width)
  }
  else{
    firth.tbl <- as.data.frame(firth.fit$CI.tbl)
    f.bias <- firth.tbl$ml.beta - beta1
    firth.bias <- rbind.data.frame(firth.bias, f.bias)
    f.covg <- c()
    for(r in 1:length(beta1)){
      f.covg[r] <- ifelse(beta1[r] >= firth.tbl$lower[r] & beta1[r]<=firth.tbl$upper[r], 1, 0)
    }
    firth.coverage <- rbind.data.frame(firth.coverage, f.covg)
    f.width <- firth.tbl$upper - firth.tbl$lower
    firth.width <- rbind.data.frame(firth.width, f.width)
  }
  
  # crr
  crr.fit=try(crr(dat$ftime,dat$status,Z,failcode=failcode), silent=TRUE)
  if(inherits(crr.fit, "try-error")){
    c.bias <- rep(NA, times=length(beta1))
    crr.bias <- rbind.data.frame(crr.bias, c.bias)
    c.covg <- rep(NA, times=length(beta1))
    crr.coverage <- rbind.data.frame(crr.coverage, c.covg)
    #temp.dif <- rep(NA, times=length(beta1))
    #beta.dif <- rbind.data.frame(beta.dif, temp.dif)
    c.width <- rep(NA, times=length(beta1))
    crr.width <- rbind.data.frame(crr.width, c.width)
  }
  else{
    c.bias <- crr.fit$coef - beta1
    crr.bias <- rbind.data.frame(crr.bias, c.bias)
    #temp.dif <- firth.tbl$ml.beta - crr.fit$coef
    #beta.dif <- rbind.data.frame(beta.dif, temp.dif)
    sum.crr.fit <- summary(crr.fit)
    c.covg <- c()
    c.width <- c()
    for(s in 1:length(beta1)){
      temp.b.exp <- exp(beta1[s])
      c.covg[s] <- ifelse(temp.b.exp >= sum.crr.fit$conf.int[s, 3] & temp.b.exp <= sum.crr.fit$conf.int[s, 4], 1, 0)
      c.width[s] <- sum.crr.fit$conf.int[s, 4] - sum.crr.fit$conf.int[s, 3]
    }
    crr.coverage <- rbind.data.frame(crr.coverage, c.covg)
    crr.width <- rbind.data.frame(crr.width, c.width)
  }
}
colnames(event.prop) <- c("0", "1", "2")
colnames(crr.coverage) <- c("x1", "x2", "x3")
colnames(firth.coverage) <- c("x1", "x2", "x3")
colnames(firth.bias) <- c("x1", "x2", "x3")
colnames(crr.bias) <- c("x1", "x2", "x3")
#colnames(beta.dif) <- c("x1", "x2", "x3")
colnames(firth.width) <- c("x1", "x2", "x3")
colnames(crr.width) <- c("x1", "x2", "x3")
sim.setting <- paste0("N=", N, ", umax=", umax, ", p=", p)

coverage.list <- list(crr.coverage, firth.coverage)
names(coverage.list) <- c("CRR", "Firth")
bias.list <- list(crr.bias, firth.bias)
names(bias.list) <- c("CRR", "Firth")
width.list <- list(crr.width, firth.width)
names(width.list) <- c("CRR", "Firth")

res <- list(data=data.list, events=event.prop, coverage=coverage.list,
            bias=bias.list, ci.width=width.list,
            sim.setting=sim.setting, sample.size=N, beta1=beta1, 
            umax=umax, event.rate=p, sim.reps=B)

# Save results
save(res, file=paste0(res.dir,"N",N,"umax",umax*100,"p",p*100, "_test.RData"))

