###########################################
# Results Sandbox
# 2024-10-07

# figure out how to get all the names of the files in the folder
# then loop through each file
res.dir <- "Y:/Anna/CRR/Results/"
res.list <- list.files(path=res.dir, pattern="RData")

# data to save:
# tables
event.prop.tab <- data.frame()
avg.bias.tab <- data.frame()
avg.cov.tab <- data.frame()
avg.width.tab <- data.frame()
mse.tab <- data.frame()
sum.tab <- data.frame()
# for figures
bias.df <- data.frame()
width.df <- data.frame()

for (i in 1:length(res.list)){
  load(paste0(res.dir, res.list[[i]])) # load file name
  # Save event proportion
  event.prop <- as.vector(colMeans(res$events))
  event.prop.name <- c(res$sample.size, res$umax, res$event.rate, round(event.prop, 3))
  event.prop.tab <- rbind.data.frame(event.prop.tab, event.prop.name)
  # Save bias table and for figures
  bias.crr.temp <- res$bias$CRR
  avg.bias.crr <- c(colMeans(bias.crr.temp, na.rm=T))
  avg.bias.crr.name <- c(res$sample.size, res$umax, res$event.rate,
                         round(avg.bias.crr, 3), 
                         length(which(is.na(bias.crr.temp$x1))), 
                         "CRR")
  bias.crr.temp$Method <- rep("CRR", times=nrow(bias.crr.temp))
  bias.crr.temp$Setting <- rep(res$sim.setting, times=nrow(bias.crr.temp))
  bias.crr.temp$SampleSize <- rep(res$sample.size, times=nrow(bias.crr.temp))
  bias.crr.temp$umax <- rep(res$umax, times=nrow(bias.crr.temp))
  bias.crr.temp$EventRate <- rep(res$event.rate, times=nrow(bias.crr.temp))
  bias.firth.temp <- res$bias$Firth
  avg.bias.firth <- c(colMeans(bias.firth.temp, na.rm=T))
  avg.bias.firth.name <- c(res$sample.size, res$umax, res$event.rate,
                           round(avg.bias.firth, 3), 
                           length(which(is.na(bias.firth.temp$x1))), 
                           "Firth")
  avg.bias.tab <- rbind.data.frame(avg.bias.tab, avg.bias.crr.name, avg.bias.firth.name)
  bias.firth.temp$Method <- rep("Firth", times=nrow(bias.firth.temp))
  bias.firth.temp$Setting <- rep(res$sim.setting, times=nrow(bias.firth.temp))
  bias.firth.temp$SampleSize <- rep(res$sample.size, times=nrow(bias.firth.temp))
  bias.firth.temp$umax <- rep(res$umax, times=nrow(bias.firth.temp))
  bias.firth.temp$EventRate <- rep(res$event.rate, times=nrow(bias.firth.temp))
  bias.df <- rbind.data.frame(bias.df, bias.crr.temp, bias.firth.temp)
  # Save coverage table
  cov.crr.temp <- res$coverage$CRR
  cov.firth.temp <- res$coverage$Firth
  avg.cov.crr <- c(colMeans(cov.crr.temp, na.rm=T))
  avg.cov.crr.name <- c(res$sample.size, res$umax, res$event.rate,
                        round(avg.cov.crr, 3), 
                        length(which(is.na(cov.crr.temp$x1))), 
                        "CRR", res$sim.setting)
  avg.cov.firth <- c(colMeans(cov.firth.temp, na.rm=T))
  avg.cov.firth.name <- c(res$sample.size, res$umax, res$event.rate,
                          round(avg.cov.firth, 3), 
                          length(which(is.na(cov.firth.temp$x1))), 
                          "Firth", res$sim.setting)
  avg.cov.tab <- rbind.data.frame(avg.cov.tab, avg.cov.crr.name, avg.cov.firth.name)
  # save the MSE
  bias.crr.temp.sq <- res$bias$CRR^2
  crr.divsor <- length(which(complete.cases(bias.crr.temp.sq)))
  mse.crr <- colSums(bias.crr.temp.sq, na.rm=T)/crr.divsor
  mse.tab.crr.name <- c(res$sample.size, res$umax, res$event.rate,
                    round(mse.crr, 3), length(which(is.na(bias.crr.temp.sq))),
                    "CRR")
  bias.firth.temp.sq <- res$bias$Firth^2
  firth.divsor <- length(which(complete.cases(bias.firth.temp.sq)))
  mse.firth <- colSums(bias.firth.temp.sq, na.rm=T)/firth.divsor
  mse.tab.firth.name <- c(res$sample.size, res$umax, res$event.rate,
                        round(mse.firth, 3), length(which(is.na(bias.firth.temp.sq))),
                        "Firth")
  mse.tab <- rbind.data.frame(mse.tab, mse.tab.crr.name, mse.tab.firth.name)
  # Create summary table
  sum.tab.name <- c(res$sample.size, res$umax, res$event.rate,
                    paste0(round(avg.bias.crr[1], 3), " (", round(mse.crr[1],3), ")"),
                    paste0(round(avg.bias.firth[1], 3), " (", round(mse.firth[1],3), ")"),
                    paste0(round(avg.bias.crr[2], 3), " (", round(mse.crr[2],3), ")"),
                    paste0(round(avg.bias.firth[2], 3), " (", round(mse.firth[2],3), ")"),
                    paste0(round(avg.bias.crr[3], 3), " (", round(mse.crr[3],3), ")"),
                    paste0(round(avg.bias.firth[3], 3), " (", round(mse.firth[3],3), ")"),
                    length(which(is.na(bias.crr.temp.sq))),
                    length(which(is.na(bias.firth.temp.sq))))
  sum.tab <- rbind.data.frame(sum.tab, sum.tab.name)
  
}
colnames(event.prop.tab) <- c("Sample Size", "umax", "Event Rate",
                              "Proportion Censored","Proportion Event 1",
                              "Proportion Event 2")
#colnames(bias.df) <- c("x1", "x2", "x3", "Method", "Setting")
colnames(avg.bias.tab) <- c("Sample Size", "umax", "Event Rate", "x1", "x2", "x3", "Missing", "Method")
colnames(avg.cov.tab) <- c("Sample Size", "umax", "Event Rate", "x1", "x2", "x3", "Missing", "Method", "Setting")
colnames(mse.tab) <- c("Sample Size", "umax", "Event Rate", "x1", "x2", "x3", "Missing", "Method")
colnames(sum.tab) <- c("Sample Size", "umax", "Event Rate",
                       "x1.crr", "x1.firth", "x2.crr", "x2.firth", "x3.crr", "x3.firth",
                       "Missing.crr", "Missing.Firth")

# Save Tables
write.csv(avg.bias.tab, file=paste0(res.dir, "Avg_Bias.csv"))
write.csv(avg.cov.tab, file=paste0(res.dir, "Avg_Covg.csv"))
write.csv(event.prop.tab, file=paste0(res.dir, "Event_Prop.csv"))
write.csv(sum.tab, file=paste0(res.dir, "Summary_Tab.csv"))
save(bias.df, file=paste0(res.dir, "Bias.RData"))

# Make Figures
bias.df.N50 <- bias.df[grep("N=50", bias.df$Setting),]

library(ggplot2)
p <- ggplot(aes(x=Setting, y=x1, fill=Method, col=Method), data=bias.df.N50)+geom_boxplot()
p + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + ylab("Bias in x1")

bias.df.N100 <- bias.df[grep("N=100", bias.df$Setting),]

library(ggplot2)
p <- ggplot(aes(x=Setting, y=x1, fill=Method, col=Method), data=bias.df.N100)+geom_boxplot()
p + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + ylab("Bias in x1")

bias.df.N200 <- bias.df[grep("N=200", bias.df$Setting),]

library(ggplot2)
p <- ggplot(aes(x=Setting, y=x1, fill=Method), data=bias.df.N200)+geom_boxplot()
p + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + ylab("Bias in x1")

# Coverage plots
temp.string <- avg.cov.tab$Setting[1]
library(stringr)
trimws(str_split_i(temp.string, pattern=",", 3))

sample.size <- c()
setting.sm <- c()
for(i in 1:nrow(avg.cov.tab)){
  temp.string <- avg.cov.tab$Setting[i]
  sample.size[i] <- str_split_i(temp.string, pattern=",", 1)
  setting.sm[i]  <- paste(trimws(str_split_i(temp.string, pattern=",", 2)), trimws(str_split_i(temp.string, pattern=",", 3)), sep=",")
}
avg.cov.tab$SampleSize <- sample.size
avg.cov.tab$setting.sm <- setting.sm

library(ggplot2)
avg.cov.N50 <- avg.cov.tab[grep("N=50", avg.cov.tab$Setting),]
avg.cov.N50$x1 <- as.numeric(avg.cov.N50$x1)
p <- ggplot(aes(x=setting.sm, y=x1, col=Method, fill=Method), data=avg.cov.N50) + 
  geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + 
  ylab("Average Coverage") + xlab("Setting") + geom_abline(slope=0, intercept=0.95) + 
  ggtitle("Average Coverage for Settings with N=50") + coord_cartesian(ylim=c(0.5, 1.00))
p

avg.cov.N100 <- avg.cov.tab[grep("N=100", avg.cov.tab$Setting),]
avg.cov.N100$x1 <- as.numeric(avg.cov.N100$x1)
p <- ggplot(aes(x=setting.sm, y=x1, col=Method, fill=Method), data=avg.cov.N100) + 
  geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + 
  ylab("Average Coverage") + xlab("Setting") + geom_abline(slope=0, intercept=0.95) + 
  ggtitle("Average Coverage for Settings with N=100")
p

avg.cov.N200 <- avg.cov.tab[grep("N=200", avg.cov.tab$Setting),]
avg.cov.N200$x1 <- as.numeric(avg.cov.N200$x1)
p <- ggplot(aes(x=setting.sm, y=x1, col=Method, fill=Method), data=avg.cov.N200) + 
  geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + 
  ylab("Average Coverage") + xlab("Setting") + geom_abline(slope=0, intercept=0.95) + 
  ggtitle("Average Coverage for Settings with N=200")
p
