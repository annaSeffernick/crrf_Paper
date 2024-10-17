#########################################################
# Application Analysis
# Clay ACT Data
# AES
# 2024-06-26
###########################################################

# Set up
library(openxlsx)
library("devtools", quietly=T)
# devtools::install_github("annaSeffernick/crrf")
library(crrf)

# Load data
# supplementary table
act.dat <- read.xlsx("Y:/Anna/CRR/Data/Clay-JCO-PO-Supplement-DS_po.19.00163-1.xlsx",  rows = 3:51)

# full dataset from Stan
full.dat <- read.xlsx("Y:/Anna/CRR/Data/Adrenocortical Tumor Master Clinical Excel 1.30.2019 no MRNs for Anna.xlsx",
                      sheet="ACT Active Excel 6.13.2018", detectDates=TRUE)

# Compare the datasets
length(which(act.dat$Sentrix_ID %in% full.dat$Sentrix_ID)) # all 48 act cases in full.dat

# Calculate time to relapse for the patients that relapsed
rel.ids <- act.dat[which(act.dat$Relapse_YN=="YES"),]$Sentrix_ID # ids of patients who relapsed
rel.dat <- full.dat[which(full.dat$Sentrix_ID %in% rel.ids),] # filter the full dataset to relapse pts
rel.dat$Relapse_Time <- difftime(as.Date(rel.dat$Relapse_Date,format="%Y-%m-%d"), as.Date(rel.dat$Collection.Date, format="%Y-%m-%d"), unit="days")

# Create variables needed for crrf package
act.dat$RFS <- act.dat$Follow_Up_Days # Define relapse-free survival variable: time until relapse or death w/o relapse or last follow-up if censored
rel.dat.or <- rel.dat[match(rel.ids, rel.dat$Sentrix_ID),] # reorder rel.dat
act.dat[which(act.dat$Relapse_YN=="YES"),]$RFS <- rel.dat.or$Relapse_Time # replase follow-up time with time to relapse for those who relapse
act.dat$RFS_status <- act.dat$Mortality
act.dat[which(act.dat$Relapse_YN=="YES"),]$RFS_status <- "Relapse"
table(act.dat$RFS_status)
act.dat$RFS_stat_num <- factor(ifelse(act.dat$RFS_status=="Relapse", 1, 
                               ifelse(act.dat$RFS_status=="Deceased", 2, 0)))
table(act.dat$RFS_stat_num)

write.csv(act.dat, file="Y:/Anna/CRR/Data/ACT_Data.csv", row.names = FALSE)
act.dat <- read.csv("Y:/Anna/CRR/Data/ACT_Data.csv")
act.dat$RFS_stat_num <- factor(act.dat$RFS_stat_num) # this part is annoying

# Fit crrf
res <- crrf(Surv(RFS, RFS_stat_num)~Methylation_Group + Gender, etype="1", 
            dset=act.dat, firth=TRUE, CI=TRUE)
res

crrf(Surv(RFS, RFS_stat_num)~Age_years, etype="1", 
     dset=act.dat, firth=TRUE, CI=TRUE)

act.dat$RFS_STATUS <- factor(act.dat$RFS_status)
crrf(Surv(RFS, RFS_STATUS)~Age_years, etype="Relapse", 
     dset=act.dat, firth=TRUE, CI=TRUE)

# Paper code
res <- crrf(Surv(RFS, RFS_status)~Methylation_Group, etype="Relapse", 
            dset=act.dat, firth=TRUE, CI=TRUE)
res
res.test <- crrftest(Surv(RFS, RFS_status)~Methylation_Group, etype="Relapse", 
         dset=act.dat, test=~1, firth=TRUE)
res.test$pvalue
table(act.dat$Wieneke_Classification, act.dat$RFS_status)

res2 <- crrf(Surv(RFS, RFS_status)~Wieneke_Classification, etype="Relapse", 
             dset=act.dat, firth=TRUE, CI=TRUE)
res2
res3 <- crrf(Surv(RFS, RFS_status)~Age_years, etype="Relapse", 
             dset=act.dat, firth=TRUE, CI=TRUE)
res3
res4 <- crrf(Surv(RFS, RFS_status)~ Wieneke_Classification + Age_years+Methylation_Group, etype="Relapse", 
             dset=act.dat, firth=TRUE, CI=TRUE)
res4

# Make Table A3
x1 <- c("Factor", "HR", "CILB", "CIUB", "P")
x2 <- c("Methylation group alone", round(exp(coef(res)),4), 
        round(exp(res$CI.tbl)[,2:3],4), round(res.test$pvalue,4))
# Adjust for age
x3 <- c("Adjusted for age", " ", " ", " ", " ")
res.age <- crrf(Surv(RFS, RFS_status)~Methylation_Group + Age_years,
                etype="Relapse", dset=act.dat, firth=TRUE, CI=TRUE)
res.age.test1 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + Age_years, 
                          etype="Relapse", dset=act.dat, 
                          test=~Methylation_Group, firth=TRUE)
res.age.test2 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + Age_years, 
                          etype="Relapse", dset=act.dat, 
                          test=~Age_years, firth=TRUE)
x4 <- c("Methylation group", round(exp(coef(res.age)[1]), 4), 
        round(exp(res.age$CI.tbl)[1,2:3], 4), round(res.age.test2$pvalue, 4))
x5 <- c("Age", round(exp(coef(res.age)[2]), 4), 
        round(exp(res.age$CI.tbl)[2,2:3], 4), round(res.age.test1$pvalue, 4))
act.tab <- rbind.data.frame(x2, x3, x4, x5)
colnames(act.tab) <- x1
# Adjust for gender
x6 <- c("Adjusted for biological sex", "", " ", " ", " ")
res.sex <- crrf(Surv(RFS, RFS_status)~Methylation_Group + Gender,
                etype="Relapse", dset=act.dat, firth=TRUE, CI=TRUE)
res.sex.test1 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + Gender, 
                          etype="Relapse", dset=act.dat, 
                          test=~Methylation_Group, firth=TRUE)
res.sex.test2 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + Gender, 
                          etype="Relapse", dset=act.dat, 
                          test=~Gender, firth=TRUE)
x7 <- c("Methylation group", round(exp(coef(res.sex)[1]), 4), 
        round(exp(res.sex$CI.tbl)[1,2:3], 4), round(res.sex.test2$pvalue, 4))
x8 <- c("Biological sex", round(exp(coef(res.sex)[2]), 4), 
        round(exp(res.sex$CI.tbl)[2,2:3], 4), round(res.sex.test1$pvalue, 4))
act.tab <- rbind.data.frame(act.tab, x6, x7, x8)
# Adjust for therapy
x9 <- c("Adjusted for therapy", " ", " ", " ", " ")
res.therapy <- crrf(Surv(RFS, RFS_status)~Methylation_Group + Initial_Treatment,
                    etype="Relapse", dset=act.dat, firth=TRUE, CI=TRUE)
res.therapy.test1 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + Initial_Treatment, 
                              etype="Relapse", dset=act.dat, 
                              test=~Methylation_Group, firth=TRUE)
res.therapy.test2 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + Initial_Treatment, 
                              etype="Relapse", dset=act.dat, 
                              test=~Initial_Treatment, firth=TRUE)
x10 <- c("Methylation group", round(exp(coef(res.therapy)[1]), 4), 
         round(exp(res.therapy$CI.tbl)[1,2:3], 4), round(res.therapy.test2$pvalue, 4))
x11 <- c("Therapy S+C:S", round(exp(coef(res.therapy)[2]), 4), 
         round(exp(res.therapy$CI.tbl)[2,2:3], 4), round(res.therapy.test1$pvalue, 4))
act.tab <- rbind.data.frame(act.tab, x9, x10, x11)
# Adjust for Clinical Presentation - too small of groups
# x12 <- c("Adjusted for clinical presentation", " ", " ", " ", " ")
# res.clin <- crrf(Surv(RFS, RFS_status)~Methylation_Group + Clinical_Presentation,
#                  etype="Relapse", dset=act.dat, firth=TRUE, CI=TRUE)
# res.clin.test1 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + Clinical_Presentation, 
#                               etype="Relapse", dset=act.dat, 
#                               test=~Methylation_Group, firth=TRUE)
# res.clin.test2 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + Clinical_Presentation, 
#                            etype="Relapse", dset=act.dat, 
#                            test=~Clinical_Presentation, firth=TRUE)
# x13 <- c("Methylation group", round(exp(coef(res.clin)[1]), 4), 
#          round(exp(res.clin$CI.tbl)[1,2:3], 4), round(res.clin.test2$pvalue, 4))
# x14 <- c("Clinical presentation", round(exp(coef(res.clin)[2]), 4), 
#          round(exp(res.clin$CI.tbl)[2,2:3], 4), round(res.clin.test1$pvalue, 4))

# Adjust for Stage
#act.dat$Stage[which(act.dat$Stage=="III ")] <- "III"
act.dat$stage <- as.factor(act.dat$Stage)
x15 <- c("Adjusted for tumor stage", " ", " ", " ", " ")
res.stage <- crrf(Surv(RFS, RFS_status)~Methylation_Group + stage,
                  etype="Relapse", dset=act.dat, firth=TRUE, CI=TRUE)
res.stage.test1 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + stage, 
                            etype="Relapse", dset=act.dat, 
                            test=~Methylation_Group, firth=TRUE)
res.stage.test2 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + stage, 
                            etype="Relapse", dset=act.dat, 
                            test=~stage, firth=TRUE)
x16 <- c("Methylation group", round(exp(coef(res.stage)[1]), 4), 
         round(exp(res.stage$CI.tbl)[1,2:3], 4), 
         round(res.stage.test2$pvalue, 4))
x17 <- c("Tumor Stage II:I*", round(exp(coef(res.stage)[2]), 4), 
         round(exp(res.stage$CI.tbl)[2,2:3], 4), 
         round(res.stage.test1$pvalue, 4))
x18 <- c("Tumor Stage III:I*", round(exp(coef(res.stage)[3]), 4), 
         round(exp(res.stage$CI.tbl)[3,2:3], 4), 
         round(res.stage.test1$pvalue, 4))
x19 <- c("Tumor Stage IV:I*", round(exp(coef(res.stage)[4]), 4), 
         round(exp(res.stage$CI.tbl)[4,2:3], 4), 
         round(res.stage.test1$pvalue, 4))
act.tab <- rbind.data.frame(act.tab, x15, x16, x17, x18, x19)
# Adjust for Wieneke Classification
x20 <- c("Adjusted for Wieneke Classification", " ", " ", " ", " ")
res.wie <- crrf(Surv(RFS, RFS_status)~Methylation_Group + Wieneke_Classification,
                  etype="Relapse", dset=act.dat, firth=TRUE, CI=TRUE)
res.wie.test1 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + Wieneke_Classification, 
                            etype="Relapse", dset=act.dat, 
                            test=~Methylation_Group, firth=TRUE)
res.wie.test2 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + Wieneke_Classification, 
                            etype="Relapse", dset=act.dat, 
                            test=~Wieneke_Classification, firth=TRUE)
x21 <- c("Methylation group", round(exp(coef(res.wie)[1]), 4), 
         round(exp(res.wie$CI.tbl)[1,2:3], 4), 
         round(res.wie.test2$pvalue, 4))
x22 <- c("Wieneke Classification Carcinoma:Adenoma*", round(exp(coef(res.wie)[2]), 4), 
         round(exp(res.wie$CI.tbl)[2,2:3], 4), 
         round(res.wie.test1$pvalue, 4))
x23 <- c("Wieneke Classification UMP:Adenoma*", round(exp(coef(res.wie)[3]), 4), 
         round(exp(res.wie$CI.tbl)[3,2:3], 4), 
         round(res.wie.test1$pvalue, 4))
act.tab <- rbind.data.frame(act.tab, x20, x21, x22, x23)
# Adjust for Tumor Weight - doesn't fit, because missing?
# x24 <- c("Adjusted for tumor weight", " ", " ", " ", " ")
# res.tw <- crrf(Surv(RFS, RFS_status)~Methylation_Group + Tumor_Weight_Grams,
#                 etype="Relapse", dset=act.dat, firth=TRUE, CI=TRUE)

# Adjust for Tumor Size
x24 <- c("Adjusted for tumor size", " ", " ", " ", " ")
res.ts <- crrf(Surv(RFS, RFS_status)~Methylation_Group + Tumor_Size_cm,
               etype="Relapse", dset=act.dat, firth=TRUE, CI=TRUE)
res.ts.test1 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + Tumor_Size_cm, 
                          etype="Relapse", dset=act.dat, 
                          test=~Methylation_Group, firth=TRUE)
res.ts.test2 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + Tumor_Size_cm, 
                         etype="Relapse", dset=act.dat, 
                         test=~Tumor_Size_cm, firth=TRUE)
x25 <- c("Methylation group", round(exp(coef(res.ts)[1]), 4), 
         round(exp(res.ts$CI.tbl)[1,2:3], 4), 
         round(res.ts.test2$pvalue, 4))
x26 <- c("Tumor Size", round(exp(coef(res.ts)[2]), 4), 
         round(exp(res.ts$CI.tbl)[2,2:3], 4), 
         round(res.ts.test1$pvalue, 4))
act.tab <- rbind.data.frame(act.tab, x24, x25, x26)
# Adjust for Ki67_Percent
act.dat$Ki67_15 <- ifelse(act.dat$Ki67_Percent <=15, "No", "Yes")
x27 <- c("Adjust for Ki67 Above 15%", " ", " ", " ", " ")
res.ki <- crrf(Surv(RFS, RFS_status)~Methylation_Group + Ki67_15,
               etype="Relapse", dset=act.dat, firth=TRUE, CI=TRUE)
res.ki.test1 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + Ki67_15, 
                         etype="Relapse", dset=act.dat, 
                         test=~Methylation_Group, firth=TRUE)
res.ki.test2 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + Ki67_15, 
                         etype="Relapse", dset=act.dat, 
                         test=~Ki67_15, firth=TRUE)
x28 <- c("Methylation group", round(exp(coef(res.ki)[1]), 4), 
         round(exp(res.ki$CI.tbl)[1,2:3], 4), 
         round(res.ki.test2$pvalue, 4))
x29 <- c("Ki67 > 15% yes:no", round(exp(coef(res.ki)[2]), 4), 
         round(exp(res.ki$CI.tbl)[2,2:3], 4), 
         round(res.ki.test1$pvalue, 4))
act.tab <- rbind.data.frame(act.tab, x27, x28, x29)
# adjust for p53_IHC - probably want to change the level?
x30 <- c("Adjust for p53 IHC", " ", " ", " ", " ")
res.p53 <- crrf(Surv(RFS, RFS_status)~Methylation_Group + P53_IHC,
                etype="Relapse", dset=act.dat, firth=TRUE, CI=TRUE)
res.p53.test1 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + P53_IHC, 
                         etype="Relapse", dset=act.dat, 
                         test=~Methylation_Group, firth=TRUE)
res.p53.test2 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + P53_IHC, 
                          etype="Relapse", dset=act.dat, 
                          test=~P53_IHC, firth=TRUE)
x31 <- c("Methylation group", round(exp(coef(res.p53)[1]), 4), 
         round(exp(res.p53$CI.tbl)[1,2:3], 4), 
         round(res.p53.test2$pvalue, 4))
x32 <- c("P53 IHC Negative: Absent*", round(exp(coef(res.p53)[2]), 4), 
         round(exp(res.p53$CI.tbl)[2,2:3], 4), 
         round(res.p53.test1$pvalue, 4))
x33 <- c("P53 IHC Positive: Absent*", round(exp(coef(res.p53)[2]), 4), 
         round(exp(res.p53$CI.tbl)[2,2:3], 4), 
         round(res.p53.test1$pvalue, 4))
act.tab <- rbind.data.frame(act.tab, x30, x31, x32, x33)
# adjust for TP53 Mutation
x34 <- c("Adjust for TP53 Mutation Category", " ", " ", " ", " ")
res.p53.mc <- crrf(Surv(RFS, RFS_status)~Methylation_Group + TP53_Mutation_Category,
                etype="Relapse", dset=act.dat, firth=TRUE, CI=TRUE)
res.p53.mc.test1 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + TP53_Mutation_Category, 
                          etype="Relapse", dset=act.dat, 
                          test=~Methylation_Group, firth=TRUE)
res.p53.mc.test2 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + TP53_Mutation_Category, 
                             etype="Relapse", dset=act.dat, 
                             test=~TP53_Mutation_Category, firth=TRUE)
x35 <- c("Methylation group", round(exp(coef(res.p53.mc)[1]), 4), 
         round(exp(res.p53.mc$CI.tbl)[1,2:3], 4), 
         round(res.p53.mc.test2$pvalue, 4))
x36 <- c("TP53 Mutation Cateogry Somatic:Germline*",
         round(exp(coef(res.p53.mc)[2]), 4), 
         round(exp(res.p53.mc$CI.tbl)[2,2:3], 4), 
         round(res.p53.mc.test1$pvalue, 4))
x37 <- c("TP53 Mutation Cateogry Somatic:Wildtype*",
         round(exp(coef(res.p53.mc)[3]), 4), 
         round(exp(res.p53.mc$CI.tbl)[3,2:3], 4), 
         round(res.p53.mc.test1$pvalue, 4))
act.tab <- rbind.data.frame(act.tab, x34, x35, x36, x37)
# Adjust for ATRX IHC
x38 <- c("Adjust for ATRX IHC", " ", " ", " ", " ")
res.atrx <- crrf(Surv(RFS, RFS_status)~Methylation_Group + ATRX_IHC,
                   etype="Relapse", dset=act.dat, firth=TRUE, CI=TRUE)
res.atrx.test1 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + ATRX_IHC, 
                           etype="Relapse", dset=act.dat, 
                           test=~Methylation_Group, firth=TRUE)
res.atrx.test2 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + ATRX_IHC, 
                           etype="Relapse", dset=act.dat, 
                           test=~ATRX_IHC, firth=TRUE)
x39 <- c("Methylation group", round(exp(coef(res.atrx)[1]), 4), 
         round(exp(res.atrx$CI.tbl)[1,2:3], 4), 
         round(res.atrx.test2$pvalue, 4))
x40 <- c("ATRX IHC Loss:Intact", round(exp(coef(res.atrx)[2]), 4), 
         round(exp(res.atrx$CI.tbl)[2,2:3], 4), 
         signif(res.atrx.test1$pvalue, digits=4))
act.tab <- rbind.data.frame(act.tab, x38, x39, x40)
# Adjust for ATRX Mutation Category
x41 <- c("Adjust for ATRX Mutation Category", " ", " ", " ", " ")
res.atrx.mc <- crrf(Surv(RFS, RFS_status)~Methylation_Group + ATRX_Mutation_Category,
                    etype="Relapse", dset=act.dat, firth=TRUE, CI=TRUE)
res.atrx.mc.test1 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + ATRX_Mutation_Category, 
                           etype="Relapse", dset=act.dat, 
                           test=~Methylation_Group, firth=TRUE)
res.atrx.mc.test2 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + ATRX_Mutation_Category, 
                              etype="Relapse", dset=act.dat, 
                              test=~ATRX_Mutation_Category, firth=TRUE)
x42 <- c("Methylation group", round(exp(coef(res.atrx.mc)[1]), 4), 
         round(exp(res.atrx.mc$CI.tbl)[1,2:3], 4), 
         round(res.atrx.mc.test2$pvalue, 4))
x43 <- c("ATRX Mutation Category Wildtype:Somatic", round(exp(coef(res.atrx.mc)[2]), 4), 
         round(exp(res.atrx.mc$CI.tbl)[2,2:3], 4), 
         round(res.atrx.mc.test1$pvalue, 4))
act.tab <- rbind.data.frame(act.tab, x41, x42, x43)
# Adjust for CTNNB1 Mutation Category
x44 <- c("Adjust for CTNNB1 Mutation Category", " ", " ", " ", " ")
res.ct <- crrf(Surv(RFS, RFS_status)~Methylation_Group + CTNNB1_Mutation_Category,
                    etype="Relapse", dset=act.dat, firth=TRUE, CI=TRUE)
res.ct.test1 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + CTNNB1_Mutation_Category, 
                              etype="Relapse", dset=act.dat, 
                              test=~Methylation_Group, firth=TRUE)
res.ct.test2 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + CTNNB1_Mutation_Category, 
                         etype="Relapse", dset=act.dat, 
                         test=~CTNNB1_Mutation_Category, firth=TRUE)
x45 <- c("Methylation group", round(exp(coef(res.ct)[1]), 4), 
         round(exp(res.ct$CI.tbl)[1,2:3], 4), 
         round(res.ct.test2$pvalue, 4))
x46 <- c("CTNNB1 Mutation Category Wildtype:Somatic", round(exp(coef(res.ct)[2]), 4), 
         round(exp(res.ct$CI.tbl)[2,2:3], 4), 
         round(res.ct.test1$pvalue, 4))
act.tab <- rbind.data.frame(act.tab, x44, x45, x46)
# Adjust for X11p15 status
x47 <- c("Adjust for x11p15 status", " ", " ", " ", "")
res.x11 <- crrf(Surv(RFS, RFS_status)~Methylation_Group + X11p15_status,
                etype="Relapse", dset=act.dat, firth=TRUE, CI=TRUE)
res.x11.test1 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + X11p15_status, 
                         etype="Relapse", dset=act.dat, 
                         test=~Methylation_Group, firth=TRUE)
res.x11.test2 <- crrftest(Surv(RFS, RFS_status)~Methylation_Group + X11p15_status, 
                          etype="Relapse", dset=act.dat, 
                          test=~X11p15_status, firth=TRUE)
x48 <- c("Methylation group", round(exp(coef(res.x11)[1]), 4), 
         round(exp(res.x11$CI.tbl)[1,2:3], 4), 
         round(res.x11.test2$pvalue, 4))
x49 <- c("x11p15 Status Wildtype:LOH", round(exp(coef(res.x11)[2]), 4), 
         round(exp(res.x11$CI.tbl)[2,2:3], 4), 
         round(res.x11.test1$pvalue, 4))
act.tab <- rbind.data.frame(act.tab, x47, x48, x49)
save(act.tab, file=paste0(res.dir, "ACT_Multi_Table.RData"))
