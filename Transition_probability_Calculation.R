# ============================================================
# UK Biobank Disease Transition Probability Analysis
# Author: JXW
# Last Modified: Feb 2025
# 
# This script analyzes transition probabilities between:
# 1. Baseline to first disease (CMD/DEP/Death)
# 2. CMD to multimorbidity or mortality
# 3. DEP to multimorbidity or mortality 
# 4. CMD-DEP multimorbidity prognosis
# 5. DEP-CMD multimorbidity prognosis
# ============================================================

### 1. Baseline to 1st disease-------------------------------
library(dplyr)
track<-read.csv("CMD_DEP.csv")
covariate<-read.csv("covariates_impute.csv")
data<-merge(track,covariate,by="n_eid")
datatrans1<-data[,c(1,2,3,4,30:47)]
datatrans2<-data[,c(1,5,6,7,30:47)]
datatrans3<-data[,c(1,8,9,10,30:47)]
datatrans1$transpattern<-1
datatrans2$transpattern<-2
datatrans3$transpattern<-3

dataall1$timeto<-dataall1$timeto-dataall1$timef
dataall1$timef<-0
library("mstate")
tmat <- transMat(x = list(c(2,3,4), c(), c(),c()),
                 names = c( "Baseline","CMD","DEP","Death"))

covs_P <- c("Mets", "Ethnicity", "Sleeptime", "Edu", "Employed",
            "Smoke", "Drink", "T2dhis", "Cvdhis","hisdepre","Dietscore", "Age", "Sex", "BMI","tdi")
msebmt_allname <- expand.covs(dataall1,  covs_P , append = TRUE,longnames = FALSE)
variables <- colnames(msebmt_allname[,24:77])
formula_str <- paste("Surv(timef, timeto, status) ~", paste(variables, collapse = " + "), "+ strata(failcode)")
formula <- as.formula(formula_str)
the.expression_allname <- coxph(formula, data = msebmt_allname, method = 'breslow')
mean_values <- colMeans(data[, c("Mets", "Ethnicity", "Sleeptime", "Edu", "Employed",
                                 "Smoke", "Drink", "T2dhis", "Cvdhis","hisdepre","Dietscore", "Age", "Sex", "BMI","tdi")], na.rm = TRUE)
newd <- data.frame(
  Age = rep(mean_values["Age"], 3),
  Sex = rep(mean_values["Sex"], 3),
  Ethnicity = rep(mean_values["Ethnicity"], 3),
  Edu = rep(mean_values["Edu"], 3),
  Employed = rep(mean_values["Employed"], 3),
  tdi = rep(mean_values["tdi"], 3),
  Smoke = rep(mean_values["Smoke"], 3),
  Drink = rep(mean_values["Drink"], 3),
  Mets = rep(mean_values["Mets"], 3),
  Dietscore = rep(mean_values["Dietscore"], 3),
  Sleeptime = rep(mean_values["Sleeptime"], 3),
  BMI = rep(mean_values["BMI"], 3),
  T2dhis = rep(mean_values["T2dhis"], 3),
  Cvdhis = rep(mean_values["Cvdhis"], 3),
  hisdepre = rep(mean_values["hisdepre"], 3),
  trans = 1:3
)

attr(newd, "trans") <- tmat
class(newd) <- c("msdata", "data.frame")
covs_P1 <- c("Mets", "Ethnicity", "Sleeptime", "Edu", "Employed",
             "Smoke", "Drink", "T2dhis", "Cvdhis","hisdepre","Dietscore", "Age", "Sex", "BMI","tdi")
newd <- expand.covs(newd, covs_P1, longnames = FALSE)
newd$strata<-1:3
msf.WW<-msfit(object = the.expression_allname ,newdata = newd,trans=tmat)
plot(msf.WW)
pt <- probtrans(msf.WW, predt = 0,)
pc<-summary(pt,from=1)

prob1_5<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$upper2),3),")")
prob1_10<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$upper2),3),")")
prob1_15<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$upper2),3),")")
prob2_5<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$pstate3),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$lower3),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$upper3),3),")")
prob2_10<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$pstate3),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$lower3),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$upper3),3),")")
prob2_15<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$pstate3),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$lower3),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$upper3),3),")")
prob3_5<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$pstate4),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$lower4),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$upper4),3),")")
prob3_10<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$pstate4),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$lower4),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$upper4),3),")")
prob3_15<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$pstate4),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$lower4),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$upper4),3),")")
se1_5<-round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$se2),3)
se1_10<-round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$se2),3)
se1_15<-round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$se2),3)
se2_5<-round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$se3),3)
se2_10<-round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$se3),3)
se2_15<-round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$se3),3)
se3_5<-round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$se4),3)
se3_10<-round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$se4),3)
se3_15<-round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$se4),3)

### 2. CMD to multimorbidity or mortality-------------------------------
track<-read.csv("CMD_Dep.csv")
covariate<-read.csv("covariates_impute.csv")
data<-merge(track,covariate,by="n_eid")
values<- data %>% filter (trans1==1|trans2==1)
mean_values <- colMeans(values[, c("Mets", "Ethnicity", "Sleeptime", "Edu", "Employed",
                                   "Smoke", "Drink", "T2dhis", "Cvdhis","hisdepre","Dietscore", "Age", "Sex", "BMI","tdi")], na.rm = TRUE)
data<-data %>% filter (trans1==1)
datatrans4<-data[,c(1,11,12,13,30:47)]
datatrans5<-data[,c(1,14,15,16,30:47)]
datatrans4$transpattern<-1
datatrans5$transpattern<-2
dataall1<-rbind(datatrans4,datatrans5)

dataall1$timeto<-dataall1$timeto-dataall1$timef
dataall1$timef<-0
tmat <- transMat(x = list(c(2,3), c(), c()),
                 names = c( "CMD","CMD-Dep","Death"))
covs_P <- c("Mets", "Ethnicity", "Sleeptime", "Edu", "Employed",
            "Smoke", "Drink", "T2dhis", "Cvdhis","hisdepre","Dietscore", "Age", "Sex", "BMI","tdi")

msebmt_allname <- expand.covs(dataall1,  covs_P , append = TRUE,longnames = FALSE)
variables <- colnames(msebmt_allname[,24:59])
formula_str <- paste("Surv(timef, timeto, status) ~", paste(variables, collapse = " + "), "+ strata(failcode)")
formula <- as.formula(formula_str)
the.expression_allname <- coxph(formula, data = msebmt_allname, method = 'breslow')
newd <- data.frame(
  Age = rep(mean_values["Age"], 3),
  Sex = rep(mean_values["Sex"], 3),
  Ethnicity = rep(mean_values["Ethnicity"], 3),
  Edu = rep(mean_values["Edu"], 3),
  Employed = rep(mean_values["Employed"], 3),
  tdi = rep(mean_values["tdi"], 3),
  Smoke = rep(mean_values["Smoke"], 3),
  Drink = rep(mean_values["Drink"], 3),
  Mets = rep(mean_values["Mets"], 3),
  Dietscore = rep(mean_values["Dietscore"], 3),
  Sleeptime = rep(mean_values["Sleeptime"], 3),
  BMI = rep(mean_values["BMI"], 3),
  T2dhis = rep(mean_values["T2dhis"], 3),
  Cvdhis = rep(mean_values["Cvdhis"], 3),
  hisdepre = rep(mean_values["hisdepre"], 3),
  trans = 1:2
)
attr(newd, "trans") <- tmat
class(newd) <- c("msdata", "data.frame")
covs_P1 <- c("Mets", "Ethnicity", "Sleeptime", "Edu", "Employed",
             "Smoke", "Drink", "T2dhis", "Cvdhis","hisdepre","Dietscore", "Age", "Sex", "BMI","tdi")
newd <- expand.covs(newd, covs_P1, longnames = FALSE)
newd$strata<-1:2
msf.WW<-msfit(object = the.expression_allname ,newdata = newd,trans=tmat)
plot(msf.WW)
pt <- probtrans(msf.WW, predt = 0,)
pc<-summary(pt,from=1)
prob4_5<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$upper2),3),")")
prob4_10<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$upper2),3),")")
prob4_15<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$upper2),3),")")
prob5_5<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$pstate3),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$lower3),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$upper3),3),")")
prob5_10<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$pstate3),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$lower3),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$upper3),3),")")
prob5_15<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$pstate3),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$lower3),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$upper3),3),")")
se4_5<-round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$se2),3)
se4_10<-round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$se2),3)
se4_15<-round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$se2),3)
se5_5<-round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$se3),3)
se5_10<-round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$se3),3)
se5_15<-round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$se3),3)

### 3. Dep to multimorbidity or mortality-------------------------------
track<-read.csv("CMD_Dep.csv")
covariate<-read.csv("covariates_impute.csv")
data<-merge(track,covariate,by="n_eid")
data<-data %>% filter (trans2==1)
datatrans6<-data[,c(1,17,18,19,30:47)]
datatrans7<-data[,c(1,20,21,22,30:47)]
datatrans6$transpattern<-1
datatrans7$transpattern<-2
dataall1<-rbind(datatrans6,datatrans7)
dataall1$timeto<-dataall1$timeto-dataall1$timef
dataall1$timef<-0
tmat <- transMat(x = list(c(2,3), c(), c()),
                 names = c( "Dep","Dep-CMD","Death"))
covs_P <- c("Mets", "Ethnicity", "Sleeptime", "Edu", "Employed",
            "Smoke", "Drink", "T2dhis", "Cvdhis","hisdepre","Dietscore", "Age", "Sex", "BMI","tdi")
msebmt_allname <- expand.covs(dataall1,  covs_P , append = TRUE,longnames = FALSE)
variables <- colnames(msebmt_allname[,24:59])
formula_str <- paste("Surv(timef, timeto, status) ~", paste(variables, collapse = " + "), "+ strata(failcode)")
formula <- as.formula(formula_str)
the.expression_allname <- coxph(formula, data = msebmt_allname, method = 'breslow')
newd <- data.frame(
  Age = rep(mean_values["Age"], 3),
  Sex = rep(mean_values["Sex"], 3),
  Ethnicity = rep(mean_values["Ethnicity"], 3),
  Edu = rep(mean_values["Edu"], 3),
  Employed = rep(mean_values["Employed"], 3),
  tdi = rep(mean_values["tdi"], 3),
  Smoke = rep(mean_values["Smoke"], 3),
  Drink = rep(mean_values["Drink"], 3),
  Mets = rep(mean_values["Mets"], 3),
  Dietscore = rep(mean_values["Dietscore"], 3),
  Sleeptime = rep(mean_values["Sleeptime"], 3),
  BMI = rep(mean_values["BMI"], 3),
  T2dhis = rep(mean_values["T2dhis"], 3),
  Cvdhis = rep(mean_values["Cvdhis"], 3),
  hisdepre = rep(mean_values["hisdepre"], 3),
  trans = 1:2
)

attr(newd, "trans") <- tmat
class(newd) <- c("msdata", "data.frame")
newd <- expand.covs(newd, covs_P1, longnames = FALSE)

newd$strata<-1:2
msf.WW<-msfit(object = the.expression_allname ,newdata = newd,trans=tmat)
plot(msf.WW)
pt <- probtrans(msf.WW, predt = 0,)
pc<-summary(pt,from=1)
prob6_5<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$upper2),3),")")
prob6_10<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$upper2),3),")")
prob6_15<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$upper2),3),")")
prob7_5<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$pstate3),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$lower3),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$upper3),3),")")
prob7_10<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$pstate3),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$lower3),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$upper3),3),")")
prob7_15<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$pstate3),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$lower3),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$upper3),3),")")
se6_5<-round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$se2),3)
se6_10<-round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$se2),3)
se6_15<-round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$se2),3)
se7_5<-round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$se3),3)
se7_10<-round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$se3),3)
se7_15<-round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$se3),3)

### 4. CMD-Dep multimorbidity prognosis-------------------------------
track<-read.csv("CMD_Dep.csv")
covariate<-read.csv("covariates_impute.csv")
data<-merge(track,covariate,by="n_eid")
values<- data %>% filter (trans4==1|trans6==1)
mean_values <- colMeans(values[, c("Mets", "Ethnicity", "Sleeptime", "Edu", "Employed",
                                   "Smoke", "Drink", "T2dhis", "Cvdhis","hisdepre","Dietscore", "Age", "Sex", "BMI","tdi")], na.rm = TRUE)
data <-data %>% filter (trans4==1)
data$trans10[data$trans8==1]=0
data$trans10[data$trans8==0]=1
datatrans8<-data[,c(1,23,24,25,30:47)]
datatrans10<-data[,c(1,58,24,25,30:47)]
datatrans8$transpattern<-1
datatrans10$transpattern<-2
dataall1<-rbind(datatrans8,datatrans10)
dataall1$timeto<-dataall1$timeto-dataall1$timef
dataall1$timef<-0
tmat <- transMat(x = list(c(2,3), c(), c()),
                 names = c( "CMD-Dep","Death","Non-Death"))
covs_P <- c("Mets", "Ethnicity", "Sleeptime", "Edu", "Employed",
            "Smoke", "Drink", "T2dhis", "Cvdhis","hisdepre","Dietscore", "Age", "Sex", "BMI","tdi")
msebmt_allname <- expand.covs(dataall1,  covs_P , append = TRUE,longnames = FALSE)
variables <- colnames(msebmt_allname[,24:59])
formula_str <- paste("Surv(timef, timeto, status) ~", paste(variables, collapse = " + "), "+ strata(failcode)")
formula <- as.formula(formula_str)
the.expression_allname <- coxph(formula, data = msebmt_allname, method = 'breslow')
newd <- data.frame(
  Age = rep(mean_values["Age"], 3),
  Sex = rep(mean_values["Sex"], 3),
  Ethnicity = rep(mean_values["Ethnicity"], 3),
  Edu = rep(mean_values["Edu"], 3),
  Employed = rep(mean_values["Employed"], 3),
  tdi = rep(mean_values["tdi"], 3),
  Smoke = rep(mean_values["Smoke"], 3),
  Drink = rep(mean_values["Drink"], 3),
  Mets = rep(mean_values["Mets"], 3),
  Dietscore = rep(mean_values["Dietscore"], 3),
  Sleeptime = rep(mean_values["Sleeptime"], 3),
  BMI = rep(mean_values["BMI"], 3),
  T2dhis = rep(mean_values["T2dhis"], 3),
  Cvdhis = rep(mean_values["Cvdhis"], 3),
  hisdepre = rep(mean_values["hisdepre"], 3),
  trans = 1:2
)

attr(newd, "trans") <- tmat
class(newd) <- c("msdata", "data.frame")
covs_P1 <- c("Mets", "Ethnicity", "Sleeptime", "Edu", "Employed",
             "Smoke", "Drink", "T2dhis", "Cvdhis","hisdepre","Dietscore", "Age", "Sex", "BMI","tdi")
newd <- expand.covs(newd, covs_P1, longnames = FALSE)
newd$strata<-1:2
msf.WW<-msfit(object = the.expression_allname ,newdata = newd,trans=tmat)
plot(msf.WW)
pt <- probtrans(msf.WW, predt = 0,)
pc<-summary(pt,from=1)
prob8_5<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$upper2),3),")")
prob8_10<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$upper2),3),")")
prob8_15<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$upper2),3),")")
se8_5<-round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$se2),3)
se8_10<-round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$se2),3)
se8_15<-round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$se2),3)

### 5. Dep-CMD multimorbidity prognosis-------------------------------
track<-read.csv("CMD_Dep.csv")
covariate<-read.csv("covariates_impute.csv")
data<-merge(track,covariate,by="n_eid")
data <-data %>% filter (trans6==1)
data$trans11[data$trans9==1]=0
data$trans11[data$trans9==0]=1
datatrans9<-data[,c(1,26,27,28,30:47)]
datatrans11<-data[,c(1,58,27,28,30:47)]
datatrans9$transpattern<-1
datatrans11$transpattern<-2
dataall1<-rbind(datatrans9,datatrans11)
dataall1$timeto<-dataall1$timeto-dataall1$timef
dataall1$timef<-0
tmat <- transMat(x = list(c(2,3), c(), c()),
                 names = c( "Dep-CMD","Death","Non-Death"))
covs_P <- c("Mets", "Ethnicity", "Sleeptime", "Edu", "Employed",
            "Smoke", "Drink", "T2dhis", "Cvdhis","hisdepre","Dietscore", "Age", "Sex", "BMI","tdi")
msebmt_allname <- expand.covs(dataall1,  covs_P , append = TRUE,longnames = FALSE)
variables <- colnames(msebmt_allname[,24:59])
formula_str <- paste("Surv(timef, timeto, status) ~", paste(variables, collapse = " + "), "+ strata(failcode)")
formula <- as.formula(formula_str)
the.expression_allname <- coxph(formula, data = msebmt_allname, method = 'breslow')
newd <- data.frame(
  Age = rep(mean_values["Age"], 3),
  Sex = rep(mean_values["Sex"], 3),
  Ethnicity = rep(mean_values["Ethnicity"], 3),
  Edu = rep(mean_values["Edu"], 3),
  Employed = rep(mean_values["Employed"], 3),
  tdi = rep(mean_values["tdi"], 3),
  Smoke = rep(mean_values["Smoke"], 3),
  Drink = rep(mean_values["Drink"], 3),
  Mets = rep(mean_values["Mets"], 3),
  Dietscore = rep(mean_values["Dietscore"], 3),
  Sleeptime = rep(mean_values["Sleeptime"], 3),
  BMI = rep(mean_values["BMI"], 3),
  T2dhis = rep(mean_values["T2dhis"], 3),
  Cvdhis = rep(mean_values["Cvdhis"], 3),
  hisdepre = rep(mean_values["hisdepre"], 3),
  trans = 1:2
)
attr(newd, "trans") <- tmat
class(newd) <- c("msdata", "data.frame")
covs_P1 <- c("Mets", "Ethnicity", "Sleeptime", "Edu", "Employed",
             "Smoke", "Drink", "T2dhis", "Cvdhis","hisdepre","Dietscore", "Age", "Sex", "BMI","tdi")
newd <- expand.covs(newd, covs_P1, longnames = FALSE)
newd$strata<-1:2
msf.WW<-msfit(object = the.expression_allname ,newdata = newd,trans=tmat)
plot(msf.WW)
pt <- probtrans(msf.WW, predt = 0,)
pc<-summary(pt,from=1)
prob9_5<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$upper2),3),")")
prob9_10<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$upper2),3),")")
prob9_15<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$upper2),3),")")
se9_5<-round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$se2),3)
se9_10<-round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$se2),3)
se9_15<-round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$se2),3)

proball<-rbind(prob1_5,prob1_10,prob1_15,prob2_5,prob2_10,prob2_15,prob3_5,prob3_10,prob3_15,
               prob4_5,prob4_10,prob4_15,prob5_5,prob5_10,prob5_15,prob6_5,prob6_10,prob6_15,
               prob7_5,prob7_10,prob7_15,prob8_5,prob8_10,prob8_15,prob9_5,prob9_10,prob9_15,
               se1_5,se1_10,se1_15,se2_5,se2_10,se2_15,
               se3_5,se3_10,se3_15,se4_5,se4_10,se4_15,
               se5_5,se5_10,se5_15,se6_5,se6_10,se6_15,
               se7_5,se7_10,se7_15,se8_5,se8_10,se8_15,
               se9_5,se9_10,se9_15
)

write.csv(proball,"all_trans.csv")