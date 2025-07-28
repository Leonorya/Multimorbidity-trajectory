# ============================================================
# Omics and CMD-Dep multimorbidity associations
# Author: Guangruiyang
# Last Modified: Feb 2025
# ============================================================
library("mstate")
library(dplyr)
library(mice)
data1 <- read.csv("/ProtePattern_CMD_DEP.csv")
data1[, 46:2965] <- as.data.frame(scale(data1[, 46:2965]))
covariate <- c("Mets", "Ethnicity", "Sleeptime", "Edu", "Employed", "Smoke",
               "Drink", "T2dhis", "Cvdhis", "Tdi", "hisdepre","Dietscore","Age","Sex","BMI")
data_to_impute <- data1[, covariate]
imputed_data <- mice(data_to_impute, method = 'pmm', m = 10, seed = 123)
completed_data <- complete(imputed_data)
data1[, covariate] <- completed_data
data1$tdi <- cut(data1$Tdi, breaks = quantile(data1$Tdi, probs = seq(0, 1, 0.2)), labels = FALSE)

trans_columns <- c("trans1", "trans2", "trans3", "trans4", "trans5","trans6","trans7","trans8","trans9")
trans_sums <- sapply(trans_columns, function(col) sum(data1[[col]] == 1))
trans_sums

datatrans1<-data1[,c(1,2,3,4,29,30,31,32,33,34,35,36,37,39,41,42,43,44,2966,46:2965)]
datatrans2<-data1[,c(1,5,6,7,29,30,31,32,33,34,35,36,37,39,41,42,43,44,2966,46:2965)]
datatrans3<-data1[,c(1,8,9,10,29,30,31,32,33,34,35,36,37,39,41,42,43,44,2966,46:2965)]
datatrans4<-data1[,c(1,11,12,13,29,30,31,32,33,34,35,36,37,39,41,42,43,44,2966,46:2965)]
datatrans5<-data1[,c(1,14,15,16,29,30,31,32,33,34,35,36,37,39,41,42,43,44,2966,46:2965)]
datatrans6<-data1[,c(1,17,18,19,29,30,31,32,33,34,35,36,37,39,41,42,43,44,2966,46:2965)]
datatrans7<-data1[,c(1,20,21,22,29,30,31,32,33,34,35,36,37,39,41,42,43,44,2966,46:2965)]
datatrans8<-data1[,c(1,23,24,25,29,30,31,32,33,34,35,36,37,39,41,42,43,44,2966,46:2965)]
datatrans9<-data1[,c(1,26,27,28,29,30,31,32,33,34,35,36,37,39,41,42,43,44,2966,46:2965)]
datatrans1$transpattern<-1
datatrans2$transpattern<-2
datatrans3$transpattern<-3
datatrans4$transpattern<-4
datatrans5$transpattern<-5
datatrans6$transpattern<-6
datatrans7$transpattern<-7
datatrans8$transpattern<-8
datatrans9$transpattern<-9
dataall1<-rbind(datatrans1,datatrans2,datatrans3,datatrans4,datatrans5,datatrans6,datatrans7,datatrans8,datatrans9)
dataall1$trans<-dataall1$failcode
dataall1 <- dataall1 %>%
  mutate_at(vars(Ethnicity, Edu, Employed,Smoke, Drink, T2dhis, Cvdhis,hisdepre,Sex,tdi), as.factor)
covs_P <-   c("Mets", "Ethnicity", "Sleeptime", "Edu", "Employed",
              "Smoke", "Drink", "T2dhis", "Cvdhis","hisdepre","Dietscore", "Age", "Sex", "BMI","tdi")
protein <- dput(names(data1[46:2965]))
covs_allname <- c(covs_P, protein)
msebmt_allname <- expand.covs(dataall1, covs_allname, append = TRUE,longnames = FALSE)

results <- data.frame(metabolite = character(),
                      HR = character(),
                      P_value = numeric(),
                      stringsAsFactors = FALSE)

count <- 0  
save_count <- 1  
for (met in metabolism) {
  tryCatch({
    count <- count + 1  
    if (count %% 100 == 0) {
      message(sprintf("Finishe %d Proteins，Undergoing %s", count, met))
    }
    formula_str <- paste("Surv(timef, timeto, status) ~", 
                         paste(met, ".1 +", met, ".2 +", met, ".3 +", met, ".4 +", met, ".5 +", met, ".6 +", met, ".7 +", met, ".8 +", met, ".9 ", sep = ""),
                         "+ Sex.1 + Sex.2 + Sex.3 + Sex.4 + Sex.5 + Sex.6 + Sex.7 + Sex.8 + Sex.9+" ,
                         "tdi1.1 + tdi1.2 + tdi1.3 + tdi1.4 + tdi1.5 + tdi1.6 + tdi1.7 + tdi1.8 + tdi1.9+",
                         "tdi2.1 + tdi2.2 + tdi2.3 + tdi2.4 + tdi2.5 + tdi2.6 + tdi2.7 + tdi2.8 + tdi2.9+",
                         "tdi3.1 + tdi3.2 + tdi3.3 + tdi3.4 + tdi3.5 + tdi3.6 + tdi3.7 + tdi3.8 + tdi3.9 +",
                         "BMI.1 + BMI.2 + BMI.3 + BMI.4 + BMI.5 + BMI.6 + BMI.7 + BMI.8 + BMI.9 +",
                         "Age.1 + Age.2 + Age.3 + Age.4 + Age.5 + Age.6 + Age.7 + Age.8 + Age.9 + ",
                         "Ethnicity.1 + Ethnicity.2 + Ethnicity.3 + Ethnicity.4 + Ethnicity.5 + Ethnicity.6 + Ethnicity.7 + Ethnicity.8 + Ethnicity.9 +",
                         "Edu1.1 + Edu1.2 + Edu1.3 + Edu1.4 + Edu1.5 + Edu1.6 + Edu1.7 + Edu1.8 + Edu1.9 + ",
                         "Edu2.1 + Edu2.2 + Edu2.3 + Edu2.4 + Edu2.5 + Edu2.6 + Edu2.7 + Edu2.8 + Edu2.9 +",
                         "Edu3.1 + Edu3.2 + Edu3.3 + Edu3.4 + Edu3.5 + Edu3.6 + Edu3.7 + Edu3.8 + Edu3.9 +", 
                         "Employed.1 + Employed.2 + Employed.3 + Employed.4 + Employed.5 + Employed.6 + Employed.7 + Employed.8 + Employed.9 +",
                         "Smoke1.1 + Smoke1.2 + Smoke1.3 + Smoke1.4 + Smoke1.5 + Smoke1.6 + Smoke1.7 + Smoke1.8 + Smoke1.9 +",
                         "Smoke2.1 + Smoke2.2 + Smoke2.3 + Smoke2.4 + Smoke2.5 + Smoke2.6 + Smoke2.7 + Smoke2.8 + Smoke2.9 +",
                         "Drink1.1 + Drink1.2 + Drink1.3 + Drink1.4 + Drink1.5 + Drink1.6 + Drink1.7 + Drink1.8 + Drink1.9 +", 
                         "Drink2.1 + Drink2.2 + Drink2.3 + Drink2.4 + Drink2.5 + Drink2.6 + Drink2.7 + Drink2.8 + Drink2.9 + ",
                         "Mets.1 + Mets.2 + Mets.3 + Mets.4 + Mets.5 + Mets.6 + Mets.7 + Mets.8 + Mets.9 +",
                         "Sleeptime.1 + Sleeptime.2 + Sleeptime.3 + Sleeptime.4 + Sleeptime.5 + Sleeptime.6 + Sleeptime.7 + Sleeptime.8 + Sleeptime.9 +",
                         "Dietscore.1 + Dietscore.2 + Dietscore.3 + Dietscore.4 + Dietscore.5 + Dietscore.6 + Dietscore.7 + Dietscore.8 + Dietscore.9 +",
                         "T2dhis.1 + T2dhis.2 + T2dhis.3 + T2dhis.4 + T2dhis.5  + T2dhis.6 + T2dhis.7 + T2dhis.8 + T2dhis.9 +", 
                         "Cvdhis.1 + Cvdhis.2 + Cvdhis.3 + Cvdhis.4 + Cvdhis.5 + Cvdhis.6 + Cvdhis.7 + Cvdhis.8 + Cvdhis.9 +",
                         "hisdepre.1 + hisdepre.2+ hisdepre.3+ hisdepre.4+ hisdepre.5+ hisdepre.6+ hisdepre.7+ hisdepre.8+ hisdepre.9+ strata(trans)", 
                         sep = "")
    
    model_formula <- as.formula(formula_str)
    model <- coxph(model_formula, data = msebmt_allname, method = 'breslow')
    
    for (i in 1:9) {
      met_var <- paste(met, ".", i, sep = "")
      if (met_var %in% rownames(summary(model)$coefficients)) {
        HR <- paste0(round(summary(model)$conf.int[met_var, "exp(coef)"], 2), " (",
                     round(summary(model)$conf.int[met_var, "lower .95"], 2), "-",
                     round(summary(model)$conf.int[met_var, "upper .95"], 2), ")")
        P_value <- sprintf("%.2E", summary(model)$coefficients[met_var, "Pr(>|z|)"])
        results <- rbind(results, data.frame(metabolite = met_var, HR = HR, P_value = P_value, stringsAsFactors = FALSE))
      }
    }
    
    if (count %% 200 == 0) {
      filename <- paste0("/dssg/home/acct-wenze.zhong/gryang/Transition/CMDDEP/ProtMice_out200_", save_count, ".csv")
      write.csv(results, file = filename, row.names = FALSE)
      save_count <- save_count + 1
    }
    
  }, error = function(e) {
    message(paste("Error in processing", met, ":", e$message))
  })
}

results$Fdr_pvalue <- sprintf("%.2E", p.adjust(results$P_value, method = "fdr"))
results$Bonferroni_pvalue <- sprintf("%.2E", p.adjust(results$P_value, method = "bonferroni"))
write.csv(results, file = "/proteomicsign_CMD_DEP1_mice.csv", row.names = FALSE)


