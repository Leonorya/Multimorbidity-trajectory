# ============================================================
# LASSO for prediction
# Author: Guangruiyang
# Last Modified: Feb 2025
# ============================================================
library(readxl)
library(dplyr)
library(tidyr)
library(glmnet)
library(caret)
library(survival)
library(parallel)
library(doParallel)
library(openxlsx)
library(CsChange)
library(survIDINRI)
data1 <-read.csv("/ProtePattern_CMD_DEP_mice.csv")
set.seed(726)

data1$outcome1<-0
data1$outcome1[data1$trans1==1&data1$trans4==0&data1$trans5==0]=1
data1$outcome2<-0
data1$outcome2[data1$trans2==1&data1$trans6==0&data1$trans7==0]=1
data1$outcome3<-0
data1$outcome3[data1$trans3==1]=1
data1$outcome4<-0
data1$outcome4[data1$trans4==1&data1$trans8==0]=1
data1$outcome5<-0
data1$outcome5[data1$trans5==1]=1
data1$outcome6<-0
data1$outcome6[data1$trans6==1&data1$trans9==0]=1
data1$outcome7<-0
data1$outcome7[data1$trans7==1]=1
data1$outcome8<-0
data1$outcome8[data1$trans8==1]=1
data1$outcome9<-0
data1$outcome9[data1$trans9==1]=1

data1[, 46:2965] <- as.data.frame(scale(data1[, 46:2965]))
proteomic_cmd_dep <- read.csv("/ProteomicCMDDEP.csv")

significant_proteins_by_outcome <- list()
for (i in 1:9) {
  fdr_col <- paste0("T", i, "FDR")
  significant_proteins_by_outcome[[i]] <- proteomic_cmd_dep %>%
    filter(!!sym(fdr_col) < 0.05) %>%
    pull(ID)
}

ids <- data1[, 1]
num_cores <- 20   
cl <- makeCluster(num_cores)
registerDoParallel(cl)
outer_folds <- createFolds(ids, k = 10, list = TRUE, returnTrain = FALSE)
all_predictions <- list()
beta_matrix <- list() 

results <- foreach(i = 1:10, .packages = c("glmnet", "survival", "caret")) %dopar% {
  tryCatch({
    test_indices <- outer_folds[[i]]
    train_indices <- setdiff(1:nrow(data1), test_indices)
    
    train.data <- data1[train_indices, ]
    test.data <- data1[test_indices, ]
    
    coef <- list()
    scores <- data.frame(ID = ids[test_indices])
    
    for (j in 1:9) {
      selected_protein_ids <- significant_proteins_by_outcome[[j]]
      selected_cols <- colnames(data1[, 46:2965]) %in% selected_protein_ids
      x <- as.matrix(train.data[, 46:2965][, selected_cols])
      y <- eval(parse(text = paste0("Surv(train.data$time", j, "to, train.data$outcome", j, ")")))
      
      set.seed(123)
      cv <- cv.glmnet(x, y, family = 'cox', alpha = 1, nfolds = 10)
      model <- glmnet(x, y, family = 'cox', alpha = 1, lambda = cv$lambda.min)
      coef[[j]] <- as.vector(coefficients(model))
      
      Proteset <- test.data[, 46:2965][, selected_cols]
      calculate_protscore <- function(coefProte, Proteset) {
        rowSums(as.data.frame(lapply(1:ncol(Proteset), function(i) Proteset[, i] * coefProte[i])))
      }
      scores[[paste0("Protescore", j)]] <- calculate_protscore(coef[[j]], Proteset)
    }
    
    list(coef = coef, scores = scores)
    
  }, error = function(e) {
    cat("Error in fold", i, ":", conditionMessage(e), "\n")
    return(NULL)
  })
}

coef_list <- vector("list", 10)
all_predictions <- vector("list", 10)
for (i in 1:10) {
  if (!is.null(results[[i]])) {
    coef_list[[i]] <- results[[i]]$coef
    all_predictions[[i]] <- results[[i]]$scores
  }
}

coef_by_outcome <- list()
for (j in 1:9) {
  outcome_coefs <- do.call(cbind, lapply(coef_list, function(l) l[[j]]))
  colnames(outcome_coefs) <- paste0("Fold", 1:10)
  coef_by_outcome[[j]] <- outcome_coefs
}

wb <- createWorkbook()
for (j in 1:9) {
  addWorksheet(wb, sheetName = paste0("Outcome", j))
  writeData(wb, sheet = j, x = coef_by_outcome[[j]], rowNames = TRUE)
  setColWidths(wb, sheet = j, cols = 1:(ncol(coef_by_outcome[[j]])+1), widths = "auto")
}
saveWorkbook(wb, file = "/coef_Protescore9.xlsx", overwrite = TRUE)

all_predictions_df <- do.call(rbind, all_predictions)
id_protescore_df <- all_predictions_df[, c("ID", grep("Protescore", colnames(all_predictions_df), value = TRUE))]
write.csv(id_protescore_df, "/Protecore.csv", row.names = FALSE)


