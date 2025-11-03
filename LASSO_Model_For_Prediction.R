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

set.seed(123)
outer_folds_by_outcome <- list()

for (j in 1:9) {
  outcome_col <- paste0("outcome", j)
  outer_folds <- createFolds(
    y = factor(data1[[outcome_col]]),
    k = 10,
    list = TRUE,
    returnTrain = FALSE
  )
  outer_folds_by_outcome[[j]] <- outer_folds
}

num_cores <- 20   
cl <- makeCluster(num_cores)
registerDoParallel(cl)
outer_folds <- createFolds(ids, k = 10, list = TRUE, returnTrain = FALSE)
all_predictions <- list()
beta_matrix <- list() 

results <- foreach(i = 1:10, .packages = c("glmnet", "survival", "dplyr", "Matrix", "caret")) %dopar% {
  tryCatch({
    cat("Processing outer fold", i, "of 10\n")
    
    coef <- list()
    scores <- data.frame(ID = data1$n_eid)
    
    for (j in 1:9) {
      folds <- outer_folds_by_outcome[[j]]
      
      test_indices <- folds[[i]]
      train_indices <- setdiff(1:nrow(data1), test_indices)
      
      train.data <- data1[train_indices, ]
      test.data  <- data1[test_indices, ]
      
      selected_protein_ids <- significant_proteins_by_outcome[[j]]
      selected_cols <- colnames(data1[, 44:2963]) %in% selected_protein_ids
      
      x <- as.matrix(train.data[, 44:2963][, selected_cols, drop = FALSE])
      y <- eval(parse(text = paste0("Surv(train.data$time", j, "to, train.data$outcome", j, ")")))
      
      set.seed(123)
      cv <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)
      
      # 集成3个lambda值
      lambda_values <- c(
        cv$lambda.min,
        exp((log(cv$lambda.min) + log(cv$lambda.1se)) / 2),
        cv$lambda.1se
      )
      
      coef_ensemble <- matrix(0, nrow = ncol(x), ncol = length(lambda_values))
      rownames(coef_ensemble) <- colnames(x)   # ✅ 保留蛋白名
      
      for (k in 1:length(lambda_values)) {
        model <- glmnet(x, y, family = "cox", alpha = 1, lambda = lambda_values[k])
        coef_temp <- as.matrix(coef(model, s = lambda_values[k]))
        if (rownames(coef_temp)[1] == "(Intercept)") {
          coef_temp <- coef_temp[-1, , drop = FALSE]
        }
        coef_ensemble[, k] <- coef_temp[, 1]
      }
      
      final_coef <- rowMeans(coef_ensemble)
      coef[[j]] <- final_coef   # ✅ 现在是带蛋白名的 named vector
      
      Proteset <- test.data[, 44:2963][, selected_cols, drop = FALSE]
      
      scores$ID[test_indices] <- test.data$n_eid
      calculate_protscore <- function(coefProte, Proteset) {
        rowSums(sweep(Proteset, 2, coefProte, `*`))
      }
      scores[test_indices, paste0("Protescore", j)] <- calculate_protscore(coef[[j]], Proteset)
      cat("Finished outcome", j, "in fold", i, "\n")
    }
    
    return(list(coef = coef, scores = scores))
    
  }, error = function(e) {
    cat("Error in fold", i, ":", conditionMessage(e), "\n")
    NULL
  })
}

all_predictions <- list()
coef_list <- list()

for (i in 1:10) {
  if (!is.null(results[[i]])) {
    all_predictions[[i]] <- results[[i]]$scores
    coef_list[[i]] <- results[[i]]$coef
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



