# ============================================================
# KEGG pathways
# Author: Guangruiyang
# Last Modified: Feb 2025
# ============================================================
library(clusterProfiler) 
library(org.Hs.eg.db)    
library(openxlsx)        

protein <- read.csv("Prot_significant.csv")
wb <- createWorkbook()
addWorksheet(wb, sheetName = "Summary")

summary_data <- data.frame(
  Variable = character(),
  Input_Proteins = integer(),
  Mapped_Proteins = integer(),
  GO_Pathways = integer(),
  KEGG_Pathways = integer(),
  stringsAsFactors = FALSE
)

variables <- c("id1", "id2", "id3", "id4", "id5", "id6", "id7", "id8", "id9")

for (variable in variables) {
  if (!(variable %in% names(protein))) {
    cat("Warning: Variable", variable, "does not exist in the dataset, skipping.\n")
    next 
  }
  protein_list <- protein[[variable]]
  if (length(protein_list[!is.na(protein_list) & protein_list != ""]) == 0) {
    cat("Warning: Variable", variable, "has no valid protein data, skipping.\n")
    next
  }
  
  results <- enrich_analysis(protein_list, variable)
  input_count <- length(protein_list[!is.na(protein_list) & protein_list != ""])
  mapped_count <- if (!is.null(results$GO) || !is.null(results$KEGG)) {
    nrow(bitr(toupper(protein_list[!is.na(protein_list) & protein_list != ""]), 
              fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  } else {
    0  
  }
  go_count <- if (!is.null(results$GO)) nrow(results$GO) else 0
  kegg_count <- if (!is.null(results$KEGG)) nrow(results$KEGG) else 0
  summary_data <- rbind(summary_data, data.frame(
    Variable = variable,
    Input_Proteins = input_count,
    Mapped_Proteins = mapped_count,
    GO_Pathways = go_count,
    KEGG_Pathways = kegg_count,
    stringsAsFactors = FALSE
  ))
  
  if (!is.null(results$GO) && nrow(results$GO) > 0) {
    sheet_name <- paste0(variable, "_GO")  
    addWorksheet(wb, sheetName = sheet_name) 
    writeData(wb, sheet = sheet_name, x = results$GO) 
    setColWidths(wb, sheet = sheet_name, cols = 1:ncol(results$GO), widths = "auto")  
    setColWidths(wb, sheet = sheet_name, cols = which(colnames(results$GO) == "Proteins"), widths = 50)  
  }
  
  if (!is.null(results$KEGG) && nrow(results$KEGG) > 0) {
    sheet_name <- paste0(variable, "_KEGG")  
    addWorksheet(wb, sheetName = sheet_name) 
    writeData(wb, sheet = sheet_name, x = results$KEGG) 
    setColWidths(wb, sheet = sheet_name, cols = 1:ncol(results$KEGG), widths = "auto") 
    setColWidths(wb, sheet = sheet_name, cols = which(colnames(results$KEGG) == "Proteins"), widths = 50)  
  }
}

writeData(wb, sheet = "Summary", x = summary_data)
setColWidths(wb, sheet = "Summary", cols = 1:ncol(summary_data), widths = "auto")
saveWorkbook(wb, "Protein_Enrichment_Results1.xlsx", overwrite = TRUE)