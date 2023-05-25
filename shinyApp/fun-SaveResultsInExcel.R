SaveResultsInExcel <- function(Results, gMCS.info, filename){
  
  wb <- createWorkbook()
  ## Add worksheets
  addWorksheet(wb, "num single KO")
  addWorksheet(wb, "ratio single KO")
  addWorksheet(wb, "gMCS single KO")
  if (!is.null(Results$num.essential.pair.gene)){
    addWorksheet(wb, "num Double KO")
    addWorksheet(wb, "ratio Double KO")
    addWorksheet(wb, "gMCS Double KO")
    addWorksheet(wb, "num combined single double")
    addWorksheet(wb, "ratio combined single double")
  }
  
  ## Add the data
  writeDataTable(wb, "num single KO", Results$num.essential.gene, tableStyle = "TableStyleLight1")
  writeDataTable(wb, "ratio single KO", Results$ratio.essential.gene, tableStyle = "TableStyleLight1")
  aux <- Results$list.gMCS.essential
  aux$gMCS_ENSEMBL <- gMCS.info$gMCSs.ENSEMBL.txt[as.numeric(as.character(aux$gMCS))]
  aux$gMCS_SYMBOL <- gMCS.info$gMCSs.SYMBOL.txt[as.numeric(as.character(aux$gMCS))]
  writeDataTable(wb, "gMCS single KO", aux, tableStyle = "TableStyleLight1", colNames = T, rowNames = F)
  
  if (!is.null(Results$num.essential.pair.gene)){
    writeDataTable(wb, "num Double KO", Results$num.essential.pair.gene, tableStyle = "TableStyleLight1")
    writeDataTable(wb, "ratio Double KO", Results$ratio.essential.pair.gene, tableStyle = "TableStyleLight1")
    aux <- Results$list.gMCS.essential.pair
    aux$gMCS_ENSEMBL <- gMCS.info$gMCSs.ENSEMBL.txt[as.numeric(as.character(aux$gMCS))]
    aux$gMCS_SYMBOL <- gMCS.info$gMCSs.SYMBOL.txt[as.numeric(as.character(aux$gMCS))]
    writeDataTable(wb, "gMCS Double KO", aux, tableStyle = "TableStyleLight1", colNames = T, rowNames = F)
    writeDataTable(wb, "num combined single double", Results$num.essential.single.pair.combined.gene, tableStyle = "TableStyleLight1")
    writeDataTable(wb, "ratio combined single double", Results$ratio.essential.single.pair.combined.gene, tableStyle = "TableStyleLight1")
  }
  
  
  saveWorkbook(wb, filename, overwrite = TRUE)
  
}
