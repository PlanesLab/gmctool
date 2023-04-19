#############################################
######                                 ######
######        Analysis of gMCSs        ######
######        and gene expression      ######
######                                 ######
#############################################
# 
# Author: Luis V Valc?rcel
# Date: 2019-05-02
# 


rm (list = ls())

# library(foreach)
# library(doParallel)
# library(parallel)
library(data.table)
# library(biomaRt)
library(Matrix)

# setwd("D:/PhD/1 - Multiple Myeloma Human-GEM NY - CoMMpass - CCLE EBI 2021-06/0 - GENERATE APP/Code_for_tool_v1")
setwd("C:/Users/lvalcarcel/OneDrive - Tecnun/Postdoc/002 - MM Human-GEM NY - CoMMpass - CCLE EBI 2021-06/00 - GENERATE APP/gMCStool-shinyApp-v5/00_Generate_gmcs_database_info")

dir.name <- "D:/PhD/1 - Multiple Myeloma Human-GEM NY - CoMMpass - CCLE EBI 2021-06/0 - GENERATE APP/Code_for_tool_v1/Datasets_gene_expression"

# Load the gene information
system.time({load(file.path(dir.name,paste0("all_genes_HumanGEM_CCLE_MM_TPM.Rdata")))})


gMCSs.ENSEMBL <- c(paste0("./Data-gMCS-Human-GEM-1.4.0/Consensus_Human-GEM-task-",1:57,"_gMCSs_ENSEMBL_20201028.txt"),
                   "./Data-gMCS-Human-GEM-1.4.0/Human-GEM-1.4.0-FullBiomass-gMCSs-ENSEMBL-2020-07-13.txt")
names(gMCSs.ENSEMBL) <- 1:58

gMCSs.ENSEMBL <- lapply(gMCSs.ENSEMBL, function(x){
  y <- matrix("", nrow = 0, ncol = 0)
  try({
    y <- as.matrix(as.data.frame(fread(file = x)))
    y <- unique(y)
    
    y <- y[apply(!is.na(y),1,sum)>0, ]
    y <- y[apply(y!="",1,sum)>0, ]
    
    y <- y[, apply(!is.na(y),2,sum)>0]
    y <- y[, apply(y!="",2,sum)>0]
  })
  return(y)
})

sapply(gMCSs.ENSEMBL, dim)

gMCSs.ENSEMBL.all <- gMCSs.ENSEMBL


# generate info for Cell Lines with culture media ####


source("fun-GenerateGmcsDatabase.R")


gMCS.db.info <- GenerateGmcsDatabase(gmcs.list.tab = gMCSs.ENSEMBL.all[1:58], table.genes.HumanGEM = table.genes.HumanGEM)
names(gMCS.db.info$gMCSs.ENSEMBL) <- 1:58


gMCS.info.raw <- list()
gMCS.info.raw[["EssentialTasks_CultureMedium"]] <- GenerateGmcsDatabase(gmcs.list.tab = gMCSs.ENSEMBL.all[1:57], table.genes.HumanGEM = table.genes.HumanGEM)
gMCS.info.raw[["EssentialTasks_FullMedium"]] <- GenerateGmcsDatabase(gmcs.list.tab = gMCSs.ENSEMBL.all[c(1:56,58)], table.genes.HumanGEM = table.genes.HumanGEM)
gMCS.info.raw[["EssentialTasks_FullMedium"]]$table.gMCSs$task[gMCS.info.raw[["EssentialTasks_FullMedium"]]$table.gMCSs$task==57] <- 58
gMCS.info.raw[["Only_CultureMedium"]] <- GenerateGmcsDatabase(gmcs.list.tab = gMCSs.ENSEMBL.all[c(57)], table.genes.HumanGEM = table.genes.HumanGEM)
gMCS.info.raw[["Only_CultureMedium"]]$table.gMCSs$task <- 57
gMCS.info.raw[["Only_FullMedium"]] <- GenerateGmcsDatabase(gmcs.list.tab = gMCSs.ENSEMBL.all[c(58)], table.genes.HumanGEM = table.genes.HumanGEM)
gMCS.info.raw[["Only_FullMedium"]]$table.gMCSs$task <- 58

gMCS.info.raw <- lapply(gMCS.info.raw, as.list)
# remove unnecessary field (gMCS.ENSEMBL.list)
# gMCS.info.raw <- lapply(gMCS.info.raw, function(x){x[c("gMCSs.ENSEMBL.txt", "table.gMCSs", #"gMCSs.ENSEMBL.length",
#                                                        "gMCSs.ENSEMBL.mat", "genes.gMCSs.ENSEMBL", "table.genes.HumanGEM",
#                                                        "gMCSs.ENSEMBL.txt.SYMBOL", "gMCSs.ENSEMBL")]})

gMCS.info.raw <- lapply(gMCS.info.raw, function(x){x[ names(x) != "gMCSs.ENSEMBL"]})


# generate field with reduced information
gMCS.info.raw[["Custom_CultureMedium"]] <- gMCS.info.raw[["EssentialTasks_CultureMedium"]][c("table.gMCSs", "genes.gMCSs.ENSEMBL", "table.genes.HumanGEM")]
gMCS.info.raw[["Custom_FullMedium"]] <- gMCS.info.raw[["EssentialTasks_FullMedium"]][c("table.gMCSs", "genes.gMCSs.ENSEMBL", "table.genes.HumanGEM")]

gMCS.info.raw[["EssentialTasks_CultureMedium"]]$fullname <- "Essential Tasks and Growth on Ham's medium"
gMCS.info.raw[["EssentialTasks_FullMedium"]]$fullname <- "Essential Tasks and Growth on unconstrained medium"
gMCS.info.raw[["Only_CultureMedium"]]$fullname <- "Only growth on Ham's medium"
gMCS.info.raw[["Only_FullMedium"]]$fullname <- "Only growth on unconstrained medium"
gMCS.info.raw[["Custom_CultureMedium"]]$fullname <- "Selected metabolic tasks and growth on Ham's medium"
gMCS.info.raw[["Custom_FullMedium"]]$fullname <- "Selected metabolic tasks and growth on unconstrained medium"

format(object.size(gMCS.info.raw), units = "auto")
sapply(gMCS.info.raw, function(x) format(object.size(x), units = "auto"))
sapply(gMCS.info.raw[[1]], function(x) format(object.size(x), units = "auto"))

save(gMCS.info.raw, file = "./gMCSs_all_cases_HumanGEMv1.4.0_ENSEMBL.Rdata", compress = "xz", compression_level = 9)
save(gMCS.info.raw, file = "../Data/gMCSs_all_cases_HumanGEMv1.4.0_ENSEMBL.Rdata", compress = "xz", compression_level = 9)
saveRDS(gMCS.db.info$gMCSs.ENSEMBL, file = "../Data/gMCSs_all_cases_HumanGEMv1.4.0_ENSEMBL_raw.RDS", compress = "xz")
