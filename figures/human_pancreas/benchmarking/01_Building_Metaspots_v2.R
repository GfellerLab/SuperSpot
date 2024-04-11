library(Seurat)
library(SuperSpot)
library(SuperCell)
library(tidyr)
library(tidyverse)
library(igraph)

## Import data
panc <- readRDS("/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/01_Data/PancreasCosMx/Seurat_Pancreas_withTranscripts.rds")
pattern_to_exclude <- "Negative[0-9]+|SystemControl[0-9]+"
mtx <- read_csv("/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/01_Data/PancreasCosMx/Pancreas_exprMat_file.csv")
columns_to_keep <- grep(pattern_to_exclude, colnames(mtx), invert = TRUE, value = TRUE)
mtx <- mtx[, columns_to_keep]
mtx.flt <- mtx[,3:ncol(mtx)] %>% t() %>% as.sparse()
metad <- read_csv("/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/01_Data/PancreasCosMx/Pancreas_metadata_file.csv")
metad <- column_to_rownames(metad,"cell_id")
meta <- readRDS("/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/01_Data/PancreasCosMx/CellType_Accessory_Data/Pancreas_celltype_InSituType.rds")
metad$cell_type <- meta[rownames(metad), "cell_types"]
metad$cell_type[is.na(metad$cell_type)] <- "QC_dropped"
colnames(mtx.flt) <- rownames(metad)
spotPosition <- GetTissueCoordinates(panc)[,c("x","y")]
colnames(spotPosition) <- c("imagerow","imagecol")
effective.gammas <- c()

for (g in 1:25){
  n.pc = 30 # number of first PC to use
  k.knn = 16 # number of neighbors to connect to each spot
  MC <- SCimplify_SpatialDLS(X = mtx.flt,
                                 spotPositions = spotPosition ,
                                 method_similarity = "1",
                                 split_not_connected = T,
                                 genes.use = NULL,
                                 gamma = g,
                                 n.pc = n.pc,
                                 method_knn = "1",
                                 k.knn = k.knn,
                                 method_normalization = "log_normalize",
                                 #cell.annotation = metad$cell_type,
                                 return.seurat.object = FALSE)
  metad[,str_c("MC_membership_",g)] <- MC$membership %>% as.character()
  MC$cell_type <- supercell_assign(clusters = metad$cell_type,
                                            supercell_membership = MC$membership,
                                            method = "absolute")
  MC.ge <- superspot_GE(MC = MC,
                        ge = mtx.flt %>% as.matrix(),
                        #ge = GetAssayData(MC$seurat.object,layer = "data"),
                        groups = as.numeric(MC$membership),
                        mode = "sum"
  )

  MC.seurat <- supercell_2_Seuratv5(
    SC.GE = MC.ge,
    SC = MC,
    fields = c("cell_type")
  )
  saveRDS(MC.seurat,file = paste0("/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/01_Data/MC.seurat_g",g,".rds"))
  #saveRDS(MC.seurat,file = paste0("/Volumes/Analyses/Thèse/2023/data/MC.seurat_g",g,".rds"))

  # print("Creating metaspots")
  # MC.spl <- SCimplify_SpatialDLS(X = mtx.flt,
  #                                spotPositions = spotPosition ,
  #                                method_similarity = "1",
  #                                split_not_connected = T,
  #                                genes.use = NULL,
  #                                gamma = g,
  #                                n.pc = n.pc,
  #                                method_knn = "1",
  #                                k.knn = k.knn,
  #                                method_normalization = "log_normalize",
  #                                cell.annotation = metad$cell_type,
  #                                return.seurat.object = FALSE)
  #
  # print("Done")
  # print("Splitting")
  # MC.full.spl <- split_unconnected(MC.spl)
  # metad[,str_c("MC_membership_spl_",g)] <- MC.full.spl$membership %>%
  #   as.character()
  #
  # MC.full.spl$cell_type <- supercell_assign(clusters = metad$cell_type,
  #                                           supercell_membership = MC.full.spl$membership,
  #                                           method = "absolute")
  #
  # effective.gammas <- c(effective.gammas,MC.full.spl$effect.gamma)
  #
  # MC.fs.ge <- superspot_GE(MC = MC.full.spl,
  #                       ge = mtx.flt %>% as.matrix(),
  #                       #ge = GetAssayData(MC.full.spl$seurat.object,layer = "data"),
  #                       groups = as.numeric(MC.full.spl$membership),
  #                       mode = "sum"
  # )
  #
  # MC.fs.seurat <- supercell_2_Seuratv5(
  #   SC.GE = MC.fs.ge,
  #   SC = MC.full.spl,
  #   fields = c("cell_type")
  # )
  # saveRDS(MC.fs.seurat,file = paste0("/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/01_Data/MC.fs.seurat_g",g,".rds"))
  #saveRDS(MC.fs.seurat,file = paste0("/Volumes/Analyses/Thèse/2023/data/MC.fs.seurat_g",g,".rds"))
}

write.csv(tibble(effective_gammas = effective.gammas),"/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/01_Data/effective_gammas.csv")
write.csv(metad,"/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/01_Data/metad.csv")
