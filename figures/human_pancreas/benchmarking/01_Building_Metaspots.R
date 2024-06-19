library(Seurat)
library(SuperSpot)
library(SuperCell)
library(tidyr)
library(tidyverse)
library(igraph)

## Import data
panc <- readRDS("./SuperSpot/01_Data/PancreasCosMx/Seurat_Pancreas_withTranscripts.rds")
pattern_to_exclude <- "Negative[0-9]+|SystemControl[0-9]+"
mtx <- read_csv("./SuperSpot/01_Data/PancreasCosMx/Pancreas_exprMat_file.csv")
columns_to_keep <- grep(pattern_to_exclude, colnames(mtx), invert = TRUE, value = TRUE)
mtx <- mtx[, columns_to_keep]
mtx.flt <- mtx[,3:ncol(mtx)] %>% t() %>% as.sparse()
metad <- read_csv("./SuperSpot/01_Data/PancreasCosMx/Pancreas_metadata_file.csv")
metad <- column_to_rownames(metad,"cell_id")
meta <- readRDS("./SuperSpot/01_Data/PancreasCosMx/CellType_Accessory_Data/Pancreas_celltype_InSituType.rds")
metad$cell_type <- meta[rownames(metad), "cell_types"]
metad$cell_type[is.na(metad$cell_type)] <- "QC_dropped"
colnames(mtx.flt) <- rownames(metad)
spotPosition <- GetTissueCoordinates(panc)[,c("x","y")]
colnames(spotPosition) <- c("imagerow","imagecol")
effective.gammas <- c()

for (g in 1:25){
  n.pc = 1:30 # number of first PC to use
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
                                 method_normalization = "SCT",
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
  saveRDS(MC.seurat,file = paste0("./SuperSpot/01_Data/MC.seurat_g",g,".rds"))
  metadata  <- read.csv("./SuperSpot/01_Data/PancreasCosMx/Pancreas_metadata_file.csv")
  MC_centroids <- supercell_spatial_centroids(MC,spotPositions = spotPosition)
  #metadata$cell <- paste(metadata$cell_ID,metadata$fov, sep = "_")

  metadata_filt <- metadata[metadata$cell_id  %in% (names(MC$membership)) == T,]

  metadata_filt$membership <- MC$membership

  rownames(metadata_filt) <- NULL
  metadata_filt <- metadata_filt %>% column_to_rownames("cell")

  #metadata_mean <- stats::aggregate(metadata_filt %>% as.matrix() ~membership,metadata_filt,mean)

  metadata_mean <- stats::aggregate(subset(metadata_filt,select = -c(cell_id,version,Run_name,Run_Tissue_name,tissue,Panel,assay_type)) %>% as.matrix() ~membership,metadata_filt,mean)

  metadata_mean$fov <- as.integer(metadata_mean$fov)

  metadata_mean$cell_ID <- metadata_mean$membership


  write_csv(x=MC_centroids[,c("imagecol","imagerow")],file = paste0("./SuperSpot/01_Data/PancreasCosMx/Centroids_metaspot_g",g,".csv"))
  write_csv(x=metadata_mean,file = paste0("./SuperSpot/01_Data/PancreasCosMx/Pancreas_metaspot_g",g,".csv"))

  counts.matrix_ms <- t(GetAssayData(MC.seurat ,assay = "RNA") %>% as.matrix()) %>% as.data.frame() %>% rownames_to_column("cell_ID")

  counts.matrix_ms$fov <- metadata_mean$fov

  counts.matrix_ms <- counts.matrix_ms[columns_to_keep]

  write_csv(x=counts.matrix_ms,file = paste0("./SuperSpot/01_Data/PancreasCosMx/count_ms_g",g,".csv"))
}

write.csv(metad,"./SuperSpot/01_Data/metad.csv")
