library(Seurat)
library(SuperSpot)
library(SuperCell)
library(tidyr)
library(tidyverse)
library(igraph)

## Import data
BCBA.so = Seurat::Load10X_Spatial("/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/BCBA/"  )
SpatialDimPlot(BCBA.so)
bcba.md <- read.csv("/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/BCBA/label.csv",header = T,sep = ";") %>% column_to_rownames("X")


BCBA.so@meta.data <- cbind(BCBA.so@meta.data,bcba.md)
spotPosition <- GetTissueCoordinates(BCBA.so)[,c("x","y")]
colnames(spotPosition) <- c("imagerow","imagecol")
effective.gammas <- c()

for (g in 1:25){
  n.pc = 1:30 # number of first PC to use
  k.knn = 6 # number of neighbors to connect to each spot
  MC <- SCimplify_SpatialDLS(X =  GetAssayData(BCBA.so,assay = "Spatial",layer = "count"),
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
  bcba.md[,str_c("MC_membership_",g)] <- MC$membership %>% as.character()
  MC$cell_type <- supercell_assign(clusters = bcba.md$cell_type,
                                   supercell_membership = MC$membership,
                                   method = "absolute")
  MC.ge <- superspot_GE(MC = MC,
                        ge =  GetAssayData(BCBA.so,assay = "Spatial",layer = "count"),
                        #ge = GetAssayData(MC$seurat.object,layer = "data"),
                        groups = as.numeric(MC$membership),
                        mode = "sum"
  )
  
  MC.seurat <- supercell_2_Seuratv5(
    SC.GE = MC.ge,
    SC = MC,
    fields = c("cell_type")
  )
  
  tissue_position  <- read.csv("/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/BCBA/spatial/tissue_positions_list.csv",header = F)
  MC_centroids <- supercell_spatial_centroids(MC,spotPositions = spotPosition)
  #metadata$cell <- paste(metadata$cell_ID,metadata$fov, sep = "_")
  
  #metadata_filt <- metadata[metadata$cell_id  %in% (names(MC$membership)) == T,]
  tissue_position_sb <- subset(tissue_position, V1 %in% names(MC$membership))
  tissue_position_sb$membership <- MC$membership
  
  rownames(tissue_position_sb) <- NULL
  #metadata_filt <- metadata_filt %>% column_to_rownames("cell")
  
  tissue_position_mean <- stats::aggregate(tissue_position_sb[,2:7] %>% as.matrix() ~membership,tissue_position_sb[,2:7],mean)
  
  #metadata_mean <- stats::aggregate(subset(metadata_filt,select = -c(cell_id,version,Run_name,Run_Tissue_name,tissue,Panel,assay_type)) %>% as.matrix() ~membership,metadata_filt,mean)
  
  #metadata_mean$fov <- as.integer(metadata_mean$fov)
  
  #metadata_mean$cell_ID <- metadata_mean$membership
  
  image <- png::readPNG(source = "/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/BCBA/spatial/tissue_lowres_image.png")
  scale.factors <- Read10X_ScaleFactors(filename = "/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/BCBA/spatial/scalefactors_json.json")
  coordinates <- tissue_position_mean[,1:6]
  coordinates <- column_to_rownames(coordinates,var = "membership")
  colnames(coordinates) <- c("tissue","row","col","imagerow","imagecol")
  coordinates$imagecol <- MC_centroids$imagecol
  coordinates$imagerow <- MC_centroids$imagerow
  #coordinates <- Read10X_Coordinates(filename = Sys.glob(file.path(image.dir, 
  #                                                                 "*tissue_positions*")), filter.matrix)
  fov <- CreateFOV(coordinates[, c("imagerow", "imagecol")], 
                   type = "centroids", radius = scale.factors[["spot"]], 
                   assay = "RNA", key = Key("slice1", quiet = TRUE))
  visium.fov <- new(Class = "VisiumV2", boundaries = fov@boundaries, 
                    molecules = fov@molecules, assay = fov@assay, key = fov@key, 
                    image = image, scale.factors = scale.factors)
  
  MC.seurat[["slice1"]] <- visium.fov
  #write_csv(x=MC_centroids[,c("imagecol","imagerow")],file = paste0("./SuperSpot/01_Data/PancreasCosMx/Centroids_metaspot_g",g,".csv"))
  #write_csv(x=metadata_mean,file = paste0("./SuperSpot/01_Data/PancreasCosMx/Pancreas_metaspot_g",g,".csv"))
  
  #counts.matrix_ms <- t(GetAssayData(MC.seurat ,assay = "RNA") %>% as.matrix()) %>% as.data.frame() %>% rownames_to_column("cell_ID")
  
  #counts.matrix_ms$fov <- metadata_mean$fov
  
  #counts.matrix_ms <- counts.matrix_ms[columns_to_keep]
  
  #write_csv(x=counts.matrix_ms,file = paste0("./SuperSpot/01_Data/PancreasCosMx/count_ms_g",g,".csv"))
  saveRDS(MC.seurat,file = paste0("/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/BCBA/Output/MC_seurat_g",g,".rds"))
}

#write.csv(metad,"./SuperSpot/01_Data/metad.csv")
