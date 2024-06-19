library(Seurat)
library(SuperSpot)
library(SuperCell)
library(tidyr)
library(tidyverse)
library(igraph)
library(peakRAM)
#devtools::install_github("RfastOfficial/Rfast2")
#install.packages("Rfast2")

#print(!requireNamespace('Rfast2'))
library(Rfast2)
print(library(Rfast2))

#starting_times <- c()
#ending_times <- c()

processing <- function(g){
  print(paste0("gamma is: ",g))
  assay <- "Nanostring"
  data <- Seurat::ReadNanostring(data.dir = "./SuperSpot/01_Data/PancreasCosMx",
                                 mtx.file = paste0("count_ms_g",g,".csv"),
                                 metadata.file = paste0("Pancreas_metaspot_g",g,".csv"),
                                 type = "centroids")
  #segs <- CreateSegmentation(data$segmentations)
  segs <- NULL
  cents <- CreateCentroids(data$centroids)
  segmentations.data <- list(centroids = cents, segmentation = segs)
  #coords <- CreateFOV(coords = segmentations.data, type = c("segmentation","centroids"), molecules = data$pixels, assay = assay)
  coords <- CreateFOV(coords = cents, type = c("segmentation","centroids"), molecules = data$pixels, assay = assay)
  obj <- CreateSeuratObject(counts = data$matrix, assay = assay)
  #cells <- intersect(Cells(x = coords, boundary = "segmentation"),
  #Cells(x = coords, boundary = "centroids"))
  cells <- Cells(x = coords, boundary = "centroids")
  cells <- intersect(Cells(obj), cells)
  coords <- subset(x = coords, cells = cells)
  obj[["pancreas"]] <- coords


  #obj <- FindVariableFeatures(obj)
  #obj <- NormalizeData(obj)
  #obj <- ScaleData(obj)
  obj <- SCTransform(obj,assay = assay)
  MC_centroids <- read.csv(file = paste0("./SuperSpot/01_Data/PancreasCosMx/Centroids_metaspot_g",g,".csv"))
  MC_centroids <- as.matrix(MC_centroids)
  colnames(MC_centroids) <- c("x","y")
  obj@images[["pancreas"]]@boundaries[["centroids"]]@coords <- t(MC_centroids)
  print("Computing SVF")
  obj <- FindSpatiallyVariableFeatures(obj, assay = "SCT", features = VariableFeatures(obj)[1:1000],selection.method = "moransi")
  svf.df_sorted <- obj@assays[["SCT"]]@meta.features[order(obj@assays[["SCT"]]@meta.features$moransi.spatially.variable.rank, decreasing = F), ]
  svf <- rownames(svf.df_sorted)
  print("Done")
  print(paste0("top 10 svf are: ",svf[1:10]))
  jpeg(file=paste0("./01_Data/svf_mc_",g,".jpeg"),width = 3440,height = 1440)
  plot(ImageFeaturePlot(obj, fov = "pancreas", features =  svf[1:6], max.cutoff = "q95"))
  dev.off()
  save(svf,file = paste0("/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/01_Data/svf_mc_",g,".rda"))
}



mem.df <- peakRAM::peakRAM(
  processing(1),
  processing(2),
  processing(3),
  processing(4),
  processing(5),
  processing(6),
  processing(7),
  processing(8),
  processing(9),
  processing(10),
  processing(11),
  processing(12),
  processing(13),
  processing(14),
  processing(15),
  processing(16),
  processing(17),
  processing(18),
  processing(19),
  processing(20),
  processing(21),
  processing(22),
  processing(23),
  processing(24),
  processing(25))


saveRDS(mem.df,"./SuperSpot/01_Data/results_memory.rds")
write.csv(mem.df,"./SuperSpot/01_Data/results_memory.csv")
