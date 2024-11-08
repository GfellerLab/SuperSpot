library(Seurat)
library(SuperSpot)
library(SuperCell)
library(tidyr)
library(tidyverse)
library(igraph)
library(peakRAM)
library(Rfast2)






processing <- function(g){
  print(paste0("gamma is: ",g))
  obj <- readRDS(paste0("./SuperSpot/VisiumMouseBrain/MC_seurat_g",g,".rds"))
  obj <- SCTransform(obj,return.only.var.genes = T)



  print("Computing SVF")
  obj <- FindSpatiallyVariableFeatures(obj,assay = "SCT" ,features = VariableFeatures(obj)[1:1000],selection.method = "moransi")
  #svf.df_sorted <- obj@assays[["SCT"]]@meta.data[order(obj@assays[["SCT"]]@meta.data$moransi.spatially.variable.rank, decreasing = F), ]
  #svf_log_norm <- svf.df_sorted$var.features
  print("Done")
  #print(paste0("top 10 svf are: ",svf_log_norm[1:10]))
  #jpeg(file=paste0("./SuperSpot/BCBA/Output/svf_mc_",g,".jpeg"),width = 3440,height = 1440)
  #jpeg(file=paste0("/Volumes/Analyses/Thèse/2023/data/umap_mc_spl_",g,".jpeg"))
  #plot(ImageFeaturePlot(obj, fov = "pancreas", features =  svf_log_norm[1:5], max.cutoff = "q95"))
  #dev.off()
  #save(svf_log_norm,file = paste0("./SuperSpot/BCBA/Output/svf_mc_log_norm_",g,".rda"))
}

#time_df <- tibble::tibble(starting_times = starting_times, ending_times = ending_times)
#write.csv(time_df,"./SuperSpot/01_Data/time_results.csv")



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


saveRDS(mem.df,"./SuperSpot/VisiumMouseBrain/results_memory_log_norm.rds")
write.csv(mem.df,"./SuperSpot/VisiumMouseBrain/results_memory_log_norm.csv")
