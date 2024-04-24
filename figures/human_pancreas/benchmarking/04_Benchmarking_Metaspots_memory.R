library(Seurat)
library(SuperSpot)
library(SuperCell)
library(tidyr)
library(tidyverse)
library(igraph)
library(peakRAM)

#starting_times <- c()
#ending_times <- c()

processing <- function(g){
  print(paste0("gamma is: ",g))
  #start_time <- Sys.time()
  #starting_times <- c(starting_times,start_time)
  so <- readRDS(paste0("/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/01_Data/MC.seurat_g",g,".rds"))
  #so <- readRDS(paste0("/Volumes/Analyses/Thèse/2023/data/MC.fs.seurat_g",g,".rds"))
  #so <- NormalizeData(so)
  # MC.seurat <- SetAssayData(
  #     object = MC.seurat,
  #     layer = "data",
  #     new.data = GetAssayData(MC.seurat,layer = "counts"),
  #     assay = "RNA"
  # )

  #so <- ScaleData(so)

  #so <- FindVariableFeatures(so)
  so <- SCTransform(so)
  so <- RunPCA(so)

  so <- RunUMAP(so, dims = 1:30)


  Idents(so) <- "cell_type"
  levels(so) <- sort(levels(so))
  jpeg(file=paste0("/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/01_Data/umap_mc_",g,".jpeg"))
  #jpeg(file=paste0("/Volumes/Analyses/Thèse/2023/data/umap_mc_spl_",g,".jpeg"))
  plot(DimPlot(so))
  dev.off()

  #VlnPlot(so,"nFeature_RNA")+NoLegend()

  #Compute marker genes
  mc.markers <-  FindAllMarkers(so,
                                only.pos = TRUE,
                                #min.pct = 0.25,
                                logfc.threshold = 0.25 ) %>%
    filter(p_val_adj < 0.05)

  mc.top.markers <- mc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log2FC)
  #end_time <- Sys.time()
  #ending_times <- c(ending_times,end_time)
}

#time_df <- tibble::tibble(starting_times = starting_times, ending_times = ending_times)
#write.csv(time_df,"/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/01_Data/time_results.csv")



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


saveRDS(mem.df,"/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/01_Data/results_memory.rds")
write.csv(mem.df,"/work/FAC/FBM/LLB/dgfeller/scrnaseq/mteleman/SuperSpot/01_Data/results_memory.csv")
