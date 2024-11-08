library(Seurat)
library(SuperSpot)
library(SuperCell)
library(tidyr)
library(tidyverse)
library(igraph)
library(arrow)
library(BPCells)

spotPositions <- arrow::read_parquet("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/binned_outputs/square_002um/spatial/tissue_positions.parquet") %>% column_to_rownames("barcode")
spotPositions.filt <- subset(spotPositions, in_tissue == 1)
rm(spotPositions)
mtx.count <- Seurat::Read10X_h5("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/binned_outputs/square_002um/filtered_feature_bc_matrix.h5")
#mtx.count <- open_matrix_10x_hdf5("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/binned_outputs/square_002um/filtered_feature_bc_matrix.h5", feature_type="Gene Expression") %>%
  #write_matrix_dir("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/binned_outputs/square_002um/filtered_feature_bc_matrix_raw")
spotPositions.filt <- spotPositions.filt[colnames(mtx.count),c("pxl_col_in_fullres","pxl_row_in_fullres")]
colnames(spotPositions.filt) <- c("imagecol","imagerow")

metadat_metaspot <- tibble(spots = rownames(spotPositions.filt))


neighbor_graph <- function(spotPositions, k.spots, countMatrix, n.pc, method_similarity, method_knn, k.knn, sigs,
                           method_normalization,split_not_connected,split_annotation,split_vector,n.cpu, plot.graph, pct){
  if (method_knn == "1"){
    print(Sys.time())
    print(paste0("Building KNN graph with nn2"))
    nn2.res <- RANN::nn2(data = spotPositions,k = k.knn+1)
    if (split_not_connected == TRUE){
      #min_dist <- quantile(nn2.res$nn.dists[,2:k.knn+1] %>% as.vector(),names = F)[3]
      min_dist <- quantile(nn2.res$nn.dists[,2:k.knn+1] %>% as.vector(),names = F,probs = pct)
      print(paste0("Neighbors with distance > ",min_dist, " are removed"))
      plot(ggplot2::ggplot(data = tibble(distances = as.vector(nn2.res$nn.dists[,2:k.knn+1]),
                                         distribution = rep(".",length(as.vector(nn2.res$nn.dists[,2:k.knn+1])))))+
             ggplot2::geom_violin(aes(x = distribution, y = distances))+
             ggplot2::geom_boxplot(aes(x = distribution, y = distances),width = 0.1)+
             ggplot2::geom_hline(yintercept=min_dist, linetype="dashed", color = "red"))

      #bad_neighbors_rows <- which(nn2.res$nn.dists > round(min_dist+0.6),arr.ind = T)[,1]
      #bad_neighbors_cols <- which(nn2.res$nn.dists > round(min_dist+0.6),arr.ind = T)[,2]
      bad_neighbors_rows <- which(nn2.res$nn.dists > min_dist,arr.ind = T)[,1]
      bad_neighbors_cols <- which(nn2.res$nn.dists > min_dist,arr.ind = T)[,2]
      new_neighbors <- nn2.res$nn.idx
      for (ind in 1:length(bad_neighbors_rows)){
        new_neighbors[bad_neighbors_rows[ind],bad_neighbors_cols[ind]] <- new_neighbors[bad_neighbors_rows[ind],1]
      }
      spot.graph <- bluster::neighborsToKNNGraph(new_neighbors[,2:k.knn+1],directed = F)

    }
    else if (split_not_connected == FALSE){
      spot.graph <- bluster::neighborsToKNNGraph(nn2.res$nn.idx[,2:k.knn+1])
    }

    if (!is.null(split_annotation)){
      cell_type.id <- which(split_vector == split_annotation)
      spot.graph <- igraph::subgraph(spot.graph,cell_type.id)
      spotPositions <- spotPositions[cell_type.id,]
      countMatrix <- countMatrix[,cell_type.id]
    }

    max_gamma <- (nrow(spotPositions)/igraph::components(spot.graph)$no)
    print(paste0("Maximum gamma is ",max_gamma))
    print(paste0("Done"))
    print(Sys.time())
  }
  else if (method_knn == "2"){
    print(paste0("Building KNN graph by DLS"))
    print(paste0("Computing distances between spots"))
    distCoord <- parallelDist::parDist(as.matrix(spotPositions),threads = n.cpu)
    min.dist <- min(distCoord)
    distCoord[distCoord > k.spots*round(min.dist)] <- 0
    distCoord[distCoord >0] <- 1
    spot.graph <- igraph::graph.adjacency(as.matrix(distCoord),mode = "undirected")
    print(paste0("Done"))
  }

  if (plot.graph == TRUE){
    igraph::V(spot.graph)$x <- spotPositions$imagecol
    igraph::V(spot.graph)$y <- -spotPositions$imagerow
    plot(spot.graph,
         layout = matrix(c(igraph::V(spot.graph)$x,igraph::V(spot.graph)$y), ncol = 2),
         vertex.label.cex = 0.1,
         vertex.size = 1,
         edge.arrow.size = 0.7,
         main = "Graph with Custom X and Y Coordinates")
  }
  library(BPCells)
  # Normalize by reads-per-cell
  #mat <- multiply_cols(countMatrix, 1/Matrix::colSums(countMatrix))

  # Log normalization
  #mat <- log1p(mat * 10000) # Log normalization
  #stats <- matrix_stats(mat, row_stats="variance")

  # To keep the example small, we'll do a very naive variable gene selection
  #variable_genes <- order(stats$row_stats["variance",], decreasing=TRUE) %>%
  #  head(1000) %>%
  #  sort()

  #mat_norm <- mat[variable_genes,]
  #mat_norm <- mat_norm %>% write_matrix_dir(tempfile("mat"))#write_matrix_memory(compress=FALSE)
  #gene_means <- stats$row_stats["mean",variable_genes]
  #gene_vars <- stats$row_stats["variance", variable_genes]
  #mat_norm <- (mat_norm - gene_means) / gene_vars
  #svd <- BPCells::svds(mat_norm, k=50)
  print("Computing PCA using BPCells")
  print(Sys.time())
  svd <- BPCells::svds(countMatrix, k=25)
  # Alternate option: irlba::irlba(mat_norm, nv=50)
  pca_matrix <- multiply_cols(svd$v, svd$d)
  print("Done")
  print(Sys.time())
  cat(sprintf("PCA dimensions: %s\n", toString(dim(pca_matrix))))
  pca_matrix[1:4,1:3]

  # so <- Seurat::CreateSeuratObject(countMatrix,assay = "RNA")
  # if (method_normalization == "log_normalize"){
  #   print(paste0("Performing Log Normalization"))
  #   so <- Seurat::NormalizeData(so)
  #   print(paste0("Done"))
  #   print(paste0("Running PCA"))
  #   so <- Seurat::FindVariableFeatures(so)
  #   so <- Seurat::ScaleData(so,do.scale = FALSE,do.center = FALSE)
  #   if (max_gamma != 1){
  #     so <- Seurat::RunPCA(so,verbose = F)
  #     #pca_matrix <- so@reductions[["pca"]]@cell.embeddings[,1:n.pc]
  #     pca_matrix <- so@reductions[["pca"]]@cell.embeddings[,n.pc]
  #   }
  #   else {
  #     pca_matrix <- matrix(0,ncol(countMatrix),n.pc)
  #   }
  # }
  # else if (method_normalization == "no"){
  #   print(paste0("Using Raw Data"))
  #   so<-Seurat::SetAssayData(so,new.data = Seurat::GetAssayData(so,assay = "RNA",layer = "counts"),assay = "RNA",layer = "data")
  #   print(paste0("Done"))
  #   print(paste0("Running PCA"))
  #   so <- Seurat::FindVariableFeatures(so)
  #   so <- Seurat::ScaleData(so,do.scale = FALSE,do.center = FALSE)
  #   if (max_gamma != 1){
  #     so <- Seurat::RunPCA(so,verbose = F)
  #     #pca_matrix <- so@reductions[["pca"]]@cell.embeddings[,1:n.pc]
  #     pca_matrix <- so@reductions[["pca"]]@cell.embeddings[,n.pc]
  #   }
  #   else {
  #     pca_matrix <- matrix(0,ncol(countMatrix),n.pc)
  #   }
  # }
  # else if (method_normalization == "SCT"){
  #   print(paste0("Performing SCT"))
  #   so <- Seurat::SCTransform(so, new.assay.name = "SCT" , assay = "RNA",verbose = F)
  #   print(paste0("Done"))
  #   print(paste0("Running PCA"))
  #   so <-Seurat::RunPCA(so, assay = "SCT",verbose = F)
  #   pca_matrix <- so@reductions[["pca"]]@cell.embeddings[,n.pc]
  # }
  #
  print(paste0("Done"))
  print(paste0("Computing PCA dist for each edge"))
  print(Sys.time())
  #pca_dist <- pbapply::pbapply(igraph::get.edgelist(spot.graph),
  #                             1 ,
  #                             function(x){parallelDist::parDist(pca_matrix[x,],threads = n.cpu)})
  pca_similarity <- pbapply::pbapply(igraph::get.edgelist(spot.graph),
                               1 ,
                               function(x){1/parallelDist::parDist(pca_matrix[x,],threads = n.cpu)})
  print(head(pca_similarity))
  pca_similarity <- pca_similarity%>% unlist()
  print(paste0("Dimensions of pca distances: ",dim(pca_similarity)))
  print(paste0("Done"))
  print(paste0("Computing similarity from PCA distances"))
  print(Sys.time())
  # if (method_similarity == "1"){
  #   pca_similarity <- pbapply::pblapply(pca_dist,FUN = function(i){1 - (i/max(pca_dist))},) %>% unlist()
  #   pca_similarity[pca_similarity == 0 ] <- 1e-100
  #   print(paste0("Done"))
  # }
  # else if (method_similarity == "2"){
  #   pca_similarity <- pbapply::pblapply(pca_dist,FUN = function(i){exp(-(i**2)/sigs)}) %>% unlist()
  #   pca_similarity[pca_similarity == 0 ]<- 1e-100
  #   print(paste0("Done"))
  # }
  #
  # else if (method_similarity == "3"){
  #   pca_similarity <- proxy::pr_dist2simil(pca_dist)
  #   pca_similarity[pca_similarity == 0 ]<- 1e-100
  #   print(paste0("Done"))
  # }
  print(paste0("Returning graph with PCA similarity as weight"))
  igraph::E(spot.graph)$weight <- pca_similarity
  return_list <- list("graph" = spot.graph, "seurat.object" = countMatrix)
  return(return_list)
}

SCimplify_SpatialDLS_HD <- function(X,
                                 spotPositions,
                                 k.spots = 1,
                                 method_similarity = "1",
                                 method_knn = "1",
                                 method_normalization = "log_normalize",
                                 split_not_connected = TRUE,
                                 split_annotation = NULL,
                                 split_vector,
                                 sigs = 1000,
                                 cell.annotation = NULL,
                                 cell.split.condition = NULL,
                                 gamma = 10,
                                 k.knn = 6,
                                 do.scale = TRUE,
                                 n.pc = 1:30,
                                 n.cpu = NULL,
                                 plot.graph = FALSE,
                                 pct = 0.6,
                                 seed = 12345,
                                 igraph.clustering = c("walktrap", "louvain"),
                                 return.singlecell.NW = TRUE,
                                 return.hierarchical.structure = TRUE,
                                 return.seurat.object = FALSE,
                                 ...){

  neighbor_graph.output <- neighbor_graph(spotPositions = spotPositions, k.spots = k.spots,split_not_connected = split_not_connected,
                                          countMatrix = X, n.pc = n.pc, method_similarity = method_similarity,
                                          method_knn = method_knn, k.knn = k.knn, sigs = sigs,method_normalization = method_normalization,
                                          split_annotation = split_annotation,split_vector =split_vector, n.cpu = n.cpu,plot.graph = plot.graph,
                                          pct = pct)

  #X <- Seurat::GetAssayData(neighbor_graph.output$seurat.object)
  X <- neighbor_graph.output$seurat.object
  N.c <- ncol(X)

  presampled.cell.ids <- colnames(X)



  sc.nw <- neighbor_graph.output$graph

  #simplify

  k   <- round(N.c/gamma)
  print(paste0("Clustering"))
  print(Sys.time())
  if(igraph.clustering[1] == "walktrap"){
    #g.s              <- igraph::cluster_walktrap(sc.nw,weights = igraph::E(sc.nw)$weight)
    #g.s              <- igraph::cluster_leiden(sc.nw,weights = igraph::E(sc.nw)$weight,resolution_parameter = 50)
    g.s              <- igraph::cluster_fast_greedy(sc.nw,weights = igraph::E(sc.nw)$weight)
    g.s$membership   <- igraph::cut_at(g.s, k)

  } else if(igraph.clustering[1] == "louvain") {
    warning(paste("igraph.clustering =", igraph.clustering, ", gamma is ignored"))
    g.s    <- igraph::cluster_louvain(sc.nw)

  } else {
    stop(paste("Unknown clustering method (", igraph.clustering, "), please use louvain or walkrtap"))
  }
  print(paste0("Done"))
  print(Sys.time())
  membership.presampled        <- g.s$membership
  names(membership.presampled) <- presampled.cell.ids


  ## Split super-cells containing cells from different annotations or conditions
  if(!is.null(cell.annotation) | !is.null(cell.split.condition)){
    if(is.null(cell.annotation)) cell.annotation <- rep("a", N.c)
    if(is.null(cell.split.condition)) cell.split.condition <- rep("s", N.c)
    names(cell.annotation) <- names(cell.split.condition) <- colnames(X)

    split.cells <- interaction(cell.annotation[presampled.cell.ids], cell.split.condition[presampled.cell.ids], drop = TRUE)

    membership.presampled.intr <- interaction(membership.presampled, split.cells, drop = TRUE)
    membership.presampled <- as.numeric(membership.presampled.intr)

    map.membership <- unique(membership.presampled)
    names(map.membership) <- unique(as.vector(membership.presampled.intr))

    names(membership.presampled) <- presampled.cell.ids
  }



  SC.NW                        <- igraph::contract(sc.nw, membership.presampled)
  SC.NW                        <- igraph::simplify(SC.NW, remove.loops = T, edge.attr.comb="sum")



  membership.all <- membership.presampled[colnames(X)]

  membership <- membership.all

  supercell_size <- as.vector(table(membership))

  igraph::E(SC.NW)$width         <- sqrt(igraph::E(SC.NW)$weight/10)
  igraph::V(SC.NW)$size          <- supercell_size
  igraph::V(SC.NW)$sizesqrt      <- sqrt(igraph::V(SC.NW)$size)

  res <- list(graph.supercells = SC.NW,
              gamma = gamma,
              N.SC = length(unique(membership)),
              membership = membership,
              supercell_size = supercell_size,
              simplification.algo = igraph.clustering[1],
              n.pc = n.pc,
              k.spots = k.spots,
              sigs = sigs,
              method_similarity = method_similarity,
              method_knn = method_knn,
              method_normalization = method_normalization,
              sc.cell.annotation. = cell.annotation,
              sc.cell.split.condition. = cell.split.condition
  )

  if(return.singlecell.NW){res$graph.singlecell <- sc.nw}
  if(!is.null(cell.annotation) | !is.null(cell.split.condition)){
    res$SC.cell.annotation. <- supercell_assign(cell.annotation, res$membership)
    res$SC.cell.split.condition. <- supercell_assign(cell.split.condition, res$membership)
  }


  if(igraph.clustering[1] == "walktrap" & return.hierarchical.structure){ res$h_membership <- g.s}
  gc()

  if(return.seurat.object == TRUE){res$seurat.object = neighbor_graph.output$seurat.object}

  return(res)
}

supercell_metaspots_shape_v2 <- function(MC, spotpositions, annotation, concavity, membership_name = "membership") {
  message("Creating polygons for visualization")

  # Load pbapply for progress bar
  if (!requireNamespace("pbapply", quietly = TRUE)) {
    install.packages("pbapply")
  }
  library(pbapply)

  unique_memberships <- unique(MC[[membership_name]])

  # Use pblapply to apply function over unique memberships with progress bar
  hull_list <- pbapply::pblapply(unique_memberships, function(memb) {
    index.tmp <- which(MC[[membership_name]] == memb)
    polygons.tmp <- concaveman::concaveman(as.matrix(spotpositions[index.tmp, ]), concavity = concavity)
    index.tmp2 <- which(names(MC[[annotation]]) == memb)
    hull_df.tmp <- data.frame(
      x = polygons.tmp[, 2],
      y = polygons.tmp[, 1],
      cell_type = rep(MC[[annotation]][index.tmp2], nrow(polygons.tmp)),
      membership = rep(memb, nrow(polygons.tmp))
    )
    return(hull_df.tmp)
  })

  # Combine the list of data frames into one data frame
  hull_df_final <- do.call(rbind, hull_list)

  message("Done")
  return(hull_df_final)
}

options(future.globals.maxSize = 1e12)
g = 16 # gamma
n.pc = 1:30 # number of first PC to use
k.knn = 8 # number of neighbors to connect to each spot

print("Creating metaspots")
MC.spl <- SCimplify_SpatialDLS_HD(X = mtx.count,
                               spotPositions = spotPositions.filt ,
                               method_similarity = "1",
                               split_not_connected = T,
                               genes.use = NULL,
                               gamma = g,
                               n.pc = n.pc,
                               method_knn = "1",
                               k.knn = k.knn,
                               method_normalization = "log_normalize",
                               n.cpu = 16,
                               return.seurat.object = FALSE,)

print("Done")

metadat_metaspot$MC_membership <- MC.spl$membership %>% as.character()
MC.spl$cell_type <- rep("unk", MC.spl$N.SC)
names(MC.spl$cell_type) <- paste0(sort(unique(MC.spl$membership)))
#names(MC.spl$cell_type) <- names(MC.spl$membership)
MC.spl$polygons <- supercell_metaspots_shape_v2(MC = MC.spl,
                                                       spotpositions = spotPositions.filt,
                                                       annotation = "cell_type",
                                                       concavity = 2,membership_name = "membership")
MC.spl$centroids <- supercell_spatial_centroids(MC = MC.spl,spotPositions = spotPositions.filt)
saveRDS(MC.spl,"./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/MC_spl_g16_fg.rds")

MC.fs.ge <- superspot_GE(MC = MC.spl,
                         ge = mtx.count, #%>% as.matrix(),
                         #ge = GetAssayData(MC.full.spl$seurat.object,layer = "data"),
                         groups = as.numeric(MC.spl$membership),
                         mode = "sum"
)

MC.fs.seurat <- supercell_2_Seuratv5(
  SC.GE = MC.fs.ge,
  SC = MC.spl,
  fields = c("cell_type")
)
#MC.fs.seurat <- NormalizeData(MC.fs.seurat)
MC.fs.seurat <- SetAssayData(MC.fs.seurat, layer = "data", GetAssayData(MC.fs.seurat,layer = "counts"))
MC.fs.seurat <- FindVariableFeatures(MC.fs.seurat)
MC.fs.seurat <- ScaleData(MC.fs.seurat,do.scale = F,do.center = F)
#MC.fs.seurat <- ScaleData(MC.fs.seurat)
MC.fs.seurat <- RunPCA(MC.fs.seurat)
MC.fs.seurat <- FindNeighbors(MC.fs.seurat)
MC.fs.seurat <- FindClusters(MC.fs.seurat,resolution = 0.5)
MC.fs.seurat <- RunUMAP(MC.fs.seurat,dims = 1:30)
jpeg(file=paste0("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/umap_mc_16.jpeg"),width = 1920,height = 1080)
plot(DimPlot(MC.fs.seurat))
dev.off()
#saveRDS(MC.fs.seurat,"/Users/admin/Downloads/square_008um/MC_fs_seurat.rds")
saveRDS(MC.fs.seurat,"./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/MC_fs_seurat_g16_fg.rds")
pca_embeddings_mc <- Embeddings(MC.fs.seurat, "pca")

# if (!file.exists("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/MC_fs_ge_raw")) {
# MC.fs.ge %>% write_matrix_dir("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/MC_fs_ge_raw")
# } else {
#   MC.fs.ge <- open_matrix_dir("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/MC_fs_ge_raw")
# }
#
# print( "Normalize by reads-per-cell")
# print(Sys.time())
# mat <- multiply_cols(MC.fs.ge, 1/Matrix::colSums(MC.fs.ge))
#
# print( "Log normalization")
# print(Sys.time())
# mat <- log1p(mat * 10000) # Log normalization
# stats <- matrix_stats(mat, row_stats="variance")
#
# print( "variable gene selection")
# print(Sys.time())
# variable_genes <- order(stats$row_stats["variance",], decreasing=TRUE) %>%
#   head(1000) %>%
#   sort()
#
# mat_norm <- mat[variable_genes,]
# mat_norm
# print( "scaling")
# print(Sys.time())
# mat_norm <- mat_norm %>% write_matrix_dir(tempfile("mat"))
# gene_means <- stats$row_stats["mean",variable_genes]
# gene_vars <- stats$row_stats["variance", variable_genes]
# mat_norm <- (mat_norm - gene_means) / gene_vars
# print( "PCA")
# print(Sys.time())
# svd <- BPCells::svds(mat_norm, k=50)
# # Alternate option: irlba::irlba(mat_norm, nv=50)
# pca_embeddings_mc <- multiply_cols(svd$v, svd$d)
#
# cat(sprintf("PCA dimensions: %s\n", toString(dim(pca_embeddings_mc))))
# pca_embeddings_mc[1:4,1:3]
# set.seed(12341512)
# umap <- uwot::umap(pca_embeddings_mc)
# umap[1:4,]
# clusts <- knn_hnsw(pca_embeddings_mc, ef=2000) %>% # Find approximate nearest neighbors
#   knn_to_snn_graph() %>% # Convert to a SNN graph
#   cluster_graph_louvain() # Perform graph-based clustering
# cat(sprintf("Clusts length: %s\n", length(clusts)))
# clusts[1:10]
# jpeg(file=paste0("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/umap_mc_9.jpeg"),width = 1920,height = 1080)
# plot(plot_embedding(clusts, umap))
# dev.off()
# mc_bpc <- list(raw.couts = MC.fs.ge,
#                norm.counts = mat_norm,
#                pca = pca_embeddings_mc,
#                umap = umap,
#                clusts = clusts)
# saveRDS(mc_bpc,"./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/MC_fs_bpc_g16_fg.rds")

# Extract the cluster assignments (assuming clustering is already done)
clusters_mc <- Idents(MC.fs.seurat)
#clusters_mc <- clusts

so_8 <- Seurat::Load10X_Spatial("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/binned_outputs/square_008um")
so_8 <- SetAssayData(so_8, layer = "data", GetAssayData(so_8,layer = "counts"))
so_8 <- FindVariableFeatures(so_8)
so_8 <- ScaleData(so_8,do.scale = F,do.center = F)
so_8 <- RunPCA(so_8)
so_8 <- FindNeighbors(so_8)
so_8 <- FindClusters(so_8,resolution = 0.5)
so_8 <- RunUMAP(so_8,dims = 1:30)
jpeg(file=paste0("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/umap_vis_8.jpeg"),width = 1920,height = 1080,quality = 100)
plot(DimPlot(so_8))
dev.off()
saveRDS(so_8,"./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/so_8.rds")
pca_embeddings_8 <- Embeddings(so_8, "pca")

# Extract the cluster assignments (assuming clustering is already done)
clusters_8 <- Idents(so_8)

library(cluster)
library(RcppArmadillo)

# Assuming your PCA-reduced data is in a data frame called pca_data
set.seed(42)  # For reproducibility

# Sample mini-batches and compute silhouette scores
num_batches <- 5
batch_size <- 10000
silhouette_scores_mc <- numeric(num_batches)

for (i in 1:num_batches) {
  ids <- sample(nrow(pca_embeddings_mc), batch_size)
  batch_pca <- pca_embeddings_mc[ids]
  batch_clust <- clusters_mc[ids]
  silhouette_scores_mc[i] <- mean(silhouette(as.numeric(batch_clust), dist(batch_pca))[, 3])  # Adjust k based on your needs
}

# Average silhouette score across batches
avg_silhouette_score_mc <- mean(silhouette_scores_mc)

print(paste("Average Silhouette Score:", avg_silhouette_score_mc))
saveRDS(avg_silhouette_score_mc,"./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/sil_score_g16_fg.rds")

set.seed(42)  # For reproducibility

# Sample mini-batches and compute silhouette scores
num_batches <- 5
batch_size <- 10000
silhouette_scores_8 <- numeric(num_batches)

for (i in 1:num_batches) {
  ids <- sample(nrow(pca_embeddings_8), batch_size)
  batch_pca <- pca_embeddings_8[ids]
  batch_clust <- clusters_8[ids]
  silhouette_scores_8[i] <- mean(silhouette(as.numeric(batch_clust), dist(batch_pca))[, 3])  # Adjust k based on your needs
}

# Average silhouette score across batches
avg_silhouette_score_8 <- mean(silhouette_scores_8)

print(paste("Average Silhouette Score:", avg_silhouette_score_8))
saveRDS(avg_silhouette_score_mc,"./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/sil_score_8um_fg.rds")
rm(so_8)

g = 64 # gamma
n.pc = 1:30 # number of first PC to use
k.knn = 8 # number of neighbors to connect to each spot

print("Creating metaspots")
MC.spl <- SCimplify_SpatialDLS_v2(X = mtx.count,
                               spotPositions = spotPositions.filt ,
                               method_similarity = "1",
                               split_not_connected = T,
                               genes.use = NULL,
                               gamma = g,
                               n.pc = n.pc,
                               method_knn = "1",
                               k.knn = k.knn,
                               method_normalization = "log_normalize",
                               n.cpu = 16,
                               return.seurat.object = FALSE,)

print("Done")

metadat_metaspot$MC_membership <- MC.spl$membership %>% as.character()
MC.spl$cell_type <- rep("unk", MC.spl$N.SC)
names(MC.spl$cell_type) <- paste0(sort(unique(MC.spl$membership)))
MC.spl$polygons <- supercell_metaspots_shape_v2(MC = MC.spl,
                                                       spotpositions = spotPositions.filt,
                                                       annotation = "cell_type",
                                                       concavity = 2,membership_name = "membership")
MC.spl$centroids <- supercell_spatial_centroids(MC = MC.spl,spotPositions = spotPositions.filt)
saveRDS(MC.spl,"./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/MC_spl_g64_fg.rds")

MC.fs.ge <- superspot_GE(MC = MC.spl,
                         ge = mtx.count, #%>% as.matrix(),
                         #ge = GetAssayData(MC.full.spl$seurat.object,layer = "data"),
                         groups = as.numeric(MC.spl$membership),
                         mode = "sum"
)

MC.fs.seurat <- supercell_2_Seuratv5(
  SC.GE = MC.fs.ge,
  SC = MC.spl,
  fields = c("cell_type")
)
MC.fs.seurat <- SetAssayData(MC.fs.seurat, layer = "data", GetAssayData(MC.fs.seurat,layer = "counts"))
MC.fs.seurat <- FindVariableFeatures(MC.fs.seurat)
MC.fs.seurat <- ScaleData(MC.fs.seurat,do.scale = F,do.center = F)
MC.fs.seurat <- RunPCA(MC.fs.seurat)
MC.fs.seurat <- FindNeighbors(MC.fs.seurat)
MC.fs.seurat <- FindClusters(MC.fs.seurat,resolution = 0.5)
MC.fs.seurat <- RunUMAP(MC.fs.seurat,dims = 1:30)
jpeg(file=paste0("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/umap_mc_64.jpeg"),width = 1920,height = 1080)
plot(DimPlot(MC.fs.seurat))
dev.off()
#saveRDS(MC.fs.seurat,"/Users/admin/Downloads/square_0016um/MC_fs_seurat.rds")
saveRDS(MC.fs.seurat,"./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/MC_fs_seurat_g64_fg.rds")
pca_embeddings_mc <- Embeddings(MC.fs.seurat, "pca")

# Extract the cluster assignments (assuming clustering is already done)
clusters_mc <- Idents(MC.fs.seurat)

so_16 <- Seurat::Load10X_Spatial("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/binned_outputs/square_016um")
so_16 <- SetAssayData(so_16, layer = "data", GetAssayData(so_16,layer = "counts"))
so_16 <- FindVariableFeatures(so_16)
so_16 <- ScaleData(so_16,do.scale = F,do.center = F)
so_16 <- RunPCA(so_16,)
so_16 <- FindNeighbors(so_16,)
so_16 <- FindClusters(so_16,resolution = 0.5)
so_16 <- RunUMAP(so_16,dims = 1:30)
jpeg(file=paste0("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/umap_vis_16.jpeg"),width = 1920,height = 1080,quality = 100)
plot(DimPlot(so_16))
dev.off()
saveRDS(so_16,"./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/so_16.rds")
pca_embeddings_16 <- Embeddings(so_16, "pca")

# Extract the cluster assignments (assuming clustering is already done)
clusters_16 <- Idents(so_16)

library(cluster)
library(RcppArmadillo)

# Assuming your PCA-reduced data is in a data frame called pca_data
set.seed(42)  # For reproducibility

# Sample mini-batches and compute silhouette scores
num_batches <- 5
batch_size <- 10000
silhouette_scores_mc <- numeric(num_batches)

for (i in 1:num_batches) {
  ids <- sample(nrow(pca_embeddings_mc), batch_size)
  batch_pca <- pca_embeddings_mc[ids]
  batch_clust <- clusters_mc[ids]
  silhouette_scores_mc[i] <- mean(silhouette(as.numeric(batch_clust), dist(batch_pca))[, 3])  # Adjust k based on your needs
}

# Average silhouette score across batches
avg_silhouette_score_mc <- mean(silhouette_scores_mc)

print(paste("Average Silhouette Score:", avg_silhouette_score_mc))
saveRDS(avg_silhouette_score_mc,"./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/sil_score_g64_fg.rds")

set.seed(42)  # For reproducibility

# Sample mini-batches and compute silhouette scores
num_batches <- 5
batch_size <- 10000
silhouette_scores_16 <- numeric(num_batches)

for (i in 1:num_batches) {
  ids <- sample(nrow(pca_embeddings_16), batch_size)
  batch_pca <- pca_embeddings_16[ids]
  batch_clust <- clusters_16[ids]
  silhouette_scores_16[i] <- mean(silhouette(as.numeric(batch_clust), dist(batch_pca))[, 3])  # Adjust k based on your needs
}

# Average silhouette score across batches
avg_silhouette_score_16 <- mean(silhouette_scores_16)

print(paste("Average Silhouette Score:", avg_silhouette_score_16))
saveRDS(avg_silhouette_score_mc,"./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/sil_score_16um_fg.rds")
