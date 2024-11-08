library(Seurat)
library(SuperSpot)
library(SuperCell)
library(tidyr)
library(tidyverse)
library(igraph)




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
  print("Computing PCA using BPCells")
  print(Sys.time())
  svd <- BPCells::svds(countMatrix, k=25)
  # Alternate option: irlba::irlba(mat_norm, nv=50)
  pca_matrix <- multiply_cols(svd$v, svd$d)
  print("Done")
  print(Sys.time())
  cat(sprintf("PCA dimensions: %s\n", toString(dim(pca_matrix))))
  pca_matrix[1:4,1:3]

  print(paste0("Done"))
  print(paste0("Computing PCA dist for each edge"))
  print(Sys.time())
  pca_similarity <- pbapply::pbapply(igraph::get.edgelist(spot.graph),
                                     1 ,
                                     function(x){1/parallelDist::parDist(pca_matrix[x,],threads = n.cpu)})
  print(head(pca_similarity))
  pca_similarity <- pca_similarity%>% unlist()
  print(paste0("Dimensions of pca distances: ",dim(pca_similarity)))
  print(paste0("Done"))
  print(paste0("Computing similarity from PCA distances"))
  print(Sys.time())
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


## Import data
spotPositions <- arrow::read_parquet("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/binned_outputs/square_002um/spatial/tissue_positions.parquet") %>% column_to_rownames("barcode")
spotPositions.filt <- subset(spotPositions, in_tissue == 1)
rm(spotPositions)
mtx.count <- Seurat::Read10X_h5("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/binned_outputs/square_002um/filtered_feature_bc_matrix.h5")
#mtx.count <- open_matrix_10x_hdf5("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/binned_outputs/square_002um/filtered_feature_bc_matrix.h5", feature_type="Gene Expression") %>%
#write_matrix_dir("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/binned_outputs/square_002um/filtered_feature_bc_matrix_raw")
spotPositions.filt <- spotPositions.filt[colnames(mtx.count),c("pxl_col_in_fullres","pxl_row_in_fullres")]
colnames(spotPositions.filt) <- c("imagecol","imagerow")

options(future.globals.maxSize = 1e12)
metadat_metaspot <- tibble(spots = rownames(spotPositions.filt))
for (g in c(64,32,16,8,4,2,1)){
  n.pc = 1:30 # number of first PC to use
  k.knn = 8 # number of neighbors to connect to each spot
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
  #MC.spl$polygons <- supercell_metaspots_shape_v2(MC = MC.spl,
  #                                                spotpositions = spotPositions.filt,
  #                                                annotation = "cell_type",
  #                                                concavity = 2,membership_name = "membership")
  MC.spl$centroids <- supercell_spatial_centroids(MC = MC.spl,spotPositions = spotPositions.filt)
  #saveRDS(MC.spl,"./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/MC_spl_g16_fg.rds")

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
  #saveRDS(MC.seurat,file = paste0("./SuperSpot/01_Data/MC.seurat_g",g,".rds"))
  #tissue_position  <- read.csv("./SuperSpot/BCBA/spatial/tissue_positions_list.csv",header = F)
  #tissue_position  <- arrow::read_parquet("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/binned_outputs/square_002um/spatial/tissue_positions.parquet")
  #MC_centroids <- supercell_spatial_centroids(MC,spotPositions = spotPosition)
  #tissue_position_sb <- subset(tissue_position, barcodes %in% names(MC.spl$membership))
  #tissue_position_sb$membership <- MC.spl$membership

  #rownames(tissue_position_sb) <- NULL
  #metadata_filt <- metadata_filt %>% column_to_rownames("cell")

  #tissue_position_mean <- stats::aggregate(tissue_position_sb[,2:7] %>% as.matrix() ~membership,tissue_position_sb[,2:7],mean)

  #metadata_mean <- stats::aggregate(subset(metadata_filt,select = -c(cell_id,version,Run_name,Run_Tissue_name,tissue,Panel,assay_type)) %>% as.matrix() ~membership,metadata_filt,mean)

  #metadata_mean$fov <- as.integer(metadata_mean$fov)

  #metadata_mean$cell_ID <- metadata_mean$membership

  image <- png::readPNG(source = "./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/binned_outputs/square_002um/spatial/tissue_lowres_image.png")
  scale.factors <- Read10X_ScaleFactors(filename = "./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/binned_outputs/square_002um/spatial/scalefactors_json.json")
  #coordinates <- tissue_position_mean[,1:6]
  #coordinates <- column_to_rownames(coordinates,var = "membership")
  #colnames(coordinates) <- c("tissue","row","col","imagerow","imagecol")
  #coordinates$imagecol <- MC.spl$centroids$imagecol
  #coordinates$imagerow <- MC.spl$centroids$imagerow
  #coordinates <- Read10X_Coordinates(filename = Sys.glob(file.path(image.dir,
  #                                                                 "*tissue_positions*")), filter.matrix)
  fov <- CreateFOV(MC.spl$centroids[, c("imagerow", "imagecol")],
                   type = "centroids", radius = scale.factors[["spot"]],
                   assay = "RNA", key = Key("slice1", quiet = TRUE))
  visium.fov <- new(Class = "VisiumV2", boundaries = fov@boundaries,
                    molecules = fov@molecules, assay = fov@assay, key = fov@key,
                    image = image, scale.factors = scale.factors)

  MC.fs.seurat[["slice1"]] <- visium.fov
  #write_csv(x=MC_centroids[,c("imagecol","imagerow")],file = paste0("./SuperSpot/01_Data/PancreasCosMx/Centroids_metaspot_g",g,".csv"))
  #write_csv(x=metadata_mean,file = paste0("./SuperSpot/01_Data/PancreasCosMx/Pancreas_metaspot_g",g,".csv"))

  #counts.matrix_ms <- t(GetAssayData(MC.seurat ,assay = "RNA") %>% as.matrix()) %>% as.data.frame() %>% rownames_to_column("cell_ID")

  #counts.matrix_ms$fov <- metadata_mean$fov

  #counts.matrix_ms <- counts.matrix_ms[columns_to_keep]

  #write_csv(x=counts.matrix_ms,file = paste0("./SuperSpot/01_Data/PancreasCosMx/count_ms_g",g,".csv"))
  saveRDS(MC.fs.seurat,file = paste0("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/benchmarking/MC_seurat_g",g,".rds"))
}

#write.csv(metad,"./SuperSpot/01_Data/metad.csv")
