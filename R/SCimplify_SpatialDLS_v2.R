#' Normalization of the transcriptomic original data
#'
#' This function normalize the transcriptomic original data based on the method asked by the user. The available methods are "log_normalize", "no" for using the raw counts and "SCT".
#'
#'
#' @param so a Seurat object
#' @param method_normalization a string indicating the method. Choose between "log_normalize", "no" and "SCT".
#'
#' @return a Seurat object
#'
#' @export

normalize.spot <- function(so,method_normalization){
  if (method_normalization == "log_normalize"){
    message(paste0("Performing Log Normalization"))
    so <- Seurat::NormalizeData(so)
    message(paste0("Done"))
    so <- Seurat::FindVariableFeatures(so)
    so <- Seurat::ScaleData(so)
    }
  else if (method_normalization == "no"){
    message(paste0("Using Raw Data"))
    so<-Seurat::SetAssayData(so,new.data = Seurat::GetAssayData(so,assay = "RNA",layer = "counts"),assay = "RNA",layer = "data")
    message(paste0("Done"))
    message(paste0("Running PCA"))
    so <- Seurat::FindVariableFeatures(so)
    so <- Seurat::ScaleData(so,do.scale = FALSE,do.center = FALSE)
    }
  else if (method_normalization == "SCT"){
    message(paste0("Performing SCT"))
    so <- Seurat::SCTransform(so, new.assay.name = "SCT" , assay = "RNA",verbose = F)
    message(paste0("Done"))
    }
  return(so)
}

#' Dimension Reduction of the transcriptomic original data
#'
#' This function performs dimensionality reduction of the transcriptomic original data based on the method asked by the user. The available methods are "PCA", "NMF" and "genes.list" for using directly the expression of selected genes instead of other methods for dimensionality reduction.
#'
#'
#' @param so a Seurat object
#' @param method_reduction a string indicating the method. Choose between "log_normalize", "no" and "SCT".
#' @param n.dim number of principal components to use from PCA
#' @param max_gamma value of the maximum reachable gamma
#' @param genes.use genes to use when "method_reduction = 'genes.list'"
#'
#' @return a Seurat object
#'
#' @export

dimreduc.spot <- function(so,method_reduction,n.dim,max_gamma,genes.use){
  if (method_reduction == "PCA"){
    message(paste0("Running PCA"))
    if (max_gamma != 1){
      so <- Seurat::RunPCA(so,assay = Seurat::DefaultAssay(so))
      #pca_matrix <- so@reductions[["pca"]]@cell.embeddings[,1:n.pc]
      dimreduc.matrix <- so@reductions[["pca"]]@cell.embeddings[,n.dim]
    }
    else {
      dimreduc.matrix <- matrix(0,ncol(so),n.dim)
    }
  }
  else if (method_reduction == "NMF"){
    message(paste0("Running NMF"))
    if (max_gamma != 1){
      so <- GeneNMF::runNMF(so,Seurat::DefaultAssay(so),k = max(n.dim))
      #pca_matrix <- so@reductions[["pca"]]@cell.embeddings[,1:n.pc]
      dimreduc.matrix <- so@reductions[["NMF"]]@cell.embeddings[,n.dim]
    }
    else {
      dimreduc.matrix <- matrix(0,ncol(so),n.dim)
    }
  }
  else if (method_reduction == "genes.list"){
    message(paste0("Selecting genes"))
    if (max_gamma != 1){
      dimreduc.matrix <- Seurat::GetAssayData(so,assay = Seurat::DefaultAssay(so))[genes.use,]
    }
    else {
      dimreduc.matrix <- matrix(0,ncol(so),n.dim)
    }
  }
  return(dimreduc.matrix)
}






#' Creation of the KNN graph to links the spots based on the spatial proximity
#'
#' This function builds the KNN graph to links the spots based on the spatial proximity
#'
#'
#' @param spotPositions a data.frame of the coordinates of the spots
#' @param col.x a character of the name of the x coordinates to get from the coordinates data frame
#' @param col.y a character of the name of the y coordinates to get from the coordinates data frame
#' @param k.spots number of neighborhoods to take into account when choosing the second method to build the KNN
#' @param countMatrix raw matrix of the gene expression for each spot
#' @param n.dim number of principal components to use from PCA
#' @param n.cpu number of cpu to use during parallelized computation of distances. By default, the maximum amount of cpu available is chosen automatically. But if your computer doesn't support it, please specify your desired number of cpu.
#' @param method_similarity method to compute transcriptional similarity from distance
#' @param plot.graph boolean of whether or not plotting the KNN graph
#' @param pct percentage of the connections to keep
#' \itemize{
#'   \item "1" - similarity = 1 - (distance/max(distance))
#'   \item "2" - similarity = exp(-(distance^2)/sigma)
#'   \item "3" - similarity = proxy::pr_dist2simil(distance)
#' }
#' @param dist.thresh fixed value to use instead of using a percentage to filter out connections
#' @param method_knn method to compute the KNN
#' \itemize{
#'   \item "1" - KNN is built with RANN::nn2() function (default)
#'   \item "2" - KNN is built by computing every distances between each spot and filtered by choosing the number of neighborhoods to take into account
#' }
#' @param k.knn numbers of neighboring spots to compute for each spot
#' @param sigs parameter when computing the transcriptional similarity with method "2"
#' @param method_normalization method for normalization of the raw count matrix
#' \itemize{
#'   \item "log_normalize" for a Seurat Log Normalization
#'   \item "SCT" for a Seurat SCT Normalization
#'   \item "no" for using raw counts
#' }
#' @param method_reduction a string indicating the method. Choose between "log_normalize", "no" and "SCT".
#' @param genes.use genes to use when "method_reduction = 'genes.list'"
#' @param split_not_connected split or not the spots in the KNN graph if their distance is too far compared to the others
#' @param split_annotation split or not the spots in the KNN graph based on their annotation
#' @param split_vector vector of the spots annotations on which the KNN will be split
#'
#' @return a list with components
#' \itemize{
#'   \item graph - an igraph graph
#'   \item so - a seurat object
#' }
#'
#' @export
#'
#'
#'

neighbor_graph_v2 <- function(spotPositions,
                              col.x,
                              col.y,
                              k.spots, countMatrix,
                              n.dim,
                              method_similarity,
                              method_knn,
                              k.knn,
                              sigs,
                           method_normalization,
                           split_not_connected,
                           split_annotation,
                           split_vector,
                           n.cpu,
                           plot.graph,
                           pct,dist.thresh,
                           method_reduction,
                           genes.use){
  if (method_knn == "1"){
    message(paste0("Building KNN graph with nn2"))
    nn2.res <- RANN::nn2(data = spotPositions,k = k.knn+1)
    if (split_not_connected == TRUE ){
      if (is.null(dist.thresh)){
        min_dist <- stats::quantile(as.vector(nn2.res$nn.dists[,2:k.knn+1]),names = F,probs = pct)
      }
      else {
        min_dist <- dist.thresh
      }

      message(paste0("Neighbors with distance > ",min_dist, " are removed"))
      plot(ggplot2::ggplot(data = tibble::tibble(distances = as.vector(nn2.res$nn.dists[,2:k.knn+1]),
                                                 distribution = rep(".",length(as.vector(nn2.res$nn.dists[,2:k.knn+1])))))+
             ggplot2::geom_violin(ggplot2::aes(x = distribution, y = distances))+
             ggplot2::geom_boxplot(ggplot2::aes(x = distribution, y = distances),width = 0.1)+
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
    message(paste0("Maximum gamma is ",max_gamma))
    message(paste0("Done"))
  }
  else if (method_knn == "2"){
    message(paste0("Building KNN graph by DLS"))
    message(paste0("Computing distances between spots"))
    distCoord <- parallelDist::parDist(as.matrix(spotPositions),threads = n.cpu)
    min.dist <- min(distCoord)
    distCoord[distCoord > k.spots*round(min.dist)] <- 0
    distCoord[distCoord >0] <- 1
    spot.graph <- igraph::graph.adjacency(as.matrix(distCoord),mode = "undirected")
    message(paste0("Done"))
  }

  if (plot.graph == TRUE){
    igraph::V(spot.graph)$x <- spotPositions[,col.x]
    igraph::V(spot.graph)$y <- -spotPositions[,col.y]
    plot(spot.graph,
         layout = matrix(c(igraph::V(spot.graph)$x,igraph::V(spot.graph)$y), ncol = 2),
         vertex.label.cex = 0.1,
         vertex.size = 1,
         edge.arrow.size = 0.7,
         main = "Graph with Custom X and Y Coordinates")
  }


  so <- Seurat::CreateSeuratObject(countMatrix,assay = "RNA")
  so <- normalize.spot(so,method_normalization)
  dimreduc.matrix <- dimreduc.spot(so,method_reduction,n.dim, max_gamma,genes.use)

  message(paste0("Done"))
  message(paste0("Computing distance for each edge"))
  dimreduc_dist <- pbapply::pbapply(igraph::as_edgelist(spot.graph),
                               1 ,
                               function(x){parallelDist::parDist(dimreduc.matrix[x,],threads = n.cpu)})
  message(paste0("Done"))
  message(paste0("Computing similarity from distances"))
  if (method_similarity == "1"){
    dimreduc_similarity <-  unlist(pbapply::pblapply(dimreduc_dist,FUN = function(i){1 - (i/max(dimreduc_dist))}))
    dimreduc_similarity[dimreduc_similarity == 0 ] <- 1e-100
    message(paste0("Done"))
  }
  else if (method_similarity == "2"){
    dimreduc_similarity <- unlist(pbapply::pblapply(dimreduc_dist,FUN = function(i){exp(-(i**2)/sigs)}))
    dimreduc_similarity[dimreduc_similarity == 0 ]<- 1e-100
    message(paste0("Done"))
  }

  else if (method_similarity == "3"){
    dimreduc_similarity <- proxy::pr_dist2simil(dimreduc_dist)
    dimreduc_similarity[dimreduc_similarity == 0 ]<- 1e-100
    message(paste0("Done"))
  }
  message(paste0("Returning graph with similarity as weight"))
  igraph::E(spot.graph)$weight <- dimreduc_similarity
  return_list <- list("graph" = spot.graph, "seurat.object" = so)
  return(return_list)
}



#' Detection of metaspots with the SuperSpot approach
#'
#' This function detects metaspots  from spatial transcriptomics data
#'
#'
#' @param X raw count matrix with rows to be genes and cols to be cells
#' @param spotPositions a data.frame of the coordinates of the spots
#' @param col.x a character of the name of the x coordinates to get from the coordinates data frame
#' @param col.y a character of the name of the y coordinates to get from the coordinates data frame
#' @param k.spots number of neighborhoods to take into account when choosing the second method to build the KNN
#' @param n.dim number of dimensions to compute for dimension reduction
#' @param method_similarity method to compute transcriptional similarity from distance
#' \itemize{
#'   \item "1" - similarity = 1 - (distance/max(distance))
#'   \item "2" - similarity = exp(-(distance^2)/sigma)
#'   \item "3" - similarity = proxy::pr_dist2simil(distance)
#' }
#' @param method_knn method to compute the KNN
#' \itemize{
#'   \item "1" - KNN is built with RANN::nn2() function (default)
#'   \item "2" - KNN is built by computing every distances between each spot and filtered by choosing the number of neighborhoods to take into account
#' }
#' @param k.knn numbers of neighboring spots to compute for each spot
#' @param sigs parameter when computing the transcriptional similarity with method "2"
#' @param method_normalization method for normalization of the raw count matrix
#' \itemize{
#'   \item "log_normalize" for a Seurat Log Normalization
#'   \item "SCT" for a Seurat SCT Normalization
#'   \item "no" for using raw counts
#' }
#' @param method_reduction method for performing dimension reduction on which the weights for the KNN will be computed
#' \itemize{
#'   \item "PCA" for using PCA
#'   \item "NMF" for using NMF
#'   \item "genes.list" for directly using a vector of genes instead of dimension reduction
#' }
#' @param genes.use genes to use when "method_reduction = 'genes.list'"
#' @param split_not_connected split or not the spots in the KNN graph if their distance is too far compared to the others
#' @param split_annotation split or not the spots in the KNN graph based on their annotation
#' @param split_vector vector of the spots annotations on which the KNN will be split
#' @param cell.annotation a vector of cell type annotation, if provided, metacells that contain single cells of different cell type annotation will be split in multiple pure metacell (may result in slightly larger numbe of metacells than expected with a given gamma)
#' @param cell.split.condition a vector of cell conditions that must not be mixed in one metacell. If provided, metacells will be split in condition-pure metacell (may result in significantly(!) larger number of metacells than expected)
#' @param gamma graining level of data (proportion of number of single cells in the initial dataset to the number of metacells in the final dataset)
#' @param k.knn parameter to compute single-cell kNN network
#' @param n.cpu number of cpu to use during parallelized computation of distances. By default, the maximum amount of cpu available is chosen automatically. But if your computer doesn't support it, please specify your desired number of cpu.
#' @param plot.graph boolean of whether or not plotting the KNN graph
#' @param pct percentage of the connections to keep
#' @param dist.thresh fixed value to use instead of using a percentage to filter out connections
#' @param seed seed to use to subsample cells for an approximate approach
#' @param igraph.clustering clustering method to identify metacells (available methods "walktrap" (default) and "fast greedy" (recommended for VisiumHD so far)).
#' @param return.singlecell.NW whether return single-cell network
#' @param return.hierarchical.structure whether return hierarchical structure of metacell
#' @param return.seurat.object whether return seurat object of spots used to computed PCA
#' @param ... other parameters of \link{build_knn_graph} function
#'
#' @return a list with components
#' \itemize{
#'   \item graph.supercells - igraph object of a simplified network (number of nodes corresponds to number of metacells)
#'   \item membership - assigmnent of each single cell to a particular metacell
#'   \item graph.singlecells - igraph object (kNN network) of single-cell data
#'   \item supercell_size - size of metacells (former super-cells)
#'   \item gamma - requested graining level
#'   \item N.SC - number of obtained metacells
#'   \item n.pc - number of principal components used for metacells construction
#'   \item k.knn - number of neighbors to build single-cell graph
#'   \item k.spots - number of neighborhood used if method_knn was "2"
#'   \item sigs - sigma parameter used if method_similarity was "2"
#'   \item method_similarity - used similarity method
#'   \item method_knn - used KNN method
#'   \item method_normalization - used normalization method
#'   \item sc.cell.annotation. - single-cell cell type annotation (if provided)
#'   \item sc.cell.split.condition. - single-cell split condition (if provided)
#'   \item SC.cell.annotation. - super-cell cell type annotation (if was provided for single cells)
#'   \item SC.cell.split.condition. - super-cell split condition (if was provided for single cells)
#'
#' }
#'
#' @export
#'
#'
#'


SCimplify_SpatialDLS_v2 <- function(X,
                                 spotPositions,
                                 col.x = "imagecol",
                                 col.y = "imagerow",
                                 k.spots = 1,
                                 method_similarity = "3",
                                 method_knn = "1",
                                 method_normalization = "log_normalize",
                                 method_reduction = "PCA",
                                 split_not_connected = TRUE,
                                 split_annotation = NULL,
                                 split_vector,
                                 sigs = 1000,
                                 cell.annotation = NULL,
                                 cell.split.condition = NULL,
                                 gamma = 10,
                                 k.knn = 6,
                                 n.dim = 1:30,
                                 n.cpu = NULL,
                                 plot.graph = FALSE,
                                 pct = 0.6,
                                 dist.thresh = NULL,
                                 seed = 12345,
                                 genes.use = NULL,
                                 igraph.clustering = c("walktrap", "fast_greedy"),
                                 return.singlecell.NW = TRUE,
                                 return.hierarchical.structure = TRUE,
                                 return.seurat.object = FALSE,
                                 ...){

  neighbor_graph.output <- neighbor_graph_v2(spotPositions = spotPositions, col.x = col.x, col.y = col.y, k.spots = k.spots,split_not_connected = split_not_connected,
                                          countMatrix = X, n.dim = n.dim, method_similarity = method_similarity,method_reduction = method_reduction,
                                          method_knn = method_knn, k.knn = k.knn, sigs = sigs,method_normalization = method_normalization,
                                          split_annotation = split_annotation,split_vector =split_vector, n.cpu = n.cpu,plot.graph = plot.graph,
                                          pct = pct,dist.thresh = dist.thresh,genes.use=genes.use)

  X <- Seurat::GetAssayData(neighbor_graph.output$seurat.object)

  N.c <- ncol(X)

  presampled.cell.ids <- colnames(X)



  sc.nw <- neighbor_graph.output$graph

  #simplify

  k   <- round(N.c/gamma)
  message(paste0("Clustering"))
  if(igraph.clustering[1] == "walktrap"){
    g.s              <- igraph::cluster_walktrap(sc.nw,weights = igraph::E(sc.nw)$weight)
    g.s$membership   <- igraph::cut_at(g.s, k)

  } else if(igraph.clustering[1] == "fast_greedy") {
    g.s              <- igraph::cluster_fast_greedy(sc.nw,weights = igraph::E(sc.nw)$weight)
    g.s$membership   <- igraph::cut_at(g.s, k)

  } else {
    stop(paste("Unknown clustering method (", igraph.clustering, "), please use fast_greedy or walkrtap"))
  }
  message(paste0("Done"))
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
              n.dim = n.dim,
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
    res$SC.cell.annotation. <- SuperCell::supercell_assign(cell.annotation, res$membership)
    res$SC.cell.split.condition. <- SuperCell::supercell_assign(cell.split.condition, res$membership)
  }


  if(return.hierarchical.structure){ res$h_membership <- g.s}
  gc()

  if(return.seurat.object == TRUE){res$seurat.object = neighbor_graph.output$seurat.object}

  return(res)
}
