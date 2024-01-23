#' Creation of the KNN graph to links the spots based on the spatial proximity
#'
#' This function builds the KNN graph to links the spots based on the spatial proximity
#'
#'
#' @param spotpositions a data.frame of the coordinates of the spots
#' @param k.spots number of neighborhoods to take into account when choosing the second method to build the KNN
#' @param countMatrix raw matrix of the gene expression for each spot
#' @param n.pc number of principal components to use from PCA
#' @param n.cpu number of cpu to use during parallelized computation of distances. By default, the maximum amount of cpu available is chosen automatically. But if your computer doesn't support it, please specify your desired number of cpu.
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
#' }
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

neighbor_graph <- function(spotPositions, k.spots, countMatrix, n.pc, method_similarity, method_knn, k.knn, sigs,
                           method_normalization,split_not_connected,split_annotation,split_vector,n.cpu){
  if (method_knn == "1"){
    print(paste0("Building KNN graph with nn2"))
    nn2.res <- RANN::nn2(data = spotPositions,k = k.knn+1)
    if (split_not_connected == TRUE){
      min_dist <- quantile(nn2.res$nn.dists[,2:k.knn+1] %>% as.vector(),names = F)[3]
      print(paste0("Neighbors with distance > ",min_dist, " are removed"))
      plot(ggplot(data = tibble(distances = as.vector(nn2.res$nn.dists[,2:k.knn+1]),
                                distribution = rep(".",length(as.vector(nn2.res$nn.dists[,2:k.knn+1])))))+
             geom_violin(aes(x = distribution, y = distances))+
             geom_boxplot(aes(x = distribution, y = distances),width = 0.1)+
             geom_hline(yintercept=min_dist, linetype="dashed", color = "red"))

      bad_neighbors_rows <- which(nn2.res$nn.dists > round(min_dist+0.6),arr.ind = T)[,1]
      bad_neighbors_cols <- which(nn2.res$nn.dists > round(min_dist+0.6),arr.ind = T)[,2]
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
      spot.graph <- subgraph(spot.graph,cell_type.id)
      spotPositions <- spotPositions[cell_type.id,]
      countMatrix <- countMatrix[,cell_type.id]
    }

    max_gamma <- (nrow(spotPositions)/igraph::components(spot.graph)$no)
    print(paste0("Maximum gamma is ",max_gamma))
    print(paste0("Done"))
    V(spot.graph)$x <- spotPositions$imagecol
    V(spot.graph)$y <- -spotPositions$imagerow
  }
  else if (method_knn == "2"){
    print(paste0("Building KNN graph by DLS"))
    print(paste0("Computing distances between spots"))
    distCoord <- parallelDist::parDist(spotPositions %>% as.matrix(),threads = n.cpu)
    min.dist <- min(distCoord)
    distCoord[distCoord > k.spots*round(min.dist)] <- 0
    distCoord[distCoord >0] <- 1
    spot.graph <- graph.adjacency(distCoord %>% as.matrix(),mode = "undirected")
    print(paste0("Done"))
  }

  so <- CreateSeuratObject(countMatrix,assay = "RNA")
  if (method_normalization == "log_normalize"){
    print(paste0("Performing Log Normalization"))
    so <- NormalizeData(so)
    print(paste0("Done"))
    print(paste0("Running PCA"))
    so <- FindVariableFeatures(so)
    so <- ScaleData(so)
    if (max_gamma != 1){
      so <- RunPCA(so,verbose = F)
      pca_matrix <- so@reductions[["pca"]]@cell.embeddings[,1:n.pc]
    }
    else {
      pca_matrix <- matrix(0,ncol(countMatrix),n.pc)
    }
  }
  else if (method_normalization == "SCT"){
    print(paste0("Performing SCT"))
    so <- SCTransform(so, new.assay.name = "SCT" , assay = "RNA",verbose = F)
    print(paste0("Done"))
    print(paste0("Running PCA"))
    so <- RunPCA(so, assay = "SCT",verbose = F)
    pca_matrix <- so@reductions[["pca"]]@cell.embeddings[,1:n.pc]
  }

  print(paste0("Done"))
  print(paste0("Computing PCA dist for each edge"))
  pca_dist <- pbapply::pbapply(igraph::get.edgelist(spot.graph),
                    1 ,
                    function(x){parallelDist::parDist(pca_matrix[x,],threads = n.cpu)})
  print(paste0("Done"))
  print(paste0("Computing similarity from PCA distances"))
  if (method_similarity == "1"){
    pca_similarity <- pbapply::pblapply(pca_dist,FUN = function(i){1 - (i/max(pca_dist))}) %>% unlist()
    pca_similarity[pca_similarity == 0 ] <- 1e-100
    print(paste0("Done"))
  }
  else if (method_similarity == "2"){
    pca_similarity <- pbapply::pblapply(pca_dist,FUN = function(i){exp(-(i**2)/sigs)}) %>% unlist()
    pca_similarity[pca_similarity == 0 ]<- 1e-100
    print(paste0("Done"))
  }

  else if (method_similarity == "3"){
    pca_similarity <- proxy::pr_dist2simil(pca_dist)
    pca_similarity[pca_similarity == 0 ]<- 1e-100
    print(paste0("Done"))
  }
  print(paste0("Returning graph with PCA similarity as weight"))
  igraph::E(spot.graph)$weight <- pca_similarity
  return_list <- list("graph" = spot.graph, "seurat.object" = so)
  return(return_list)
}

#' Detection of metaspots with the SuperSpot approach
#'
#' This function detects metaspots  from spatial transcriptomics data
#'
#'
#' @param X raw count matrix with rows to be genes and cols to be cells
#' @param spotpositions a data.frame of the coordinates of the spots
#' @param k.spots number of neighborhoods to take into account when choosing the second method to build the KNN
#' @param n.pc number of principal components to use from PCA
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
#' }
#' @param split_not_connected split or not the spots in the KNN graph if their distance is too far compared to the others
#' @param split_annotation split or not the spots in the KNN graph based on their annotation
#' @param split_vector vector of the spots annotations on which the KNN will be split
#' @param genes.use a vector of genes used to compute PCA
#' @param genes.exclude a vector of genes to be excluded when computing PCA
#' @param cell.annotation a vector of cell type annotation, if provided, metacells that contain single cells of different cell type annotation will be split in multiple pure metacell (may result in slightly larger numbe of metacells than expected with a given gamma)
#' @param cell.split.condition a vector of cell conditions that must not be mixed in one metacell. If provided, metacells will be split in condition-pure metacell (may result in significantly(!) larger number of metacells than expected)
#' @param n.var.genes if \code{"genes.use"} is not provided, \code{"n.var.genes"} genes with the largest variation are used
#' @param gamma graining level of data (proportion of number of single cells in the initial dataset to the number of metacells in the final dataset)
#' @param k.knn parameter to compute single-cell kNN network
#' @param do.scale whether to scale gene expression matrix when computing PCA
#' @param n.pc number of principal components to use for construction of single-cell kNN network
#' @param n.cpu number of cpu to use during parallelized computation of distances. By default, the maximum amount of cpu available is chosen automatically. But if your computer doesn't support it, please specify your desired number of cpu.
#' @param fast.pca use \link[irlba]{irlba} as a faster version of prcomp (one used in Seurat package)
#' @param do.approx compute approximate kNN in case of a large dataset (>50'000)
#' @param approx.N number of cells to subsample for an approximate approach
#' @param block.size number of cells to map to the nearest metacell at the time (for approx coarse-graining)
#' @param seed seed to use to subsample cells for an approximate approach
#' @param igraph.clustering clustering method to identify metacells (available methods "walktrap" (default) and "louvain" (not recommended, gamma is ignored)).
#' @param return.singlecell.NW whether return single-cell network (which consists of approx.N if \code{"do.approx"} or all cells otherwise)
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
#'   \item genes.use - used genes
#'   \item do.approx - whether approximate coarse-graining was perfirmed
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


SCimplify_SpatialDLS <- function(X,
                                spotPositions,
                                k.spots = 1,
                                method_similarity = "3",
                                method_knn = "1",
                                method_normalization = "log_normalize",
                                split_not_connected,
                                split_annotation = NULL,
                                split_vector,
                                sigs = 1000,
                                genes.use = NULL,
                                genes.exclude = NULL,
                                cell.annotation = NULL,
                                cell.split.condition = NULL,
                                n.var.genes = min(1000, nrow(X)),
                                gamma = 10,
                                k.knn = 6,
                                do.scale = TRUE,
                                n.pc = 5,
                                n.cpu = NULL,
                                fast.pca = TRUE,
                                do.approx = FALSE,
                                approx.N = 20000,
                                block.size = 10000,
                                seed = 12345,
                                igraph.clustering = c("walktrap", "louvain"),
                                return.singlecell.NW = TRUE,
                                return.hierarchical.structure = TRUE,
                                return.seurat.object = FALSE,
                                ...){

  neighbor_graph.output <- neighbor_graph(spotPositions = spotPositions, k.spots = k.spots,split_not_connected = split_not_connected,
                                          countMatrix = X, n.pc = n.pc, method_similarity = method_similarity,
                                          method_knn = method_knn, k.knn = k.knn, sigs = sigs,method_normalization = method_normalization,
                                          split_annotation = split_annotation,split_vector =split_vector, n.cpu = n.cpu)

  X <- GetAssayData(neighbor_graph.output$seurat.object)

  N.c <- ncol(X)

  presampled.cell.ids <- colnames(X)



  sc.nw <- neighbor_graph.output$graph

  #simplify

  k   <- round(N.c/gamma)
  print(paste0("Clustering"))
  if(igraph.clustering[1] == "walktrap"){
    g.s              <- igraph::cluster_walktrap(sc.nw,weights = igraph::E(sc.nw)$weight)
    g.s$membership   <- igraph::cut_at(g.s, k)

  } else if(igraph.clustering[1] == "louvain") {
    warning(paste("igraph.clustering =", igraph.clustering, ", gamma is ignored"))
    g.s    <- igraph::cluster_louvain(sc.nw)

  } else {
    stop(paste("Unknown clustering method (", igraph.clustering, "), please use louvain or walkrtap"))
  }
  print(paste0("Done"))
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
  if(!do.approx){
    SC.NW                        <- igraph::simplify(SC.NW, remove.loops = T, edge.attr.comb="sum")
  }


  if(do.approx){

    PCA.averaged.SC      <- as.matrix(Matrix::t(supercell_GE(t(PCA.presampled$x[,n.pc]), groups = membership.presampled)))
    X.for.roration       <- Matrix::t(X[genes.use, rest.cell.ids])



    if(do.scale){ X.for.roration <- scale(X.for.roration) }
    X.for.roration[is.na(X.for.roration)] <- 0


    membership.omitted   <- c()
    if(is.null(block.size) | is.na(block.size)) block.size <- 10000

    N.blocks <- length(rest.cell.ids)%/%block.size
    if(length(rest.cell.ids)%%block.size > 0) N.blocks <- N.blocks+1


    if(N.blocks>0){
      for(i in 1:N.blocks){ # compute knn by blocks
        idx.begin <- (i-1)*block.size + 1
        idx.end   <- min(i*block.size,  length(rest.cell.ids))

        cur.rest.cell.ids    <- rest.cell.ids[idx.begin:idx.end]

        PCA.ommited          <- X.for.roration[cur.rest.cell.ids,] %*% PCA.presampled$rotation[, n.pc] ###

        D.omitted.subsampled <- proxy::dist(PCA.ommited, PCA.averaged.SC) ###

        membership.omitted.cur        <- apply(D.omitted.subsampled, 1, which.min) ###
        names(membership.omitted.cur) <- cur.rest.cell.ids ###

        membership.omitted   <- c(membership.omitted, membership.omitted.cur)
      }
    }

    membership.all       <- c(membership.presampled, membership.omitted)
    #membership.all       <- membership.all[colnames(X)]


    names_membership.all <- names(membership.all)
    ## again split super-cells containing cells from different annotation or split conditions
    if(!is.null(cell.annotation) | !is.null(cell.split.condition)){

      split.cells <- interaction(cell.annotation[names_membership.all],
                                 cell.split.condition[names_membership.all], drop = TRUE)


      membership.all.intr <- interaction(membership.all, split.cells, drop = TRUE)

      membership.all.intr.v <- as.vector(membership.all.intr)
      membership.all.intr.v.u <- unique(membership.all.intr.v)

      ## add new nodes to SC.NW
      adj <- igraph::get.adjlist(SC.NW, mode = "all")

      add.node.id <- igraph::vcount(SC.NW) + 1
      membership.all.const <- membership.all

      for(i in sort(unique(membership.all.const))){

        cur.sc <- which(membership.all == i)
        cur.main.node <- membership.all.intr.v[cur.sc[1]]
        n.add.nodes <- length(unique(membership.all.intr.v[cur.sc])) - 1

        additional.nodes <- setdiff(unique(membership.all.intr.v[cur.sc]), cur.main.node)

        a.n <- 1
        if(n.add.nodes > 0){
          f.node.id <- add.node.id
          l.node.id <- add.node.id + n.add.nodes -1

          for(j in f.node.id:l.node.id){


            membership.all[membership.all.intr.v == additional.nodes[a.n]] <- j
            a.n <- a.n+1
            adj[[j]] <- c(as.numeric(adj[[i]]), i, f.node.id:l.node.id) # split super-cell node by adding additional node and connecting it to the same neighbours
            add.node.id <- add.node.id + 1
          }

          adj[[i]] <- c(as.numeric(adj[[i]]), f.node.id:l.node.id)
        }
      }


      SC.NW                        <- igraph::graph_from_adj_list(adj, duplicate = F)
      SC.NW                        <- igraph::as.undirected(SC.NW)


    }
    SC.NW                        <- igraph::simplify(SC.NW, remove.loops = T, edge.attr.comb="sum")
    names(membership.all) <- names_membership.all
    membership.all <- membership.all[colnames(X)]

  } else {
    membership.all       <- membership.presampled[colnames(X)]
  }
  membership       <- membership.all

  supercell_size   <- as.vector(table(membership))

  igraph::E(SC.NW)$width         <- sqrt(igraph::E(SC.NW)$weight/10)
  igraph::V(SC.NW)$size          <- supercell_size
  igraph::V(SC.NW)$sizesqrt      <- sqrt(igraph::V(SC.NW)$size)

  res <- list(graph.supercells = SC.NW,
              gamma = gamma,
              N.SC = length(unique(membership)),
              membership = membership,
              supercell_size = supercell_size,
              genes.use = genes.use,
              simplification.algo = igraph.clustering[1],
              do.approx = do.approx,
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

  if(return.seurat.object == TRUE){seurat.object = neighbor_graph.output$seurat.object}

  return(res)
}

#' Computation of the spatial centroids about metaspots
#'
#' This function computes the spatial centroids of metaspots from the coordinates of the original spots
#'
#'
#' @param spotpositions a data.frame of the coordinates of the spots
#' @param MC Metaspot object obtained from SCimplify_SpatialDLS function
#'
#' @return a dataframe of the centroids coordinates for each metaspot
#'
#' @export
#'
#'
#'


supercell_spatial_centroids <- function(MC,spotPositions){
  membership <-MC$membership
  seuratCoordMetacell <-  cbind(spotPositions,membership)

  centroids <- stats::aggregate(spotPositions %>% as.matrix() ~membership,spotPositions,mean) #should be taken from object slot
  centroids[["supercell_size"]] <- MC[["supercell_size"]]
  return(centroids)
}

#' Update of the memberships of unconnected components
#'
#' This function checks if the spots assigned as a metaspots are connected in the KNN graph and split them by updating their membership
#'
#'
#' @param m given membership
#' @param MC Metaspot object obtained from SCimplify_SpatialDLS function
#'
#' @return a vector of membership with updated split memberships
#'
#' @export
#'
#'
#'

split_membership <- function(m,MC){
  memberships <- MC$membership
  vertex.id <- which(memberships == m)
  sub_network <- subgraph(MC$graph.singlecell, vertex.id)
  if (!igraph::is.connected(sub_network)) {
    #plot(sub_network)
    #print(m)
    new_memberships <- clusters(sub_network)$membership
    old_memberships <- rep(m, length(vertex.id))
    new_memberships <- paste(old_memberships, new_memberships, sep = "_") #%>% as.numeric()
    memberships[vertex.id] <- new_memberships
    #print(memberships[vertex.id])

  }
  return(memberships[vertex.id])
}

#' Split metaspots that are unconnect in KNN graph
#'
#' This function splits metaspots that are unconnect in KNN graph
#'
#'
#' @param MC Metaspot object obtained from SCimplify_SpatialDLS function
#'
#' @return a vector of membership with updated split memberships
#'
#' @export
#'
#'
#'

split_unconnected <- function(MC){
  memberships <- MC$membership
  split.memb <- pbapply::pbsapply(unique(memberships),function(m) split_membership(m,MC = MC))
  split.memb <- unlist(split.memb)
  split.memb <- split.memb[names(MC$membership)]
  split.memb <- as.character(split.memb)
  rv <- 1:length(unique(split.memb))
  names(rv) <- unique(split.memb)
  final_memb <- plyr::revalue(split.memb,rv)
  MC$membership <- final_memb %>% as.integer()
  names(MC$membership) <- names(memberships)
  #MC$supercell_size <- table(MC$membership) %>% unname()
  MC$supercell_size <- as.vector(table(MC$membership))
  MC$N.SC <- length(MC$supercell_size)
  return(MC)
}

#' Simplification of spatial transcriptomic dataset
#'
#' This function converts (i.e., averages or sums up) gene-expression matrix of spatial transcriptomic data into a gene expression
#' matrix of metacells
#'
#' @param MC Metaspot object obtained from SCimplify_SpatialDLS function
#' @param ge gene expression matrix (or any coordinate matrix) with genes as rows and cells as cols
#' @param groups vector of membership (assignment of single-cell to metacells)
#' @param mode string indicating whether to average or sum up `ge` within metacells
#' @param weights vector of a cell weight (NULL by default), used for computing average gene expression withing cluster of metaells
#' @param do.median.norm whether to normalize by median value (FALSE by default)
#'
#' @return a matrix of simplified (averaged withing groups) data with ncol equal to number of groups and nrows as in the initial dataset
#' @export

superspot_GE <- function (MC, ge, groups, mode = c("average", "sum"), weights = NULL,
                         do.median.norm = FALSE)
{
  if (ncol(ge) != length(groups)) {
    stop("Length of the vector groups has to be equal to the number of cols in matrix ge")
  }
  mode <- mode[1]
  if (!(mode %in% c("average", "sum"))) {
    stop(paste("mode", mode, "is unknown. Available values are 'average' and 'sum'."))
  }
  N.SC <- MC$N.SC
  supercell_size <- as.vector(MC$supercell_size)
  j <- rep(1:N.SC, supercell_size)
  goups.idx <- plyr::split_indices(groups %>% as.integer())
  i <- unlist(goups.idx)
  if (is.null(weights)) {
    M.AV <- Matrix::sparseMatrix(i = i, j = j)
    GE.SC <- ge %*% M.AV
    if (mode == "average") {
      GE.SC <- sweep(GE.SC, 2, supercell_size, "/")
    }
  }
  else {
    if (length(weights) != length(groups)) {
      stop("weights must be the same length as groups or NULL in case of unweighted averaging")
    }
    if (mode != "average") {
      stop(paste("weighted averaging is supposted only for mode = 'average', not for",
                 mode))
    }
    M.AV <- Matrix::sparseMatrix(i = i, j = j, x = weights[i])
    GE.SC <- ge %*% M.AV
    weighted_supercell_size <- unlist(lapply(goups.idx,
                                             FUN = function(x) {
                                               sum(weights[x])
                                             }))
    GE.SC <- t(t(GE.SC) / weighted_supercell_size)
  }
  if (do.median.norm) {
    GE.SC <- (GE.SC + 0.01)/apply(GE.SC + 0.01, 1, stats::median)
  }
  return(GE.SC)
}

#' Computation of the geometric values of metaspots
#'
#' This function computes the geometric values of metaspots based on the previously computed polygons using pliman R package
#'
#' @param MC Metaspot object obtained from SCimplify_SpatialDLS function
#' @param polygon_col character of the name where the polygons representing the metaspots are stored in the MC object
#' @param membership_col character of the name where the memberships of the spots are stored in the MC object
#'
#' @return a dataframe with columns
#' \itemize{
#'   \item memberships - memberships of the metaspots
#'   \item elongations - elongations of the metaspots
#'   \item circularities - circularities of the metaspots
#'   \item convexities - convexities of the metaspots
#' }
#'
#' @export

shape_metrics <- function(MC,polygon_col,membership_col){
  full.df <- MC[[polygon_col]]
  names(table(MC[[membership_col]]))[table(MC[[membership_col]]) > 2] -> big.memb
  #print(dupli.memb)
  #full.df[full.df$membership %in% big.memb,] -> full.df_final
  memberships <- c()
  elongations <- c()
  circularities <- c()
  convexities <- c()
  #for (m in unique(full.df_final$membership)){
  for (m in unique(full.df$membership)){
    #print(length(unique(full.df_final$membership)))
    #df.tmp <- subset(full.df_final, membership == m)[,1:2] %>% as.matrix()
    df.tmp <- subset(full.df, membership == m)[,1:2] %>% as.matrix()
    #if (nrow(df.tmp)> 2)  {
    if (m %in% big.memb)  {
      #plot_polygon(df.tmp,points = T)
      aligned.tmp <- pliman::poly_align(df.tmp,plot = F)
      #plot_polygon(aligned.tmp)
      memberships <- c(memberships,m)
      elongations <- c(elongations,pliman::poly_elongation(aligned.tmp) %>% abs())
      circularities <- c(circularities,pliman::poly_circularity_norm(aligned.tmp))
      convexities <- c(convexities,pliman::poly_convexity(aligned.tmp))
    }
    else {
      memberships <- c(memberships,m)
      elongations <- c(elongations,0)
      circularities <- c(circularities,0)
      convexities <- c(convexities,0)
    }

  }
  metrics <- tibble(memberships = memberships, elongations = elongations, circularities = circularities, convexities = convexities)
  return(metrics)
}

#' Superspots to Seurat object
#'
#' This function transforms superspot gene expression and superspot partition into \link[Seurat]{Seurat} object
#'
#'
#' @param SC.GE gene expression matrix with genes as rows and cells as columns
#' @param SC super-cell (output of \link{SCimplify} function)
#' @param fields which fields of \code{SC} to use as cell metadata
#'
#'
#'
#' @return \link[Seurat]{Seurat} object
#'
#' @export

supercell_2_Seuratv5 <- function (SC.GE, SC, fields = c())
{
  N.c <- ncol(SC.GE)
  if (is.null(SC$supercell_size)) {
    warning(paste0("supercell_size field of SC is missing, size of all super-cells set to 1"))
    supercell_size <- rep(1, N.c)
  }
  else {
    supercell_size <- SC$supercell_size
  }
  if (length(supercell_size) != N.c) {
    stop(paste0("length of SC$supercell_size has to be the same as number of super-cells ",
                N.c))
  }
  if (is.null(colnames(SC.GE))) {
    colnames(SC.GE) <- as.character(1:N.c)
  }
  counts <- SC.GE
  if (is.numeric(fields)) {
    fields <- names(SC)[fields]
  }
  fields <- intersect(fields, names(SC))
  if (length(fields) > 0) {
    SC.fields <- SC[fields]
  }
  else {
    SC.fields <- NULL
  }
  SC.field.length <- lapply(SC.fields, length)
  SC.fields <- SC.fields[which(SC.field.length == N.c)]
  meta <- data.frame(size = supercell_size, row.names = colnames(SC.GE),
                     stringsAsFactors = FALSE)
  if (length(SC.fields) > 0) {
    meta <- cbind(meta, SC.fields)
  }
  m.seurat <- Seurat::CreateSeuratObject(counts = SC.GE, meta.data = meta)
  return(m.seurat)
}

