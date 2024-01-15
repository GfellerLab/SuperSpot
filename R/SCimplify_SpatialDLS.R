neighbor_graph <- function(spotPositions, k.spots, countMatrix, n.pc, method_similarity, method_knn, k.knn, sigs,
                           method_normalization,split_not_connected,split_annotation = NULL,split_vector){
  if (method_knn == "1"){
    #coord_assay <- Seurat::CreateAssayObject(counts = t(spotPositions))
    #coord_neighbors <- Seurat::FindNeighbors(coord_assay,k.param = k.knn)
    #coord_knn_am <- coord_neighbors$nn
    #spot.graph <- graph.adjacency(coord_knn_am,mode = "undirected")
    print(paste0("Building KNN graph with nn2"))
    nn2.res <- RANN::nn2(data = spotPositions,k = k.knn+1)
    if (split_not_connected == TRUE){
      
      #min_dist <- min(nn2.res$nn.dists[,2:k.knn])
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
    
    #grah.matrix <- as_adj(spot.graph,sparse = F)
    #grah.matrix[grah.matrix < 4 ] <- 0
    #spot.graph <- graph_from_adjacency_matrix(grah.matrix,mode = "undirected")
    max_gamma <- (nrow(spotPositions)/igraph::components(spot.graph)$no)
    print(paste0("Maximum gamma is ",max_gamma))
    print(paste0("Done"))
    V(spot.graph)$x <- spotPositions$imagecol
    V(spot.graph)$y <- -spotPositions$imagerow
    #plot(spot.graph,
    #     layout = matrix(c(V(spot.graph)$x, V(spot.graph)$y), ncol = 2),
    #     vertex.label.cex = 0.1,
    #     vertex.size = 1,
    #     edge.arrow.size = 0.7,
    #     main = "Graph with Custom X and Y Coordinates")
  }
  else if (method_knn == "2"){
    print(paste0("Building KNN graph by DLS"))
    print(paste0("Computing distances between spots"))
    distCoord <- parallelDist::parDist(spotPositions %>% as.matrix())
    min.dist <- min(distCoord)
    #distCoord2 <- distCoord
    distCoord[distCoord > k.spots*round(min.dist)] <- 0
    distCoord[distCoord >0] <- 1
    #print(rownames(distCoord))
    spot.graph <- graph.adjacency(distCoord %>% as.matrix(),mode = "undirected")
    #print(E(spot.graph))
    #print(V(spot.graph))
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
    #if (max_gamma != 1){
    so <- RunPCA(so, assay = "SCT",verbose = F)
    pca_matrix <- so@reductions[["pca"]]@cell.embeddings[,1:n.pc]
    #}
    #else {
    #  pca_matrix <- matrix(0,ncol(countMatrix),n.pc)
    #}
  }
  
  print(paste0("Done"))
  print(paste0("Computing PCA dist for each edge"))
  #pb <- 
  pca_dist <- pbapply::pbapply(igraph::get.edgelist(spot.graph),
                    1 ,
                    function(x){parallelDist::parDist(pca_matrix[x,])})
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
  #print(length(igraph::E(spot.graph)$weight))
  #print(length(pca_similarity))
  igraph::E(spot.graph)$weight <- pca_similarity
  #print(igraph::E(spot.graph)$weight)
  #hist(igraph::E(spot.graph)$weight)
  #hist(pca_similarity)
  return_list <- list("graph" = spot.graph, "seurat.object" = so)
  return(return_list)
}

SCimplify_SpatialDLS <- function(X,
                                spotPositions,
                                k.spots = 1,
                                method_similarity,
                                method_knn,
                                method_normalization,
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
                                n.pc = 30,
                                fast.pca = TRUE,
                                do.approx = FALSE,
                                approx.N = 20000,
                                block.size = 10000,
                                seed = 12345,
                                igraph.clustering = c("walktrap", "louvain"),
                                return.singlecell.NW = TRUE,
                                return.hierarchical.structure = TRUE,
                                ...){
  
  neighbor_graph.output <- neighbor_graph(spotPositions = spotPositions, k.spots = k.spots,split_not_connected = split_not_connected,
                                          countMatrix = X, n.pc = n.pc, method_similarity = method_similarity,
                                          method_knn = method_knn, k.knn = k.knn, sigs = sigs,method_normalization = method_normalization,
                                          split_annotation = split_annotation,split_vector =split_vector)
  
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
  
  res <- list(seurat.object = neighbor_graph.output$seurat.object,
              graph.supercells = SC.NW,
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
  return(res)
}


supercell_spatial_centroids <- function(MC,spotPositions){
  membership <-MC$membership
  seuratCoordMetacell <-  cbind(spotPositions,membership)
  
  centroids <- stats::aggregate(spotPositions %>% as.matrix() ~membership,spotPositions,mean) #should be taken from object slot
  centroids[["supercell_size"]] <- MC[["supercell_size"]]
  return(centroids)
}


split_unconnected <- function(MC){
  memberships <- MC$membership
  memberships_tmp <- memberships #%>% as.character()
  for (m in unique(memberships)[11:20]){
    #print(m)
    vertex.id <- which(memberships %in% m)
    sub_network <- subgraph(MC$graph.singlecell, vertex.id)
    
    if (is.connected(sub_network) == FALSE){
      #plot(sub_network)
      new_memberships <- components(sub_network,mode = "strong")$membership
      #dg <- decompose(sub_network)
      old_memberships <- rep(m,length(vertex.id))
      updated_memberships <- paste(old_memberships,new_memberships,sep = ".") %>% as.numeric()
      updated_memberships.chr <- paste(old_memberships,new_memberships,sep = ".")
      updated_memberships.chr2 <- as.character(updated_memberships)
      memberships_tmp[vertex.id] <- updated_memberships
      plot(ggplot(spotPosition.sb)+geom_point(aes(x = imagecol,y = imagerow),alpha = 0.1)+geom_point(data = spotPosition.sb[vertex.id,],mapping = aes(x = imagecol,y = imagerow),color='red')+NoLegend())
      plot(ggplot(spotPosition.sb)+geom_point(aes(x = imagecol,y = imagerow),alpha = 0.1)+geom_point(data = cbind(spotPosition.sb[vertex.id,],updated_memberships.chr),mapping = aes(x = imagecol,y = imagerow,color = updated_memberships.chr)) + 
  ggrepel::geom_label_repel(data = cbind(spotPosition.sb[vertex.id,],updated_memberships.chr),aes(x = imagecol,y = imagerow,label = updated_memberships.chr),
                              box.padding   = 0.35, 
                              point.padding = 0.5,
                              segment.color = 'grey50') +NoLegend())
      plot(ggplot(spotPosition.sb)+geom_point(aes(x = imagecol,y = imagerow),alpha = 0.1)+geom_point(data = cbind(spotPosition.sb[vertex.id,],updated_memberships.chr2),mapping = aes(x = imagecol,y = imagerow,color = updated_memberships.chr2)) + 
             ggrepel::geom_label_repel(data = cbind(spotPosition.sb[vertex.id,],updated_memberships.chr2),aes(x = imagecol,y = imagerow,label = updated_memberships.chr2),
                                       box.padding   = 0.35, 
                                       point.padding = 0.5,
                                       segment.color = 'grey50') +NoLegend())
    }
    
  }
  #MC$membership <- as.integer(memberships_tmp*1000)
  MC$membership <- memberships_tmp
  #names(MC$membership) <- names(memberships)
  MC$supercell_size <- as.vector(table(MC$membership))
  MC$N.SC <- length(MC$supercell_size)
  return(MC)
}

split_unconnected3 <- function(MC){
  memberships <- MC$membership
  memberships_tmp <- memberships #%>% as.character()
  for (m in unique(memberships)){
    #print(m)
    vertex.id <- which(memberships %in% m)
    sub_network <- subgraph(MC$graph.singlecell, vertex.id)
    
    if (is.connected(sub_network) == FALSE){
      #plot(sub_network)
      new_memberships <- igraph::components(sub_network,mode = "strong")$membership
      max_membership <- max(memberships_tmp)
      #print(max_membership)
      old_memberships <- rep(m,length(vertex.id))
      updated_memberships <- paste(old_memberships,new_memberships,sep = "") %>% as.numeric()
      updated_memberships_shifted <- updated_memberships+max_membership
      updated_memberships.chr2 <- as.character(updated_memberships_shifted)
      memberships_tmp[vertex.id] <- updated_memberships_shifted
      #plot(ggplot(spotPosition.sb)+geom_point(aes(x = imagecol,y = imagerow),alpha = 0.1)+geom_point(data = spotPosition.sb[vertex.id,],mapping = aes(x = imagecol,y = imagerow),color='red')+NoLegend())
      #plot(ggplot(spotPosition.sb)+geom_point(aes(x = imagecol,y = imagerow),alpha = 0.1)+geom_point(data = cbind(spotPosition.sb[vertex.id,],updated_memberships.chr),mapping = aes(x = imagecol,y = imagerow,color = updated_memberships.chr)) + 
             #ggrepel::geom_label_repel(data = cbind(spotPosition.sb[vertex.id,],updated_memberships.chr),aes(x = imagecol,y = imagerow,label = updated_memberships.chr),
              #                         box.padding   = 0.35, 
              #                         point.padding = 0.5,
              #                         segment.color = 'grey50') +NoLegend())
      #plot(ggplot(spotPosition.sb)+geom_point(aes(x = imagecol,y = imagerow),alpha = 0.1)+geom_point(data = cbind(spotPosition.sb[vertex.id,],updated_memberships.chr2),mapping = aes(x = imagecol,y = imagerow,color = updated_memberships.chr2)) + 
      #       ggrepel::geom_label_repel(data = cbind(spotPosition.sb[vertex.id,],updated_memberships.chr2),aes(x = imagecol,y = imagerow,label = updated_memberships.chr2),
                                       #box.padding   = 0.35, 
                                       #point.padding = 0.5,
      #                                 segment.color = 'grey50') +NoLegend())
    }
    
  }
  #MC$membership <- as.integer(memberships_tmp*1000)
  MC$membership <- memberships_tmp
  names(MC$membership) <- names(memberships)
  MC$supercell_size <- as.vector(table(MC$membership))
  MC$N.SC <- length(MC$supercell_size)
  return(MC)
}


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

split_unconnected2 <- function(MC){
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

supercell_GEv2 <- function (MC, ge, groups, mode = c("average", "sum"), weights = NULL, 
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
    #GE.SC <- sweep(GE.SC, 2, weighted_supercell_size, "/")
    GE.SC <- t(t(GE.SC) / weighted_supercell_size)
  }
  if (do.median.norm) {
    GE.SC <- (GE.SC + 0.01)/apply(GE.SC + 0.01, 1, stats::median)
  }
  return(GE.SC)
}

supercell_metaspots_shape <- function(MC, spotpositions,annotation,concavity,membership_name){
  hull_df_final <- data.frame(x = c(NULL), y = c(NULL), cell_type = c(NULL))
  for (memb in unique(MC[[membership_name]])){
    index.tmp <- which(MC[[membership_name]] == memb)
    polygons.tmp <- concaveman::concaveman(as.matrix(spotpositions[index.tmp,]),concavity = concavity)
    index.tmp2 <- which(names(MC[[annotation]]) == memb)
    hull_df.tmp <- data.frame(x = polygons.tmp[, 2], y = polygons.tmp[, 1],
                              cell_type = rep(MC[[annotation]][index.tmp2],nrow(polygons.tmp)),
                              membership = rep(memb,nrow(polygons.tmp)))
    hull_df_final <- rbind(hull_df_final,hull_df.tmp)
  }
  return(hull_df_final)
}

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



supercell_2_Seuratv5 <- function (SC.GE, SC, fields = c(), var.genes = NULL, is.log.normalized = TRUE, 
          do.center = TRUE, do.scale = TRUE, N.comp = NULL) 
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
  #m.seurat <- Seurat::NormalizeData(m.seurat)
  #print("Done: NormalizeData")
  #if (is.log.normalized) {
  #  print("Doing: data to normalized data")
  #  m.seurat@assays$RNA@layers$data <- m.seurat@assays$RNA@layers$counts
  #}
  #if (length(unique(meta$size)) > 1) {
  #  print("Doing: weighted scaling")
  #  m.seurat@assays$RNA@layers$scale.data <- t(as.matrix(corpcor::wt.scale(Matrix::t((m.seurat@assays$RNA@layers$data)), 
  #                                                                  w = meta$size, center = do.center, scale = do.scale)))
  #  print("Done: weighted scaling")
  #}
  #else {
  #  print("Doing: unweighted scaling")
  #  m.seurat <- Seurat::ScaleData(m.seurat)
  #  print("Done: unweighted scaling")
  #}
  #m.seurat@assays$RNA@misc[["scale.data.weighted"]] <- m.seurat@assays$RNA@layers$scale.data
  #if (is.null(var.genes)) {
  #  var.genes <- sort(SC$genes.use)
  #}
  #if (is.null(N.comp)) {
  #  N.comp <- min(50, ncol(m.seurat@assays$RNA@layers$counts) - 
  #                  1)
  #}
  #Seurat::VariableFeatures(m.seurat) <- var.genes
  #m.seurat <- Seurat::RunPCA(m.seurat, verbose = F, npcs = max(N.comp),features = var.genes,)
  #m.seurat@reductions$pca_seurat <- m.seurat@reductions$pca
  #my_pca <- supercell_prcomp(X = Matrix::t(SC.GE[var.genes, 
  #]), genes.use = var.genes, fast.pca = TRUE, supercell_size = meta$supercell_size, 
  #k = dim(m.seurat@reductions$pca_seurat)[2], do.scale = do.scale, 
  #do.center = do.center)
  #dimnames(my_pca$x) <- dimnames(m.seurat@reductions$pca_seurat)
  #m.seurat@reductions$pca@cell.embeddings <- my_pca$x
  #m.seurat@reductions$pca@feature.loadings <- my_pca$rotation
  #m.seurat@reductions$pca@stdev <- my_pca$sdev
  #m.seurat@reductions$pca_weighted <- m.seurat@reductions$pca
  #m.seurat <- Seurat::FindNeighbors(m.seurat, compute.SNN = TRUE, 
  #                                  verbose = TRUE)
  #if (!is.null(SC$graph.supercells)) {
  #  adj.mtx <- igraph::get.adjacency(SC$graph.supercells, 
  #                                   attr = "weight")
  #  m.seurat@graphs$RNA_nn@i <- adj.mtx@i
  #  m.seurat@graphs$RNA_nn@p <- adj.mtx@p
  #  m.seurat@graphs$RNA_nn@Dim <- adj.mtx@Dim
  #  m.seurat@graphs$RNA_nn@x <- adj.mtx@x
  #  m.seurat@graphs$RNA_nn@factors <- adj.mtx@factors
  #  m.seurat@graphs$RNA_super_cells <- m.seurat@graphs$RNA_nn
  #}
  #else {
  #  warning("Super-cell graph was not found in SC object, no super-cell graph was added to Seurat object")
  #}
  return(m.seurat)
}

