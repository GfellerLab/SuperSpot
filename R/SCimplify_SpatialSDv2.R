knn_from_spotPositions <- function(spotPositions,k.knn){
  coord_assay <- Seurat::CreateAssayObject(counts = t(spotPositions))
  coord_neighbors <- Seurat::FindNeighbors(coord_assay,k.param = k.knn)
  coord_knn_am <- coord_neighbors$nn
  return(coord_knn_am)
}

#test <- knn_from_spotPositions(spotPositions,30)

knn_from_transcriptomics <- function(countMatrix,k.knn,n.pc){
  so <- CreateSeuratObject(countMatrix,assay = "RNA")
  so <- SCTransform(so, new.assay.name = "SCT" , assay = "RNA")
  so <- RunPCA(so, assay = "SCT")
  so <- FindNeighbors(so, dims = 1:n.pc, reduction = "pca",k.param = k.knn)
  
  transcr_knn_am <- so@graphs$SCT_nn
  return(transcr_knn_am)
}

#test2 <- knn_from_transcriptomics(countMatrix = mtx,n.pc = 30,k = 30,genes.use = hvg)


combining_sp_tr_knn <- function(spotPositions,countMatrix,k.knn_sp, k.knn_tr, n.pc, space_weight){
  am.sp <- knn_from_spotPositions2(spotPositions = spotPositions,k.knn = k.knn_sp)
  am.transc <- knn_from_transcriptomics2(countMatrix = countMatrix, k.knn = k.knn_tr, n.pc = n.pc)
  #am.sp2 <- am.sp*((max(am.transc)*space_weight)/(max(am.sp)))
  am.sp.transc <- (am.sp*space_weight)+am.transc
  knn.sp.transc <- graph_from_adjacency_matrix(am.sp.transc %>% as.matrix(),mode = "undirected",weighted = TRUE)
  return(knn.sp.transc)
}

SCimplify_SpatialSDv2 <- function(X,
                                spotPositions,
                                space_weight,
                                genes.use = NULL,
                                genes.exclude = NULL,
                                cell.annotation = NULL,
                                cell.split.condition = NULL,
                                n.var.genes = min(1000, nrow(X)),
                                gamma = 10,
                                k.knn_sp = 30,
                                k.knn_tr = 30,
                                do.scale = TRUE,
                                n.pc = 10,
                                fast.pca = TRUE,
                                do.approx = FALSE,
                                approx.N = 20000,
                                block.size = 10000,
                                seed = 12345,
                                igraph.clustering = c("walktrap", "louvain"),
                                return.singlecell.NW = TRUE,
                                return.hierarchical.structure = TRUE,
                                ...){
  
  N.c <- ncol(X)
  
  
  presampled.cell.ids <- colnames(X)
  
  

  sc.nw <- combining_sp_tr_knn(spotPositions = spotPositions,countMatrix = X,k.knn_sp = k.knn_sp, k.knn_tr = k.knn_tr,
                               n.pc = n.pc, space_weight = space_weight)
    
  #simplify
  
  k   <- round(N.c/gamma)
  
  if(igraph.clustering[1] == "walktrap"){
    g.s              <- igraph::cluster_walktrap(sc.nw)
    g.s$membership   <- igraph::cut_at(g.s, k)
    
  } else if(igraph.clustering[1] == "louvain") {
    warning(paste("igraph.clustering =", igraph.clustering, ", gamma is ignored"))
    g.s    <- igraph::cluster_louvain(sc.nw)
    
  } else {
    stop(paste("Unknown clustering method (", igraph.clustering, "), please use louvain or walkrtap"))
  }
  
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
              k.knn_sp = k.knn_sp,
              k.knn_tr = k.knn_tr,
              sc.cell.annotation. = cell.annotation,
              sc.cell.split.condition. = cell.split.condition
  )
  
  if(return.singlecell.NW){res$graph.singlecell <- sc.nw}
  if(!is.null(cell.annotation) | !is.null(cell.split.condition)){
    res$SC.cell.annotation. <- supercell_assign(cell.annotation, res$membership)
    res$SC.cell.split.condition. <- supercell_assign(cell.split.condition, res$membership)
  }
  
  if(igraph.clustering[1] == "walktrap" & return.hierarchical.structure)  res$h_membership <- g.s
  
  return(res)
}

###############################################################################
## V2
##############################################################################





knn_from_spotPositions2 <- function(spotPositions,k.knn){
  coord_assay <- Seurat::CreateAssayObject(counts = t(spotPositions))
  coord_neighbors <- Seurat::FindNeighbors(coord_assay,k.param = k.knn)
  coord_knn_am <- coord_neighbors$nn
  coord_knn <- graph.adjacency(as.matrix(coord_knn_am), mode = "undirected", weighted = T)
  return(coord_knn)
}

#test <- knn_from_spotPositions2(spotPositions,30)

knn_from_transcriptomics2 <- function(countMatrix,k.knn,n.pc){
  so <- CreateSeuratObject(countMatrix,assay = "RNA")
  so <- SCTransform(so, new.assay.name = "SCT" , assay = "RNA")
  so <- RunPCA(so, assay = "SCT")
  so <- FindNeighbors(so, dims = 1:n.pc, reduction = "pca",k.param = k.knn)
  transcr_knn_am <- so@graphs$SCT_nn
  transcr_knn <- graph.adjacency(as.matrix(transcr_knn_am), mode = "undirected", weighted = T)
  pca_dist <- apply(igraph::get.edgelist(transcr_knn),
                    1 ,
                    function(i){stats::dist(so@reductions[["pca"]]@cell.embeddings[,1:n.pc][i,])})
  pca_similarity <- lapply(pca_dist,FUN = function(i){1 - (i/max(pca_dist))}) %>% unlist()
  igraph::E(transcr_knn)$weight <- pca_similarity
  return(transcr_knn)
}

#test2 <- knn_from_transcriptomics(countMatrix = mtx,n.pc = 30,k = 30,genes.use = hvg)


combining_sp_tr_knn2 <- function(spotPositions,countMatrix,k.knn_sp, k.knn_tr, n.pc, space_weight){
  knn.sp <- knn_from_spotPositions2(spotPositions = spotPositions,k.knn = k.knn_sp)
  knn.transc <- knn_from_transcriptomics2(countMatrix = countMatrix, k.knn = k.knn_tr, n.pc = n.pc)
  knn.sp.transc <- igraph::intersection(knn.transc,knn.sp)
  return(knn.sp.transc)
}

SCimplify_SpatialSDv3 <- function(X,
                                  spotPositions,
                                  space_weight,
                                  genes.use = NULL,
                                  genes.exclude = NULL,
                                  cell.annotation = NULL,
                                  cell.split.condition = NULL,
                                  n.var.genes = min(1000, nrow(X)),
                                  gamma = 20,
                                  k.knn_sp = 100,
                                  k.knn_tr = 30,
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
  
  N.c <- ncol(X)
 
  presampled.cell.ids <- colnames(X) 
 
  sc.nw <- combining_sp_tr_knn2(spotPositions = spotPositions,countMatrix = X,k.knn_sp = k.knn_sp, k.knn_tr = k.knn_tr,
                               n.pc = n.pc,space_weight = space_weight)
  
  #simplify
  
  k   <- round(N.c/gamma)
  
  if(igraph.clustering[1] == "walktrap"){
    g.s              <- igraph::cluster_walktrap(sc.nw,weights = igraph::E(sc.nw)$weight_1)
    g.s$membership   <- igraph::cut_at(g.s, k)
    
  } else if(igraph.clustering[1] == "louvain") {
    warning(paste("igraph.clustering =", igraph.clustering, ", gamma is ignored"))
    g.s    <- igraph::cluster_louvain(sc.nw)
    
  } else {
    stop(paste("Unknown clustering method (", igraph.clustering, "), please use louvain or walkrtap"))
  }
  
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
              k.knn_sp = k.knn_sp,
              k.knn_tr = k.knn_tr,
              sc.cell.annotation. = cell.annotation,
              sc.cell.split.condition. = cell.split.condition
  )
  
  if(return.singlecell.NW){res$graph.singlecell <- sc.nw}
  if(!is.null(cell.annotation) | !is.null(cell.split.condition)){
    res$SC.cell.annotation. <- supercell_assign(cell.annotation, res$membership)
    res$SC.cell.split.condition. <- supercell_assign(cell.split.condition, res$membership)
  }
  
  if(igraph.clustering[1] == "walktrap" & return.hierarchical.structure)  res$h_membership <- g.s
  
  return(res)
}

