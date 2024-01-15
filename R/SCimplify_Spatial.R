SCimplify_SpatialSD <- function(X,
                      genes.use = NULL,
                      genes.exclude = NULL,
                      cell.annotation = NULL,
                      cell.split.condition = NULL,
                      n.var.genes = min(1000, nrow(X)),
                      gamma = 10,
                      k.knn = 5,
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
  
  if(gamma > 100 & N.c < 100000){
    warning(paste0("Graining level (gamma = ", gamma, ") seems to be very large! Please, consider using smaller gamma, the suggested range is 10-50."))
  }
  
  if(is.null(rownames(X))){
    if(!(is.null(genes.use) | is.null(genes.exclude))){
      stop("rownames(X) is Null \nGene expression matrix X is expected to have genes as rownames")
    } else {
      warning("colnames(X) is Null, \nGene expression matrix X is expected to have genes as rownames! \ngenes will be created automatically in a form 'gene_i' ")
      #rownames(X) <- paste("gene", 1:nrow(X), sep = "_")
    }
  }
  
  if(is.null(colnames(X))){
    warning("colnames(X) is Null, \nGene expression matrix X is expected to have cellIDs as colnames! \nCellIDs will be created automatically in a form 'cell_i' ")
    #colnames(X) <- paste("cell", 1:N.c, sep = "_")
  }
  
  #keep.genes    <- setdiff(rownames(X), genes.exclude)
  #X             <- X[keep.genes,]
  
  
  if(is.null(genes.use)){
    n.var.genes <- min(n.var.genes, nrow(X))
    if(N.c > 50000){
      set.seed(seed)
      idx         <- sample(N.c, 50000)
      gene.var    <- apply(X[,idx], 1, stats::var)
    } else {
      gene.var    <- apply(X, 1, stats::var)
    }
    
    genes.use   <- names(sort(gene.var, decreasing = TRUE))[1:n.var.genes]
  }
  
  if(length(intersect(genes.use, genes.exclude)) > 0){
    stop("Sets of genes.use and genes.exclude have non-empty intersection")
  }
  
  #genes.use <- genes.use[genes.use %in% rownames(X)]
  #X <- X[genes.use,]
  
  if(do.approx & approx.N >= N.c){
    do.approx <- FALSE
    warning("N.approx is larger or equal to the number of single cells, thus, an exact simplification will be performed")
  }
  
  if(do.approx & (approx.N < round(N.c/gamma))){
    approx.N <- round(N.c/gamma)
    warning(paste("N.approx is set to N.SC", approx.N))
  }
  
  if(do.approx & ((N.c/gamma) > (approx.N/3))){
    warning("N.approx is not much larger than desired number of super-cells, so an approximate simplification may take londer than an exact one!")
  }
  
  if(do.approx){
    set.seed(seed)
    approx.N            <- min(approx.N, N.c)
    presample           <- sample(1:N.c, size = approx.N, replace = FALSE)
    presampled.cell.ids <- colnames(X)[sort(presample)]
    rest.cell.ids       <- setdiff(colnames(X), presampled.cell.ids)
  } else {
    presampled.cell.ids <- colnames(X)
    rest.cell.ids       <- c()
  }
  
  #X.for.pca            <- Matrix::t(X[genes.use, presampled.cell.ids])
  #if(do.scale){ X.for.pca            <- scale(X.for.pca) }
  #X.for.pca[is.na(X.for.pca)] <- 0
  
  #if(is.null(n.pc[1]) | min(n.pc) < 1){stop("Please, provide a range or a number of components to use: n.pc")}
  #if(length(n.pc)==1) n.pc <- 1:n.pc
  
  #if(fast.pca & (N.c < 1000)){
  #  warning("Normal pca is computed because number of cell is low for irlba::irlba()")
  #  fast.pca <- FALSE
  #}
  
  #if(!fast.pca){
  #  PCA.presampled          <- stats::prcomp(X.for.pca, rank. = max(n.pc), scale. = F, center = F)
  #} else {
  #  set.seed(seed)
  #  PCA.presampled          <- irlba::irlba(X.for.pca, nv = max(n.pc, 25))
  #  PCA.presampled$x        <- PCA.presampled$u %*% diag(PCA.presampled$d)
  #  PCA.presampled$rotation <- PCA.presampled$v
  #}
  
  
  sc.nw <- build_knn_graph(
    X = X,
    #k = k.knn, 
    from = "dist",
    use.nn2 =F,
    #dist_method = "euclidean",
    #directed = directed,
    #DoSNN = DoSNN,
    #pruning = pruning,
    #which.snn = which.snn,
    #kmin = kmin,
    ...
  )
  
  
  #simplify
  
  k   <- round(N.c/gamma)
  
  if(igraph.clustering[1] == "walktrap"){
    g.s              <- igraph::cluster_walktrap(sc.nw$graph.knn)
    g.s$membership   <- igraph::cut_at(g.s, k)
    
  } else if(igraph.clustering[1] == "louvain") {
    warning(paste("igraph.clustering =", igraph.clustering, ", gamma is ignored"))
    g.s    <- igraph::cluster_louvain(sc.nw$graph.knn)
    
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
  
  
  
  SC.NW                        <- igraph::contract(sc.nw$graph.knn, membership.presampled)
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
              k.knn = k.knn,
              sc.cell.annotation. = cell.annotation,
              sc.cell.split.condition. = cell.split.condition
  )
  
  if(return.singlecell.NW){res$graph.singlecell <- sc.nw$graph.knn}
  if(!is.null(cell.annotation) | !is.null(cell.split.condition)){
    res$SC.cell.annotation. <- supercell_assign(cell.annotation, res$membership)
    res$SC.cell.split.condition. <- supercell_assign(cell.split.condition, res$membership)
  }
  
  if(igraph.clustering[1] == "walktrap" & return.hierarchical.structure)  res$h_membership <- g.s
  
  return(res)
}
