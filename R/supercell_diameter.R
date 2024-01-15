supercell_diameter <- function(MC){
  clusters_from_graph <- igraph::groups(MC$h_membership)
  len_sh_paths <- c()
  diam <- c()
  utils::txtProgressBar()
  for (i in 1:length(clusters_from_graph)){
    #print(paste0(round((i/length(clusters_from_graph))*100,digits = 0),"%"))
    
    sub_network <- subgraph(MC$graph.singlecell,clusters_from_graph[[i]])
    
    #sh_paths <- shortest.paths(sub_network)
    #hist(sh_paths)
    
    #shpaths_tmp <- sh_paths
    #sh_paths_tmp[sh_paths_tmp == 0 ] <- max(sh_paths_tmp)
    #len_sh_paths <- c(len_sh_paths, min(sh_paths_tmp))
    #print(len_sh_paths)
    
    diam <- c(diam,igraph::diameter(sub_network,weights = rep(1,length(E(sub_network)))))
    #print(diam[length(diam)])
    
  }
  return(diam)
}


supercell_int_clust_simil <- function(MC){
  clusters_from_graph <- igraph::groups(MC$h_membership)
  int_clust_simils <- c()
  for (i in 1:length(clusters_from_graph)){
    int_clust_simils <- c(int_clust_simils,igraph::E(subgraph(MC$graph.singlecell,clusters_from_graph[[i]]))$weight)
  }
  return(int_clust_simils)
}

supercell_clust_simils <- function(MC,md){
  so <- MC$seurat.object
  clusters_from_graph <- igraph::groups(MC$h_membership)
  bad_purities <- which(MC$purity_annot <= 0.8)
  intra_simi <- c()
  membership_intra <- c()
  annot_intra <- c()
  extra_simi <- c()
  membership_extra <- c()
  annot_extra <- c()
  len_mix <- c()
  for (i in bad_purities){
    #print(i)
    md_tmp <- md[clusters_from_graph[[i]],"Classification"]
    len_mix <- c(len_mix, length(unique(md_tmp)))
    for (j in unique(md_tmp)){
      #print(unique(md_tmp))
      samples <- rownames(subset(md[clusters_from_graph[[i]],], Classification == j))
      from = c()
      to = c()
      samples2 <- rownames(subset(md[clusters_from_graph[[i]],], Classification != j))
      from2 = c()
      to2 = c()
      
      for (s in samples){
        from <- c(from,rep(s,length(samples)-1))
        to <- c(to,samples[samples!=s])
        from2 <- c(from2,rep(s,length(samples2)))
        to2 <- c(to2,samples2)
      }
      
      g <- graph_from_data_frame(data.frame(from = from, to = to), directed=FALSE,)
                                 #vertices=data.frame(name = c(samples,samples)))
      g2 <- graph_from_data_frame(data.frame(from = from2, to = to2), directed=FALSE,)
      #plot(g)
      #plot(g2)
      pca_dist <- pbapply::pbapply(igraph::get.edgelist(g),
                                   1 ,
                                   function(x){parallelDist::parDist(so@reductions[["pca"]]@cell.embeddings[,1:MC$n.pc][x,])})
      pca_dist2 <- pbapply::pbapply(igraph::get.edgelist(g2),
                                   1 ,
                                   function(x){parallelDist::parDist(so@reductions[["pca"]]@cell.embeddings[,1:MC$n.pc][x,])})
      print(paste0("Done"))
      print(paste0("Computing similarity from PCA distances"))
      if (MC$method_similarity == "1"){
        pca_similarity <- pbapply::pblapply(pca_dist,FUN = function(i){1 - (i/max(pca_dist))}) %>% unlist()
        pca_similarity[pca_similarity == 0 ]<- 1e-100
        pca_similarity2 <- pbapply::pblapply(pca_dist2,FUN = function(i){1 - (i/max(pca_dist2))}) %>% unlist()
        pca_similarity2[pca_similarity2 == 0 ]<- 1e-100
        print(paste0("Done"))
      }
      else if (MC$method_similarity == "2"){
        pca_similarity <- pbapply::pblapply(pca_dist,FUN = function(i){exp(-(i**2)/MC$sigs)}) %>% unlist()
        pca_similarity[pca_similarity == 0 ]<- 1e-100
        pca_similarity2 <- pbapply::pblapply(pca_dist2,FUN = function(i){exp(-(i**2)/MC$sigs)}) %>% unlist()
        pca_similarity2[pca_similarity2 == 0 ]<- 1e-100
        print(paste0("Done"))
      }
      
      else if (MC$method_similarity == "3"){
        pca_similarity <- proxy::as.simil(pca_dist)
        pca_similarity[pca_similarity == 0 ]<- 1e-100
        pca_similarity2 <- proxy::as.simil(pca_dist2)
        pca_similarity2[pca_similarity2 == 0 ]<- 1e-100
        print(paste0("Done"))
      }
      intra_simi <- c(intra_simi,pca_similarity)
      membership_intra <- c(membership_intra,rep(i,length(pca_similarity)))
      annot_intra <- c(annot_intra,rep("intra",length(pca_similarity)))
      extra_simi <- c(extra_simi,pca_similarity2)
      membership_extra <- c(membership_extra,rep(i,length(pca_similarity2)))
      annot_extra <- c(annot_extra,rep("inter",length(pca_similarity2)))
    }
    
  }
  df <- tibble(intra_simi = c(intra_simi,extra_simi), membership_intra = c(membership_intra,membership_extra), annot_intra = c(annot_intra,annot_extra))
  return(df)
}

supercell_ext_clust_simil <- function(seuratCoord,MC,countmatrix = NULL){
  membership <- MC$membership
  
  #seuratCoord <- original_coord
  
  seuratCoordMetacell <-  cbind(seuratCoord,membership)
  
  centroids <- stats::aggregate(seuratCoord %>% as.matrix() ~membership,seuratCoord,mean)
  
  total_ext_clust_simil <- c()
  
  if (is.null(MC$seurat.object)){
  print(paste0("Performing SCT"))
  so <- CreateSeuratObject(countmatrix,assay = "RNA")
  so <- SCTransform(so, new.assay.name = "SCT" , assay = "RNA",verbose = F)
  print(paste0("Done"))
  print(paste0("Running PCA"))
  so <- RunPCA(so, assay = "SCT",verbose = F)
  print(paste0("Done"))}
  
  else {
    so <- MC$seurat.object
  }
  
  for (m in centroids$membership){
    
    print(paste0(m,"/",length(centroids$membership)))
    dist_centroid <- proxy::dist(x = select(centroids, c("imagerow","imagecol"))[m,],
                                 y = select(centroids, c("imagerow","imagecol")))
    dist_centroid[dist_centroid == 0] <- max(dist_centroid)
    closest_ms <- which.min(dist_centroid)
    
    samples_ms1 <- rownames(subset(seuratCoordMetacell, membership == m))
    samples_ms2 <- rownames(subset(seuratCoordMetacell, membership == closest_ms))
    
    from = c()
    to = c()
    
    for (s in samples_ms1){
      from <- c(from,rep(s,length(samples_ms2)))
      to <- c(to,samples_ms2)
    }
    
    g <- graph_from_data_frame(data.frame(from = from, to = to), directed=FALSE,
                               vertices=data.frame(name = c(samples_ms1,samples_ms2)))
    
    
    print(paste0("Computing PCA dist for each edge"))
    
    pca_dist <- pbapply::pbapply(igraph::get.edgelist(g),
                                 1 ,
                                 function(x){parallelDist::parDist(so@reductions[["pca"]]@cell.embeddings[,1:MC$n.pc][x,])})
    print(paste0("Done"))
    print(paste0("Computing similarity from PCA distances"))
    if (MC$method_similarity == "1"){
      pca_similarity <- pbapply::pblapply(pca_dist,FUN = function(i){1 - (i/max(pca_dist))}) %>% unlist()
      pca_similarity[pca_similarity == 0 ]<- 1e-100
      print(paste0("Done"))
    }
    else if (MC$method_similarity == "2"){
      pca_similarity <- pbapply::pblapply(pca_dist,FUN = function(i){exp(-(i**2)/MC$sigs)}) %>% unlist()
      pca_similarity[pca_similarity == 0 ]<- 1e-100
      print(paste0("Done"))
    }
    
    else if (MC$method_similarity == "3"){
      pca_similarity <- proxy::as.simil(pca_dist)
      pca_similarity[pca_similarity == 0 ]<- 1e-100
      print(paste0("Done"))
    }
    #hist(pca_similarity)
    total_ext_clust_simil <- c(total_ext_clust_simil, pca_similarity)
  }
  return(total_ext_clust_simil)
}
