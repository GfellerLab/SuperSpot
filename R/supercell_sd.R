supercell_sd <- function(spotPositions, MC){

  membership <- MC$membership

  seuratCoord <- spotPositions
  seuratCoordMetacell <-  cbind(seuratCoord,membership)

  std_var <- stats::aggregate(seuratCoord %>% as.matrix() ~membership,seuratCoord,sd)

  simpl_sd <- c()
  for (i in 1:dim(std_var)[1]){
    sd_comb <- sqrt( (std_var[i,2])**2 + (std_var[i,3])**2)
    simpl_sd <- c(simpl_sd,sd_comb)
  }
  
  return(simpl_sd)
}

supercell_maxdist <- function(spotPositions, MC){
  membership <- MC$membership
  
  seuratCoord <- spotPositions
  seuratCoordMetacell <-  cbind(seuratCoord,membership)
  
  centroid <- stats::aggregate(seuratCoord %>% as.matrix() ~membership,seuratCoord,mean)
  
  max_dists <- c()
  for (i in 1:dim(centroid)[1]){
    merge_centr_spots <- rbind(centroid[i,2:3], subset(seuratCoordMetacell, membership == centroid[i,1]) %>% select(c("imagerow","imagecol")))
    whole_dists <- dist(merge_centr_spots)
    dist_from_centroid <- usedist::dist_get(whole_dists, paste0(i),rownames(merge_centr_spots))
    m <- max(dist_from_centroid)
    max_dists <- c(max_dists,m)
  }
  return(max_dists)
}


supercell_sd2 <- function(spotPositions, MC){
  membership <- MC$membership
  
  seuratCoord <- spotPositions
  seuratCoordMetacell <-  cbind(seuratCoord,membership)
  
  centroid <- stats::aggregate(seuratCoord %>% as.matrix() ~membership,seuratCoord,mean)
  
  sds <- c()
  for (i in 1:dim(centroid)[1]){
    merge_centr_spots <- rbind(centroid[i,2:3], subset(seuratCoordMetacell, membership == centroid[i,1]) %>% select(c("imagerow","imagecol")))
    whole_dists <- dist(merge_centr_spots)
    dist_from_centroid <- usedist::dist_get(whole_dists, paste0(i),rownames(merge_centr_spots))
    s <- sd(dist_from_centroid)
    sds <- c(sds,s)
  }
  return(sds)
}
