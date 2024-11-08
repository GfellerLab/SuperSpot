max.dist <- c()
gammas <- c()
for (g in 2:20) {
  print(g)
  n.pc = 1:5 # number of first PC to use
  k.knn = 15 # number of neighbors to connect to each spot

  print("Creating metaspots")
  # By default, SCimplify_SpatialDLS computes distances in a parallalized way. By default, all the available pus are used. If your computer doesn't suppont, you can change the
  MC.well7.tmp <- SCimplify_SpatialDLS_v2(X = well7.mtx ,
  spotPositions = spotPosition,
  col.x = "X",
  col.y = "Y",
  method_similarity = "1", split_not_connected = TRUE, method_reduction = "PCA", gamma = g,
n.dim = n.pc,
  method_knn = "1",
    k.knn = k.knn,
  method_normalization = "log_normalize",
  #cell. annotation = well7.md$Main_molecular_cell_type
  )
  l <- 0
for (m in unique(MC.well7.tmp$membership)) {
    cell.id <- names(MC.well7.tmp$membership[MC.well7.tmp$membership== m])
    #print(cell.id)
    #print(spotPosition[cell.id,])
    if (length(cell.id) > 1){
      max.dist.tmp <- max(dist(spotPosition[cell.id,]))
      max.dist <- c(max.dist,max.dist.tmp)
      l <- l+1
    }
}
  #plot(boxplot(max.dist))
  gammas <- c(gammas,rep(g,l))
}

dist_vs_g_df <- tibble(max_dists = max.dist, gammas = gammas)

ggplot() + geom_boxplot(aes(x = gammas, y = max.dist, group = gammas))


max.dist.t <- c()
thresh.dist <- c()
nb.sing <- c()
thresh.sing <- c()
for (t in (2:10)/10) {
  print(t)
  g = 10
  n.pc = 1:5 # number of first PC to use
  k.knn = 15 # number of neighbors to connect to each spot

  print("Creating metaspots")
  # By default, SCimplify_SpatialDLS computes distances in a parallalized way. By default, all the available pus are used. If your computer doesn't suppont, you can change the
  MC.well7.tmp <- SCimplify_SpatialDLS_v2(X = well7.mtx ,
                                          spotPositions = spotPosition,
                                          col.x = "X",
                                          col.y = "Y",
                                          method_similarity = "1", split_not_connected = TRUE, method_reduction = "PCA", gamma = g,
                                          n.dim = n.pc,
                                          method_knn = "1",
                                          k.knn = k.knn,
                                          method_normalization = "log_normalize",
                                          pct = t
                                          #cell. annotation = well7.md$Main_molecular_cell_type
  )
  l <- 0
  for (m in unique(MC.well7.tmp$membership)) {
    cell.id <- names(MC.well7.tmp$membership[MC.well7.tmp$membership== m])
    #print(cell.id)
    #print(spotPosition[cell.id,])
    if (length(cell.id) > 1){
      max.dist.tmp <- max(dist(spotPosition[cell.id,]))
      max.dist.t <- c(max.dist.t,max.dist.tmp)
      l <- l+1
    }
  }
  #plot(boxplot(max.dist))
  thresh.dist <- c(thresh.dist,rep(t,l))
  nb.sing <- c(nb.sing,length(MC.well7.tmp$supercell_size[MC.well7.tmp$supercell_size==1]))
  thresh.sing <- c(thresh.sing,t)
}

dist_vs_thresh_df <- tibble(max_dists = max.dist.t, pct = 1-thresh.dist)
sing_vs_thresh_df <- tibble(nb_singletons = nb.sing, pct = 1-thresh.sing)

ggplot(sing_vs_thresh_df) + geom_boxplot(aes(x = pct, y = max.dist.t, group = pct),data = dist_vs_thresh_df) + geom_point(aes(x=pct, y = nb_singletons),data = sing_vs_thresh_df, color = "red") + geom_line(aes(x=pct, y = nb_singletons),data = sing_vs_thresh_df, color = "red")+scale_y_continuous(trans='log10')


max.dist.k <- c()
k.dist <- c()
nb.sing.k <- c()
k.sing <- c()
for (k in (3:20)) {
  print(k)
  g = 10
  n.pc = 1:5 # number of first PC to use
  k.knn = k # number of neighbors to connect to each spot
  pct = 0.6

  print("Creating metaspots")
  # By default, SCimplify_SpatialDLS computes distances in a parallalized way. By default, all the available pus are used. If your computer doesn't suppont, you can change the
  MC.well7.tmp <- SCimplify_SpatialDLS_v2(X = well7.mtx ,
                                          spotPositions = spotPosition,
                                          col.x = "X",
                                          col.y = "Y",
                                          method_similarity = "1", split_not_connected = TRUE, method_reduction = "PCA", gamma = g,
                                          n.dim = n.pc,
                                          method_knn = "1",
                                          k.knn = k.knn,
                                          method_normalization = "log_normalize",
                                          pct = pct
                                          #cell. annotation = well7.md$Main_molecular_cell_type
  )
  l <- 0
  for (m in unique(MC.well7.tmp$membership)) {
    cell.id <- names(MC.well7.tmp$membership[MC.well7.tmp$membership== m])
    #print(cell.id)
    #print(spotPosition[cell.id,])
    if (length(cell.id) > 1){
      max.dist.tmp <- max(dist(spotPosition[cell.id,]))
      max.dist.k <- c(max.dist.k,max.dist.tmp)
      l <- l+1
    }
  }
  #plot(boxplot(max.dist))
  k.dist <- c(k.dist,rep(k,l))
  nb.sing.k <- c(nb.sing.k,length(MC.well7.tmp$supercell_size[MC.well7.tmp$supercell_size==1]))
  k.sing <- c(k.sing,k)
}

dist_vs_k_df <- tibble(max_dists = max.dist.k, k = k.dist)
sing_vs_k_df <- tibble(nb_singletons = nb.sing.k, k = k.sing)

ggplot(sing_vs_k_df) + geom_boxplot(aes(x = k, y = max.dist.k, group = k),data = dist_vs_k_df) + geom_point(aes(x=k, y = nb_singletons),data = sing_vs_k_df, color = "red") + geom_line(aes(x=k, y = nb_singletons),data = sing_vs_k_df, color = "red")+scale_y_continuous(trans='log10')
