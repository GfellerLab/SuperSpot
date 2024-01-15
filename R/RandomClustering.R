random_clustering <- function(cells, n.clusters){
  clusts <- paste0(1:n.clusters)
  md.c <- c()
  md.cl <- c()
  for (c in cells){
    random.cl <- sample(clusts,1)
    md.c <- c(md.c,c)
    md.cl <- c(md.cl,random.cl)
  }
  random.clusters <- tibble(cell = md.c, random_clusters = md.cl)
  return(random.clusters)
}

big.random.clusters.purities <- tibble(.rows = 12)
target <-  select(ST.so@meta.data,"Classification")[sort(rownames(select(ST.so@meta.data,"Classification"))), ]
for (i in 1:10000){
  print(i)
  rdn.clusts <- random_clustering(cells = colnames(countMatrix),12) %>% column_to_rownames("cell")
  rdn.clusts.md <- rdn.clusts[sort(rownames(rdn.clusts)), ]
  links <- data.frame(
    source=rdn.clusts.md, 
    target=target, 
    value=rep(1,length(target))
  )
  method_purity <- c("max_proportion", "entropy")[1]
  purity <- supercell_purity(
    clusters = links$target,
    supercell_membership = links$source, 
    method = method_purity
  )
  big.random.clusters.purities <- cbind(big.random.clusters.purities,purity)
}

big.random.clusters.purities <- cbind(big.random.clusters.purities,rowMeans(big.random.clusters.purities))
summary(big.random.clusters.purities[,ncol(big.random.clusters.purities)])



big.random.clusters.purities2 <- tibble(.rows = 18)
target2 <-  select(nano.obj.sb@meta.data,"predicted.annotation.l1")[sort(rownames(select(nano.obj.sb@meta.data,"predicted.annotation.l1"))), ]
for (i in 1:10000){
  print(i)
  rdn.clusts <- clusteval::random_clustering(mtx %>% t(),18)
  rdn.clusts.df <- as.data.frame(rdn.clusts)
  rownames(rdn.clusts.df) <- colnames(mtx)
  rdn.clusts.md <- rdn.clusts.df[sort(rownames(rdn.clusts.df)), ]
  links <- data.frame(
    source=rdn.clusts.md, 
    target=target2, 
    value=rep(1,length(target2))
  )
  method_purity <- c("max_proportion", "entropy")[1]
  purity <- supercell_purity(
    clusters = links$target,
    supercell_membership = links$source, 
    method = method_purity
  )
  big.random.clusters.purities2 <- cbind(big.random.clusters.purities2,purity)
}

big.random.clusters.purities2 <- cbind(big.random.clusters.purities2,rowMeans(big.random.clusters.purities2))
summary(big.random.clusters.purities2[,ncol(big.random.clusters.purities2)])
