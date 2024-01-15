basicST <- setClass("basicST", slots = c(count = "data.frame",positions = "data.frame", meta_data = "data.frame", color = "data.frame", for_plot = "data.frame"))

#couls <- scales::hue_pal()(length(md$Classification %>% unique()))
#couls_df <- tibble(md$Classification %>% unique(), couls = couls)

CreateBasicSTobject <- function(countmatrix, positions, meta_data){
  couls <- scales::hue_pal()(length(meta_data[,1] %>% unique()))
  couls_df <- tibble(meta_data[,1] %>% unique(), couls = couls)
  positions <- positions[sort(rownames(positions)),]
  meta_data <- meta_data[sort(rownames(meta_data)),]
  if ( is.data.frame(meta_data) == F){
    annotation_df <- tibble(Classification = meta_data)
  }
  else {
    #annotation_df <- tibble(Classification = annotation[,1])
    annotation_df <- meta_data
    #colnames(annotation_df) <- c("Classification")
  }
  rownames(annotation_df) <- sort(rownames(meta_data))
  for_plot <- cbind(positions,annotation_df)
  bST <- basicST(count =countmatrix, positions = positions, meta_data = annotation_df, color = couls_df,for_plot = for_plot)
}

##bST <- CreateBasicSTobject(countmatrix = mtx %>% as.data.frame(), positions = spotPositions, annotation = select(md,"Classification"))


basicSTplot <- function( bST.object, techno , polys = F, md.coul, md.poly) {
  if (techno == "10x"){
    if (polys == F){
    ggplot(data = bST.object@for_plot) + geom_point(mapping = aes(imagecol %>% as.numeric(),imagerow %>% as.numeric(), colour = bST.object@for_plot[,md.coul] ))+scale_y_reverse()
    }
    else {
      ggplot(data = bST.object@for_plot) + geom_point(mapping = aes(imagecol %>% as.numeric(),imagerow %>% as.numeric(), colour = bST.object@for_plot[,md.coul]))+
           ggforce::geom_mark_hull( aes(group = bST.object@for_plot[,md.poly],x = imagecol %>% as.numeric(), y = imagerow %>% as.numeric()),
                                    expand = 0.0001,radius = 0.0001)+scale_y_reverse() + NoLegend()
      #table(ST.so@meta.data$grid_clusts)
    }
  }
  
  else if (techno == "Nano"){
    if (polys == F){
      ggplot(data = bST.object@for_plot) + geom_point(mapping = aes(imagecol %>% as.numeric(),imagerow %>% as.numeric(),
                                                                  colour =bST.object@for_plot[,md.coul]), size = 0.3, pch = ".")+scale_y_reverse()
    }
    else {
      ggplot(data = bST.object@for_plot) + geom_point(mapping = aes(imagecol %>% as.numeric(),imagerow %>% as.numeric(),
                                                                    colour =bST.object@for_plot[,md.coul]), size = 0.3, pch = ".") +
        ggforce::geom_mark_hull( aes(group =bST.object@for_plot[,md.poly],x = imagecol %>% as.numeric(), y = imagerow %>% as.numeric()),expand = 0.01,radius = 0.01)+scale_y_reverse() + NoLegend()
    }
  }
}

#basicSTplot(bST,"10x")


#bST2 <- CreateBasicSTobject(countmatrix = mtx <- nano.obj@assays$Nanostring@counts %>% as.data.frame(), positions =  GetTissueCoordinates(nano.obj ) %>% column_to_rownames("cell"), annotation = select(nano.obj@meta.data,"predicted.annotation.l1"))

#basicSTplot(bST2,"Nano")                          

SpatialDimPlotSC <- function(original_coord,
                      MC,
                      metacell.col = NULL,
                      dims = c(1,2),
                      reduction = "umap",
                      #color = NULL,
                      bST,techno) {
 
  annot <- c(bST@annotation$Classification %>% unique(),"Metacell")
  annot <- sort(annot)
  couleurs <- c(scales::hue_pal()(length(annot)))
  indice <- which(annot == "Metacell")
  couleurs[indice] <- "black"
  
  membership <- MC$membership
  

  seuratCoord <- original_coord
  seuratCoordMetacell <-  cbind(seuratCoord,membership)
  
  centroids <- aggregate(seuratCoord %>% as.matrix() ~ membership,seuratCoord,mean) #should be taken from object slot
  
  centroids$size <- MC$supercell_size
  
  if(is.null(metacell.col)) {
    metacell.col <- "SuperCell_MC"
    centroids[[metacell.col]] <- rep("Metacell",length(MC$supercell_size))
  } else {
    centroids[[metacell.col]] <- MC[[metacell.col]][,1]
  }
  
  p <- basicSTplot(bST,techno)+
    geom_point(data=centroids,aes_string(colnames(centroids)[1+dims[2]],
                                         colnames(centroids)[1+dims[1]],
                                         color = metacell.col,size = "size")) 
  #if (!is.null(color)) {
  p <- p + scale_color_manual(values = couleurs)
  #}
  
  return(p)
  #plot(ggplot(data=centroids,size=0) + geom_density2d_filled(aes(x= UMAP_1,y=UMAP_2)) + NoLegend())
}


SpatialDimPlotSC2 <- function(original_coord, 
                      seurat.mc, 
                      metacell.col = NULL,
                      sc.col = NULL,
                      sc.col2 = NULL,
                      dims = c(1,2),
                      #reduction = "wnn.umap",
                      mc.color = NULL,
                      sc.color = NULL,
                      polygons_col,
                      alpha = 1,
                      pt_size = 0,
                      metric = "supercell_size",
                      continuous_metric = F,
                      meta_data,
                      concavity = 1) {
  
  membership <- seurat.mc$membership
  
  hull_df_final <- seurat.mc[[polygons_col]]
  
  #print(head(hull_df_final))
  
  #seurat$Metacell <- membership
  # p <- PCAPlot(sc,group.by = "Metacell") + NoLegend()
  # centroids <- aggregate(cbind(PC_1,PC_2)~Metacell,p$data,var)
  # boxplot(centroids$PC_1)
  seuratCoord <- original_coord
  
  seuratCoordMetacell <-  cbind(seuratCoord,membership)
  
  centroids <- stats::aggregate(seuratCoord %>% as.matrix() ~membership,seuratCoord,mean) #should be taken from object slot
  centroids[[metric]] <- seurat.mc[[metric]]
  
  #print(head(centroids))
  
  if(is.null(metacell.col)) {
    metacell.col <- "SuperCell_MC"
    #centroids[[metacell.col]] <- rep("red",length(seurat.mc$supercell_size))
    centroids[[metacell.col]] <- as.character(1:length(seurat.mc$supercell_size))
  } else {
    centroids[[metacell.col]] <- seurat.mc[[metacell.col]]#[,1]
  }
  
  if(!is.null(sc.col)){
    if (!is.null(sc.col2)){
    seuratCoord <- data.frame(seuratCoord)
    seuratCoord[[sc.col]] <- meta_data[[sc.col]]
    seuratCoord[[sc.col2]] <- meta_data[[sc.col2]]
    #print(seuratCoord[[sc.col2]][1:5])
    
    #p <- ggplot2::ggplot(seuratCoord,
    #                     aes_string(colnames(seuratCoord)[dims[2]],
    #                                colnames(seuratCoord)[dims[1]],
    #                                color=sc.col)) +
    #  ggplot2::geom_point(size=pt_size, alpha = alpha) + 
    #  ggforce::geom_mark_hull(aes_string(group = sc.col2, ),expand = 0.0001,radius = 0.0001)
    p <- ggplot2::ggplot(seuratCoord) +
      ggplot2::geom_point( aes_string(colnames(seuratCoord)[dims[2]],
                                      colnames(seuratCoord)[dims[1]],
                                      colour=sc.col),size=pt_size, alpha = alpha) + 
      geom_polygon(data = hull_df_final, aes(x = x, y = y, #fill = cell_type, 
                                             group = membership),
                   alpha = 0, color = "black", linetype = "solid")
    }
    else{
      seuratCoord <- data.frame(seuratCoord)
      seuratCoord[[sc.col]] <- meta_data[[sc.col]]
      p <- ggplot2::ggplot(seuratCoord,
                           aes_string(colnames(seuratCoord)[dims[2]],
                                      colnames(seuratCoord)[dims[1]],
                                      color=sc.col)) +
        ggplot2::geom_point(size=pt_size, alpha = alpha)#+scale_y_reverse()
    }
  
  }
  else{
    p <- ggplot2::ggplot(data.frame(seuratCoord),
                         aes_string(colnames(seuratCoord)[dims[2]],
                                    colnames(seuratCoord)[dims[1]])) +
      ggplot2::geom_point(size=pt_size,color = "grey",alpha = 1)#+scale_y_reverse()
  } 
  
  if(!continuous_metric){
    p <- p + ggplot2::geom_point(data=centroids,aes_string(colnames(centroids)[1+dims[2]],
                                                           colnames(centroids)[1+dims[1]],
                                                           fill = metacell.col,size = metric),colour="black",pch=21) +scale_y_reverse()
  }else{
    p <- p + ggplot2::geom_point(data=centroids,aes_string(colnames(centroids)[1+dims[2]],
                                                           colnames(centroids)[1+dims[1]],
                                                           fill = metric,label = "membership"),colour="black",pch=21,size=2) +scale_y_reverse()
  } 
  
  if(!is.null(metacell.col) & !is.null(sc.col) & !is.null(mc.color)){
    if(metacell.col == sc.col & !is.null(mc.color)){
      sc.color = mc.color
    } 
  } 
  
  if (!is.null(mc.color) & !continuous_metric) {
    p <- p + ggplot2::scale_fill_manual(values = mc.color) +  ggplot2::theme_classic() +scale_y_reverse()
  }
  #if (!is.null(sc.color)) {
   # p <- p + ggplot2::scale_color_manual(values = sc.color) +  ggplot2::theme_classic()
  #}
  return(p)
  #plot(ggplot(data=centroids,size=0) + geom_density2d_filled(aes(x= UMAP_1,y=UMAP_2)) + NoLegend())
}

SpatialDimPlotSC3 <- function(original_coord, 
                              seurat.mc, 
                              metacell.col = NULL,
                              sc.col = NULL,
                              sc.col2 = NULL,
                              dims = c(1,2),
                              #reduction = "wnn.umap",
                              mc.color = NULL,
                              sc.color = NULL,
                              alpha = 1,
                              pt_size = 0,
                              metric = "supercell_size",
                              continuous_metric = F,
                              meta_data,
                              concavity = 1,
                              expand = 0.0001,
                              radius = 0.0001) {
  
  membership <- seurat.mc$membership
  seuratCoord <- original_coord
  
  seuratCoordMetacell <-  cbind(seuratCoord,membership)
  #print(head(seuratCoordMetacell))
  
  centroids <- stats::aggregate(seuratCoord %>% as.matrix() ~membership,seuratCoord,mean) #should be taken from object slot
  centroids[[metric]] <- seurat.mc[[metric]]
  
  #print(head(centroids))
  #print(dim(centroids))
  
  if(is.null(metacell.col)) {
    metacell.col <- "SuperCell_MC"
    #centroids[[metacell.col]] <- rep("red",length(seurat.mc$supercell_size))
    centroids[[metacell.col]] <- as.character(1:length(seurat.mc$supercell_size))
  } else {
    centroids[[metacell.col]] <- seurat.mc[[metacell.col]]#[,1]
  }
  table(seuratCoordMetacell$membership) %>% unname() -> centroids$supercell_size
  
  #print(head(seuratCoord))
  seuratCoordMetacell$membership[duplicated(seuratCoordMetacell$membership) == TRUE] -> dupli.memb
  #print(dupli.memb)
  seuratCoordMetacell[! seuratCoordMetacell$membership %in% dupli.memb,] -> seuratCoord.uni
  
  two_size_mc <- subset(centroids, supercell_size == 2)$membership
  seuratCoordMetacell[seuratCoordMetacell$membership %in% two_size_mc,] -> seuratCoord.double
  
  seuratCoordMetacell <- data.frame(seuratCoordMetacell)
  seuratCoordMetacell[[sc.col]] <- meta_data[[sc.col]]
  seuratCoordMetacell[[sc.col2]] <- meta_data[[sc.col2]]
  seuratCoordMetacell[["cell_type"]] <- seuratCoordMetacell[[sc.col]]
  seuratCoordMetacell[["MC_membership"]] <- seuratCoordMetacell[[sc.col2]]
  head(seuratCoordMetacell) %>% print()
      p <- ggplot2::ggplot(seuratCoordMetacell) + 
        ggplot2::geom_point(aes(x = imagecol,y = imagerow,color = cell_type)) + 
        ggforce::geom_mark_hull(data = seuratCoordMetacell ,
                                mapping = aes(group = MC_membership, x = imagecol, y=imagerow, fill = cell_type),
                                expand = expand,radius = radius,concavity = concavity)+
        geom_point(data = seuratCoord.uni,mapping = aes(x = imagecol, y=imagerow),size = 0.5)+geom_line(data = seuratCoord.double,mapping = aes(x = imagecol, y = imagerow, group = membership))+
        scale_y_reverse()
    
    
  return(p)
}

SpatialDimPlotSC4 <- function(original_coord, 
                              seurat.mc, 
                              metacell.col = NULL,
                              sc.col = NULL,
                              sc.col2 = NULL,
                              dims = c(1,2),
                              #reduction = "wnn.umap",
                              mc.color = NULL,
                              sc.color = NULL,
                              alpha = 1,
                              pt_size = 0,
                              metric = "supercell_size",
                              continuous_metric = F,
                              meta_data,
                              concavity = 1,
                              expand = 0.0001,
                              radius = 0.0001) {
  
  membership <- seurat.mc$membership
  seuratCoord <- original_coord
  
  seuratCoordMetacell <-  cbind(seuratCoord,membership)
  #print(head(seuratCoordMetacell))
  
  centroids <- stats::aggregate(seuratCoord %>% as.matrix() ~membership,seuratCoord,mean) #should be taken from object slot
  centroids[[metric]] <- seurat.mc[[metric]]
  
  #print(head(centroids))
  #print(dim(centroids))
  
  if(is.null(metacell.col)) {
    metacell.col <- "SuperCell_MC"
    #centroids[[metacell.col]] <- rep("red",length(seurat.mc$supercell_size))
    centroids[[metacell.col]] <- as.character(1:length(seurat.mc$supercell_size))
  } else {
    centroids[[metacell.col]] <- seurat.mc[[metacell.col]]#[,1]
  }
  table(seuratCoordMetacell$membership) %>% unname() -> centroids$supercell_size
  
  #print(head(seuratCoord))
  seuratCoordMetacell$membership[duplicated(seuratCoordMetacell$membership) == TRUE] -> dupli.memb
  #print(dupli.memb)
  seuratCoordMetacell[! seuratCoordMetacell$membership %in% dupli.memb,] -> seuratCoord.uni
  
  two_size_mc <- subset(centroids, supercell_size == 2)$membership
  seuratCoordMetacell[seuratCoordMetacell$membership %in% two_size_mc,] -> seuratCoord.double
  
  seuratCoordMetacell <- data.frame(seuratCoordMetacell)
  seuratCoordMetacell[[sc.col]] <- meta_data[[sc.col]]
  seuratCoordMetacell[[sc.col2]] <- meta_data[[sc.col2]]
  seuratCoordMetacell[["cell_type"]] <- seuratCoordMetacell[[sc.col]]
  seuratCoordMetacell[["MC_membership"]] <- seuratCoordMetacell[[sc.col2]]
  #head(seuratCoordMetacell) %>% print()
  #head(seuratCoordMetacell[!(seuratCoordMetacell$membership %in% two_size_mc)]) %>% print()
  singleton_doublets <- union(rownames(seuratCoord.uni),rownames(seuratCoord.double))
  colors_vect <- scales::hue_pal()(length(unique(seuratCoordMetacell[["cell_type"]])))
  colors_palette <- setNames(colors_vect,sort(unique(seuratCoordMetacell[["cell_type"]])))
  p <- ggplot2::ggplot(seuratCoordMetacell) + 
    ggplot2::geom_point(aes(x = imagecol,y = imagerow,color = cell_type)) + 
    ggforce::geom_mark_hull(data = seuratCoordMetacell[!(rownames(seuratCoordMetacell) %in% singleton_doublets),] ,
                            mapping = aes(group = MC_membership, x = imagecol, y=imagerow, fill = cell_type),
                            expand = expand,radius = radius,concavity = concavity)+
    scale_fill_manual(values = colors_palette)+
    geom_point(data = seuratCoord.uni,mapping = aes(x = imagecol, y=imagerow),size = 0.5)+geom_line(data = seuratCoord.double,mapping = aes(x = imagecol, y = imagerow, group = membership))+
    scale_y_reverse()
  
  
  return(p)
}

SpatialDimPlotSC5 <- function(original_coord, 
                              seurat.mc, 
                              metacell.col = NULL,
                              sc.col = NULL,
                              sc.col2 = NULL,
                              dims = c(1,2),
                              #reduction = "wnn.umap",
                              mc.color = NULL,
                              sc.color = NULL,
                              alpha = 1,
                              pt_size = 0,
                              metric = "supercell_size",
                              continuous_metric = F,
                              meta_data,
                              concavity = 1,
                              expand = 0.0001,
                              radius = 0.0001) {
  
  membership <- seurat.mc$membership
  seuratCoord <- original_coord
  
  seuratCoordMetacell <-  cbind(seuratCoord,membership)
  #print(head(seuratCoordMetacell))
  
  centroids <- stats::aggregate(seuratCoord %>% as.matrix() ~membership,seuratCoord,mean) #should be taken from object slot
  centroids[[metric]] <- seurat.mc[[metric]]
  
  #print(head(centroids))
  #print(dim(centroids))
  
  if(is.null(metacell.col)) {
    metacell.col <- "SuperCell_MC"
    #centroids[[metacell.col]] <- rep("red",length(seurat.mc$supercell_size))
    centroids[[metacell.col]] <- as.character(1:length(seurat.mc$supercell_size))
  } else {
    centroids[[metacell.col]] <- seurat.mc[[metacell.col]]#[,1]
  }
  table(seuratCoordMetacell$membership) %>% unname() -> centroids$supercell_size
  
  #print(head(seuratCoord))
  seuratCoordMetacell$membership[duplicated(seuratCoordMetacell$membership) == TRUE] -> dupli.memb
  #print(dupli.memb)
  seuratCoordMetacell[! seuratCoordMetacell$membership %in% dupli.memb,] -> seuratCoord.uni
  
  two_size_mc <- subset(centroids, supercell_size == 2)$membership
  seuratCoordMetacell[seuratCoordMetacell$membership %in% two_size_mc,] -> seuratCoord.double
  three_size_mc <- subset(centroids, supercell_size == 3)$membership
  seuratCoordMetacell[seuratCoordMetacell$membership %in% three_size_mc,] -> seuratCoord.triple
  
  seuratCoordMetacell <- data.frame(seuratCoordMetacell)
  seuratCoordMetacell[[sc.col]] <- meta_data[[sc.col]]
  seuratCoordMetacell[[sc.col2]] <- meta_data[[sc.col2]]
  seuratCoordMetacell[["cell_type"]] <- seuratCoordMetacell[[sc.col]]
  seuratCoordMetacell[["MC_membership"]] <- seuratCoordMetacell[[sc.col2]]
  #head(seuratCoordMetacell) %>% print()
  #head(seuratCoordMetacell[!(seuratCoordMetacell$membership %in% two_size_mc)]) %>% print()
  singleton_doublets <- union(rownames(seuratCoord.uni),rownames(seuratCoord.double))
  p <- ggplot2::ggplot(seuratCoordMetacell) + 
    ggplot2::geom_point(aes(x = imagecol,y = imagerow,color = cell_type)) + 
    ggpubr::stat_chull(data = seuratCoordMetacell,aes(group = MC_membership, x = imagecol, y=imagerow, fill = cell_type),alpha=0.5,geom ="polygon",color = "black")+
    geom_point(data = seuratCoord.uni,mapping = aes(x = imagecol, y=imagerow),size = 0.5)+
    #geom_line(data = seuratCoord.double,mapping = aes(x = imagecol, y = imagerow, group = membership))+
    #geom_line(data = seuratCoord.triple,mapping = aes(x = imagecol, y = imagerow, group = membership))+
    scale_y_reverse()
  
  
  return(p)
}

SpatialDimPlotSC6 <- function(original_coord, 
                              seurat.mc, 
                              metacell.col = NULL,
                              sc.col = NULL,
                              sc.col2 = NULL,
                              polygons_col,
                              dims = c(1,2),
                              alpha = 1,
                              pt_size = 0,
                              metric = "supercell_size",
                              meta_data,
                              concavity = 2,
                              alpha_hull = 0.5) {
  
  membership <- seurat.mc$membership
  seuratCoord <- original_coord
  
  seuratCoordMetacell <-  cbind(seuratCoord,membership)
  #print(head(seuratCoordMetacell))
  
  #hull_df_final <- data.frame(x = c(NULL), y = c(NULL), cell_type = c(NULL))
  #for (memb in unique(seurat.mc$membership)){
  #  index.tmp <- which(seurat.mc$membership == memb)
  #  polygons.tmp <- concaveman::concaveman(as.matrix(seuratCoord[index.tmp,]),concavity = concavity)
  #  index.tmp2 <- which(names(seurat.mc[[metacell.col]]) == memb)
  #  hull_df.tmp <- data.frame(x = polygons.tmp[, 2], y = polygons.tmp[, 1],
  #                            cell_type = rep(seurat.mc[[metacell.col]][index.tmp2],nrow(polygons.tmp)),
  #                            membership = rep(memb,nrow(polygons.tmp)))
  #  hull_df_final <- rbind(hull_df_final,hull_df.tmp)
  #}
  hull_df_final <- seurat.mc[[polygons_col]]
  
  #print(head(seuratCoord))
  #seuratCoordMetacell$membership[duplicated(seuratCoordMetacell$membership) == TRUE] -> dupli.memb
  #print(dupli.memb)
  #seuratCoordMetacell[! seuratCoordMetacell$membership %in% dupli.memb,] -> seuratCoord.uni
  
  seuratCoordMetacell <- data.frame(seuratCoordMetacell)
  seuratCoordMetacell[[sc.col]] <- meta_data[[sc.col]]
  seuratCoordMetacell[[sc.col2]] <- meta_data[[sc.col2]]
  seuratCoordMetacell[["cell_type"]] <- seuratCoordMetacell[[sc.col]]
  seuratCoordMetacell[["MC_membership"]] <- seuratCoordMetacell[[sc.col2]]
  #head(seuratCoordMetacell) %>% print()
  #head(seuratCoordMetacell[!(seuratCoordMetacell$membership %in% two_size_mc)]) %>% print()
  
  seuratCoordMetacell[["MC_membership"]][duplicated(seuratCoordMetacell[["MC_membership"]]) == TRUE] -> dupli.memb
  #print(dupli.memb)
  seuratCoordMetacell[! seuratCoordMetacell[["MC_membership"]] %in% dupli.memb,] -> seuratCoord.uni
  
  
  p <- ggplot2::ggplot(seuratCoordMetacell) + 
    ggplot2::geom_point(aes(x = imagecol,y = imagerow,color = cell_type),alpha = alpha) + 
    geom_polygon(data = hull_df_final, aes(x = x, y = y, fill = cell_type, group = membership), alpha = alpha_hull, color = "black", linetype = "solid")+
    geom_point(data = seuratCoord.uni,mapping = aes(x = imagecol, y=imagerow),size = 0.5)+
    scale_y_reverse()
  
  
  return(p)
}
