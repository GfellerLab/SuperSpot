#' Creation of polygons around the metaspots
#'
#' This function detects creates of around the metaspots based of the memberships of spots
#'
#'
#' @param MC Metaspot object obtained from SCimplify_SpatialDLS function
#' @param spotpositions a data.frame of the coordinates of the spots
#' @param annotation a vector of the annotation for each spot
#' @param concavity a number to define the concavity of the polygons
#' @param membership_name a character of the name where the memberships have been stored in the MC object
#'
#' @return a dataframe with columns
#' \itemize{
#'   \item x - x coordinates of the polygons
#'   \item y - y coordinates of the polygons
#'   \item cell_type - annotation of the polygons
#'   \item membership - membership of the polygons
#' }
#'
#' @export
#'
#'
#'


supercell_metaspots_shape <- function(MC, spotpositions,annotation,concavity,membership_name = "membership"){
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

#' Plotting metaspots on the original spots
#'
#' This function projects the metaspots on the original spots based on the spatial coordinates
#'
#'
#' @param original_coord dataframe of the coordinates of the original spots
#' @param MC Metaspot object obtained from SCimplify_SpatialDLS function
#' @param sc.col a character of the name of the annotation of the spots to get in the meta.data
#' @param sc.col2 a character of the name of the corresponding membership of the spots to get in the meta.data
#' @param polygons_col a character of the name where the polygons have been stored in the MC object
#' @param alpha parameter of the alpha for the spots
#' @param meta_data dataframe containing both the annotations of the spots and their membership
#' @param alpha_hull parameter of the alpha for the polygons of the metaspots
#'
#' @return a ggplot2 plot
#'
#' @export
#'
#'
#'

SpatialDimPlotSC <- function(original_coord, 
                              MC,
                              sc.col = NULL,
                              sc.col2 = NULL,
                              polygons_col,
                              alpha = 1,
                              meta_data,
                              alpha_hull = 0.5) {
  
  membership <- MC$membership
  seuratCoord <- original_coord
  
  seuratCoordMetacell <-  cbind(seuratCoord,membership)

  hull_df_final <- MC[[polygons_col]]
  
  seuratCoordMetacell <- data.frame(seuratCoordMetacell)
  seuratCoordMetacell[[sc.col]] <- meta_data[[sc.col]]
  seuratCoordMetacell[[sc.col2]] <- meta_data[[sc.col2]]
  seuratCoordMetacell[["cell_type"]] <- seuratCoordMetacell[[sc.col]]
  seuratCoordMetacell[["MC_membership"]] <- seuratCoordMetacell[[sc.col2]]
  
  seuratCoordMetacell[["MC_membership"]][duplicated(seuratCoordMetacell[["MC_membership"]]) == TRUE] -> dupli.memb
  seuratCoordMetacell[! seuratCoordMetacell[["MC_membership"]] %in% dupli.memb,] -> seuratCoord.uni
  
  
  p <- ggplot2::ggplot(seuratCoordMetacell) + 
    ggplot2::geom_point(aes(x = imagecol,y = imagerow,color = cell_type),alpha = alpha) + 
    geom_polygon(data = hull_df_final, aes(x = x, y = y, fill = cell_type, group = membership), alpha = alpha_hull, color = "black", linetype = "solid")+
    geom_point(data = seuratCoord.uni,mapping = aes(x = imagecol, y=imagerow),size = 0.5)+
    scale_y_reverse()
  
  
  return(p)
}
