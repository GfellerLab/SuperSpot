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
  message("Creating polygons for visualization")
  hull_df_final <- data.frame(x = c(NULL), y = c(NULL), cell_type = c(NULL))
  i <- 0
  pb <- progress::progress_bar$new(
    format = "  Progress [:bar] :percent ETA: :eta",
    total = length(unique(MC[[membership_name]])),    # Total number of iterations
    clear = FALSE,               # If TRUE, clears the progress bar when done
    width = 60                   # Width of the progress bar
  )
  for (memb in unique(MC[[membership_name]])){
    index.tmp <- which(MC[[membership_name]] == memb)
    polygons.tmp <- concaveman::concaveman(as.matrix(spotpositions[index.tmp,]),concavity = concavity)
    index.tmp2 <- which(names(MC[[annotation]]) == memb)
    hull_df.tmp <- data.frame(x = polygons.tmp[, 2], y = polygons.tmp[, 1],
                              cell_type = rep(MC[[annotation]][index.tmp2],nrow(polygons.tmp)),
                              membership = rep(memb,nrow(polygons.tmp)))
    hull_df_final <- rbind(hull_df_final,hull_df.tmp)
    i <- i+1
    #print( paste0((i/length(unique(MC[[membership_name]])))*100,"%") )
    pb$tick()
  }
  message("Done")
  return(hull_df_final)
}

supercell_metaspots_shape_v2 <- function(MC, spotpositions, annotation, concavity, membership_name = "membership") {
  message("Creating polygons for visualization")

  # Load pbapply for progress bar
  if (!requireNamespace("pbapply", quietly = TRUE)) {
    install.packages("pbapply")
  }
  library(pbapply)

  unique_memberships <- unique(MC[[membership_name]])

  # Use pblapply to apply function over unique memberships with progress bar
  hull_list <- pbapply::pblapply(unique_memberships, function(memb) {
    index.tmp <- which(MC[[membership_name]] == memb)
    polygons.tmp <- concaveman::concaveman(as.matrix(spotpositions[index.tmp, ]), concavity = concavity)
    index.tmp2 <- which(names(MC[[annotation]]) == memb)
    hull_df.tmp <- data.frame(
      x = polygons.tmp[, 2],
      y = polygons.tmp[, 1],
      cell_type = rep(MC[[annotation]][index.tmp2], nrow(polygons.tmp)),
      membership = rep(memb, nrow(polygons.tmp))
    )
    return(hull_df.tmp)
  })

  # Combine the list of data frames into one data frame
  hull_df_final <- do.call(rbind, hull_list)

  message("Done")
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
#' @param col.x a character of the name of the x coordinates to get from the coordinates data frame
#' @param col.y a character of the name of the y coordinates to get from the coordinates data frame
#' @param polygons_col a character of the name where the polygons have been stored in the MC object
#' @param alpha parameter of the alpha for the spots
#' @param meta_data dataframe containing both the annotations of the spots and their membership
#' @param alpha_hull parameter of the alpha for the polygons of the metaspots
#' @param spot.color vector for the spot colors
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
                             col.x = "imagecol",
                             col.y = "imagerow",
                              polygons_col,
                              alpha = 1,
                              meta_data,
                              alpha_hull = 0.5,
                             spot.color = NULL) {

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

  if (is.null(spot.color)){
    p <- ggplot2::ggplot(seuratCoordMetacell) +
      ggplot2::geom_point(ggplot2::aes(x = .data[[col.x]],y = .data[[col.y]],color = cell_type),
                          alpha = alpha) +
      ggplot2::geom_polygon(data = hull_df_final,
                            ggplot2::aes(x = x, y = y, fill = cell_type, group = membership),
                   alpha = alpha_hull,
                   color = "black",
                   linetype = "solid")+
      ggplot2::geom_point(data = seuratCoord.uni,
                 mapping = ggplot2::aes(x = .data[[col.x]], y=.data[[col.y]]),
                 size = 0.5)+
      ggplot2::scale_y_reverse()
  }
  else{
    p <- ggplot2::ggplot(seuratCoordMetacell) +
      ggplot2::geom_point(ggplot2::aes(x = .data[[col.x]],y = .data[[col.y]],color = cell_type),
                          alpha = alpha) +
      ggplot2::geom_polygon(data = hull_df_final,
                            ggplot2::aes(x = x, y = y, fill = cell_type, group = membership),
                   alpha = alpha_hull,
                   color = "black",
                   linetype = "solid") +
      ggplot2::scale_fill_manual(values = spot.color) +
      ggplot2::scale_color_manual(values = spot.color) +
      ggplot2::geom_point(data = seuratCoord.uni,
                 mapping = ggplot2::aes(x = .data[[col.x]], y=.data[[col.y]]),
                 size = 0.5)+
      ggplot2::scale_y_reverse()
  }

  return(p)
}
