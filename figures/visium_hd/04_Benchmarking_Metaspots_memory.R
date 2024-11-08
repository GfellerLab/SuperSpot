library(Seurat)
library(SuperSpot)
library(SuperCell)
library(tidyr)
library(tidyverse)
library(igraph)
library(peakRAM)
#devtools::install_github("RfastOfficial/Rfast2")
#install.packages("Rfast2")

#print(!requireNamespace('Rfast2'))
library(Rfast2)
#print(library(Rfast2))

#starting_times <- c()
#ending_times <- c()


FindSpatiallyVariableFeatures.Seurat_EDITED <- function (object, assay = NULL, slot = "scale.data", features = NULL,
                                                         image = NULL, selection.method = c("markvariogram", "moransi"),
                                                         r.metric = 5, x.cuts = NULL, y.cuts = NULL, nfeatures = 2000,
                                                         verbose = TRUE, ...)
{
  assay <- assay %||% DefaultAssay(object = object)
  features <- features %||% rownames(x = object[[assay]])
  image <- image %||% Seurat:::DefaultImage(object = object) # ADDING ::: for workaround
  tc <- GetTissueCoordinates(object = object[[image]])[,1:2]


  # MODIFICED FOR LOCAL COPY
  object[[assay]] <- FindSpatiallyVariableFeatures.StdAssay_EDITED(object = object[[assay]],
                                                                   slot = slot, features = features, spatial.location = tc,
                                                                   selection.method = selection.method, r.metric = r.metric,
                                                                   x.cuts = x.cuts, y.cuts = y.cuts, nfeatures = nfeatures,
                                                                   verbose = verbose, ...)
  object <- SeuratObject:::LogSeuratCommand(object = object)
}






FindSpatiallyVariableFeatures.StdAssay_EDITED <- function (object, layer = "scale.data", spatial.location, selection.method = c("markvariogram",
                                                                                                                                "moransi"), features = NULL, r.metric = 5, x.cuts = NULL,
                                                           y.cuts = NULL, nfeatures = nfeatures, verbose = TRUE, ...)
{
  features <- features %||% rownames(x = object)
  if (selection.method == "markvariogram" && "markvariogram" %in%
      names(x = Misc(object = object))) {
    features.computed <- names(x = Misc(object = object,
                                        slot = "markvariogram"))
    features <- features[!features %in% features.computed]
  }
  data <- GetAssayData(object = object, layer = layer)
  data <- as.matrix(x = data[features, ])
  data <- data[Seurat:::RowVar(x = data) > 0, ]
  if (nrow(x = data) != 0) {
    svf.info <- FindSpatiallyVariableFeatures(object = data,
                                              spatial.location = spatial.location, selection.method = selection.method,
                                              r.metric = r.metric, x.cuts = x.cuts, y.cuts = y.cuts,
                                              verbose = verbose, ...)
  }
  else {
    svf.info <- c()
  }
  if (selection.method == "markvariogram") {
    if ("markvariogram" %in% names(x = Misc(object = object))) {
      svf.info <- c(svf.info, Misc(object = object, slot = "markvariogram"))
    }
    suppressWarnings(expr = Misc(object = object, slot = "markvariogram") <- svf.info)
    svf.info <- ComputeRMetric(mv = svf.info, r.metric)
    svf.info <- svf.info[order(svf.info[, 1]), , drop = FALSE]
  }
  if (selection.method == "moransi") {
    colnames(x = svf.info) <- paste0("MoransI_", colnames(x = svf.info))
    svf.info <- svf.info[order(svf.info[, 2], -abs(svf.info[,
                                                            1])), , drop = FALSE]
  }
  var.name <- paste0(selection.method, ".spatially.variable")
  var.name.rank <- paste0(var.name, ".rank")
  svf.info[[var.name]] <- FALSE
  svf.info[[var.name]][1:(min(nrow(x = svf.info), nfeatures))] <- TRUE
  svf.info[[var.name.rank]] <- 1:nrow(x = svf.info)

  ############################
  ## This line, attempting to store metadata,  produces error:

  print(">>>> USING EDITED FUNCTION!!!! <<<")

  # Error in `LayerData<-`:
  # ! 'layer' must be a single non-empty string
  #https://github.com/satijalab/seurat/issues/8226
  #object[names(x = svf.info)] <- svf.info

  ## Substituting:
  object <- AddMetaData(object, metadata = svf.info) # Works

  ############################


  return(object)
}

library(Matrix)
dist_to_sparse <- function(d) {
  n <- attr(d, "Size")
  labels <- attr(d, "Labels")
  if (is.null(labels)) {
    labels <- as.character(1:n)
  }

  # Total number of elements in the dist object
  m <- n * (n - 1) / 2

  # Generate indices for the lower triangle
  # The positions correspond to combinations of indices (i, j) where i > j
  get_lower_triangle_indices <- function(n) {
    i <- rep(2:n, times = 1:(n - 1))
    j <- unlist(sapply(2:n, function(x) 1:(x - 1)))
    return(list(i = i, j = j))
  }

  indices <- get_lower_triangle_indices(n)
  i_lower <- indices$i
  j_lower <- indices$j

  # Extract the distances
  x_lower <- as.vector(d)

  # Create indices for both lower and upper triangles to make the matrix symmetric
  i_all <- c(i_lower, j_lower)
  j_all <- c(j_lower, i_lower)
  x_all <- c(x_lower, x_lower)

  # Include the diagonal elements (optional, usually zeros)
  i_diag <- 1:n
  j_diag <- 1:n
  x_diag <- rep(0, n)

  # Combine all indices and values
  i <- c(i_all, i_diag)
  j <- c(j_all, j_diag)
  x <- c(x_all, x_diag)

  # Create the sparse symmetric matrix
  sparse_dist <- sparseMatrix(
    i = i,
    j = j,
    x = x,
    dims = c(n, n),
    symmetric = TRUE,
    dimnames = list(labels, labels)
  )

  return(sparse_dist)
}

RunMoransII <- function (data, pos, verbose = TRUE)
{
  mysapply <- sapply
  if (verbose) {
    message("Computing Moran's I")
    mysapply <- pbapply::pbsapply
  }
  Rfast2.installed <- PackageCheck("Rfast2", error = FALSE)
  if (Rfast2.installed) {
    MyMoran <- Rfast2::moranI
  }
  else if (!PackageCheck("ape", error = FALSE)) {
    stop("'RunMoransI' requires either Rfast2 or ape to be installed",
         call. = FALSE)
  }
  else {
    MyMoran <- ape::Moran.I
    if (getOption("Seurat.Rfast2.msg", TRUE)) {
      message("For a more efficient implementation of the Morans I calculation,",
              "\n(selection.method = 'moransi') please install the Rfast2 package",
              "\n--------------------------------------------",
              "\ninstall.packages('Rfast2')", "\n--------------------------------------------",
              "\nAfter installation of Rfast2, Seurat will automatically use the more ",
              "\nefficient implementation (no further action necessary).",
              "\nThis message will be shown once per session")
      options(Seurat.Rfast2.msg = FALSE)
    }
  }
  print("pos is:")
  print(head(pos))
  print(anyNA(pos))
  print(str(pos))
  #pos$x <- as.numeric(as.character(pos$x))
  #pos$y <- as.numeric(as.character(pos$y))
  #na_rows <- which(is.na(pos$x) | is.na(pos$y))
  # Remove cells with missing spatial data
  #cells_to_keep <- setdiff(Cells(obj), rownames(pos)[na_rows])
  #obj <- subset(obj, cells = cells_to_keep)
  pos.dist <- parallelDist::parDist(x = as.matrix(pos))
  #pos.dist.mat <-Seurat::as.sparse(x = pos.dist)
  pos.dist.mat <-dist_to_sparse(pos.dist)
  print("done distance matrix")
  weights <- 1/pos.dist.mat^2
  diag(x = weights) <- 0
  results <- mysapply(X = 1:nrow(x = data), FUN = function(x) {
    tryCatch(expr = MyMoran(data[x, ], weights), error = function(x) c(1,
                                                                       1, 1, 1))
  })
  pcol <- ifelse(test = Rfast2.installed, yes = 2, no = 4)
  results <- data.frame(observed = unlist(x = results[1, ]),
                        p.value = unlist(x = results[pcol, ]))
  rownames(x = results) <- rownames(x = data)
  return(results)
}

processing <- function(g){
  print(paste0("gamma is: ",g))
  obj <- readRDS(paste0("./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/benchmarking/MC_seurat_g",g,".rds"))
  pos <- Seurat::GetTissueCoordinates(obj)
  cells.to.keep <- subset(pos, x < 22695.42/30, y < 65252.03/30)$cell
  print(length(cells.to.keep))
  obj <- subset(obj, cells = cells.to.keep )
  svf.df <- Seurat::RunMoransI(GetAssayData(obj,assay = "RNA",layer = "scale.data")[VariableFeatures(obj)[1:1000],],GetTissueCoordinates(obj)[,1:2])
  #obj=obj[,unname(which(colSums(GetAssayData(obj))!=0))]
  #obj <- SCTransform(obj)
  #obj <- FindSpatiallyVariableFeatures.Seurat_EDITED(obj, assay = "RNA", features = VariableFeatures(obj)[1:1000],selection.method = "moransi")
  #obj <- FindSpatiallyVariableFeatures(obj, assay = "SCT", features = VariableFeatures(obj)[1:1000],selection.method = "moransi")
  #svf.df_sorted <- obj@assays[["SCT"]]@meta.features[order(obj@assays[["SCT"]]@meta.features$moransi.spatially.variable.rank, decreasing = F), ]
  #svf <- rownames(svf.df_sorted)
  print("Done")
  #print(paste0("top 10 svf are: ",svf[1:10]))
  #jpeg(file=paste0("./01_Data/svf_mc_",g,".jpeg"),width = 3440,height = 1440)
  #plot(ImageFeaturePlot(obj, fov = "pancreas", features =  svf[1:6], max.cutoff = "q95"))
  #dev.off()
  #save(svf,file = paste0("./SuperSpot/01_Data/svf_mc_",g,".rda"))
}



mem.df <- peakRAM::peakRAM(
  processing(64),
  processing(32),
  processing(16),
  processing(8),
  processing(4),
  processing(2),
  processing(1))


saveRDS(mem.df,"./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/benchmarking/results_memory.rds")
write.csv(mem.df,"./SuperSpot/VisiumHD/01_Data/VisiumHDColonCancer/benchmarking/results_memory.csv")

benchmark.results <- read_csv("/./SuperSpot/figures/visium_hd/results_memory.csv")
benchmark.results$Peak_RAM_Used_GB <- (benchmark.results$Peak_RAM_Used_MiB)/953.7
benchmark.results$Elapsed_Time_hour <- (benchmark.results$Elapsed_Time_sec)/3600
plot(c(64,32,16,8,4,2,1),benchmark.results$Peak_RAM_Used_GB,type = "b")
plot(c(64,32,16,8,4,2,1),benchmark.results$Elapsed_Time_hour,type = "b")
