# SuperSpot

Before installing SuperSpot, packages bluster and Supercell are required. You can install them with the following commands :
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("bluster")

if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("GfellerLab/SuperCell")

SuperSpot can be installed with the command:
if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("GfellerLab/SuperSpot")
