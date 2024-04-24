# SuperSpot
## Description
SuperSpot is an R package bringing metacells's concept ([*Baran et al., 2019*](https://doi.org/10.1186/s13059-019-1812-2), [*Ben-Kiki et al. 2022*](https://doi.org/10.1186/s13059-022-02667-1), [*Bilous et al., 2022*](https://doi.org/10.1186/s12859-022-04861-1) and [*Persad et al., 2022*](https://doi.org/10.1038/s41587-023-01716-9)) and extending [SuperCell](https://github.com/GfellerLab/SuperCell) to spatial transcriptomic data.

SuperSpot combines adjacent and transcriptionally similar spots into "metaspots". The process involves representing spots as nodes in a graph with edges connecting spots in spatial proximity and edge weights representing transcriptional similarity. Hierarchical clustering is used to aggregate spots into metaspots at a user-defined resolution.



## Installation
Before installing SuperSpot, packages bluster and SuperCell are required. You can install them with the following commands:
``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("bluster")

if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("GfellerLab/SuperCell")
```

SuperSpot can be installed with the command:
``` r
if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("GfellerLab/SuperSpot")
```

## License
SuperSpot can be used freely by academic groups for non-commercial purposes.
The product is provided free of charge, and, therefore, on an "as is" basis, without warranty of any kind. 
Commercial users are requested to an obtain a separate license for SuperCell. To do so, please contact Nadette Bulgin (nbulgin@lcr.org) at the Ludwig Institute for  Cancer Research Ltd.

For scientific questions, please contact Matei Teleman ([matei.teleman\@unil.ch](mailto:matei.teleman@unil.ch)) or David Gfeller ([David.Gfeller\@unil.ch](mailto:David.Gfeller@unil.ch)).
