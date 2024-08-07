% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCimplify_SpatialDLS.R
\name{split_membership}
\alias{split_membership}
\title{Update of the memberships of unconnected components}
\usage{
split_membership(m, MC)
}
\arguments{
\item{m}{given membership}

\item{MC}{Metaspot object obtained from SCimplify_SpatialDLS function}
}
\value{
a vector of membership with updated split memberships
}
\description{
This function checks if the spots assigned as a metaspots are connected in the KNN graph and split them by updating their membership
}
