# Author: Jacob Carlson

#################################
#                               #
#                               #
#      "dcr": source code       #
#                               #
#                               #
#################################

##############################################
### Robust variance estimation function(s) ###
##############################################

robust.se.nodfc <- function(model, cluster){
  M <- length(unique(cluster))
  N <- length(cluster)
  dfc <- 1
  uj <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum))
  rcse.cov <- dfc * sandwich(model, meat = crossprod(uj)/N)
  rcse.se <- coeftest(model, rcse.cov)
  return(rcse.cov)
}

robust.se.nodfc.ef <- function(model, ef, cluster) {
  M <- length(unique(cluster))
  N <- length(cluster)
  dfc <- 1
  uj <- apply(ef, 2, function(x) tapply(x, cluster, sum))
  rcse.cov <- dfc * sandwich(model, meat = crossprod(uj)/N)
  rcse.se <- coeftest(model, rcse.cov)
  return(rcse.cov)
}

hc0.robust <- function(model){
  dfc <- 1
  uj <- estfun(model)
  cov.mat <- dfc * sandwich(model, meat = Matrix::crossprod(uj)/nrow(uj))
  return(cov.mat)
}

hc0.robust.ef <- function(model, ef){
  dfc <- 1
  uj <- ef
  cov.mat <- dfc * sandwich(model, meat = Matrix::crossprod(uj)/nrow(uj))
  return(cov.mat)
}

###############
### SOURCES ###
###############

# Samii, Cyrus, 2015, "Cluster-Robust Variance Estimation for Dyadic Data",
# https://doi.org/10.7910/DVN/OMJYE5, Harvard Dataverse, V2, UNF:6:WJJ3ZmDS7COvpy1kwztcMQ==[fileUNF]
