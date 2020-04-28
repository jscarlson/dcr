# Author: Jacob Carlson

#################################
#                               #
#                               #
#      "dcr": source code       #
#                               #
#                               #
#################################

##############################################################
### Compute cluster–robust standard errors for dyadic data ###
##############################################################

#' @export
dyad_robust <- function(model, dyad_id, dyad_mem1, dyad_mem2, data) {

  # set up starting params
  data[[dyad_mem1]] <- as.character(data[[dyad_mem1]]) # convert to character
  data[[dyad_mem2]] <- as.character(data[[dyad_mem2]]) # convert to character
  data[[dyad_id]] <- as.character(data[[dyad_id]]) # convert to character
  gp.tag <- -99 # arbitrary tag
  unique.dyad.mem <- na.omit(unique(c(data[[dyad_mem1]], data[[dyad_mem2]]))) # unique members
  unique.dyad.mem <- unique.dyad.mem[order(unique.dyad.mem)] # order members
  N_dyad <- length(unique.dyad.mem)
  dyad.by.obs <- data[,c(dyad_mem1, dyad_mem2)] # create dyad matrix

  # sum variance estimators for clustering on all dyads containing member i
  for(i in 1:N_dyad) {
    dyad.mem.i <- unique.dyad.mem[i] # set member i
    print(paste0("Loop status: ", as.character(round(100*(i/N_dyad),2)), "%")) # completion status
    dyad.with.i <- apply(dyad.by.obs, 1, function(x) as.numeric(dyad.mem.i %in% x)) # identify all dyads with member i
    # print(sum(dyad.with.i == 1))
    dyad.category <- dyad.with.i*gp.tag + (1-dyad.with.i)*1:nrow(dyad.by.obs) # give unique group tag to dyads containing i
    if (i==1) {
      cov.mat.sum <- vcovCL(model, dyad.category, type = "HC0", multi0 = TRUE)
    } else if (i!=1) {
      cov.mat.sum <- cov.mat.sum + vcovCL(model, dyad.category, type = "HC0", multi0 = TRUE)
    }
  }
  print(cov.mat.sum)

  # substract repeated variance estimator for repeated dyads
  cov.mat.sum.intermed <- cov.mat.sum - vcovCL(model, data[, dyad_id], type = "HC0", multi0 = TRUE)
  print(cov.mat.sum.intermed)

  # substract HC variance estimator
  V.hat <- cov.mat.sum.intermed - (N_dyad-2)*vcovHC(model, type="HC0")
  print(V.hat)

  # return standard errors
  coef.var <- diag(V.hat)
  coef.se <- sqrt(coef.var)
  outputlst <- list(coef.se, coef.var)
  return(outputlst[[1]])

}

# parallel version
#' @export
dyad_robust_parallel <- function(model, dyad_id, dyad_mem1, dyad_mem2, data) {

  # set up starting params
  data[[dyad_mem1]] <- as.character(data[[dyad_mem1]]) # convert to character
  data[[dyad_mem2]] <- as.character(data[[dyad_mem2]]) # convert to character
  data[[dyad_id]] <- as.character(data[[dyad_id]]) # convert to character
  gp.tag <- -99
  unique.dyad.mem <- na.omit(unique(c(data[[dyad_mem1]], data[[dyad_mem2]]))) # unique members
  unique.dyad.mem <- unique.dyad.mem[order(unique.dyad.mem)] # order members
  N_dyad <- length(unique.dyad.mem)
  dyad.by.obs <- data[,c(dyad_mem1, dyad_mem2)] # create dyad matrix
  N_cores <- detectCores() # number of cores available for computation
  id_iter <- 1:N_dyad # list of indices for all dyads

  # sum variance estimators for clustering on all dyads containing member i

  dyad_cont_mem_clust <- function(i, s, unique.dyads, dyad.pairs) {
    dyad.mem.i <- unique.dyads[i] # set member i
    dyad.with.i <- apply(dyad.pairs, 1, function(x) as.numeric(dyad.mem.i %in% x)) # identify all dyads with member i
    dyad.category <- dyad.with.i*gp.tag + (1-dyad.with.i)*1:nrow(dyad.pairs) # give unique group tag to dyads containing i
    if (length(unique(dyad.category)) == 1) { # two few cluster protection
      stop("Stop! Only one cluster!")
    } else {
      cov.mat <- vcovCL(model, dyad.category, type = "HC0", multi0 = TRUE)
    }
  }
  cov.mat.sum.mc <- mclapply(id_iter, dyad_cont_mem_clust, mc.cores = ceiling(N_cores/2), s = cov.mat.sum, unique.dyads = unique.dyad.mem, dyad.pairs = dyad.by.obs)
  cov.mat.sum <- Reduce("+", cov.mat.sum.mc)

  # substract repeated variance estimator for repeated dyads
  cov.mat.sum.intermed <- cov.mat.sum - vcovCL(model, data[, dyad_id], type = "HC0", multi0 = TRUE)

  # substract HC variance estimator
  V.hat <- cov.mat.sum.intermed - (N_dyad-2)*vcovHC(model, type="HC0")

  # return standard errors
  coef.var <- diag(V.hat)
  coef.se <- sqrt(coef.var)
  outputlst <- list(coef.se, coef.var)
  return(outputlst[[1]])

}

# all purpose function; slower than dyad_robust but handles more model types more flexibly
#' @export
dyad_robust_any <- function(model, dyad_id, dyad_mem1, dyad_mem2, spec_vars, data, vcvHC = NULL, ef = NULL) {

  # set up starting params
  data[[dyad_mem1]] <- as.character(data[[dyad_mem1]]) # convert to character
  data[[dyad_mem2]] <- as.character(data[[dyad_mem2]]) # convert to character
  data[[dyad_id]] <- as.character(data[[dyad_id]]) # convert to character
  data <- na_fix_data(data, spec_vars, c(dyad_id, dyad_mem1, dyad_mem2))
  gp.tag <- -99
  unique.dyad.mem <- na.omit(unique(c(data[[dyad_mem1]], data[[dyad_mem2]]))) # unique members
  unique.dyad.mem <- unique.dyad.mem[order(unique.dyad.mem)] # order members
  N_dyad <- length(unique.dyad.mem)
  print(N_dyad)
  dyad.by.obs <- data[,c(dyad_mem1, dyad_mem2)] # create dyad matrix

  # substract repeated variance estimator for repeated dyads
  rc <- robust.se.nodfc(model, data[, dyad_id])
  print(sqrt(diag(rc)))
  cov.mat.sum.intermed <- -1*rc

  # substract HC variance estimator
  hc <- hc0.robust(model)
  print(sqrt(diag(hc)))
  cov.mat.sum <- cov.mat.sum.intermed - (N_dyad-2)*hc

  # sum variance estimators for clustering on all dyads containing member i
  for(i in 1:N_dyad) {
    dyad.mem.i <- unique.dyad.mem[i] # set member i
    print(paste0("Loop status: ", as.character(round(100*(i/N_dyad),2)), "%")) # completion status
    dyad.with.i <- apply(dyad.by.obs, 1, function(x) as.numeric(dyad.mem.i %in% x)) # identify all dyads with member i
    dyad.category <- dyad.with.i*gp.tag + (1-dyad.with.i)*1:nrow(dyad.by.obs) # give unique group tag to dyads containing i
    print(sqrt(diag(robust.se.nodfc(model, dyad.category)))[1:5])
    cov.mat.sum <- cov.mat.sum + robust.se.nodfc(model, dyad.category)
  }

  # return standard errors
  coef.var <- diag(cov.mat.sum)
  coef.se <- sqrt(coef.var)
  outputlst <- list(coef.se, coef.var)
  return(outputlst[[1]])

}

###############
### SOURCES ###
###############

# Samii, Cyrus, 2015, "Cluster-Robust Variance Estimation for Dyadic Data",
# https://doi.org/10.7910/DVN/OMJYE5, Harvard Dataverse, V2, UNF:6:WJJ3ZmDS7COvpy1kwztcMQ==[fileUNF]
