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

#' Compute dyadic clustering robust (DCR) variance estimates
#'
#' This function computes DCR standard error and variance estimates for
#' model parameters.
#'
#' @param model Model object (note: must be compatible with `r sandwich` package functions)
#' @param dyad_id String for name of dyad identifier variable
#' @param dyad_mem1 String for name of first dyad member identifier variable
#' @param dyad_mem2 String for name of second dyad member identifier variable
#' @param data `r data.frame` object containing dyadic data, including dyad and dyad member identifier variables
#' @param posdef Logical value indicating whether or not DCR variance-covariance matrix is forced to be positive semi-definite
#' @return A list containing DCR standard errors and variances for model parameters
#' @export
dcr <- function(model, dyad_id, dyad_mem1, dyad_mem2, data, posdef = TRUE) {

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
    dyad.category <- dyad.with.i*gp.tag + (1-dyad.with.i)*1:nrow(dyad.by.obs) # give unique group tag to dyads containing i
    if (i==1) {
      cov.mat.sum <- sandwich::vcovCL(model, dyad.category, type = "HC0", multi0 = TRUE)
    } else if (i!=1) {
      cov.mat.sum <- cov.mat.sum + sandwich::vcovCL(model, dyad.category, type = "HC0", multi0 = TRUE)
    }
  }

  # substract repeated variance estimator for repeated dyads
  cov.mat.sum.intermed <- cov.mat.sum - sandwich::vcovCL(model, data[, dyad_id], type = "HC0", multi0 = TRUE)

  # substract HC variance estimator
  V.hat <- cov.mat.sum.intermed - (N_dyad-2)*sandwich::vcovHC(model, type="HC0")

  # force posdef
  param.names <- colnames(V.hat)
  if (sum(diag(V.hat) < 0) >= 1 & posdef == TRUE) {
    decomp <- eigen(V.hat, symmetric = TRUE)
    pos_eigens <- pmax(decomp$values, rep.int(0, length(decomp$values)))
    V.hat <- decomp$vectors %*% diag(pos_eigens) %*% t(decomp$vectors)
  }

  # return standard errors
  coef.var <- diag(V.hat)
  coef.se <- sqrt(coef.var)
  outputlst <- list(coef.se, coef.var)
  names(outputlst) <- c("SE", "VAR")
  names(outputlst$SE) <- param.names
  names(outputlst$VAR) <- param.names
  return(outputlst)

}

# parallel version

#' Compute dyadic clustering robust (DCR) variance estimates (parallelized)
#'
#' This function enables the parallel computation of DCR standard error
#' and variance estimates.
#'
#' The default number of cores is set to the ceiling of half the total number of cores
#' available.
#'
#' @param model Model object (note: must be compatible with `r sandwich` package functions)
#' @param dyad_id String for name of dyad identifier variable
#' @param dyad_mem1 String for name of first dyad member identifier variable
#' @param dyad_mem2 String for name of second dyad member identifier variable
#' @param ncore Integer specifying the number of cores to be used for parallel computation
#' @param data `r data.frame` object containing dyadic data, including dyad and dyad member identifier variables
#' @param posdef Logical value indicating whether or not DCR variance-covariance matrix is forced to be positive semi-definite
#' @return A list containing DCR standard errors and variances for model parameters
#' @export
dcr_parallel <- function(model, dyad_id, dyad_mem1, dyad_mem2, ncore = ceiling(parallel::detectCores()/2), data, posdef=TRUE) {

  # set up starting params
  data[[dyad_mem1]] <- as.character(data[[dyad_mem1]]) # convert to character
  data[[dyad_mem2]] <- as.character(data[[dyad_mem2]]) # convert to character
  data[[dyad_id]] <- as.character(data[[dyad_id]]) # convert to character
  gp.tag <- -99
  unique.dyad.mem <- na.omit(unique(c(data[[dyad_mem1]], data[[dyad_mem2]]))) # unique members
  unique.dyad.mem <- unique.dyad.mem[order(unique.dyad.mem)] # order members
  N_dyad <- length(unique.dyad.mem)
  dyad.by.obs <- data[,c(dyad_mem1, dyad_mem2)] # create dyad matrix
  N_cores <- parallel::detectCores() # number of cores available for computation
  id_iter <- 1:N_dyad # list of indices for all dyads

  # sum variance estimators for clustering on all dyads containing member i

  dyad_cont_mem_clust <- function(i, s, unique.dyads, dyad.pairs) {
    dyad.mem.i <- unique.dyads[i] # set member i
    dyad.with.i <- apply(dyad.pairs, 1, function(x) as.numeric(dyad.mem.i %in% x)) # identify all dyads with member i
    dyad.category <- dyad.with.i*gp.tag + (1-dyad.with.i)*1:nrow(dyad.pairs) # give unique group tag to dyads containing i
    if (length(unique(dyad.category)) == 1) { # two few cluster protection
      stop("Stop! Only one cluster!")
    } else {
      cov.mat <- sandwich::vcovCL(model, dyad.category, type = "HC0", multi0 = TRUE)
    }
  }
  cov.mat.sum.mc <- parallel::mclapply(id_iter, dyad_cont_mem_clust, mc.cores = ncore, s = cov.mat.sum, unique.dyads = unique.dyad.mem, dyad.pairs = dyad.by.obs)
  cov.mat.sum <- Reduce("+", cov.mat.sum.mc)

  # substract repeated variance estimator for repeated dyads
  cov.mat.sum.intermed <- cov.mat.sum - sandwich::vcovCL(model, data[, dyad_id], type = "HC0", multi0 = TRUE)

  # substract HC variance estimator
  V.hat <- cov.mat.sum.intermed - (N_dyad-2)*sandwich::vcovHC(model, type="HC0")

  # force posdef
  param.names <- colnames(V.hat)
  if (sum(diag(V.hat) < 0) >= 1 & posdef == TRUE) {
    decomp <- eigen(V.hat, symmetric = TRUE)
    pos_eigens <- pmax(decomp$values, rep.int(0, length(decomp$values)))
    V.hat <- decomp$vectors %*% diag(pos_eigens) %*% t(decomp$vectors)
  }

  # return standard errors
  coef.var <- diag(V.hat)
  coef.se <- sqrt(coef.var)
  outputlst <- list(coef.se, coef.var)
  names(outputlst) <- c("SE", "VAR")
  names(outputlst$SE) <- param.names
  names(outputlst$VAR) <- param.names
  return(outputlst)

}

# all purpose function; slower than dyad_robust but handles more model types more flexibly

#' Compute dyadic clustering robust (DCR) variance estimates (customizable)
#'
#' This function computes DCR standard error and variance estimates for
#' model parameters with more model object handling flexibility than `r dcr`.
#'
#' @param model Model object
#' @param dyad_id String for name of dyad identifier variable
#' @param dyad_mem1 String for name of first dyad member identifier variable
#' @param dyad_mem2 String for name of second dyad member identifier variable
#' @param spec_vars String vector of variable names in model specification
#' @param data `r data.frame` object containing dyadic data, including dyad and dyad member identifier variables
#' @param ef Optional matrix containing empirical estimating functions
#' @param posdef Logical value indicating whether or not DCR variance-covariance matrix is forced to be positive semi-definite
#' @return A list containing DCR standard errors and variances for model parameters
#' @export
dcr_custom <- function(model, dyad_id, dyad_mem1, dyad_mem2, spec_vars, data, ef = NULL, posdef = TRUE) {

  # set up starting params
  data[[dyad_mem1]] <- as.character(data[[dyad_mem1]]) # convert to character
  data[[dyad_mem2]] <- as.character(data[[dyad_mem2]]) # convert to character
  data[[dyad_id]] <- as.character(data[[dyad_id]]) # convert to character
  data <- na_fix_data(data, spec_vars, c(dyad_id, dyad_mem1, dyad_mem2))
  gp.tag <- -99
  unique.dyad.mem <- na.omit(unique(c(data[[dyad_mem1]], data[[dyad_mem2]]))) # unique members
  unique.dyad.mem <- unique.dyad.mem[order(unique.dyad.mem)] # order members
  N_dyad <- length(unique.dyad.mem)
  dyad.by.obs <- data[,c(dyad_mem1, dyad_mem2)] # create dyad matrix

  # substract repeated variance estimator for repeated dyads
  if (is.null(ef)) {
    rc <- robust.se.nodfc(model, data[, dyad_id])
    cov.mat.sum.intermed <- -1*rc
  } else {
    rc <- robust.se.nodfc.ef(model, ef, data[, dyad_id])
    cov.mat.sum.intermed <- -1*rc
  }

  # substract HC variance estimator
  if (is.null(ef)) {
    hc <- hc0.robust(model)
    cov.mat.sum <- cov.mat.sum.intermed - (N_dyad-2)*hc
  } else {
    hc <- hc0.robust.ef(model, ef)
    cov.mat.sum <- cov.mat.sum.intermed - (N_dyad-2)*hc
  }

  # sum variance estimators for clustering on all dyads containing member i
  for(i in 1:N_dyad) {
    dyad.mem.i <- unique.dyad.mem[i] # set member i
    print(paste0("Loop status: ", as.character(round(100*(i/N_dyad),2)), "%")) # completion status
    dyad.with.i <- apply(dyad.by.obs, 1, function(x) as.numeric(dyad.mem.i %in% x)) # identify all dyads with member i
    dyad.category <- dyad.with.i*gp.tag + (1-dyad.with.i)*1:nrow(dyad.by.obs) # give unique group tag to dyads containing i
    if (is.null(ef)) {
      cov.mat.sum <- cov.mat.sum + robust.se.nodfc(model, dyad.category)
    } else {
      cov.mat.sum <- cov.mat.sum + robust.se.nodfc.ef(model, ef, dyad.category)
    }
  }

  # force posdef
  param.names <- colnames(cov.mat.sum)
  if (sum(diag(cov.mat.sum) < 0) >= 1 & posdef == TRUE) {
    decomp <- eigen(cov.mat.sum, symmetric = TRUE)
    pos_eigens <- pmax(decomp$values, rep.int(0, length(decomp$values)))
    cov.mat.sum <- decomp$vectors %*% diag(pos_eigens) %*% t(decomp$vectors)
  }

  # return standard errors
  coef.var <- diag(cov.mat.sum)
  coef.se <- sqrt(coef.var)
  outputlst <- list(coef.se, coef.var)
  names(outputlst) <- c("SE", "VAR")
  names(outputlst$SE) <- param.names
  names(outputlst$VAR) <- param.names
  return(outputlst)

}

###############
### SOURCES ###
###############

# Samii, Cyrus, 2015, "Cluster-Robust Variance Estimation for Dyadic Data",
# https://doi.org/10.7910/DVN/OMJYE5, Harvard Dataverse, V2, UNF:6:WJJ3ZmDS7COvpy1kwztcMQ==[fileUNF]
