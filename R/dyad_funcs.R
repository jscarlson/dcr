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
#' This function computes DCR variance estimates and standard errors
#' (DCRSEs) for model parameters.
#'
#' @param model A model object.
#' @param dyad_mem1 A string for the name of first dyad member identifier variable.
#' @param dyad_mem2 A string for the name of second dyad member identifier variable.
#' @param data A `data.frame` object containing dyadic data.
#' @param posdef A logical value indicating whether or not DCR variance-covariance matrix is forced to be positive semi-definite.
#' @param dofcorr A logical value indicating whether or not to apply a small sample correction to the final DCRSE estimates. This correction is equivalent to multiplying DCRSEs by sqrt(N/(N-1)), where N is the number of unique dyad members in the analytic sample. Correspondingly, when computing p-values, the test statistic should be compared to a t-distribution with DOF = N - 1.
#' @return A list containing DCRSEs and DCR variances for model parameters, as well as the number of unique dyad members in the analytic sample.
#' @export
dcr <- function(model, dyad_mem1, dyad_mem2, data, posdef = TRUE, dofcorr = FALSE) {

  create_dyadid <- function(data, name, dyad_mem1, dyad_mem2, directed = F) {
    dyadframe <- data[,c(dyad_mem1, dyad_mem2)]
    if (directed == T) {
      f <- function(x, y) c(as.character(x), as.character(y))
    } else {
      f <- function(x, y) sort(c(as.character(x), as.character(y)))
    }
    sorteddyads <- mapply(f, dyadframe[[dyad_mem1]], dyadframe[[dyad_mem2]], SIMPLIFY = F)
    dyad_id <- unlist(lapply(sorteddyads, function(x) as.factor(paste0(x, collapse = "_"))), use.names = F)
    data[[name]] <- as.character(dyad_id)
    return(data)
  }

  data[[dyad_mem1]] <- as.character(data[[dyad_mem1]])
  data[[dyad_mem2]] <- as.character(data[[dyad_mem2]])
  data <- create_dyadid(data, "_dcr_dyad_id", dyad_mem1, dyad_mem2, directed = F)

  gp.tag <- -99
  unique.dyad.mem <- na.omit(unique(c(data[[dyad_mem1]], data[[dyad_mem2]])))
  unique.dyad.mem <- unique.dyad.mem[order(unique.dyad.mem)]
  N_dyad <- length(unique.dyad.mem)
  dyad.by.obs <- data[, c(dyad_mem1, dyad_mem2)]

  progress_bar <- txtProgressBar(min = 1, max = N_dyad, style = 3, char = "=")

  print("Computing Dyad Member CRSEs")

  for (i in 1:N_dyad) {

    dyad.mem.i <- unique.dyad.mem[i]
    dyad.with.i <- apply(dyad.by.obs, 1, function(x) as.numeric(dyad.mem.i %in% x))
    dyad.category <- dyad.with.i * gp.tag + (1 - dyad.with.i) * 1:nrow(dyad.by.obs)

    if (i == 1) {
      cov.mat.sum <- multiwayvcov::cluster.vcov(model, dyad.category, df_correction = FALSE)
    } else {
      cov.mat.sum <- cov.mat.sum + multiwayvcov::cluster.vcov(model, dyad.category, df_correction = FALSE)
    }

    setTxtProgressBar(progress_bar, value = i)

  }

  close(progress_bar)

  cov.mat.sum.intermed <- cov.mat.sum - multiwayvcov::cluster.vcov(model, data[, "_dcr_dyad_id"], df_correction = FALSE)
  V.hat <- cov.mat.sum.intermed - (N_dyad - 2) * multiwayvcov::cluster.vcov(model, 1:nrow(data), df_correction = FALSE)

  param.names <- colnames(V.hat)

  if (sum(diag(V.hat) < 0) >= 1 & posdef == TRUE) {
    decomp <- eigen(V.hat, symmetric = TRUE)
    pos_eigens <- pmax(decomp$values, rep.int(0, length(decomp$values)))
    V.hat <- decomp$vectors %*% diag(pos_eigens) %*% t(decomp$vectors)
  }

  if (dofcorr == TRUE) {
    coef.var <- (N_dyad / (N_dyad - 1)) * diag(V.hat)
  } else {
    coef.var <- diag(V.hat)
  }
  coef.se <- sqrt(coef.var)

  outputlst <- list(coef.se, coef.var, N_dyad)
  names(outputlst) <- c("dcrse", "dcrvar", "N_udm")
  names(outputlst$dcrse) <- param.names
  names(outputlst$dcrvar) <- param.names
  return(outputlst)

}

# sandwich version

#' Compute dyadic clustering robust (DCR) variance estimates (`sandwich` version)
#'
#' This function computes DCR variance estimates and standard errors
#' (DCRSEs) for model parameters.
#'
#' @param model A model object.
#' @param dyad_mem1 A string for the name of first dyad member identifier variable.
#' @param dyad_mem2 A string for the name of second dyad member identifier variable.
#' @param data A `data.frame` object containing dyadic data.
#' @param posdef A logical value indicating whether or not DCR variance-covariance matrix is forced to be positive semi-definite.
#' @param dofcorr A logical value indicating whether or not to apply a small sample correction to the final DCRSE estimates. This correction is equivalent to multiplying DCRSEs by sqrt(N/(N-1)), where N is the number of unique dyad members in the analytic sample. Correspondingly, when computing p-values, the test statistic should be compared to a t-distribution with DOF = N - 1.
#' @return A list containing DCRSEs and DCR variances for model parameters, as well as the number of unique dyad members in the analytic sample.
#' @export
dcr_sandwich <- function(model, dyad_id, dyad_mem1, dyad_mem2, data, posdef = TRUE, dofcorr = FALSE) {

  # set up starting params
  data[[dyad_mem1]] <- as.character(data[[dyad_mem1]]) # convert to character
  data[[dyad_mem2]] <- as.character(data[[dyad_mem2]]) # convert to character
  data[[dyad_id]] <- as.character(data[[dyad_id]]) # convert to character
  gp.tag <- -99 # arbitrary tag
  unique.dyad.mem <- na.omit(unique(c(data[[dyad_mem1]], data[[dyad_mem2]]))) # unique members
  unique.dyad.mem <- unique.dyad.mem[order(unique.dyad.mem)] # order members
  N_dyad <- length(unique.dyad.mem)
  dyad.by.obs <- data[,c(dyad_mem1, dyad_mem2)] # create dyad matrix

  # progress bar
  progress_bar <- txtProgressBar(min=1, max=N_dyad, style = 3, char="=")

  # sum variance estimators for clustering on all dyads containing member i
  print("Computing Dyad Member CRSEs")
  for(i in 1:N_dyad) {
    dyad.mem.i <- unique.dyad.mem[i] # set member i
    dyad.with.i <- apply(dyad.by.obs, 1, function(x) as.numeric(dyad.mem.i %in% x)) # identify all dyads with member i
    dyad.category <- dyad.with.i*gp.tag + (1-dyad.with.i)*1:nrow(dyad.by.obs) # give unique group tag to dyads containing i
    if (i==1) {
      cov.mat.sum <- sandwich::vcovCL(model, as.character(dyad.category[-as.numeric(model$na.action)]), type = "HC0", multi0 = TRUE, cadjust = FALSE, fix = FALSE)
    } else if (i!=1) {
      cov.mat.sum <- cov.mat.sum + sandwich::vcovCL(model, as.character(dyad.category[-as.numeric(model$na.action)]), type = "HC0", multi0 = TRUE, cadjust = FALSE, fix = FALSE)
    }
    setTxtProgressBar(progress_bar, value = i)
  }

  # close prog bar
  close(progress_bar)

  # substract repeated variance estimator for repeated dyads
  cov.mat.sum.intermed <- cov.mat.sum - sandwich::vcovCL(model, as.character(data[, dyad_id][[1]][-as.numeric(model$na.action)]), type = "HC0", multi0 = TRUE, cadjust = FALSE, fix = FALSE)

  # substract HC variance estimator
  V.hat <- cov.mat.sum.intermed - (N_dyad - 2) * sandwich::vcovHC(model, type = "HC0")

  # force posdef
  param.names <- colnames(V.hat)
  if (sum(diag(V.hat) < 0) >= 1 & posdef == TRUE) {
    decomp <- eigen(V.hat, symmetric = TRUE)
    pos_eigens <- pmax(decomp$values, rep.int(0, length(decomp$values)))
    V.hat <- decomp$vectors %*% diag(pos_eigens) %*% t(decomp$vectors)
  }

  # return standard errors
  if (dofcorr == TRUE) {
    coef.var <- (N_dyad / (N_dyad - 1)) * diag(V.hat)
  } else {
    coef.var <- diag(V.hat)
  }
  coef.se <- sqrt(coef.var)

  outputlst <- list(coef.se, coef.var, N_dyad)
  names(outputlst) <- c("dcrse", "dcrvar", "N_udm")
  names(outputlst$dcrse) <- param.names
  names(outputlst$dcrvar) <- param.names
  return(outputlst)

}

###############
### SOURCES ###
###############

# Samii, Cyrus, 2015, "Cluster-Robust Variance Estimation for Dyadic Data",
# https://doi.org/10.7910/DVN/OMJYE5, Harvard Dataverse, V2, UNF:6:WJJ3ZmDS7COvpy1kwztcMQ==[fileUNF]
