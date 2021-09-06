# Author: Jacob Carlson

#################################
#                               #
#                               #
#      "dcr": source code       #
#                               #
#                               #
#################################

########################
### Helper functions ###
########################

#' Create a dyad identifier variable
#'
#' This function returns a `data.frame` with a new dyad ID variable based
#' on a provided `data.frame` with individual dyad member ID variables.
#'
#' @param data A `data.frame` object containing individual dyad member ID variables.
#' @param name A string for the name of new dyad ID variable to be created.
#' @param dyad_mem1 A string for the name of first dyad member ID variable.
#' @param dyad_mem2 A string for the name of second dyad member ID variable.
#' @param directed A logical value indicating whether the dyad ID should be for the directed or undirected dyad.
#' @return The input `data.frame` with a new dyad ID variable appended to it.
#' @export
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
