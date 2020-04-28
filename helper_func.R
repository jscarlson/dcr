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

# create cluster vector of appropriate length, per NA omissions
#' @export
na_fix_data <- function(specdata, clustvar, regvars) {
  return(na.omit(specdata[,c(regvars, clustvar)]))
}

# create categorical variables for dyads containing member i
#' @export
dyad_categ_create <- function(dyad_mem1, dyad_mem2, data) {

  # start up params
  unique.dyad.mem <- unique(c(data[,dyad_mem1], data[,dyad_mem2])) # unique members
  unique.dyad.mem <- unique.dyad.mem[order(unique.dyad.mem)] # order members
  N_dyad <- length(unique.dyad.mem)
  dyad.by.obs <- data[,c(dyad_mem1, dyad_mem2)] # create dyad matrix

  # create categorical variables for dyads containing i
  for(i in 1:N_dyad) {
    dyad.mem.i <- unique.dyad.mem[i] # set member i
    dyad.with.i <- apply(dyad.by.obs, 1, function(x) as.numeric(dyad.mem.i %in% x)) # identify all dyads with member i
    dyad.category <- dyad.with.i*gp.tag + (1-dyad.with.i)*1:nrow(dyad.by.obs) # give unique group tag to dyads containing i
    namelencurr <- length(names(new_data))
    new_data <- cbind(data, dyad.category)
    names(new_data)[namelencurr + 1] <- paste0("con_", dyad.mem.i)
  }

  # return updated dataframe
  return(new_data)

}

###############
### SOURCES ###
###############

# Samii, Cyrus, 2015, "Cluster-Robust Variance Estimation for Dyadic Data",
# https://doi.org/10.7910/DVN/OMJYE5, Harvard Dataverse, V2, UNF:6:WJJ3ZmDS7COvpy1kwztcMQ==[fileUNF]
