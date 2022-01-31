#' @title Resource filtering assembly
#'
#' @param sp.names vector, names of the species in the meta-community.
#' @param metaweb adjacency matrix of the meta food web (metaweb).
#' @param keep.n.basal TRUE/FALSE, if to keep the number of basal species as
#'   originally passed in sp.names.
#'
#' @return A vector with the species names.
resource_filtering <- function(sp.names, metaweb, keep.n.basal = FALSE) {
  # find isolated species
  isolated <- .find_isolated(sp.names, metaweb)
  # find potential replacements
  # search until no isolated found or reached maximum iterations
  new_sp <- sp.names
  i <- 0
  while(length(isolated) > 0 & i < 1e6) {
    # get potential replacements
    replacements <- .find_replacements(new_sp, isolated, metaweb, keep.n.basal)
    # create new pool
    new_sp <- union(setdiff(new_sp, isolated), replacements)
    isolated <- .find_isolated(new_sp, metaweb)
    i <- i + 1
  }

  if (length(sp.names) != length(new_sp)) {
    stop("New species are not the same number as original pool")
  }

  return(new_sp)
}
