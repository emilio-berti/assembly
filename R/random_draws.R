#' @title Draw random species from a meta-community
#'
#' @param n integer, number of species to draw.
#' @param sp.names vector, names of the species in the meta-community (default =
#'   NULL).
#'
#' @return A vector with the species names.
#'
#' @details Either one of \code{sp.names} or \code{S} need to be specified. If
#'   both are specified, \code{sp.names} overrules \code{S}, which will not be
#'   evaluated.
draw_random_species <- function(n, sp.names = NULL) {
  if (is.null(sp.names)) {
    stop("Specify total species richness if sp.names is not specified")
  }
  sp <- sample(sp.names, n)
  return(sp)
}
