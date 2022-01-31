#' @title Draw random species from a meta-community
#'
#' @param n integer, number of species to draw.
#' @param sp.names vector, names of the species in the meta-community (default =
#'   NULL).
#' @param S integer, meta-community species richness (default = NULL).
#'
#' @return A vector with the species names.
#'
#' @details Either one of \code{sp.names} or \code{S} need to be specified. If
#'   both are specified, \code{sp.names} overrules \code{S}, which will not be
#'   evaluated.
#'
#' @examples
#' S <- 10
#' species <- LETTERS[seq_len(S)]
#' draw_random_species(3, sp.names = species)
#' draw_random_species(3, S = S)
draw_random_species <- function(n, sp.names = NULL, S = NULL) {
  if (is.null(sp.names) & is.null(S)) {
    stop("Specify total species richness if sp.names is not specified")
  }
  if (is.null(sp.names)) sp.names <- seq_len(S)
  sp <- sample(sp.names, n)
  return(sp)
}
