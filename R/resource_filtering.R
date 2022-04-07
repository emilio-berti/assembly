#' @title Resource filtering assembly
#'
#' @param sp.names vector, names of the species in the meta-community.
#' @param metaweb adjacency matrix of the meta food web (metaweb).
#' @param keep.n.basal TRUE/FALSE, if to keep the number of basal species as
#'   originally passed in sp.names.
#'
#' @return A vector with the species names.
#'
#' @example
#' set.seed(1234)
#' data(adirondack)
#' patch <- draw_random_species(50, colnames(adirondack))
#' show_fw(patch, adirondack)
#' patch_filtered <- resource_filtering(patch, adirondack, keep.n.basal = TRUE)
#' show_fw(patch_filtered, adirondack)
#' setdiff(patch, patch_filtered)
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
    stop("Number of species changed")
  }
  if (.components(new_sp, metaweb) > 1) stop("Isolated component detected")
  return(new_sp)
}
