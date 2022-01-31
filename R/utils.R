#' @title Plot food web
#'
#' @param adj adjancency matrix of the food web.
#' @param metweb adjacency matrix of the metaweb, to use to order basal species first.
#' @param title string, title of the plot.
#'
#' @return NULL
show_fw <- function(sp.names, metaweb, title = NULL) {
  adj <- metaweb[sp.names, sp.names]
  sp_basals <- .basals(metaweb)
  sp_basals <- which(colnames(adj) %in% sp_basals)
  sp_consumers <- .consumers(metaweb)
  sp_consumers <- which(colnames(adj) %in% sp_consumers)
  adj <- adj[c(sp_basals, sp_consumers), c(sp_basals, sp_consumers)]
  S <- nrow(adj)
  adj <- adj[nrow(adj):1, ]
  adj <- t(adj)
  adj[1:length(sp_basals), ] <- -1
  image(adj, col = c("brown", "goldenrod", "steelblue"),
        frame = FALSE, axes = FALSE)
  title(title)
  grid(nx = S, ny = S, lty = 1, col = adjustcolor("grey20", alpha.f = .1))
}

#' @title Get food web adjacency matrix
#'
#' @param sp.names vector, names of the species in the meta-community.
#' @param metaweb adjacency matrix of the meta food web (metaweb).
#'
#' @return The adjacency matrix of the local community
.get_adj <- function(sp.names, metaweb) {
  adj <- metaweb[sp.names, sp.names]
  if (any((colSums(adj) + rowSums(adj)) == 0)) {
    message("Local community contains isolated species")
  }
  return(adj)
}

#' @title Consumers in the metaweb
#'
#' @param metaweb adjacency matrix of the meta food web (metaweb).
#'
#' @return vectors of basal species names.
.consumers <- function(metaweb) {
  sp <- colnames(metaweb)[colSums(metaweb) > 1]
  return(sp)
}

#' @title Basal species in the metaweb
#'
#' @param metaweb adjacency matrix of the meta food web (metaweb).
#'
#' @return vectors of basal species names.
.basals <- function(metaweb) {
  sp <- colnames(metaweb)[colSums(metaweb) == 0]
  return(sp)
}

#' @title Find isolated species
#'
#' @param sp.names vector of species names.
#' @param metaweb adjacency matrix of the metaweb.
#'
#' @return a vector with the names of isolated basal species and isolated consumers.
.find_isolated <- function(sp.names, metaweb) {
  # get isolated species in the local pool
  basal_sp <- .basals(metaweb)
  consumer_sp <- .consumers(metaweb)
  adj <- metaweb[sp.names, sp.names]
  # get species index
  basal_id <- which(colnames(adj) %in% basal_sp)
  consumer_id <- which(!colnames(adj) %in% basal_sp)
  # get isolated
  isolated_basals <- which(rowSums(adj[basal_id, consumer_id]) == 0)
  isolated_consumers <- which(colSums(adj[basal_id, consumer_id]) == 0)
  isolated <- union(names(isolated_basals), names(isolated_consumers))
  isolated <- intersect(sp.names, isolated)
  return(isolated)
}

#' @title Find replacement species
#'
#' @param sp.names vector of species names.
#' @param isolated vector of species names that are isolated.
#' @param metaweb adjacency matrix of the metaweb.
#'
#' @return a vector or replacement names.
.find_replacements <- function(sp.names, isolated, metaweb, keep.n.basal = FALSE) {
  if (keep.n.basal) {
    n_basals <- length(intersect(isolated, .basals(metaweb)))
    available <- setdiff(.consumers(metaweb), sp.names)
  } else {
    n_basals <- 0
    available <- setdiff(colnames(metaweb), sp.names)
  }
  replacements <- draw_random_species(n = length(isolated) - n_basals,
                                      sp.names = available)
  if (n_basals > 0) {
    replacements <- c(replacements, sample(setdiff(.basals(metaweb), sp.names),
                                           n_basals))
  }
  return(replacements)
}
