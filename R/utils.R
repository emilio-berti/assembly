#' @title Name metaweb species
#'
#' @details Gives the metaweb matrix rownames and columnames drawing randomly
#'   from a sequene of letters.
#' @param metaweb  adjacency matrix of the metaweb, to use to order basal species first.
#'
#' @return the named metaweb.
name_metaweb <- function(metaweb) {
  if (!is.null(rownames(metaweb)) | !is.null(rownames(metaweb))) {
    stop("Names are already present")
  }
  if (dim(metaweb)[1] != dim(metaweb)[2]) {
    stop("metaweb is not symmetric")
  }
  l <- dim(metaweb)[1]
  n <- sapply(seq_len(l), function(x) {
    paste(sample(letters, 10, replace = TRUE), collapse = "")
  })
  if (any(table(n) > 1)) name_metaweb(metaweb)
  rownames(metaweb) <- n
  colnames(metaweb) <- n
  return(metaweb)
}

#' @title Plot food web
#'
#' @param sp.names adjancency matrix of the food web.
#' @param metaweb adjacency matrix of the metaweb, to use to order basal species first.
#' @param title string, title of the plot.
#'
#' @return NULL
show_fw <- function(sp.names, metaweb, title = NULL) {
  # reorder species names to match order of metweb
  sp.names <- intersect(colnames(metaweb), sp.names)
  # get adjacency matrix
  adj <- metaweb[sp.names, sp.names]
  # get basals and consumers
  sp_basals <- .basals(metaweb)
  sp_basals <- which(colnames(adj) %in% sp_basals)
  sp_consumers <- .consumers(metaweb)
  sp_consumers <- which(colnames(adj) %in% sp_consumers)
  # reorder adjacency matrix to have basals first
  adj <- adj[c(sp_basals, sp_consumers), c(sp_basals, sp_consumers)]
  S <- nrow(adj)
  adj <- adj[S:1, ]
  adj <- t(adj)
  adj[1:length(sp_basals), ] <- -1
  # plot
  oldpar <- par(no.readonly = TRUE) #to restore default
  par(mar = c(1, 1, 2, 2))
  image(adj, col = c("brown", "goldenrod", "steelblue"),
        frame = FALSE, axes = FALSE)
  title(title)
  grid(nx = S, ny = S, lty = 1, col = adjustcolor("grey20", alpha.f = .1))
  par(oldpar) #restore default
}

#' @title Plot graph
#'
#' @param sp.names adjancency matrix of the food web.
#' @param metaweb adjacency matrix of the metaweb, to use to order basal species first.
#' @param title string, title of the plot.
#'
#' @return NULL
show_graph <- function(sp.names, metaweb, title = NULL) {
  g <- graph_from_adjacency_matrix(metaweb[sp.names, sp.names])
  # add role
  V(g)$role <- NA
  sp_basals <- .basals(metaweb)
  V(g)$role[V(g)$name %in% sp_basals] <- "basal"
  sp_consumers <- .consumers(metaweb)
  V(g)$role[V(g)$name %in% sp_consumers] <- "consumer"
  # add colors
  V(g)$col <- ifelse(V(g)$role == "basal", "forestgreen", "goldenrod")
  # plot
  oldpar <- par(no.readonly = TRUE)
  par(mar = c(0, 0, 0, 0))
  plot(g, vertex.label = NA,
       vertex.color = V(g)$col,
       vertex.size = 10,
       edge.arrow.size = .5)
  par(oldpar)
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
  sp <- colnames(metaweb)[colSums(metaweb) > 0]
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

#' @title Top consumer species in the metaweb
#'
#' @param metaweb adjacency matrix of the meta food web (metaweb).
#'
#' @return vectors of top consumer species names.
.top <- function(metaweb) {
  sp <- colnames(metaweb)[colSums(metaweb) > 0 & rowSums(metaweb) == 0]
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
  basal_sp <- intersect(basal_sp, sp.names)
  consumer_sp <- .consumers(metaweb)
  consumer_sp <- intersect(consumer_sp, sp.names)
  adj <- metaweb[sp.names, sp.names]
  # get isolated
  isolated_basals <- rowSums(adj[basal_sp, consumer_sp, drop = FALSE]) == 0
  isolated_basals <- names(isolated_basals[which(isolated_basals)])
  isolated_consumers <- colSums(adj[union(basal_sp, consumer_sp), #need to consider also present species that are not basals
                                    consumer_sp, drop = FALSE]) == 0
  isolated_consumers <- names(isolated_consumers[which(isolated_consumers)])
  isolated <- union(isolated_basals, isolated_consumers)
  isolated <- intersect(sp.names, isolated)
  return(isolated)
}

#' @title Find replacement species
#'
#' @param sp.names vector of species names.
#' @param isolated vector of species names that are isolated.
#' @param metaweb adjacency matrix of the metaweb.
#' @param keep.n.basal logical, if to keep the constant number of basal species.
#'
#' @return a vector or replacement names.
.find_replacements <- function(sp.names, isolated, metaweb, keep.n.basal = FALSE) {
  if (keep.n.basal) {
    # search available only in consumers
    n_basals <- length(intersect(isolated, .basals(metaweb)))
    available <- setdiff(.consumers(metaweb), sp.names)
  } else {
    n_basals <- 0
    available <- setdiff(colnames(metaweb), sp.names)
  }
  replacements <- draw_random_species(n = length(isolated) - n_basals,
                                      sp.names = available)
  if (keep.n.basal) {
    replacements <- c(replacements,
                      sample(setdiff(.basals(metaweb), sp.names), n_basals))
    new_basals <- length(intersect(replacements, .basals(metaweb)))
    if (new_basals != n_basals) stop("Number of basal species changed")
  } else {
    new_basals <- length(intersect(replacements, .basals(metaweb)))
    if (new_basals == 0) stop("No basal species")
  }
  return(replacements)
}

#' @title Number of connected components
#'
#' @details Return the number of connected components in the community
#'
#' @param sp.names vector of species names.
#' @param metaweb adjacency matrix of the metaweb.
.components <- function(sp.names, metaweb) {
  g <- graph_from_adjacency_matrix(metaweb[sp.names, sp.names])
  return (components(g)$no)
}

#' @title Number of modules of the food web
#'
#' @details Return the number of modules in the network. Modules are obtained
#'   using a fast greedy algorithm on the undirected food web.
#'
#' @param sp.names vector of species names.
#' @param metaweb adjacency matrix of the metaweb.
.modules <- function(sp.names, metaweb) {
  g <- graph_from_adjacency_matrix(metaweb[sp.names, sp.names])
  cl <- cluster_fast_greedy(as.undirected(g))
  return (length(cl))
}

#' @title Trophic niches
#'
#' @details Calculate the trophic niche of species as the mean and variance of
#'   the trophic level of their resources.
#'
#' @param sp.names vector of species names.
#' @param metaweb adjacency matrix of the metaweb.
#'
#' @return matrix with mean and variance of the trophic level of species'
#'   resources. These are both zero if the species is a basal resource.
trophic_niche <- function(sp.names, metaweb) {
  web <- metaweb[sp.names, sp.names]
  tl <- trophic_levels(sp.names, metaweb)[, 1]
  niche <- matrix(NA, ncol = 2, nrow = length(tl))
  rownames(niche) <- names(tl)
  colnames(niche) <- c("mean", "var")
  for (x in names(tl)) {
    diet_breadth <- tl[names(which(web[, colnames(web) == x] > 0))]
    niche[x, ] <- c(mean(diet_breadth), var(diet_breadth))
  }
  niche[is.nan(niche[, "mean"]), "mean"] <- 0
  niche[is.na(niche[, "var"]), "var"] <- 0
  return(niche)
}
