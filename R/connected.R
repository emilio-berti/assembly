#' @title Create random pool of connected species
#' @param web adjacency matrix.
#' @param S number of species.
#' @param verbose TRUE/FALSE.
#' @param max_iter maximum number of iterations.
#' @param plot TRUE/FALSE if to plot the community.
#' @examples
#' a <- connected_pool(web, 13, verbose = FALSE)
#' g <- graph_from_adjacency_matrix(web[a, a], "undirected")
#' while (igraph::components(g, "strong")$no > 1) {
#'   message("Not strongly connected - rerunning")
#'   a <- connected_pool(web, 13, verbose = FALSE)
#'   g <- graph_from_adjacency_matrix(web[a, a], "undirected")
#' }
#' plot(g)
#' g <- graph_from_adjacency_matrix(web[a, a], "directed")
#' all(as.numeric(names(as.numeric(names((degree(g, mode = "in") == 0))))) %in% basal)
connected_pool <- function(
  web,
  S,
  verbose = FALSE,
  max_iter = 1000,
  plot = TRUE,
  iter = 0 #do not modify this
) {
  pool <- rep(NA, S)
  basal <- as.vector(which(colSums(web) == 0)) # get basal species
  diff <- S - sum(!is.na(pool))
  i <- 0 #counter for while loop exit condition
  while (diff != 0 & i < 100) {
    # if (verbose) {
    #   if (i %% 10 == 0) {
    #     message("Missing ", diff, " species")
    #   }
    # }
    available <- setdiff(seq_len(nrow(web)), pool) #get available species
    sp <- sample(available, diff) #get random species
    basal_present <- intersect(basal, union(pool, sp))
    if (length(basal_present) == 0) {
      i <- i + 1
      next
    }
    pool <- replace_species(pool, basal_present)
    # remove species that have not connected to present basal species
    has_resource <- colSums(web[basal_present, union(pool, sp), drop = FALSE]) > 0
    no_resource <- which(!has_resource)
    no_resource <- setdiff(no_resource, basal)
    sp <- setdiff(sp, no_resource)
    if (length(sp) == 0) {
      i <- i + 1
      next
    }
    pool <- replace_species(pool, sp)
    pool <- sort(unique(pool))
    isolated <- which(colSums(web[pool, pool, drop = FALSE]) == 0)
    isolated <- setdiff(pool[isolated], basal)
    pool[which(pool %in% isolated)] <- NA
    diff <- S - sum(!is.na(pool))
    i <- i + 1
  }
  if (!check_all_species_have_resources(pool, web)) {
    if (iter < max_iter) {
      if (verbose) {
        message("Unconnected species detected - starting another iteration")
      }
      pool <- connected_pool(web, S, verbose, max_iter, iter = iter + 1)
    } else {
      pool <- NULL
      warning("Unconnected species detected")
    }
  } else {
    if (verbose) {
      message("Connected species pool found")
    }
    if (plot) {
      plot(graph_from_adjacency_matrix(web[pool, pool]))
    }
    return(pool)
  }
}

#' @title Replace species
#' @param species_list vector of species with some NAs to be replaced.
#' @param replacements candidates for replacement.
#' @return the original vector \code{species_list} with replaced species.
replace_species <- function(
  species_list,
  replacements
) {
  replacements <- replacements[!replacements %in% species_list]
  to_replace <- which(is.na(species_list))
  if (length(replacements) < length(to_replace)) {
    to_replace <- to_replace[seq_len(length(replacements))]
  } else if (length(replacements) > length(to_replace)) {
    replacements <- sample(replacements, length(to_replace))
  }
  species_list[to_replace] <- replacements
  return(species_list)
}

#' @param pool species list of the local community.
#' @param web the metaweb.
#' @return TRUE/FALSE
check_all_species_have_resources <- function(
  pool,
  web,
  verbose = FALSE
) {
  basal_species <- which(colSums(web) == 0)
  if (!all(which(colSums(web[pool, pool, drop = FALSE]) == 0) %in% basal_species)) {
    if (verbose) {
      message("Unconnected species detected")
    }
    return(FALSE)
  } else {
    if (verbose) {
      message("All species connected")
    }
    return(TRUE)
  }
}
