#' @title Metropolis-Hastings algorithm
#'
#' @param old.value is the similarity value of the starting species pool.
#' @param new.value is the similarity value of the changed species pool.
#' @param t is the 'temperature' of the system.
#'
#' @return boolean, if to accept the replacement or not.
metropolis.hastings <- function(old.value, new.value, t = 0){
  accept <- FALSE
  if (new.value < old.value){
    # good move, we always accept
    accept <- TRUE
  } else {
    # bad move case
    # it will be accepted given a probability defined by the temperature 't'
    delta.energy <- (old.value - new.value) / old.value
    t <- ifelse(t == 0 & delta.energy == 0, 1e-6, t) #to avoid 0 / 0
    if (exp(delta.energy / t) > runif(1)) {
      accept <- TRUE
    }
  }
  return(accept)
}

#' @title One move of the limiting similarity algorithm
#'
#' @param sp.names vector, names of the species in the meta-community.
#' @param metaweb adjacency matrix of the meta food web (metaweb).
#' @param t is the 'temperature' of the system.
#' @param method character, same as in similarity (igraph). Options are
#'   'jaccard', 'dice', and 'invlogweighted'.
#' @param stat character, statistic used to summarize similarity. Currently
#'   available are c("mean", "sum", "max").
#' @param mode character, which edges are used to compute similarity. Currently
#'   available are c("all", "in", "out").
#' @param return.similarity logical, if to return the difference in similarity
#'   score.
#'
#' @details Similarities are calculated using igraph as matrices. To summarize
#'   these into species-level metrics, the argument "stat" is needed. When stat
#'   = "mean", probability of removal of species is proportional to the average
#'   of their similarities, etc. Global similarities, i.e. of the whole food
#'   web, are also summarized depending on the "stat" argument in a similar way.
#'
#'   When mode = 'in', similarity is computed using only 'in' links, i.e.
#'   species are considered similar if they share similar resources, but not
#'   necessarily have similar consumers. The opposite is treu when mode = 'out'.
#'   When mode = 'all' (default) all edges are considered.
#'
#' @return A vector with the species names.
.move <- function(
    sp.names,
    metaweb,
    t = 0,
    method = "jaccard",
    stat = "mean",
    mode = "all",
    return.similarity = FALSE
) {
  g <- graph_from_adjacency_matrix(metaweb[sp.names, sp.names])
  consumers <- intersect(sp.names, .consumers(metaweb))
  simil <- similarity(g,
                      vids = which(sp.names %in% consumers),
                      method = method,
                      mode = mode)
  diag(simil) <- NA
  prob_removed <- apply(simil, MARGIN = 2, stat, na.rm = TRUE)
  remove <- sample(consumers, size = 1, prob = prob_removed)
  # replace and check new similarity
  repl <- .find_replacements(sp.names, remove, metaweb,
                             keep.n.basal = TRUE) #avoid pick a basal
  new.sp <- union(setdiff(sp.names, remove), repl)
  # if isolated detected, skip move
  if (length(.find_isolated(new.sp, metaweb) > 0)) {
    if (return.similarity) {
      return (list(species = sp.names, similarity = NA))
    } else{
      return(sp.names)
    }
  }
  # if more than one components, skip move
  if (.components(new.sp, metaweb) > 1) {
    if (return.similarity) {
      return (list(species = sp.names, similarity = NA))
    } else{
      return(sp.names)
    }
  }
  # else get new similarity
  new.g <- graph_from_adjacency_matrix(metaweb[new.sp, new.sp])
  consumers <- intersect(new.sp, .consumers(metaweb)) #basal species are filtered this way
  new.simil <- similarity(new.g,
                          vids = which(new.sp %in% consumers),
                          method = method,
                          mode = mode)
  diag(new.simil) <- NA
  # summarize similarities
  if (stat == "mean") {
    simil <- mean(simil, na.rm = TRUE)
    new.simil <- mean(new.simil, na.rm = TRUE)
  } else if (stat == "sum") {
    simil <- sum(simil, na.rm = TRUE)
    new.simil <- sum(new.simil, na.rm = TRUE)
  } else if (stat == "max") {
    simil <- max(simil, na.rm = TRUE)
    new.simil <- max(new.simil, na.rm = TRUE)
  } else {
    stop("'stat' must be one of c('mean', 'sum', 'max')")
  }
  if (return.similarity) {
    # compare old and new similarity
    if (metropolis.hastings(simil, new.simil, t = t)) {
      return (list(species = new.sp, similarity = new.simil))
    } else {
      return (list(species = sp.names, similarity = simil))
    }
  } else {
    # compare old and new similarity
    if (metropolis.hastings(simil, new.simil, t = t)) {
      return (new.sp)
    } else {
      return (sp.names)
    }

  }
}

#' @title Limiting similarity filtering assembly
#'
#' @param sp.names vector, names of the species in the meta-community.
#' @param metaweb adjacency matrix of the meta food web (metaweb).
#' @param t is the 'temperature' of the system.
#' @param method character, same as in similarity (igraph). Options are
#'   'jaccard', 'dice', and 'invlogweighted'.
#' @param max.iter is the maximum number of moves computed.
#' @param stat character, statistic used to summarize similarity. Currently
#'   available are c("mean", "sum", "max").
#' @param mode character, which edges are used to compute similarity. Currently
#'   available are c("all", "in", "out").
#' @param output.verbose logical, if to return also: 1) the average for species'
#'   trophic niches, defined as mean and variance of the trophic level of their
#'   resources. 2) the similiarity score. If output.verbose == TRUE, these are
#'   returned for each assembly move. Note that this will slow down computations
#'   significantly.
#'
#' @details Similarities are calculated using igraph as matrices. To summarize
#'   these into species-level metrics, the argument "stat" is needed. When stat
#'   = "mean", probability of removal of species is proportional to the average
#'   of their similarities, etc. Global similarities, i.e. of the whole food
#'   web, are also summarized depending on the "stat" argument in a similar way.
#'
#'   When mode = 'in', similarity is computed using only 'in' links, i.e.
#'   species are considered similar if they share similar resources, but not
#'   necessarily have similar consumers. The opposite is treu when mode = 'out'.
#'   When mode = 'all' (default) all edges are considered.
#'
#'   The temperature parameter 't' specifies the degree of stochasticity of the
#'   algorithm. For t > 0, an unfavourable move can be accepted, if it passes a
#'   probabilistic acceptance criterion. For t = 0, stochasticity is removed and
#'   the algorithm becomes purely deterministic. This may be result in some
#'   unwarranted behavior, e.g. strongly modular food webs.
#'
#' @return A vector with the species names.
similarity_filtering <- function(
    sp.names,
    metaweb,
    t = 0,
    method = "jaccard",
    stat = "mean",
    mode = "all",
    max.iter = 1e3,
    output.verbose = FALSE
) {
  if (t == 0) message("Temperature 't' = 0; this is a purely deterministic filtering")
  # check for isolated species and stop if detected
  isolated <- .find_isolated(sp.names, metaweb)
  if (length(isolated) > 0) stop("Isolated species detected in input")
  if (.components(sp.names, metaweb) > 1) stop("Isolated components detected in input")
  # start iterations
  new_sp <- sp.names
  if (output.verbose) {
    mean_niche <- rep(NA, max.iter)
    var_niche <- rep(NA, max.iter)
    simil = rep(NA, max.iter)
  }
  for (i in seq_len(max.iter)) {
    new_sp <- .move(new_sp, metaweb, t, method, stat, return.similarity = output.verbose)
    if (output.verbose) {
      old_niche <- trophic_niche(sp.names, metaweb)
      new_niche <- trophic_niche(new_sp$species, metaweb)
      mean_niche[i] <- mean(new_niche[, "mean"])
      var_niche[i] <- mean(new_niche[, "var"])
      simil[i] <- new_sp$similarity
      new_sp <- new_sp$species #retain only vector of species names
    }
  }
  if (length(sp.names) != length(new_sp)) {
    stop("Number of species changed")
  }
  # check for isolated species and stop if detected
  isolated <- .find_isolated(new_sp, metaweb)
  if (length(isolated) > 0) stop("Isolated species detected in output")
  if (.components(new_sp, metaweb) > 1) warning("Isolated component detected in output")
  # return types depends if niche is calculated
  if (output.verbose) {
    return(list(species = new_sp,
                mean.niche = mean_niche,
                var.niche = var_niche,
                similarity = simil))
  } else {
    return(new_sp)
  }
}
