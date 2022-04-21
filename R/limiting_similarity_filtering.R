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
#'
#' @return A vector with the species names.
.move <- function(
  sp.names,
  metaweb,
  t = 0,
  method = "jaccard"
) {
  # find isolated species
  isolated <- .find_isolated(sp.names, metaweb)
  if (length(isolated) > 0) stop("Isolated species detected in input")
  g <- graph_from_adjacency_matrix(metaweb[sp.names, sp.names])
  consumers <- intersect(sp.names, .consumers(metaweb))
  # BUG ---------------
  # SIMILARITY IS NOW ONLY FOR FIRST CONSMUER.
  simil <- similarity(g,
                      vids = which(sp.names %in% consumers),
                      method = method)
  diag(simil) <- NA
  simil <- colSums(simil, na.rm = TRUE)
  remove <- sample(consumers, size = 1, prob = simil)
  # replace and check new similarity
  repl <- .find_replacements(sp.names, remove, metaweb,
                             keep.n.basal = TRUE) #avoid pick a basal
  new.sp <- union(setdiff(sp.names, remove), repl)
  # if isolated detected, skip
  if(length(.find_isolated(new.sp, metaweb) > 0)) return (sp.names)
  # else get new similarity
  new.g <- graph_from_adjacency_matrix(metaweb[new.sp, new.sp])
  consumers <- intersect(new.sp, .consumers(metaweb)) #basal species are filtered this way
  new.simil <- similarity(new.g,
                          vids = which(new.sp %in% consumers),
                          method = method)
  diag(new.simil) <- NA
  new.simil <- colSums(new.simil, na.rm = TRUE)
  # compare old and new similarity
  if (metropolis.hastings(sum(simil), sum(new.simil), t = t)) {
    return (new.sp)
  } else {
    return (sp.names)
  }
}

#' @title Limiting similarity filtering assembly
#'
#' @param sp.names vector, names of the species in the meta-community.
#' @param metaweb adjacency matrix of the meta food web (metaweb).
#' @param t is the 'temperature' of the system.
#' @param method character, same as in similarity (igraph). Options are
#'   'jaccard', 'dice', and 'invlogweighted'.
#' @param max.iter is the number of maximum iterations.
#'
#' @return A vector with the species names.
similarity_filtering <- function(
  sp.names,
  metaweb,
  t = 0,
  method = "jaccard",
  max.iter = 1e3
) {
  new_sp <- sp.names
  for (i in seq_len(max.iter)) new_sp <- .move(new_sp, metaweb, t = t)
  if (length(sp.names) != length(new_sp)) {
    stop("Number of species changed")
  }
  if (.components(new_sp, metaweb) > 1) stop("Isolated component detected")
  return(new_sp)
}
