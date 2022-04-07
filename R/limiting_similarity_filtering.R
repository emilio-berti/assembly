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
    if (t > 0){
      delta.energy <- (old.value - new.value) / old.value
      if (exp(delta.energy / t) > runif(1)) {
        accept <- TRUE
      }
    }
  }
  return(accept)
}

#' @title One move of the limiting similarity algorithm
#'
#' @param sp.names vector, names of the species in the meta-community.
#' @param metaweb adjacency matrix of the meta food web (metaweb).
#' @param t is the 'temperature' of the system.
#'
#' @return A vector with the species names.
.move <- function(sp.names, metaweb, t = 0) {
  # find isolated species
  isolated <- .find_isolated(sp.names, metaweb)
  if (length(isolated) > 0) stop("Isolated species detected in input")
  g <- graph_from_adjacency_matrix(metaweb[sp.names, sp.names])
  consumers <- setdiff(sp.names, .basals(metaweb)) #basal species are filtered this way
  simil <- similarity(g,
                      vids = which(sp.names %in% consumers),
                      method = "jaccard")[, 1]
  remove <- sample(consumers, size = 1, prob = simil)
  # replace and check new similarity
  repl <- .find_replacements(sp.names, remove, metaweb,
                                        keep.n.basal = TRUE) #avoid pick a basal
  new.sp <- union(setdiff(sp.names, remove), repl)
  new.g <- graph_from_adjacency_matrix(metaweb[new.sp, new.sp])
  consumers <- setdiff(new.sp, .basals(metaweb)) #basal species are filtered this way
  new.simil <- similarity(g,
                          vids = which(new.sp %in% consumers),
                          method = "jaccard")[, 1]
  # if isolated detected, skip
  if(length(.find_isolated(new.sp, metaweb) > 0)) return (sp.names)
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
#' @param max.iter is the number of maximum iterations.
#'
#' @return A vector with the species names.
similarity_filtering <- function(sp.names, metaweb, t = 0, max.iter = 1e3) {
  new_sp <- sp.names
  for (i in seq_len(max.iter)) new_sp <- .move(new_sp, metaweb, t = t)
  if (length(sp.names) != length(new_sp)) {
    stop("Number of species changed")
  }
  if (.components(new_sp, metaweb) > 1) stop("Isolated component detected")
  return(new_sp)
}
