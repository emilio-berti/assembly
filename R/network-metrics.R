#' @title Calculate trophic level
#'
#' @param sp.names vector, names of the species in the meta-community.
#' @param metaweb adjacency matrix of the meta food web (metaweb).
#'
#' @return numeric vector of trophic levels.
trophic_levels <- function(sp.names, metaweb) {
  fw <- metaweb[sp.names, sp.names]
  fw <- t(fw)
  nn <- rowSums(fw)
  nn[nn == 0] <- 1
  ww <- diag(1 / nn)
  L1 <- ww %*% fw
  L2 <- L1 - diag(rep(1, length(nn)))
  b <- -1 * rep(1, length(nn))
  Tro.lev <- solve(L2) %*% b
  return(Tro.lev)
}
#' @title Compute network metrics
#'
#' @param sp.names vector, names of the species in the meta-community.
#' @param metaweb adjacency matrix of the meta food web (metaweb).
#'
#' @return a data.frame with network metrics.
network <- function(sp.names, metaweb) {
  fw <- metaweb[sp.names, sp.names]
  S <- ncol(fw) #number of species
  L <- sum(fw) #number of links
  C <- L / (S ^ 2) #connectance
  bas <- sum(colSums(fw) == 0) #number of basals
  frB = bas / S #fraction of basals
  top <- sum(rowSums(fw) == 0) #number of top
  frT <- top / S #fraction of top
  frI <- 1 - frT - frB #fraction of intermediate
  frCB <- sum(diag(fw) == 1) / S#fraction cannibals
  # standard deviation of normalized generality
  genk <- (S / L) *apply(fw, 2, sum)
  sdGen <- sd(genk)
  # standard deviation of normalized vulnerability
  vulk <- (S / L) *apply(fw, 1, sum)
  sdVul <- sd(vulk)
  # max. of TL
  # if empty matrix (0 size in both dimensions), or there is no basal species
  # then just create NA
  TL <- tryCatch(trophic_levels(sp.names, metaweb),
                 error = function(e) rep(NA, S))
  maxTL <- max(TL)
  meanTL <- mean(TL)
  sdTL <- sd(TL)
  omn <- sum(TL %% 1 != 0) #number omnivores
  frOmn <- omn / S  #fraction of omnivores
  g <- graph_from_adjacency_matrix(fw)
  Sim <- similarity(g, mode = "all") #Jaccard similarity
  diag(Sim) <- NA
  meanSim <- apply(Sim, 2, mean, na.rm = TRUE)
  MMSim <- mean(meanSim) #mean across all species
  clust <- transitivity(g) #clustering ceofficient
  members <- membership(cluster_fast_greedy(as.undirected(g))) #for modularity
  mod <-  modularity(as.undirected(g), members) #modularity
  CPL <- average.path.length(g, #characteristic path length
                             directed = FALSE,
                             unconnected = TRUE)
  Motifs <- graph.motifs(g, size = 3)
  Motifsm <- matrix(Motifs)
  SumMotifs <- sum(Motifs, na.rm = T)
  PercentMotifs <- Motifs[1:16] / SumMotifs
  PercentMotifsFrame <- as.data.frame(t(PercentMotifs))
  # output of function
  ans <- data.frame(
    connectance = C,
    clustering = clust,
    modularity = mod,
    avg.jaccard.similarity = MMSim,
    sd.generality = sdGen,
    sd.vulnerability = sdVul,
    characteristic.path.length = CPL,
    mean.trophic.level = meanTL,
    max.trophic.level = maxTL,
    sd.trophic.level = sdTL,
    fraction.basal.species = frB,
    fraction.intermediate.species = frI,
    fraction.top.species = frT,
    fraction.omnivore.species = frOmn,
    fraction.collider.motif = PercentMotifsFrame$V3,
    fraction.chain.motif = PercentMotifsFrame$V5,
    fraction.fork.motif = PercentMotifsFrame$V7,
    fraction.IGP = PercentMotifsFrame$V8
  )
  return(ans)
}
