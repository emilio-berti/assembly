# Code used to create the assembly::adirondack dataset.
#
# The original source was the GATEWAy 1.0 database:
# https://idata.idiv.de/ddm/Data/ShowData/283?version=3
#
# library(igraph)
#
# .create_d <- function() {
#   gw <- read.csv("~/Documents/databases/GATEWAy.csv")
#   gw <- gw[gw$study.site == "Adirondack lakes", c("res.taxonomy", "con.taxonomy")]
#   g <- igraph::graph_from_data_frame(gw)
#   g <- igraph::simplify(g)
#   adj <- igraph::as_adj(g, sparse = FALSE)
#   any(colSums(adj) == 0 & rowSums(adj) == 0)
#   adirondack <- adj
#   usethis::use_data(adirondack)
# }
