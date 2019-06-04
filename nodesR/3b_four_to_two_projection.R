rm(list = ls())
load(file = "data/extreme_points_four_node_graph.Rdata")
#------------------------------------------------------------------------------#
# do projection from three nodes down to two nodes
# keep mobius with meanings
# '4k1', 'co-diamond'
#------------------------------------------------------------------------------#
vmat.proj.q <- cbind(vmat.prob_mobius.q[, c(1,2)], vmat.prob_mobius.q[, c('4k1', 'co-diamond')])

## remove redundant constrains
mat_after_redundant = redundant(vmat.proj.q, representation = "V")
vmat.proj.minimal.q = mat_after_redundant$output
colnames(vmat.proj.minimal.q)[c(3, 4)] = c("2k1", "k2")
rownames(vmat.proj.minimal.q) = c("2k1", "k2")

## save
save(vmat.proj.minimal.q, file = 'data/four_nodes_to_two_nodes_projection.Rdata')

