rm(list = ls())

load(file = "data/extreme_points_three_node_graph.Rdata")


#------------------------------------------------------------------------------#
# do projection from three nodes down to two nodes
# keep mobius with meanings
# "3k1", "p3.bar"
#------------------------------------------------------------------------------#

vmat.proj.q <- cbind(vmat.prob_mobius.q[, c(1,2)], vmat.prob_mobius.q[, c("3k1", "p3.bar")])

## remove redundant constrains
vmat.proj.minimal.q = redundant(vmat.proj.q, representation = "V")$output
colnames(vmat.proj.minimal.q)[c(3, 4)] = c("2k1", "k2")
rownames(vmat.proj.minimal.q) = c("2k1", "k2")

## save
save(vmat.proj.minimal.q, file = 'data/three_nodes_to_two_nodes_projection.Rdata')
