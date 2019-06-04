rm(list = ls())

load(file = "data/extreme_points_four_node_graph.Rdata")

#------------------------------------------------------------------------------#
# do projection from three nodes down to two nodes
# keep mobius with meanings
# '4k1', 'co-diamond', 'co-paw', 'co-claw'
#------------------------------------------------------------------------------#

vmat.proj.q <- cbind(vmat.prob_mobius.q[, c(1,2)], vmat.prob_mobius.q[, c('4k1', 'co-diamond', 'co-paw', 'co-claw')])


## remove redundant constrains
mat_after_redundant = redundant(vmat.proj.q, representation = "V")
vmat.proj.minimal.q = mat_after_redundant$output
colnames(vmat.proj.minimal.q) = c(NA, NA, "3k1", "p3.bar", 'p3', 'k3')
rownames(vmat.proj.minimal.q) = c("4k1", "c4.bar", "c4", "k4")

Hmat = scdd(vmat.proj.minimal.q)$output
# number of ineq
length(which(Hmat[,1] == 0))

## save
save(vmat.proj.minimal.q, file = 'data/four_nodes_to_three_nodes_projection.Rdata')
write.csv(vmat.proj.minimal.q, file = 'data/four_nodes_to_three_nodes_projection.csv' )

