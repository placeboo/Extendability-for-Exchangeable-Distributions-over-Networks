rm(list = ls())

# load extreme points of six-node graph
load(file = "data/extreme_points_six_node_graph.Rdata")

#------------------------------------------------------------------------------#
# do projection from six nodes down to two nodes
# keep mobius with meanings
# first two cols + "k6.bar", "k2_4k1"
#------------------------------------------------------------------------------#

vmat.proj.q <- cbind(vmat.prob_mobius.q[, c(1,2)], vmat.prob_mobius.q[, c("k6.bar", "k2_4k1")])

## remove redundant constrains
vmat.proj.minimal = redundant(vmat.proj.q, representation = "V")
vmat.proj.minimal.q = vmat.proj.minimal$output

## save
save(vmat.proj.minimal, file = 'data/six_nodes_to_two_nodes_projection.Rdata')
write.csv(vmat.proj.minimal.q[, -c(1,2)], "data/six_nodes_to_two_nodes_projection.csv")
