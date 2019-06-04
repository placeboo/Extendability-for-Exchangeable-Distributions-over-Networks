rm(list = ls())

load(file = "data/extreme_points_seven_node_graph.Rdata")
#------------------------------------------------------------------------------#
# do projection from seven nodes down to six nodes
# keep mobius with meanings
# 'class1', 'class2', 'class3', 'class4', ... "class34"
#------------------------------------------------------------------------------#
name_tmp = paste("class", 1:34, sep = '')

vmat.proj.q <- cbind(vmat.prob_mobius.q[, c(1,2)], vmat.prob_mobius.q[, name_tmp])

## remove redundant constrains
mat_after_redundant = redundant(vmat.proj.q, representation = "V")
vmat.proj.minimal.q = mat_after_redundant$output

colnames(vmat.proj.minimal.q) = c(NA, NA, name_tmp)
# mark the original index kept
orig.index = c(1: nrow(vmat.proj.q))[-mat_after_redundant$redundant]
## save
save(orig.index, vmat.proj.minimal.q, file = 'data/seven_nodes_to_five_nodes_projection.Rdata')
