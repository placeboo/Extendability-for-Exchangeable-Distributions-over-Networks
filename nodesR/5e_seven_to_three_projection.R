lrm(list = ls())

load(file = "data/extreme_points_seven_node_graph.Rdata")
#------------------------------------------------------------------------------#
# do projection from seven nodes down to three nodes
# keep mobius with meanings
# 'class1', 'class2', 'class3', 'class4'
#------------------------------------------------------------------------------#
vmat.proj.q = cbind(vmat.prob_mobius.q[, c(1,2)], vmat.prob_mobius.q[, c('class1', 'class2', 'class3', 'class4')])

## remove redundant constrains
mat_after_redundant = redundant(vmat.proj.q, representation = "V")
vmat.proj.minimal.q = mat_after_redundant$output
colnames(vmat.proj.minimal.q) = c(NA, NA, 'class1', 'class2', 'class3', 'class4')
# mark the original index kept
orig.index = c(1: nrow(vmat.proj.q))[-mat_after_redundant$redundant]

#colnames(vmat.proj.minimal.q) = c(NA, NA, "3k1", "p3.bar", "p3", "k3")
rownames(vmat.proj.minimal.q) = rownames(vmat.proj.q)[orig.index]
# what graph they are
load(file = "data/seven_node_subgraph_class_unique.Rdata")
View(seven_node_graphclass_unique[rownames(vmat.proj.minimal.q),])
save(seven_node_graphclass_unique, file = "data/edge_graph_seven_to_three.Rdata")
## save
save(orig.index, vmat.proj.minimal.q, file = 'data/seven_nodes_to_three_nodes_projection.Rdata')
write.csv(vmat.proj.minimal.q, file = 'data/seven_nodes_to_three_nodes_projection.csv')


Hmat = scdd(vmat.proj.minimal.q)$output
# number of ineq
length(which(Hmat[,1] == 0))



