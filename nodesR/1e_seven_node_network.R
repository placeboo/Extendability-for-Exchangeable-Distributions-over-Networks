rm(list = ls())

load(file = 'data/seven_node_mobius_unique.Rdata')
load(file = 'data/subgraph_class_table.Rdata')

#------------------------------------------------------------------------------#
# Make hyperplan matrix
#------------------------------------------------------------------------------#
# test whether the mobius unique matrix is correct
# mobius_inv = round(solve(prob_moubius_uniq.mat))
# moubius_inv_class1 = mobius_inv['class1', ]
# diff = moubius_inv_class1 - as.vector(count.mat)
# diff[diff != 0] 

tmp_moubius = prob_moubius_uniq.mat * as.vector(count.mat)
sum_to_one = -1 * apply(tmp_moubius, 2, sum)
prob_mobius.mat2 = rbind(prob_moubius_uniq.mat, sum_to_one)
n.rows = dim(prob_mobius.mat2)[1]

prob_mobius.mat2 = cbind(c(rep(0, n.rows - 1), 1), c(rep(0, n.rows - 1), 1), prob_mobius.mat2)
colnames(prob_mobius.mat2)[c(1,2)] = c("in/eq", "const")
prob_mobius.mat2.q = d2q(prob_mobius.mat2)

# Convert from H-rep to V-rep
vmat <- scdd(prob_mobius.mat2.q, representation="H")
vmat.prob_mobius.q  <- vmat$output
colnames(vmat.prob_mobius.q) <- colnames(prob_mobius.mat2)

colnames(vmat.prob_mobius.q)[1] <- ""
colnames(vmat.prob_mobius.q)[2] <- ""
save(vmat.prob_mobius.q, file = "data/extreme_points_seven_node_graph.Rdata")
