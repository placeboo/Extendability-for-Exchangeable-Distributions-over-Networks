rm(list = ls())
#------------------------------------------------------------------------------#
#   Title: Create binary probability space, and Mobius parameters based on three nodes.
#   Date: 10/18/2017
#   Author: Jiaqi Yin
#   Purpose: Mobius parameterized
#------------------------------------------------------------------------------#

N = 3 # number of nodes
NN = 1: N # nodes
comb = combn(N, 2) # possible combination
L = ncol(comb) # number of lines
line.name = paste(comb[1,],'-', comb[2,], sep = "")

prob_space.ls = list()

#------------------------------------------------------------------------------#
# 3k1, k3
#------------------------------------------------------------------------------#
prob_space.ls[[1]] = matrix(rep(c(0,1), each = L), 2, byrow = T)
colnames(prob_space.ls[[1]]) = line.name
rownames(prob_space.ls[[1]]) = c("3k1", "k3")

#------------------------------------------------------------------------------#
# p3.bar, p3
#------------------------------------------------------------------------------#
mat = p3(NN)
prob_space.ls[[2]] = rbind(1 - fillprob(NN, mat), fillprob(NN, mat))
rownames(prob_space.ls[[2]]) = rep(c("p3.bar", "p3"), each = nrow(prob_space.ls[[2]])/2)


#------------------------------------------------------------------------------#
# binary probability space
#------------------------------------------------------------------------------#
prob_space.mat = do.call(rbind, prob_space.ls)

# save
save(prob_space.mat, file = "data/three_nodes_graphs.Rdata")
write.csv(prob_space.mat, file = "data/three_nodes_graphs.csv")

# delete the repeated rows in prob_space.mat
prob_space.mat2 = prob_space.mat[!duplicated(rownames(prob_space.mat)),]
write.csv(prob_space.mat2, file = 'data/three_nodes_unique_graphs.csv')
#------------------------------------------------------------------------------#
# Mobius parameterized
#------------------------------------------------------------------------------#
prob_mobius.mat = prob_mobius(prob_space.mat)
# save
save(prob_mobius.mat, file = 'data/three_nodes_mobius.Rdata')
write.csv(prob_mobius.mat, file = 'data/three_nodes_mobius.csv')

#------------------------------------------------------------------------------#
# Make hyperplan matrix
#------------------------------------------------------------------------------#
sum_to_one = -1 * apply(prob_mobius.mat, 2, sum)
prob_mobius.mat2 = rbind(prob_mobius.mat, sum_to_one)
n.rows = dim(prob_mobius.mat2)[1]

prob_mobius.mat2 = cbind(c(rep(0, n.rows - 1), 1), c(rep(0, n.rows - 1), 1), prob_mobius.mat2)
colnames(prob_mobius.mat2)[c(1,2)] = c("in/eq", "const")
prob_mobius.mat2.q = d2q(prob_mobius.mat2)

# Convert from H-rep to V-rep
vmat <- scdd(prob_mobius.mat2.q, representation="H")
vmat.prob_mobius.q  <- vmat$output
colnames(vmat.prob_mobius.q) <- colnames(prob_mobius.mat2)
rownames(vmat.prob_mobius.q) = m2g(vmat.prob_mobius.q[, -c(1,2)])

colnames(vmat.prob_mobius.q)[1] <- ""
colnames(vmat.prob_mobius.q)[2] <- ""

write.csv(vmat.prob_mobius.q, "data/extreme_points_three_node_graph.csv")
save(vmat.prob_mobius.q, file = "data/extreme_points_three_node_graph.Rdata")

Hmat = scdd(vmat.prob_mobius.q)$output
# number of ineq
length(which(Hmat[,1] == 0))
