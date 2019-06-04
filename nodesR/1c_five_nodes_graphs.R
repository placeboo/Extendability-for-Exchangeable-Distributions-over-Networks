rm(list = ls())

#------------------------------------------------------------------------------#
#   Title: Create binary probability space, and Mobius parameters based on five nodes.
#   Date: 10/18/2017
#   Author: Jiaqi Yin
#   Purpose: Mobius parameterized
#------------------------------------------------------------------------------#

N = 5 # number of nodes
NN = 1: N # nodes
comb = combn(N, 2) # possible combination
L = ncol(comb) # number of lines
line.name = paste(comb[1,],'-', comb[2,], sep = "")

prob_space.ls = list()

#------------------------------------------------------------------------------#
# 5k1, k5
#------------------------------------------------------------------------------#
prob_space.ls[[1]] = matrix(rep(c(0,1), each = L), 2, byrow = T)
colnames(prob_space.ls[[1]]) = line.name
rownames(prob_space.ls[[1]]) = c("5k1", "k5")

#------------------------------------------------------------------------------#
# "k2_3k1", "k2_3k1.bar"
#------------------------------------------------------------------------------#
k = 2
prob_space.ls[[k]] = rbind(diag(1, L), 1 - diag(1, L))
colnames(prob_space.ls[[k]]) = line.name
rownames(prob_space.ls[[k]]) =  rep(c("k2_3k1", "k2_3k1.bar"), each = nrow(prob_space.ls[[k]])/2)

#------------------------------------------------------------------------------#
# "p3_2k1", "p3_2k1.bar"
#------------------------------------------------------------------------------#
k = 3
mat = p3(NN)
prob_space.ls[[k]] =rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) =  rep(c("p3_2k1", "p3_2k1.bar"), each = nrow(prob_space.ls[[k]])/2)

#------------------------------------------------------------------------------#
# "w4.bar", "w4"
#------------------------------------------------------------------------------#
k = 4
mat = twoK2(NN)
prob_space.ls[[k]] =rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) =  rep(c("w4.bar", "w4"), each = nrow(prob_space.ls[[k]])/2)

#------------------------------------------------------------------------------#
# "claw_k1", "claw_k1.bar"
#------------------------------------------------------------------------------#
k = 5
mat = claw(NN)
prob_space.ls[[k]] =rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) =  rep(c("claw_k1", "claw_k1.bar"), each = nrow(prob_space.ls[[k]])/2)

#------------------------------------------------------------------------------#
# "p2_p3", "p2_p3.bar"
#------------------------------------------------------------------------------#
k = 6
mat = p2_p3(NN)
prob_space.ls[[k]] =rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) =  rep(c("p2_p3", "p2_p3.bar"), each = nrow(prob_space.ls[[k]])/2)

#------------------------------------------------------------------------------#
# "co-gem", "gem"
#------------------------------------------------------------------------------#
k = 7
mat = gem(NN)
prob_space.ls[[k]] =rbind(1 - fillprob(NN, mat), fillprob(NN, mat))
rownames(prob_space.ls[[k]]) =  rep(c("co-gem", "gem"), each = nrow(prob_space.ls[[k]])/2)

#------------------------------------------------------------------------------#
# "k3_2k1", "k3_2k1.bar"
#------------------------------------------------------------------------------#
k = 8
mat = k3_2k1.bar(NN)
prob_space.ls[[k]] =rbind(1 - fillprob(NN, mat), fillprob(NN, mat))
rownames(prob_space.ls[[k]]) =  rep(c("k3_2k1", "k3_2k1.bar"), each = nrow(prob_space.ls[[k]])/2)

#------------------------------------------------------------------------------#
#  "k14", "k14.bar"
#------------------------------------------------------------------------------#
k = 9
mat = k14(NN)
prob_space.ls[[k]] =rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) =  rep(c("k14", "k14.bar"), each = nrow(prob_space.ls[[k]])/2)

#------------------------------------------------------------------------------#
# "co-butterfly", "butterfly"
#------------------------------------------------------------------------------#
k = 10
mat = butterfly(NN)
prob_space.ls[[k]] =rbind(1 - fillprob(NN, mat), fillprob(NN, mat))
rownames(prob_space.ls[[k]]) =  rep(c("co-butterfly", "butterfly"), each = nrow(prob_space.ls[[k]])/2)

#------------------------------------------------------------------------------#
# "fork", "co-fork"
#------------------------------------------------------------------------------#
k = 11
mat = chair(NN)
prob_space.ls[[k]] =rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) =  rep(c("fork", "co-fork"), each = nrow(prob_space.ls[[k]])/2)

#------------------------------------------------------------------------------#
# "co-dart", "dart"
#------------------------------------------------------------------------------#
k = 12
mat = co.dart(NN)
prob_space.ls[[k]] =rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) =  rep(c("co-dart", "dart"), each = nrow(prob_space.ls[[k]])/2)
#------------------------------------------------------------------------------#
# "p5", "p5.bar"
#------------------------------------------------------------------------------#
k = 13
mat = p5(NN)
prob_space.ls[[k]] =rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) =  rep(c("p5", "p5.bar"), each = nrow(prob_space.ls[[k]])/2)

#------------------------------------------------------------------------------#
# "k23.bar", "k23"
#------------------------------------------------------------------------------#
k = 14
mat = k2_k3(NN)
prob_space.ls[[k]] =rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) =  rep(c("k23.bar", "k23"), each = nrow(prob_space.ls[[k]])/2)

#------------------------------------------------------------------------------#
# "P", "P.bar"
#------------------------------------------------------------------------------#
k = 15
mat = P(NN)
prob_space.ls[[k]] =rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) =  rep(c("P", "P.bar"), each = nrow(prob_space.ls[[k]])/2)

#------------------------------------------------------------------------------#
# "bull"
#------------------------------------------------------------------------------#
k = 16
mat = bull(NN)
prob_space.ls[[k]] = fillprob(NN, mat)
rownames(prob_space.ls[[k]]) =  rep(c("bull"), each = nrow(prob_space.ls[[k]]))

#------------------------------------------------------------------------------#
# "cricket", 'co-cricket'
#------------------------------------------------------------------------------#
k = 17
mat = cricket(NN)
prob_space.ls[[k]] =rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) =  rep(c("cricket", 'co-cricket'), each = nrow(prob_space.ls[[k]])/2)

#------------------------------------------------------------------------------#
# 'c5'
#------------------------------------------------------------------------------#
k = 18
mat = c5(NN)
prob_space.ls[[k]] = fillprob(NN, mat)
rownames(prob_space.ls[[k]]) =  rep(c("c5"), each = nrow(prob_space.ls[[k]]))

#------------------------------------------------------------------------------#
# binary probability space
#------------------------------------------------------------------------------#
prob_space.mat = do.call(rbind, prob_space.ls)

# save
save(prob_space.mat, file = "data/five_nodes_graphs.Rdata")
write.csv(prob_space.mat, file = "data/five_nodes_graphs.csv")

# delete the repeated rows in prob_space.mat
prob_space.mat2 = prob_space.mat[!duplicated(rownames(prob_space.mat)),]
write.csv(prob_space.mat2, file = 'data/five_nodes_unique_graphs.csv')
#------------------------------------------------------------------------------#
# Mobius parameterized
#------------------------------------------------------------------------------#
prob_mobius.mat = prob_mobius(prob_space.mat)
# save
save(prob_mobius.mat, file = 'data/five_nodes_mobius.Rdata')
write.csv(prob_mobius.mat, file = 'data/five_nodes_mobius.csv')

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

colnames(vmat.prob_mobius.q)[1] <- ""
colnames(vmat.prob_mobius.q)[2] <- ""

rownames(vmat.prob_mobius.q) = m2g(vmat.prob_mobius.q[, -c(1,2)])

write.csv(vmat.prob_mobius.q, "data/extreme_points_five_node_graph.csv")
save(vmat.prob_mobius.q, file = "data/extreme_points_five_node_graph.Rdata")
