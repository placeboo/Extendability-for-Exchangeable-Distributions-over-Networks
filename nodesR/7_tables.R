#----------------------------------------------------------------#
# hyperplan table
#----------------------------------------------------------------#


#----------------------------------------------------------------#
# the vextex of multiple hyperplanes' intersections
#----------------------------------------------------------------#
# 5->3
load(file = "data/hyperplane_vertex_five_to_three.Rdata")
print(xtable(five_to_three_hypermat, caption = "Polytopes of five-node network projecting to three-node network. 1: the corners of hyperplan in Polytopes; otherwise 0.", digits = 0, align = rep("c", ncol(five_to_three_hypermat)+1)), include.colnames = F)

# 6->3
load(file = "data/hyperplane_vertex_six_to_three.Rdata")
print(xtable(six_to_three_hypermat, caption = "Polytopes of six-node network projecting to three-node network. 1: the corners of hyperplan in Polytopes; otherwise 0.", digits = 0, align = rep("c", ncol(six_to_three_hypermat)+1)))

#----------------------------------------------------------------#
# projection, number of inequality
#----------------------------------------------------------------#
hyperplane.mat = matrix(NA, 5, 5)
colnames(hyperplane.mat) = c(2:6)
rownames(hyperplane.mat) = c(2:6)

hyperplane.mat[,1] = 1
hyperplane.mat[2:5,2] = c(4, 4, 8, 7)
hyperplane.mat[3:4,3] = c(11, 472)

print(xtable(hyperplane.mat, caption = "number of inequality", digits = 0, align = rep("c", ncol(hyperplane.mat)+1)))
