library(rcdd)
# run in server
rm(list = ls())
h.mat = matrix(NA, 2, 3)
rownames(h.mat) = c(5,6)
colnames(h.mat) = c(4:6)


six_to_four.vmat = get(load(file = 'six_nodes_to_four_nodes_projection.Rdata'))
six_to_five.vmat = get(load( file = 'six_nodes_to_five_nodes_projection.Rdata'))
six.vmat = get(load(file = "six_nodes_mobius.Rdata"))
five.vmat = get(load(file = "five_nodes_mobius.Rdata"))

# v->H
six_to_four.hmat = scdd(six_to_four.vmat)$output
save(six_to_four.hmat, file = "six_to_four_hmat.Rdata")
h.mat[2, 1] = length(which(six_to_four.hmat[, 1] == 0))


six_to_five.hmat = scdd(six_to_five.vmat)$output
h.mat[2, 2] = length(which(six_to_five.hmat[, 1] == 0))
save(six_to_five.hmat, file = "six_to_five_hmat.Rdata")

six.hmat = scdd(six.vmat)$output
h.mat[2, 3] = length(which(six.hmat[, 1] == 0))
save(six.hmat, file = "six_hmat.Rdata")

five.mat = scdd(five.mat)$output
h.mat[1, 2] = length(which(five.hmat[, 1] == 0))
save(five.hmat, file = "five_hmat.Rdata")

save(h.mat, file = "number_equ_large.Rdata")
