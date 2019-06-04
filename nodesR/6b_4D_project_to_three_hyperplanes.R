cord4.mat = diag(1, 4, 4)
colnames(cord4.mat) = c("3k1", "p3.bar", "p3", "k3")
rownames(cord4.mat) = c("x", "y", "z", "t")

pnt5to3.4mat = cord4.mat %*% Lambdas5to3
pnt6to3.4mat = cord4.mat %*% Lambdas6to3

# how those points symmetric
# 5-nodes

pair.index = cbind(c(1, 3, 5, 7, 8), c(9, 4, 6, 7, 2))
left.mat = pnt5to3.4mat[, pair.index[,1]]
right.mat = pnt5to3.4mat[, pair.index[,2]]

# middle
midde.mat = 1/2 * (left.mat + right.mat)
v.mat = cbind(rep(0, 5), rep(0,5), t(midde.mat))
v.mat = d2q(v.mat)
redundant(v.mat, representation = "V")


#----------------------------------------------------------------#
# find the hyperplan of each projectation, using the joint probability
# M = (0, b, -A)   Ax \leq b
# M = (1, b, -A)   Ax = b
#----------------------------------------------------------------#
# 5->3 
five_to_three.4vmat = makeV(t(pnt5to3.4mat))
five_to_three.4hmat = scdd(five_to_three.4vmat)$out
five_to_threeIneq.4hmat = five_to_three.4hmat[-which(five_to_three.4hmat[,1] == 1), ]

five_to_three_check.4mat = round(-five_to_threeIneq.4hmat[, -c(1,2)] %*% pnt5to3.4mat - five_to_threeIneq.4hmat[, 2],4)
five_to_three_hyper4mat = matrix(0, nrow = nrow(five_to_three_check.4mat), ncol = ncol(five_to_three_check.4mat))
for(i in 1: nrow(five_to_three_hyper4mat)){
        five_to_three_hyper4mat[i,which(five_to_three_check.4mat[i,]==0)] = 1 
}

# colname: graphs's name and its location 
five_to_three_colname4 = apply(pnt5to3.4mat, 2, function(x) paste(x[1], x[2], x[3],sep = ","))
five_to_three_colname4 = rbind(colnames(five_to_three_check.4mat), paste("(", five_to_three_colname4, ")", sep = ""))
five_to_three_hyper4mat = rbind(five_to_three_colname4, five_to_three_hyper4mat)

# 6->3
six_to_three.4vmat = makeV(t(pnt6to3.4mat))
six_to_three.4hmat = scdd(six_to_three.4vmat)$out
six_to_threeIneq.4hmat = six_to_three.4hmat[-which(six_to_three.4hmat[,1] == 1), ]
round(six_to_threeIneq.4hmat, 3)
six_to_three_check.4mat = round(-six_to_threeIneq.4hmat[, -c(1,2)] %*% pnt6to3.4mat - six_to_threeIneq.4hmat[, 2],4)
six_to_three_hyper4mat = matrix(0, nrow = nrow(six_to_three_check.4mat), ncol = ncol(six_to_three_check.4mat))
for(i in 1: nrow(six_to_three_hyper4mat)){
        six_to_three_hyper4mat[i,which(six_to_three_check.4mat[i,]==0)] = 1 
}

# colname: graphs's name and its location 
six_to_three_colname4 = apply(pnt6to3.4mat, 2, function(x) paste(x[1], x[2], x[3],sep = ","))
six_to_three_colname4 = rbind(colnames(six_to_three_check.4mat), paste("(", six_to_three_colname4, ")", sep = ""))
six_to_three_hyper4mat = rbind(six_to_three_colname4, six_to_three_hyper4mat)

