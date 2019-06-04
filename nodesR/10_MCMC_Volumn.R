#------------------------------------------------------------------------------#
# calculate the volunm of four, five, six, seven node-netwrok projected to three-node network
#------------------------------------------------------------------------------#
rm(list = ls())
set.seed(100)

#------------------------------------------------------------------------------#
# five to three projection, volumn
#------------------------------------------------------------------------------#
# import the H-mat (describe the surface of the projection shape)
load(file = "data/five_to_three_hmat.Rdata")
five_to_threeIneq.hmat = five_to_three.hmat[which(five_to_three.hmat[,1] == 0), ]

# location of the sample
n = 5000
rnd = 1000
five_to_three_prop_vec = rep(0, rnd)
for(iter in 1: rnd){
        sample = rdirichlet(n, c(1,1,1,1))[, -1]
        # check each sample inside or outside of 5->3 polytopes. All negtive, inside; otherwise, outside
        sample_loc = -five_to_threeIneq.hmat[, -c(1,2)] %*% t(sample) - five_to_threeIneq.hmat[, 2]
        inside_prop = sum(apply(sample_loc, 2, function(x) all(x < 0)))/n
        print(inside_prop)
        five_to_three_prop_vec[iter] = inside_prop
        iter = iter + 1
}
mean(five_to_three_prop_vec) # 0.9719494 

#------------------------------------------------------------------------------#
# six to three projection, volumn
#------------------------------------------------------------------------------#
load(file = "data/six_to_three_hmat.Rdata")
six_to_threeIneq.hmat = six_to_three.hmat[which(six_to_three.hmat[,1] == 0), ]
# location of the sample
n = 5000
rnd = 1000
six_to_three_prop_vec = rep(0, rnd)
for(iter in 1: rnd){
        sample = rdirichlet(n, c(1,1,1,1))[, -1]
        # check each sample inside or outside of 5->3 polytopes. All negtive, inside; otherwise, outside
        sample_loc = -six_to_threeIneq.hmat[, -c(1,2)] %*% t(sample) - six_to_threeIneq.hmat[, 2]
        inside_prop = sum(apply(sample_loc, 2, function(x) all(x < 0)))/n
        print(inside_prop)
        six_to_three_prop_vec[iter] = inside_prop
        iter = iter + 1
}
mean(six_to_three_prop_vec) # 0.9539292

#------------------------------------------------------------------------------#
# seven to three projection, volumn
#------------------------------------------------------------------------------#
load(file = "data/seven_to_three_hmat.Rdata")
seven_to_threeIneq.hmat = seven_to_three.hmat[which(seven_to_three.hmat[,1] == 0), ]
# location of the sample
n = 5000
rnd = 1000
seven_to_three_prop_vec = rep(0, rnd)
for(iter in 1: rnd){
        sample = rdirichlet(n, c(1,1,1,1))[, -1]
        # check each sample inside or outside of 5->3 polytopes. All negtive, inside; otherwise, outside
        sample_loc = -seven_to_threeIneq.hmat[, -c(1,2)] %*% t(sample) - seven_to_threeIneq.hmat[, 2]
        inside_prop = sum(apply(sample_loc, 2, function(x) all(x < 0)))/n
        print(inside_prop)
        seven_to_three_prop_vec[iter] = inside_prop
        iter = iter + 1
}
mean(seven_to_three_prop_vec) # 0.9285504

#------------------------------------------------------------------------------#
# five to four projection, volumn
#------------------------------------------------------------------------------#
four_extreme = get(load("data/extreme_points_four_node_graph.Rdata"))
five_to_four = get(load("data/five_nodes_to_four_nodes_projection.Rdata"))
# re-orgnized. The colnumn names should be the same
four_extreme1 = four_extreme[, -c(1,2)]
five_to_four1 = five_to_four[, -c(1,2)]
five_to_four1 = five_to_four1[, colnames(four_extreme1)]

four_extreme.d = q2d(four_extreme1)
five_to_four.d = q2d(five_to_four1)

Lambdas5to4 = solve(t(four_extreme.d)) %*% t(five_to_four.d)

cord.mat = diag(1, 11)
#cord.mat = cord.mat[,-1]
#cord.mat = t(cord.mat)
colnames(cord.mat) = rownames(Lambdas5to4)

# project Lambda to new cordination
pnt5to4.mat = cord.mat %*% Lambdas5to4

# find the hmat
five_to_four.vmat = makeV(t(pnt5to4.mat))
five_to_four.hmat = scdd(five_to_four.vmat)$out
colnames(five_to_four.hmat)=  c("","", rownames(Lambdas5to4))
save(five_to_four.hmat, file = "data/five_to_four_hmat.Rdata")

five_to_fourIneq.hmat = five_to_four.hmat[which(five_to_four.hmat[,1] == 0), ]
# location of the sample
n = 5000
rnd = 1000
five_to_four_prop_vec = rep(0, rnd)
for(iter in 1: rnd){
        #sample = rdirichlet(n, rep(1,11))[, -1]
        sample = rdirichlet(n, rep(1,11))
        # check each sample inside or outside of 5->4 polytopes. All negtive, inside; otherwise, outside
        sample_loc = -five_to_fourIneq.hmat[, -c(1,2)] %*% t(sample) - five_to_fourIneq.hmat[, 2]
        inside_prop = sum(apply(sample_loc, 2, function(x) all(x < 0)))/n
        #print(inside_prop)
        five_to_four_prop_vec[iter] = inside_prop
        iter = iter + 1
}
mean(five_to_four_prop_vec) # 0.6073876

# volume of n-simple
vol_n_simplex= function(n){
        return(sqrt(n+1)/(prod(1:n) * sqrt(2^n)))
}
vol_n_simplex(4) # 0.1178511
vol_n_simplex(11) # 1.917653e-09
