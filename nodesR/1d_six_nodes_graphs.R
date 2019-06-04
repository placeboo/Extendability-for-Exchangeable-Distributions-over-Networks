rm(list = ls())
#------------------------------------------------------------------------------#
#   Title: Create binary probability space, and Mobius parameters based on six nodes.
#   Date: 10/22/2017
#   Author: Jiaqi Yin
#   Purpose: Mobius parameterized
#------------------------------------------------------------------------------#


N = 6
NN = 1:N
comb = combn(N, 2) # possible combination
L = ncol(comb) # number of lines
line.name = paste(comb[1,],'-', comb[2,], sep = "")
pair3.mat = combn(NN, 3)
pair4.mat = combn(NN, 4)
pair5.mat = combn(NN, 5)

# build the isomorphism, each entry is the totally number of each isomorphism categary. name is the same as http://graphclasses.org/smallgraphs.html#nodes6
prob_space.ls = list()

################
## k6.bar, k6 ##
################
prob_space.ls[[1]] = matrix(rep(c(0,1), each = L), 2, byrow = T)
colnames(prob_space.ls[[1]]) = line.name
rownames(prob_space.ls[[1]]) = c("k6.bar", "k6")

############################################################
######################## one edges  ########################
############################################################

############################
## "k2_4k1", "k2_4k1.bar" ##
############################
prob_space.ls[[2]] = rbind(diag(1, L), 1 - diag(1, L))
rownames(prob_space.ls[[2]]) = rep(c("k2_4k1", "k2_4k1.bar"), each = L)
colnames(prob_space.ls[[2]]) = line.name

############################################################
######################## two edges  ########################
############################################################

############################
## "p3_3k1", "p3_3k1.bar" ##
############################
prob_space.ls[[3]] = rbind(fillprob(NN, p3(NN)), 1 - fillprob(NN, p3(NN)))
rownames(prob_space.ls[[3]]) = rep(c("p3_3k1", "p3_3k1.bar"), each = nrow(prob_space.ls[[3]])/2)
colnames(prob_space.ls[[3]]) = line.name

###################################
## ""twoK2_2k1", "twoK2_2k1.bar" ##
###################################
prob_space.ls[[4]] = rbind(fillprob(NN, twoK2(NN)), 1 - fillprob(NN, twoK2(NN)))
rownames(prob_space.ls[[4]]) = rep(c("twoK2_2k1", "twoK2_2k1.bar"), each = nrow(prob_space.ls[[4]])/2)
colnames(prob_space.ls[[4]]) = line.name

############################################################
######################## three edges  ######################
############################################################

###########################
## "p4_2k1", "p4_2k1.bar"##
###########################
temp.mat = matrix(NA, ncol = L)
for(i in 1: ncol(pair4.mat)){
      temp.mat = rbind(temp.mat, rbind(fillprob(NN, p4(pair4.mat[, i]))))
}
temp.mat = temp.mat[-1, ]
prob_space.ls[[5]] = rbind(temp.mat, 1 - temp.mat)
rownames(prob_space.ls[[5]]) = rep(c("p4_2k1", "p4_2k1.bar"), each = nrow(prob_space.ls[[5]])/2)
colnames(prob_space.ls[[5]]) = line.name

#############################
## "threeK2", "threeK2.bar"##
#############################
prob_space.ls[[6]] = rbind(fillprob(NN, threeK2(NN)), 1 - fillprob(NN, threeK2(NN)))
rownames(prob_space.ls[[6]]) = rep(c("threeK2", "threeK2.bar"), each = nrow(prob_space.ls[[6]])/2)
colnames(prob_space.ls[[6]]) = line.name

#############################
## "claw_2k1", "claw_2k1.bar"##
#############################
prob_space.ls[[7]] = rbind(fillprob(NN, claw(NN)), 1 - fillprob(NN, claw(NN)))
rownames(prob_space.ls[[7]]) = rep(c("claw_2k1", "claw_2k1.bar"), each = nrow(prob_space.ls[[7]])/2)
colnames(prob_space.ls[[7]]) = line.name

#######################
## "x197", "x197.bar"##
#######################
# p3 combined with p3.bar
x197.mat = matrix(NA, ncol = 3)
for(i in 1: ncol(pair3.mat)){
      pair3.vec = sort(pair3.mat[, i])
      rst = sort(NN[!NN%in%pair3.vec])
      p3.mat = p3(pair3.vec)
      p3.mat.bar = probToEdges(1 - fillprob(rst, p3(rst)))
      x197.temp = cbind(matrix(rep(p3.mat, each = 3), ncol = ncol(p3.mat), byrow = F), rep(p3.mat.bar, 3))
      x197.mat = rbind(x197.mat, x197.temp)
}
prob_space.ls[[8]] = rbind(fillprob(NN, x197.mat[-1, ]), 1 - fillprob(NN, x197.mat[-1, ]))
rownames(prob_space.ls[[8]]) = rep(c("x197", "x197.bar"), each = nrow(prob_space.ls[[8]])/2)
colnames(prob_space.ls[[8]]) = line.name

############################
## "k3_3k1", "k3_3k1.bar" ##
############################
prob_space.ls[[9]] = rbind(fillprob(NN, k3(NN)), 1 - fillprob(NN, k3(NN)))
rownames(prob_space.ls[[9]]) = rep(c("k3_3k1", "k3_3k1.bar"), each = nrow(prob_space.ls[[9]])/2)
colnames(prob_space.ls[[9]]) = line.name


############################################################
######################## four edges  #######################
############################################################

##########################
## "p2_p4", "p2_p4.bar" ##
##########################
# modify from "p4_2k1"
pair4.mat = combn(NN, 4)
temp.mat = matrix(NA , ncol = L)
for(i in 1: ncol(pair4.mat)){
      p4.mat = p4(pair4.mat[, i])
      rst = sort(NN[! NN %in% pair4.mat[, i]])
      p2.vec = rep(paste(rst[1], '-', rst[2], sep = ""), nrow(p4.mat))

      p4_p2.mat = cbind(p4.mat, p2.vec)
      temp.mat = rbind(temp.mat, rbind(fillprob(NN, p4_p2.mat)))
}
temp.mat = temp.mat[-1, ]
prob_space.ls[[10]] = rbind(temp.mat, 1 - temp.mat)
rownames(prob_space.ls[[10]]) = rep(c("p2_p4", "p2_p4.bar"), each = nrow(prob_space.ls[[10]])/2)
colnames(prob_space.ls[[10]]) = line.name

###########################
## ""twoP3", "twoP3.bar" ##
###########################
mat = twoP3(NN)
prob_space.ls[[11]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[11]]) = rep(c("twoP3", "twoP3.bar"), each = nrow(prob_space.ls[[11]])/2)
colnames(prob_space.ls[[11]]) = line.name

############################
## "c4_2k1", "c4_2k1.bar" ##
############################
c4_2k1.mat = matrix(NA, ncol = 4)
for(i in 1: ncol(pair4.mat)){
      temp.mat = c4(pair4.mat[, i])
      c4_2k1.mat = rbind(c4_2k1.mat, temp.mat)
}
c4_2k1.mat = c4_2k1.mat[-1, ]

temp.mat = fillprob(NN, c4_2k1.mat)
prob_space.ls[[12]] = rbind(temp.mat, 1 - temp.mat)
rownames(prob_space.ls[[12]]) = rep(c("c4_2k1", "c4_2k1.bar" ), each = nrow(prob_space.ls[[12]])/2)
colnames(prob_space.ls[[12]]) = line.name

##############################
## "k2_claw", "k2_claw.bar" ##
##############################
rst.mat = apply(pair4.mat, 2, function(y) NN[! NN %in% y])
k2_claw.mat = matrix(NA, ncol = 4)
for(i in 1: ncol(pair4.mat)){
      claw.temp.mat = claw(pair4.mat[, i])
      k2_claw.temp.mat = cbind(claw.temp.mat, rep(paste(rst.mat[1, i], "-", rst.mat[2, i], sep = ""), nrow(claw.temp.mat)))
      k2_claw.mat = rbind(k2_claw.mat, k2_claw.temp.mat)
}
k2_claw.mat = k2_claw.mat[-1, ]

temp.mat = fillprob(NN, k2_claw.mat)
prob_space.ls[[13]] = rbind(temp.mat, 1 - temp.mat)
rownames(prob_space.ls[[13]]) = rep(c("k2_claw", "k2_claw.bar"), each = nrow(prob_space.ls[[13]])/2)
colnames(prob_space.ls[[13]]) = line.name

############################
## "k14_k1", "k14_k1.bar" ##
############################
pair5.mat = combn(NN, 5)
k14_k1.mat = matrix(NA, ncol = 4)
for(i in 1: ncol(pair5.mat)){
      k14_k1.mat = rbind(k14_k1.mat, k14(pair5.mat[, i]))
}
k14_k1.mat = k14_k1.mat[-1, ]

temp.mat = fillprob(NN, k14_k1.mat)
prob_space.ls[[14]] = rbind(temp.mat, 1 - temp.mat)
rownames(prob_space.ls[[14]]) = rep(c("k14_k1", "k14_k1.bar"), each = nrow(prob_space.ls[[14]])/2)
colnames(prob_space.ls[[14]]) = line.name

################################
## "chair_k1", "chair_k1.bar" ##
################################
chair_k1.mat = matrix(NA, ncol = 4)
for(i in 1: ncol(pair5.mat)){
      chair_k1.mat = rbind(chair_k1.mat, chair(pair5.mat[, i]))
}
chair_k1.mat = chair_k1.mat[-1, ]

temp.mat = fillprob(NN, chair_k1.mat)
prob_space.ls[[15]] = rbind(temp.mat, 1 - temp.mat)
rownames(prob_space.ls[[15]]) = rep(c("chair_k1", "chair_k1.bar"), each = nrow(prob_space.ls[[15]])/2)
colnames(prob_space.ls[[15]]) = line.name

######################################
## "co-dart_k1", "co-dart_k1.bar" ##
######################################
co.dart_k1.mat = matrix(NA,ncol = 4)
for(i in 1: ncol(pair4.mat)){
      co.dart_k1.mat = rbind(co.dart_k1.mat, paw(pair4.mat[, i]))
}
co.dart_k1.mat = co.dart_k1.mat[-1, ]

temp.mat = fillprob(NN, co.dart_k1.mat)
prob_space.ls[[16]] = rbind(temp.mat, 1 - temp.mat)
rownames(prob_space.ls[[16]]) = rep(c("co-dart_k1", "co-dart_k1.bar"), each = nrow(prob_space.ls[[16]])/2)
colnames(prob_space.ls[[16]]) = line.name

##########################
## "p5_k1", "p5_k1.bar" ##
##########################
p5_k1.mat = matrix(NA, ncol = 4)
for(i in 1: ncol(pair5.mat)){
      p5_k1.mat = rbind(p5_k1.mat, p5(pair5.mat[, i]))
}
p5_k1.mat = p5_k1.mat[-1, ]

temp.mat = fillprob(NN, p5_k1.mat)
prob_space.ls[[17]] = rbind(temp.mat, 1 - temp.mat)
rownames(prob_space.ls[[17]]) = rep(c("p5_k1", "p5_k1.bar"), each = nrow(prob_space.ls[[17]])/2)
colnames(prob_space.ls[[17]]) = line.name

###############################
## "k1_k2_k3", "k1_k2_k3.bar"##
###############################
# first, build k3
k3.mat = matrix(NA,ncol = 3)
k2.mat = matrix(NA, ncol = 1)
rst.mat = apply(pair3.mat, 2, function(y) NN[! NN %in% y])
for(i in 1: (ncol(pair3.mat))){
      k3.mat = rbind(k3.mat, k3(pair3.mat[,i]))
      k2.mat = rbind(k2.mat, k2(rst.mat[,i]))
}

k3.mat = k3.mat[-1, ]
k2.mat = k2.mat[-1, ]

# rep k3.mat for matching k2
k3.mat = matrix(rep(k3.mat, each = 3), ncol = 3, byrow = F)
k1_k2_k3.mat = cbind(k3.mat, k2.mat)

temp.mat = fillprob(NN, k1_k2_k3.mat)
prob_space.ls[[18]] = rbind(temp.mat, 1 - temp.mat)
rownames(prob_space.ls[[18]]) = rep(c("k1_k2_k3", "k1_k2_k3.bar"), each = nrow(prob_space.ls[[18]])/2)
colnames(prob_space.ls[[18]]) = line.name


############################################################
######################## five edges  #######################
############################################################

########################
## "cross", "co-cross"##
########################
prob_space.ls[[19]] = rbind(fillprob(NN, cross(NN)), 1 - fillprob(NN, cross(NN)))
rownames(prob_space.ls[[19]]) = rep(c("cross", "co-cross"), each = nrow(prob_space.ls[[19]])/2)
colnames(prob_space.ls[[19]]) = line.name

#################
## "H", "H.bar"##
#################
prob_space.ls[[20]] = rbind(fillprob(NN, H(NN)), 1 - fillprob(NN, H(NN)))
rownames(prob_space.ls[[20]]) = rep(c("H", "H.bar"), each = nrow(prob_space.ls[[20]])/2)
colnames(prob_space.ls[[20]]) = line.name

#########################
## "c4_p2", "c4_p2.bar"##
#########################
c4_p2.mat = matrix(NA, ncol = 5)
for(i in 1: ncol(pair4.mat)){
      pair4.vec = pair4.mat[, i]
      rst.vec = sort(NN[!NN %in% pair4.vec])
      c4_p2.mat = rbind(c4_p2.mat, cbind(c4(pair4.vec), p2.temp = paste(rst.vec[1], '-', rst.vec[2],sep = '')))
}
c4_p2.mat = c4_p2.mat[-1,]
prob_space.ls[[21]] = rbind(fillprob(NN, c4_p2.mat), 1 - fillprob(NN, c4_p2.mat))
rownames(prob_space.ls[[21]]) = rep(c("c4_p2", "c4_p2.bar"), each = nrow(prob_space.ls[[21]])/2)
colnames(prob_space.ls[[21]]) = line.name

#################
## "E", "E.bar"##
#################
prob_space.ls[[22]] = rbind(fillprob(NN, E(NN)), 1 - fillprob(NN, E(NN)))
rownames(prob_space.ls[[22]]) = rep(c("E", "E.bar"), each = nrow(prob_space.ls[[22]])/2)
colnames(prob_space.ls[[22]]) = line.name

#########################
## "k3_p3", "k3_p3.bar"##
#########################
k3_p3.mat = matrix(NA,ncol = 5)
for(i in 1:ncol(pair3.mat)){
      pair3.vec = pair3.mat[,i]
      rst.vec = NN[!NN %in% pair3.vec]
      k3.mat = matrix(rep(k3(pair3.vec), 3), byrow = T, ncol = 3)
      k3_p3.temp = cbind(k3.mat, p3(rst.vec))
      k3_p3.mat = rbind(k3_p3.mat, k3_p3.temp)
}
k3_p3.mat = k3_p3.mat[-1, ]

prob_space.ls[[23]] = rbind(fillprob(NN, k3_p3.mat), 1 - fillprob(NN, k3_p3.mat))
rownames(prob_space.ls[[23]]) = rep(c("k3_p3", "k3_p3.bar"), each = nrow(prob_space.ls[[23]])/2)
colnames(prob_space.ls[[23]]) = line.name

#######################
## "x198", "x198.bar"##
#######################
x198.mat = matrix(NA, ncol = 5)
for(i in 1:ncol(pair5.mat)){
      pair5.vec = pair5.mat[, i]
      x198.mat = rbind(x198.mat, P.bar(pair5.vec))
}
x198.mat = x198.mat[-1, ]

prob_space.ls[[24]] = rbind(fillprob(NN, x198.mat), 1 - fillprob(NN, x198.mat))
rownames(prob_space.ls[[24]]) = rep(c("x198", "x198.bar"), each = nrow(prob_space.ls[[24]])/2)
colnames(prob_space.ls[[24]]) = line.name

###################
## "p6", "p6.bar"##
###################
prob_space.ls[[25]] = rbind(fillprob(NN, p6(NN)), 1 - fillprob(NN, p6(NN)))
rownames(prob_space.ls[[25]]) = rep(c("p6", "p6.bar"), each = nrow(prob_space.ls[[25]])/2)
colnames(prob_space.ls[[25]]) = line.name

###################
## "w5.bar", "w5"##
###################
w5_bar.mat = matrix(NA, ncol = 5)
for(i in 1:ncol(pair5.mat)){
      pair5.vec = pair5.mat[, i]
      w5_bar.mat = rbind(w5_bar.mat, c5(pair5.vec))
}
w5_bar.mat  = w5_bar.mat[-1, ]
prob_space.ls[[26]] = rbind(fillprob(NN, w5_bar.mat), 1 - fillprob(NN, w5_bar.mat))
rownames(prob_space.ls[[26]]) = rep(c("w5.bar", "w5"), each = nrow(prob_space.ls[[26]])/2)
colnames(prob_space.ls[[26]]) = line.name

#######################
## "x172", "x172.bar"##
#######################
prob_space.ls[[27]] = rbind(fillprob(NN, x172(NN)), 1 - fillprob(NN, x172(NN)))
rownames(prob_space.ls[[27]]) = rep(c("x172", "x172.bar"), each = nrow(prob_space.ls[[27]])/2)
colnames(prob_space.ls[[27]]) = line.name

#############################
## "bull_k1", "bull_k1.bar"##
#############################
bull_k1.mat = matrix(NA, ncol = 5)
for(i in 1:ncol(pair5.mat)){
      pair5.vec = pair5.mat[, i]
      bull_k1.mat = rbind(bull_k1.mat, bull(pair5.vec))
}
bull_k1.mat = bull_k1.mat[-1, ]
prob_space.ls[[28]] = rbind(fillprob(NN, bull_k1.mat), 1 - fillprob(NN, bull_k1.mat))
rownames(prob_space.ls[[28]]) = rep(c("bull_k1", "bull_k1.bar"), each = nrow(prob_space.ls[[28]])/2)
colnames(prob_space.ls[[28]]) = line.name

#####################################
## "cricket_k1", "cricket_k1.bar"##
#####################################
cricket_k1.mat = matrix(NA, ncol = 5)
for(i in 1:ncol(pair5.mat)){
      pair5.vec = pair5.mat[, i]
      cricket_k1.mat = rbind(cricket_k1.mat, cricket(pair5.vec))
}
cricket_k1.mat = cricket_k1.mat[-1, ]

prob_space.ls[[29]] = rbind(fillprob(NN, cricket_k1.mat), 1 - fillprob(NN, cricket_k1.mat))
rownames(prob_space.ls[[29]]) = rep(c("cricket_k1", "cricket_k1.bar"), each = nrow(prob_space.ls[[29]])/2)
colnames(prob_space.ls[[29]]) = line.name

#######################
## "P_k1", "P_k1.bar"##
#######################
p_k1.mat = matrix(NA, ncol = 5)
for(i in 1:ncol(pair5.mat)){
      pair5.vec = pair5.mat[, i]
      p_k1.mat = rbind(p_k1.mat, P(pair5.vec))
}
p_k1.mat = p_k1.mat[-1, ]

prob_space.ls[[30]] = rbind(fillprob(NN, p_k1.mat), 1 - fillprob(NN, p_k1.mat))
rownames(prob_space.ls[[30]]) = rep(c("P_k1", "P_k1.bar"), each = nrow(prob_space.ls[[30]])/2)
colnames(prob_space.ls[[30]]) = line.name

###########################
## "k5_k1.bar", "k5_k1"  ##
###########################
prob_space.ls[[31]] = rbind(fillprob(NN, k5_k1.bar(NN)), 1 - fillprob(NN, k5_k1.bar(NN)))
rownames(prob_space.ls[[31]]) = rep(c("k5_k1.bar", "k5_k1"), each = nrow(prob_space.ls[[31]])/2)
colnames(prob_space.ls[[31]]) = line.name

####################################
## "diamond_2k1", "diamond_2k1.bar",##
####################################
diamond_2k1.mat  = matrix(NA, ncol = 5)
for(i in 1: ncol(pair4.mat)){
      pair4.vec = pair4.mat[, i]
      diamond_2k1.mat = rbind(diamond_2k1.mat, diamond(pair4.vec))
}
diamond_2k1.mat = diamond_2k1.mat[-1, ]
prob_space.ls[[32]] = rbind(fillprob(NN, diamond_2k1.mat), 1 - fillprob(NN, diamond_2k1.mat))
rownames(prob_space.ls[[32]]) = rep(c("diamond_2k1", "diamond_2k1.bar"), each = nrow(prob_space.ls[[32]])/2)
colnames(prob_space.ls[[32]]) = line.name

############################
## "paw_p2", "paw_p2.bar" ##
############################
paw_p2.mat = matrix(NA, ncol = 5)
for(i in 1: ncol(pair4.mat)){
      pair4.vec = pair4.mat[, i]
      rst.vec = sort(NN[!NN %in% pair4.vec])
      paw_p2.mat = rbind(paw_p2.mat, cbind(paw(pair4.vec), p2.temp = paste(rst.vec[1], '-', rst.vec[2],sep = '')))
}
paw_p2.mat = paw_p2.mat[-1, ]
prob_space.ls[[33]] = rbind(fillprob(NN, paw_p2.mat), 1 - fillprob(NN, paw_p2.mat))
rownames(prob_space.ls[[33]]) = rep(c("paw_p2", "paw_p2.bar"), each = nrow(prob_space.ls[[33]])/2)
colnames(prob_space.ls[[33]]) = line.name

############################################################
######################## six edges  ########################
############################################################

####################################
## "co-fork_k1", "co-fork_k1.bar" ##
####################################
co.fort.mat = matrix(NA,ncol = 6)
for(i in 1:ncol(pair5.mat)){
      pair5.vec = pair5.mat[, i]
      co.fort.mat = rbind(co.fort.mat, co.fork(pair5.vec))
}
co.fort.mat = co.fort.mat[-1, ]
prob_space.ls[[34]] = rbind(fillprob(NN, co.fort.mat), 1 - fillprob(NN, co.fort.mat))
rownames(prob_space.ls[[34]]) = rep(c("co-fork_k1", "co-fork_k1.bar"), each = nrow(prob_space.ls[[34]])/2)
colnames(prob_space.ls[[34]]) = line.name

#######################################
## "butterfly_k1", "butterfly_k1.bar"##
#######################################
butterfly_k1.mat = matrix(NA,ncol = 6)
for(i in 1:ncol(pair5.mat)){
      pair5.vec = pair5.mat[, i]
      butterfly_k1.mat = rbind(butterfly_k1.mat, butterfly(pair5.vec))
}
butterfly_k1.mat = butterfly_k1.mat[-1, ]
prob_space.ls[[35]] = rbind(fillprob(NN, butterfly_k1.mat), 1 - fillprob(NN, butterfly_k1.mat))
rownames(prob_space.ls[[35]]) = rep(c("butterfly_k1", "butterfly_k1.bar"), each = nrow(prob_space.ls[[35]])/2)
colnames(prob_space.ls[[35]]) = line.name

###########################
## "co-4-fan", "four-fan"##
###########################
# house and a single node
co.4.fan.mat = matrix(NA,ncol = 6)
for(i in 1:ncol(pair5.mat)){
      pair5.vec = pair5.mat[, i]
      co.4.fan.mat = rbind(co.4.fan.mat, house(pair5.vec))
}
co.4.fan.mat = co.4.fan.mat[-1, ]
prob_space.ls[[36]] = rbind(fillprob(NN, co.4.fan.mat), 1 - fillprob(NN, co.4.fan.mat))
rownames(prob_space.ls[[36]]) = rep(c("co-4-fan", "four-fan"), each = nrow(prob_space.ls[[36]])/2)
colnames(prob_space.ls[[36]]) = line.name

#################
## "A", "A.bar"##
#################
A.mat = A(NN)
prob_space.ls[[37]] = rbind(fillprob(NN, A.mat), 1 - fillprob(NN, A.mat))
rownames(prob_space.ls[[37]]) = rep(c("A", "A.bar"), each = nrow(prob_space.ls[[37]])/2)
colnames(prob_space.ls[[37]]) = line.name

#################
## "R", "R.bar"##
#################
R.mat = R(NN)
prob_space.ls[[38]] = rbind(fillprob(NN, R.mat), 1 - fillprob(NN, R.mat))
rownames(prob_space.ls[[38]]) = rep(c("R", "R.bar"), each = nrow(prob_space.ls[[38]])/2)
colnames(prob_space.ls[[38]]) = line.name

#########################
## "Twok3", "Twok3.bar"##
#########################
twoK3.mat = twoK3(NN)
prob_space.ls[[39]] = rbind(fillprob(NN, twoK3.mat), 1 - fillprob(NN, twoK3.mat))
rownames(prob_space.ls[[39]]) = rep(c("Twok3", "Twok3.bar"), each = nrow(prob_space.ls[[39]])/2)
colnames(prob_space.ls[[39]]) = line.name

###################
## "c6", "c6.bar"##
###################
c6.mat = c6(NN)
prob_space.ls[[40]] = rbind(fillprob(NN, c6.mat), 1 - fillprob(NN, c6.mat))
rownames(prob_space.ls[[40]]) = rep(c("c6", "c6.bar"), each = nrow(prob_space.ls[[40]])/2)
colnames(prob_space.ls[[40]]) = line.name

#####################
##"x98.bar", "x98" ##
#####################
k = 41
x98.bar.mat = x98.bar(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, x98.bar.mat), 1 - fillprob(NN, x98.bar.mat))
rownames(prob_space.ls[[k]]) = rep(c("x98.bar", "x98"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name


###################
## "s3.bar", "s3"##
###################
k = 42
s3.bar.mat = s3.bar(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, s3.bar.mat), 1 - fillprob(NN, s3.bar.mat))
rownames(prob_space.ls[[k]]) = rep(c("s3.bar", "s3"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name


#####################
## "x18", "x18.bar"##
#####################
k = 43
x18.mat = x18(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, x18.mat), 1 - fillprob(NN, x18.mat))
rownames(prob_space.ls[[k]]) = rep(c("x18", "x18.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

###############################
## "five-pan", "five-pan.bar"##
###############################
k = 44
five.pan.mat = five_pan(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, five.pan.mat), 1 - fillprob(NN, five.pan.mat))
rownames(prob_space.ls[[k]]) = rep(c("five-pan", "five-pan.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name


#######################
## "x166", "x166.bar"##
#######################
k = 45
x166.mat = x166(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, x166.mat), 1 - fillprob(NN, x166.mat))
rownames(prob_space.ls[[k]]) = rep(c("x166", "x166.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

#######################
## "x169", "x169.bar"##
#######################
k = 46
x169.mat = x169(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, x169.mat), 1 - fillprob(NN, x169.mat))
rownames(prob_space.ls[[k]]) = rep(c("x169", "x169.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

#####################
## "x84", "x84.bar"##
#####################
k = 47
x84.mat = x84(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, x84.mat), 1 - fillprob(NN, x84.mat))
rownames(prob_space.ls[[k]]) = rep(c("x84", "x84.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

#####################
## "x95", "x95.bar"##
#####################
k = 48
x95.mat = x95(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, x95.mat), 1 - fillprob(NN, x95.mat))
rownames(prob_space.ls[[k]]) = rep(c("x95", "x95.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

###########################
## "k4_2k1", "k4_2k1.bar"##
###########################
k = 49
k4_2k1.mat = matrix(NA, ncol = 6)
for(i in 1:ncol(pair4.mat)){
      pair4.vec = pair4.mat[, i]
      k4.mat = k4(pair4.vec)
      k4_2k1.mat = rbind(k4_2k1.mat, k4.mat)
}
k4_2k1.mat = k4_2k1.mat[-1, ]
prob_space.ls[[k]] = rbind(fillprob(NN, k4_2k1.mat), 1 - fillprob(NN, k4_2k1.mat))
rownames(prob_space.ls[[k]]) = rep(c("k4_2k1", "k4_2k1.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

#############################
## "dart_k1", "dart_k1.bar"##
#############################
k = 50
dart_k1.mat = matrix(NA, ncol = 6)
for(i in 1:ncol(pair5.mat)){
      pair5.vec = pair5.mat[, i]
      dart.mat = dart(pair5.vec)
      dart_k1.mat = rbind(dart_k1.mat, dart.mat)
}
dart_k1.mat = dart_k1.mat[-1, ]
prob_space.ls[[k]] = rbind(fillprob(NN, dart_k1.mat), 1 - fillprob(NN, dart_k1.mat))
rownames(prob_space.ls[[k]]) = rep(c("dart_k1", "dart_k1.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

###########################
## "k23_k1", "k23_k1.bar"##
###########################
k = 51
k23_k1.mat = matrix(NA, ncol = 6)
for(i in 1:ncol(pair5.mat)){
      pair5.vec = pair5.mat[, i]
      k23.mat = k23(pair5.vec)
      k23_k1.mat = rbind(k23_k1.mat, k23.mat)
}
k23_k1.mat = k23_k1.mat[-1, ]
prob_space.ls[[k]] = rbind(fillprob(NN, k23_k1.mat), 1 - fillprob(NN, k23_k1.mat))
rownames(prob_space.ls[[k]]) = rep(c("k23_k1", "k23_k1.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

#################################
## "dimand_p2", "dimand_p2.bar"##
#################################
k = 52
dimand_p2.mat = matrix(NA, ncol = 6)
for(i in 1:ncol(pair4.mat)){
      pair4.vec = pair4.mat[, i]
      rst.vec = NN[!NN %in% pair4.vec]
      diamond.mat = diamond(pair4.vec)
      p2.vec = rep(deToIn(rst.vec[1], rst.vec[2]), nrow(diamond.mat))
      dimand_p2.mat = rbind(dimand_p2.mat, cbind(diamond.mat, p2.vec))
}
dimand_p2.mat = dimand_p2.mat[-1, ]
prob_space.ls[[k]] = rbind(fillprob(NN, dimand_p2.mat), 1 - fillprob(NN, dimand_p2.mat))
rownames(prob_space.ls[[k]]) = rep(c("dimand_p2", "dimand_p2.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

############################
## "k5-e_k1.bar", "k5-e_k1"##
############################
k = 53
mat = y1(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("k5-e_k1.bar", "k5-e_k1"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

###################
## "y2", "y2.bar"##
###################
k = 54
mat = y2(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("y2", "y2.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name


############################################################
######################## seven edges  ######################
############################################################

############################################
## "claw_k1.bar_k1", "claw_k1.bar_k1.bar" ##
############################################
k = 55
mat = matrix(NA, ncol = 7)
for(i in 1:ncol(pair5.mat)){
      pair5.vec = pair5.mat[, i]
      claw_k1.bar.mat = claw_k1.bar(pair5.vec)
      mat = rbind(mat, claw_k1.bar.mat)
}
mat = mat[-1, ]
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("claw_k1.bar_k1", "claw_k1.bar_k1.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

########################################
## "p2_p3.bar_k1", "p2_p3.bar_k1.bar" ##
########################################
k = 56
mat = matrix(NA, ncol = 7)
for(i in 1:ncol(pair5.mat)){
      pair5.vec = pair5.mat[, i]
      p2_p3.bar.mat = p2_p3.bar(pair5.vec)
      mat = rbind(mat, p2_p3.bar.mat)
}
mat = mat[-1, ]
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("p2_p3.bar_k1", "p2_p3.bar_k1.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name


############################
## "gem_k1", "gem_k1.bar" ##
############################
k = 57
mat = matrix(NA, ncol = 7)
for(i in 1:ncol(pair5.mat)){
      pair5.vec = pair5.mat[, i]
      gem.mat = gem(pair5.vec)
      mat = rbind(mat, gem.mat)
}
mat = mat[-1, ]
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("gem_k1", "gem_k1.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name


##########################################
## "k3_2k1.bar_k1", "k3_2k1.bar_k1.bar" ##
##########################################
k = 58
mat = matrix(NA, ncol = 7)
for(i in 1:ncol(pair5.mat)){
      pair5.vec = pair5.mat[, i]
      k3_2k1.bar.mat = k3_2k1.bar(pair5.vec)
      mat = rbind(mat, k3_2k1.bar.mat)
}
mat = mat[-1, ]
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("k3_2k1.bar_k1", "k3_2k1.bar_k1.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

##########################
## "k4_p2", "k4_p2.bar" ##
##########################
k = 59
mat = matrix(NA, ncol = 7)
for(i in 1:ncol(pair4.mat)){
      pair4.vec = pair4.mat[, i]
      rst = NN[!NN %in% pair4.vec]
      k4.vec = k4(pair4.vec)
      k4_p2.vec = c(k4.vec, deToIn(rst[1], rst[2]))
      mat = rbind(mat, k4_p2.vec)
}
mat = mat[-1, ]
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("k4_p2", "k4_p2.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

##########################
## "w4_k1.bar", "w4_k1" ##
##########################
k = 60
mat = w4_k1.bar(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("w4_k1.bar", "w4_k1"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

######################
## "x37", "x37.bar" ##
######################
k = 61
mat = x37(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("x37", "x37.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

#######################
## "fish", "co-fish" ##
######################
k = 62
mat = fish(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("fish", "co-fish"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

###########################
## "domino", "co-domino" ##
###########################
k = 63
mat = domino(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("domino", "co-domino"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

##############################
## "twin-c5", "twin-c5.bar" ##
##############################
k = 64
mat = twin.c5(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("twin-c5", "twin-c5.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

######################
## "x58.bar", "x58" ##
######################
k = 65
mat = x58.bar(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("x58.bar", "x58"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

##########################
## "k33-e.bar", "k33-e" ##
##########################
k = 66
mat = twoK3_e(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("k33-e.bar", "k33-e"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

####################
## "x5.bar", "x5" ##
####################
k = 67
mat = x5.bar(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("x5.bar", "x5"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

#############################
## "antenna", "co-antenna" ##
#############################
k = 68
mat = antenna(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("antenna", "co-antenna"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name


#############################
## "x45", "x45.bar" #########
#############################
k = 69
mat = x45(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("x45", "x45.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

########################################
## "twin-house.bar", "twin-house" ######
########################################
k = 70
mat = twin.house.bar(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("twin-house.bar", "twin-house"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

#############################
## "x167", "x167.bar" #######
#############################
k = 71
mat = x167(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("x167", "x167.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

#############################
## "x168", "x168.bar" #######
#############################
k = 72
mat = x168(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("x168", "x168.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

#############################
## "x170", "x170.bar" #######
#############################
k = 73
mat = x170(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("x170", "x170.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

#############################
## "x171", "x171.bar" #######
#############################
k = 74
mat = x171(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("x171", "x171.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

###########################
## "x96", "x96.bar" #######
###########################
k = 75
mat = x96(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("x96", "x96.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

#############################
## "x163", "x163.bar" #######
#############################
k = 76
mat = x163(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("x163", "x163.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

################################################
## "p3_2k1.bar_k1.bar", "p3_2k1.bar_k1" #######
################################################
k = 77
mat = y3(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("p3_2k1.bar_k1.bar", "p3_2k1.bar_k1"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name

#########################
## "y4", "y4.bar" #######
#########################
k = 78
mat = y4(NN)
prob_space.ls[[k]] = rbind(fillprob(NN, mat), 1 - fillprob(NN, mat))
rownames(prob_space.ls[[k]]) = rep(c("y4", "y4.bar"), each = nrow(prob_space.ls[[k]])/2)
colnames(prob_space.ls[[k]]) = line.name
# checking each graphs
test = prob_space.ls[[k]]
dim(test)
apply(test, 2, sum)
checkUnique(test)


################################
#
#
#
# save the probability space
prob_space.mat = do.call(rbind, prob_space.ls)
write.csv(prob_space.mat, "data/six_nodes_graphs.csv")
save(prob_space.mat, file = "data/six_nodes_graphs.Rdata")

# delete the repeated rows in prob_space.mat
prob_space.mat2 = prob_space.mat[!duplicated(rownames(prob_space.mat)),]
write.csv(prob_space.mat2, file = 'data/six_nodes_unique_graphs.csv')



# Mobuis parameter represents probability space
prob_mobius.mat = prob_mobius(prob_space.mat)

# delete the repeated rows in prob_mobius.mtx
prob_mobius.mat2 = prob_mobius.mat[!duplicated(rownames(prob_mobius.mat)),]

# save
save(prob_mobius.mat, file = 'data/six_nodes_mobius.Rdata')
write.csv(prob_mobius.mat, file = 'data/six_nodes_mobius.csv')


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

save(vmat.prob_mobius.q, file = "data/extreme_points_six_node_graph.Rdata")
write.csv(vmat.prob_mobius.q, "data/extreme_points_six_node_graph.csv")



