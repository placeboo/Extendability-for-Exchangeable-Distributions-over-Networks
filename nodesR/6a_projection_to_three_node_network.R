rm(list = ls())

#----------------------------------------------------------------#
# function, build data frame showing the points after projection,
# and how the convex combination
#----------------------------------------------------------------#
convexComb = function(Lambdas, mat1, mat2_proj){
        combina.vec = rep(0,  ncol(Lambdas))
        number_cmbn = rep(0,  ncol(Lambdas))
        is_kept = rep(0,  ncol(Lambdas)) # binary, 1: yes, 0: no
        graph_kept = rep(0,  ncol(Lambdas)) # what graph 
        
        rowname = rownames(mat1)
        graphKept = rowname
        for(i in 1: ncol(Lambdas)){
                iLambda = Lambdas[, i]
                graphs = rowname[iLambda != 0]
                is_kept[i] = ifelse(any(graphs %in% graphKept), 1, 0)
                graph_kept[i] = paste(graphs[graphs %in% graphKept], collapse = ", ")
                vec = paste(iLambda[iLambda != 0], graphs)
                combina.vec[i] = paste(vec, collapse = ", ")
                number_cmbn[i] = length(vec)
        }
        return(data.frame(number = number_cmbn, combination = combina.vec, isKept = is_kept, graphKept = graph_kept, mat2_proj, check.names = F))
}

#----------------------------------------------------------------#
# 08/13/2018
# function, output table represent hyperplane after projections
# see the bottom
# input: hmat
#----------------------------------------------------------------#
out_hmat_tab = function(hmat){
        ineq_mat = hmat[which(hmat[,1]==0),]
        tab = ineq_mat[, -c(1,2)] * -1
        
        tab = cbind(tab, paste("<=", ineq_mat[,2], sep = ""))
        return(tab)
}


#----------------------------------------------------------------#
# load data
# 3-node network, extreme points
# 4-node network projects to 3-node network, extreme points
# 5-node network projects to 3-node network, extreme points
# 6-node network projects to 3-node network, extreme points
# 7-node network projects to 3-node network, extreme points
#----------------------------------------------------------------#

three_extreme = get(load("data/extreme_points_three_node_graph.Rdata"))
four_to_three = get(load("data/four_nodes_to_three_nodes_projection.Rdata"))
five_to_three = get(load("data/five_nodes_to_three_nodes_projection.Rdata"))
six_to_three = get(load("data/six_nodes_to_three_nodes_projection.Rdata"))
load("data/seven_nodes_to_three_nodes_projection.Rdata")
assign('seven_to_three', vmat.proj.minimal.q)


three_extreme1 = three_extreme[, -c(1,2)]
four_to_three1 = four_to_three[, -c(1,2)]
five_to_three1 = five_to_three[, -c(1,2)]
six_to_three1 = six_to_three[, -c(1,2)]
seven_to_three1 = seven_to_three[, -c(1,2)]

# re-orgnized. The colnumn names should be the same
four_to_three1 = four_to_three1[, colnames(three_extreme1)] #"3k1","k3","p3.bar","p3"  
five_to_three1 = five_to_three1[, colnames(three_extreme1)]
six_to_three1 = six_to_three1[, colnames(three_extreme1)]
seven_to_three1 = seven_to_three1[, c(1,4,2,3)]

three_extreme.d = q2d(three_extreme1)
four_to_three.d = q2d(four_to_three1)
five_to_three.d = q2d(five_to_three1)
six_to_three.d = q2d(six_to_three1)
seven_to_three.d = q2d(seven_to_three1)

Lambdas5to3 = solve(t(three_extreme.d)) %*% t(five_to_three.d)
Lambdas5to3_r = round(Lambdas5to3, 3)

Lambdas6to3 = solve(t(three_extreme.d)) %*% t(six_to_three.d)
Lambdas6to3_r = round(Lambdas6to3, 3)

Lambdas7to3 = solve(t(three_extreme.d)) %*% t(seven_to_three.d)
Lambdas7to3_r = round(Lambdas7to3, 3)


five_to_three_combn = convexComb(Lambdas5to3_r, mat1 = three_extreme.d, mat2_proj = five_to_three1)

six_to_three_combn = convexComb(Lambdas6to3_r, mat1 = three_extreme.d, mat2_proj = six_to_three1)

seven_to_three_combn = convexComb(Lambdas7to3_r, mat1 = three_extreme.d, mat2_proj = seven_to_three1)


write.csv(five_to_three_combn, file = "data/five_to_three_combination.csv")
save(five_to_three_combn, file = "data/five_to_three_combination.Rdata")
write.csv(six_to_three_combn, file = "data/six_to_three_combination.csv")
save(six_to_three_combn, file = "data/six_to_three_combination.Rdata")
write.csv(seven_to_three_combn, file = "data/seven_to_three_combination.csv")
save(seven_to_three_combn, file = "data/seven_to_three_combination.Rdata")

# test = rbind(five_to_three, six_to_three)
# nonredundant.test = redundant(test, representation = "V")
# index =  c(1: nrow(test))[-nonredundant.test$redundant]
# 
# rownames(nonredundant.test$output) = rownames(test)[index]

#----------------------------------------------------------------#
# manage the location of data, prepare for graph
# empty.vec should be (0, 0, 0, 1) in 4-dim. However, we will know the value in t-cord since sum to 1
#----------------------------------------------------------------#
empty.vec = c(0, 0, 0)
p3_bar.vec = c(1, 0, 0)
p3.vec = c(0, 1, 0)
k3.vec = c(0, 0, 1)

cord.mat = cbind(p3.vec, empty.vec, p3_bar.vec, k3.vec)
rownames(cord.mat) = c("x", "y", "z")

pnt5to3.mat = cord.mat %*% Lambdas5to3
pnt6to3.mat = cord.mat %*% Lambdas6to3
pnt7to3.mat = cord.mat %*% Lambdas7to3
#----------------------------------------------------------------#
# find the hyperplan of each projectation, using joint probability
# M = (0, b, -A)   Ax \leq b
# M = (1, b, -A)   Ax = b
#----------------------------------------------------------------#
# 5->3 
five_to_three.vmat = makeV(t(pnt5to3.mat))
five_to_three.hmat = scdd(five_to_three.vmat)$out
colnames(five_to_three.hmat)=  c("","", "p3.bar", "p3", "k3")
save(five_to_three.hmat, file = "data/five_to_three_hmat.Rdata")

five_to_threeIneq.hmat = five_to_three.hmat[which(five_to_three.hmat[,1] == 0), ]

out_hmat_tab(round(five_to_three.hmat,2))

# show the ineq
five_to_threeIneq.hmat.r = round(five_to_threeIneq.hmat, 2)
hyperplane = t(apply(-five_to_threeIneq.hmat.r, 1, function(x) paste(x[-(1:2)], c("p3.bar", "p3", "k3"), sep = "")))
hyperplane = apply(hyperplane, 1, function(x) paste(x[1], x[2], x[3], sep = "+"))
hyperplane = paste(hyperplane, "<=", five_to_threeIneq.hmat.r[,2], sep = "")

five_to_three_check.mat = round(-five_to_threeIneq.hmat[, -c(1,2)] %*% pnt5to3.mat - five_to_threeIneq.hmat[, 2],4)
five_to_three_hypermat = matrix(0, nrow = nrow(five_to_three_check.mat), ncol = ncol(five_to_three_check.mat))
for(i in 1: nrow(five_to_three_hypermat)){
        five_to_three_hypermat[i,which(five_to_three_check.mat[i,]==0)] = 1 
}

# colname: graphs's name and its location 
five_to_three_colname = apply(pnt5to3.mat, 2, function(x) paste(x[1], x[2], x[3],sep = ","))
five_to_three_colname = rbind(colnames(five_to_three_check.mat), paste("(", five_to_three_colname, ")", sep = ""))
colnames(five_to_three_hypermat) = colnames(five_to_three_check.mat)
five_to_three_hypermat2 = rbind(five_to_three_colname, five_to_three_hypermat)

# output
print(xtable(five_to_three_hypermat2, caption = "Polytopes of five-node network projecting to three-node network. 1: the corners of hyperplan in Polytopes; otherwise 0.", digits = 0, align = rep("c", ncol(five_to_three_hypermat)+1)), include.colnames = F)

# save
save(five_to_three_hypermat, file = "data/hyperplane_vertex_five_to_three.Rdata")

# 6->3
six_to_three.vmat = makeV(t(pnt6to3.mat))
six_to_three.hmat = scdd(six_to_three.vmat)$out

out_hmat_tab(round(six_to_three.hmat,2))

save(six_to_three.hmat, file = "data/six_to_three_hmat.Rdata")
six_to_threeIneq.hmat = six_to_three.hmat[which(six_to_three.hmat[,1] == 0), ]
round(six_to_threeIneq.hmat, 3)

# show the ineq
six_to_threeIneq.hmat.r = round(six_to_threeIneq.hmat, 2)
six_to_three_hyperplane = t(apply(-six_to_threeIneq.hmat.r, 1, function(x) paste(x[-(1:2)], c("x", "y", "z"), sep = "")))
six_to_three_hyperplane = apply(six_to_three_hyperplane, 1, function(x) paste(x[1], x[2], x[3], sep = "+"))
six_to_three_hyperplane = paste(six_to_three_hyperplane, "<=", six_to_threeIneq.hmat.r[,2], sep = "")

six_to_three_check.mat = round(-six_to_threeIneq.hmat[, -c(1,2)] %*% pnt6to3.mat - six_to_threeIneq.hmat[, 2],4)
six_to_three_hypermat = matrix(0, nrow = nrow(six_to_three_check.mat), ncol = ncol(six_to_three_check.mat))
for(i in 1: nrow(six_to_three_hypermat)){
        six_to_three_hypermat[i,which(six_to_three_check.mat[i,]==0)] = 1 
}

# colname: graphs's name and its location 
colnames(six_to_three_hypermat) = colnames(six_to_three_check.mat)
six_to_three_colname = apply(pnt6to3.mat, 2, function(x) paste(x[1], x[2], x[3],sep = ","))
six_to_three_colname = rbind(colnames(six_to_three_check.mat), paste("(", six_to_three_colname, ")", sep = ""))

six_to_three_hypermat2 = rbind(six_to_three_colname, six_to_three_hypermat)

# output
print(xtable(six_to_three_hypermat2, caption = "Polytopes of six-node network projecting to three-node network. 1: the corners of hyperplan in Polytopes; otherwise 0.", digits = 0, align = rep("c", ncol(six_to_three_hypermat)+1)), include.colnames = F)
save(six_to_three_hypermat2, file = "data/hyperplane_vertex_six_to_three.Rdata")


# 7 -> 3
seven_to_three.vmat = makeV(t(pnt7to3.mat))
seven_to_three.hmat = scdd(seven_to_three.vmat)$out
save(seven_to_three.hmat, file = "data/seven_to_three_hmat.Rdata")

seven_to_threeIneq.hmat = seven_to_three.hmat[which(seven_to_three.hmat[,1] == 0), ]
round(seven_to_threeIneq.hmat, 3)

# show the ineq
seven_to_threeIneq.hmat.r = round(seven_to_threeIneq.hmat, 2)
seven_to_three_hyperplane = t(apply(-seven_to_threeIneq.hmat.r, 1, function(x) paste(x[-(1:2)], c("p3_bar", "p3", "k3"), sep = "")))
seven_to_three_hyperplane = apply(seven_to_three_hyperplane, 1, function(x) paste(x[1], x[2], x[3], sep = "+"))
seven_to_three_hyperplane = paste(seven_to_three_hyperplane, "<=", seven_to_threeIneq.hmat.r[,2], sep = "")

seven_to_three_check.mat = round(-seven_to_threeIneq.hmat[, -c(1,2)] %*% pnt7to3.mat - seven_to_threeIneq.hmat[, 2],4)
seven_to_three_hypermat = matrix(0, nrow = nrow(seven_to_three_check.mat), ncol = ncol(seven_to_three_check.mat))
for(i in 1: nrow(seven_to_three_hypermat)){
        seven_to_three_hypermat[i,which(seven_to_three_check.mat[i,]==0)] = 1 
}

# colname: graphs's name and its location 
colnames(seven_to_three_hypermat) = colnames(seven_to_three_check.mat)
seven_to_three_colname = apply(pnt6to3.mat, 2, function(x) paste(x[1], x[2], x[3],sep = ","))
seven_to_three_colname = rbind(colnames(seven_to_three_check.mat), paste("(", seven_to_three_colname, ")", sep = ""))

seven_to_three_hypermat2 = rbind(seven_to_three_colname, seven_to_three_hypermat)
#----------------------------------------------------------------#
# drawing
#----------------------------------------------------------------#
#setwd("/Users/thomas/madrid/stat/566/examples")
#source("http://www.stat.washington.edu/tsr/s566/labs/y0y1polytopenew-rgl-col.R")
## Set up interactive plot
# open 3d window
open3d()

# resize window

par3d(windowRect = c(200, 200, 800, 800))

# plot 3-simplx
plot3d(x = cord.mat[1, ], y = cord.mat[2, ], z = cord.mat[3,], size = 8, add = T)
# add 3-node network extreme points
text3d(cord.mat[, 1], text = "p3", adj = 1.2)
text3d(cord.mat[, 2], text = "3k1", adj = 1.2)
text3d(cord.mat[, 3], text = "p3.bar", adj = 1.2)
text3d(cord.mat[, 4], text = "k3", adj = 1.2)


## drow 3-simple
lines3d(x = c(0, 0, 1, 0, 0), y = c(0, 0, 0, 1, 0), z = c(1, 0, 0, 0, 0), col="black",lty=3)
lines3d(x = c(0, 0), y = c(0, 1), z = c(1, 0), col="black",lty=3)
lines3d(x = c(0, 1), y = c(0, 0), z = c(1, 0), col="black",lty=3)


points3d(x = pnt5to3.mat[1, ], y = pnt5to3.mat[2, ], z = pnt5to3.mat[3,], col = "blue", size = 8)

points3d(x = pnt6to3.mat[1, ], y = pnt6to3.mat[2, ], z = pnt6to3.mat[3,], col = "red", size = 15)

# add legend
legend3d("topright", legend = c("5-node network to 3-node", "6-node network to 3-node"), pch = 16, col = c("blue", "red"), cex=2, inset=c(0.02))
#legend3d("topright", legend = c("5-node network to 3-node", "6-node network to 3-node", "7-node network to 3-node"), pch = 16, col = c("blue", "red", "purple"), cex=2, inset=c(0.02))

# what are those points?
for(i in 1: ncol(pnt5to3.mat)){
        text3d(as.numeric(pnt5to3.mat[,i]), texts = colnames(pnt5to3.mat)[i], adj = 1.2)
}

for(i in 1: ncol(pnt6to3.mat)){
        text3d(as.numeric(pnt6to3.mat[,i]), texts = colnames(pnt6to3.mat)[i], adj = -0.2)
}

# six_to_three_combn['c6', 'combination'] 
# five_to_three_combn[c('c5', 'p2_p3'), 'combination'] # c6 is in the middle of c5 and p2_p3


# polytopes of five-node projection
# non-trial faces, the number of vertexes of the face =< 3
for(i in 1: nrow(five_to_three_hypermat)){
        iface = five_to_three_hypermat[i, ]
        if(sum(iface) < 4){
                iloc.mat = pnt5to3.mat[, ifelse(iface == 1, T, F)]
                # draw
                lines3d(x = c(iloc.mat[1, ], iloc.mat[1, 1]), y =c(iloc.mat[2, ], iloc.mat[2, 1]), z = c(iloc.mat[3, ], iloc.mat[3, 1]), col="blue",lty=3)
        }
}

# polytopes of six-node projection
# non-trial faces, the number of vertexes of the face =< 3
for(i in 1: nrow(six_to_three_hypermat)){
        iface = six_to_three_hypermat[i, ]
        if(sum(iface) < 4){
                iloc.mat = pnt6to3.mat[, ifelse(iface == 1, T, F)]
                # draw
                lines3d(x = c(iloc.mat[1, ], iloc.mat[1, 1]), y =c(iloc.mat[2, ], iloc.mat[2, 1]), z = c(iloc.mat[3, ], iloc.mat[3, 1]), col="red",lty=3)
        }
}

# one special case in six projection, the face is made by 2k3, c6, c6.bar, and 2k3.bar
lines3d(x = pnt6to3.mat["x", c("c6", "Twok3.bar")], y = pnt6to3.mat["y", c("c6", "Twok3.bar")], z = pnt6to3.mat["z", c("c6", "Twok3.bar")], col="red",lty=3)
lines3d(x = pnt6to3.mat["x", c("c6.bar", "Twok3")], y = pnt6to3.mat["y", c("c6.bar", "Twok3")], z = pnt6to3.mat["z", c("c6.bar", "Twok3")], col="red",lty=3)

#rgl.postscript("persp3dd.pdf","pdf")

# add the symmetric line
lines3d(x = c(0.5, 0), y = c(0.5, 0), z = c(0, 0.5), col="green",lty=3)
# add Erdos Renyi line
p_vec = seq(0, 1, by = 0.01)
p_mat = edge_p(p_vec)
curve_cord = cord.mat %*% t(p_mat)
# add edge
lines3d(x = curve_cord[1,], y = curve_cord[2,], z = curve_cord[3,], col=8,lty=3)
browseURL(paste("file://", writeWebGL(dir=file.path("plots/", "webGL"), width=700), sep=""))
# close window
rgl.close()

#----------------------------------------------------------#
# 7-> 3
#----------------------------------------------------------#
open3d()

# resize window

par3d(windowRect = c(200, 200, 800, 800))

# plot 3-simplx
plot3d(x = cord.mat[1, ], y = cord.mat[2, ], z = cord.mat[3,], size = 8, add = T)
# add 3-node network extreme points
text3d(cord.mat[, 1], text = "p3", adj = 1.2)
text3d(cord.mat[, 2], text = "3k1", adj = 1.2)
text3d(cord.mat[, 3], text = "p3.bar", adj = 1.2)
text3d(cord.mat[, 4], text = "k3", adj = 1.2)


## drow 3-simple
lines3d(x = c(0, 0, 1, 0, 0), y = c(0, 0, 0, 1, 0), z = c(1, 0, 0, 0, 0), col="black",lty=3)
lines3d(x = c(0, 0), y = c(0, 1), z = c(1, 0), col="black",lty=3)
lines3d(x = c(0, 1), y = c(0, 0), z = c(1, 0), col="black",lty=3)

# nodes
points3d(x = pnt7to3.mat[1, ], y = pnt7to3.mat[2, ], z = pnt7to3.mat[3,], col = "purple", size = 8)
# connect face with 3 nodes
for(i in 1: nrow(seven_to_three_hypermat)){
        iface = seven_to_three_hypermat[i, ]
        if(sum(iface) < 4){
                iloc.mat = pnt7to3.mat[, ifelse(iface == 1, T, F)]
                # draw
                lines3d(x = c(iloc.mat[1, ], iloc.mat[1, 1]), y =c(iloc.mat[2, ], iloc.mat[2, 1]), z = c(iloc.mat[3, ], iloc.mat[3, 1]), col="purple",lty=3)
        }
}
#for(i in 1: ncol(pnt7to3.mat)){
#        text3d(as.numeric(pnt7to3.mat[,i]), texts = colnames(pnt7to3.mat)[i], adj = 1.2)
#}

# adding untrivial face
# face with four faces
tmp_idx = which(apply(seven_to_three_hypermat, 1, sum) == 4)
seven_to_three_hypermat[tmp_idx,]
for(i in 1:length(tmp_idx)){
        i_face = seven_to_three_hypermat[tmp_idx[i], ]
        print(which(i_face == 1))
}
lines3d(x = pnt7to3.mat[1, c("class175", "class627", "class571")], y = pnt7to3.mat[2, c("class175", "class627", "class571")], z = pnt7to3.mat[3, c("class175", "class627", "class571")],col="purple",lty=3)

lines3d(x = pnt7to3.mat[1, c("class903", "class968", "class807","class780")], y = pnt7to3.mat[2, c("class903", "class968", "class807","class780")], z = pnt7to3.mat[3, c("class903", "class968", "class807","class780")],col="purple",lty=3)

lines3d(x = pnt7to3.mat[1, c("class1038", "class1029", "class1017")], y = pnt7to3.mat[2, c("class1038", "class1029", "class1017")], z = pnt7to3.mat[3, c("class1038", "class1029", "class1017")],col="purple",lty=3)

# face with five faces
tmp_idx = which(apply(seven_to_three_hypermat, 1, sum) == 5)
seven_to_three_hypermat[tmp_idx,]
for(i in 1:length(tmp_idx)){
        i_face = seven_to_three_hypermat[tmp_idx[i], ]
        print(which(i_face == 1))
}
lines3d(x = pnt7to3.mat[1, c("class1029", "class968")], y = pnt7to3.mat[2, c("class1029", "class968")], z = pnt7to3.mat[3, c("class1029", "class968")],col="purple",lty=3)
lines3d(x = pnt7to3.mat[1, c("class780", "class627")], y = pnt7to3.mat[2, c("class780", "class627")], z = pnt7to3.mat[3, c("class780", "class627")],col="purple",lty=3)

# face with six nodes
tmp_idx = which(apply(seven_to_three_hypermat, 1, sum) == 6)
seven_to_three_hypermat[tmp_idx,]
for(i in 1:length(tmp_idx)){
        i_face = seven_to_three_hypermat[tmp_idx[i], ]
        print(which(i_face == 1))
}

lines3d(x = curve_cord[1,], y = curve_cord[2,], z = curve_cord[3,], col="dark blue",lty=3)
# add the symmetric line
lines3d(x = c(0.5, 0), y = c(0.5, 0), z = c(0, 0.5), col="red",lty=3)
# save
browseURL(paste("file://", writeWebGL(dir=file.path("plots/", "webGL"), width=700), sep=""))
# close window
rgl.close()


# add the Erdos Renyi, suppose every edge occurs independently with probability p
edge_p = function(p_vec){
        # one egde probability p, vec
        p_mat = matrix(0, ncol = 4, nrow = length(p_vec))
        for (i in 1: length(p_vec)){
                p_i = p_vec[i]
                p_mat[i, ] = c(3*p_i^2*(1-p_i), (1-p_i)^3, 3* p_i*(1-p_i)^2, p_i^3)
        }
        colnames(p_mat) = c("p3", "3k1", "p3.bar", "k3")
        return(p_mat)
}

remove_edge_p = function(p_vec){
        # one egde probability p, vec
        p_mat = matrix(0, ncol = 4, nrow = length(p_vec))
        for (i in 1: length(p_vec)){
                p_i = p_vec[i]
                p_mat[i, ] = c(3* p_i*(1-p_i)^2, p_i^3, 3*p_i^2*(1-p_i),  (1-p_i)^3)
        }
        colnames(p_mat) = c("p3", "3k1", "p3.bar", "k3")
        return(p_mat)
}
remove_p_mat = remove_edge_p(p_vec)
remove_curve = cord.mat %*% t(remove_p_mat)


for(j in 1: ncol(remove_curve)){
        lines3d(x = c(curve_cord[1,j], remove_curve[1,j]), y = c(curve_cord[2,j], remove_curve[2,j]), z = c(curve_cord[3,j], remove_curve[3,j]),col="purple",lty=3)
}

for(j in 1: ncol(remove_curve)){
        lines3d(x = c(0, remove_curve[1,j]), y = c(0, remove_curve[2,j]), z = c(0, remove_curve[3,j]),col="purple",lty=3)
}
for(j in 1: ncol(remove_curve)){
        lines3d(x = c(0, remove_curve[1,j]), y = c(0, remove_curve[2,j]), z = c(1, remove_curve[3,j]),col="purple",lty=3)
}
