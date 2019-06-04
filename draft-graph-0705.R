library(rgl)
library(rcdd)

par3d()$userMatrix
user_mat =par3d()$userMatrix

user_mat = rbind(c(0.6078451, -0.7939834, -0.01071101,    0),
                 c(0.1489580,  0.1007667,  0.98369592,    0),
                 c(-0.7799590, -0.5995303,  0.17952067,    0),
                 c( 0.0000000,  0.0000000,  0.00000000,    1))

open3d(zoom = 1, userMatrix = user_mat)
# resize window
par3d(windowRect = c(200, 200, 800, 800))
plot3d(c(1,0,0,1),c(0,1,0,0),c(0,0,1,0), ylab="%P3_bar",xlab="%P3",zlab="%K3",par3d(FOV=1)) 

lines3d(x = c(1,0,0,1,0,0),y = c(0,1,0,0,0,0),z = c(0,0,1,0,0,1), col="dark grey",lwd=5)
lines3d(x=c(0,0), y = c(0,1) , z = c(0,0), col="dark grey",lwd=5)
points3d(x = c(1,0,0,0),y = c(0,1,0,0),z = c(0,0,1,0), col="red", size = 12)
rgl.snapshot("plots/threeD-three-node.png")


empty.vec = c(0, 0, 0)
p3_bar.vec = c(1, 0, 0)
p3.vec = c(0, 1, 0)
k3.vec = c(0, 0, 1)

cord.mat = cbind(p3.vec, empty.vec, p3_bar.vec, k3.vec)
rownames(cord.mat) = c("x", "y", "z")

three_extreme = get(load("data/extreme_points_three_node_graph.Rdata"))
four_to_three = get(load("data/four_nodes_to_three_nodes_projection.Rdata"))
five_to_three = get(load("data/five_nodes_to_three_nodes_projection.Rdata"))
six_to_three = get(load("data/six_nodes_to_three_nodes_projection.Rdata"))
load("data/seven_nodes_to_three_nodes_projection.Rdata")
assign('seven_to_three', vmat.proj.minimal.q)

five_to_three_hypermat = get(load(file = "data/hyperplane_vertex_five_to_three.Rdata"))
six_to_three_hypermat = get(load(file = "data/hyperplane_vertex_six_to_three.Rdata"))
seven_to_three_hypermat = get(load(file = "data/hyperplane_vertex_five_to_three.Rdata"))

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

pnt5to3.mat = cord.mat %*% Lambdas5to3
pnt6to3.mat = cord.mat %*% Lambdas6to3
pnt7to3.mat = cord.mat %*% Lambdas7to3
#----------------------------------------------------#
# 5->3
#----------------------------------------------------#
open3d(zoom = 1, userMatrix = user_mat)
# resize window
par3d(windowRect = c(200, 200, 800, 800))
plot3d(c(1,0,0,1),c(0,1,0,0),c(0,0,1,0), ylab="%P3_bar",xlab="%P3",zlab="%K3",par3d(FOV=1)) 

lines3d(x = c(1,0,0,1,0,0),y = c(0,1,0,0,0,0),z = c(0,0,1,0,0,1), col="dark grey",lwd=3)
lines3d(x=c(0,0), y = c(0,1) , z = c(0,0), col="dark grey",lwd=3)
#points3d(x = c(1,0,0,0),y = c(0,1,0,0),z = c(0,0,1,0), col="red", size = 8)
points3d(x = pnt5to3.mat[1, ], y = pnt5to3.mat[2, ], z = pnt5to3.mat[3,], col = "blue", size = 8)
# polytopes of five-node projection
# non-trial faces, the number of vertexes of the face =< 3
for(i in 1: nrow(five_to_three_hypermat)){
        iface = five_to_three_hypermat[i, ]
        if(sum(iface) < 4){
                iloc.mat = pnt5to3.mat[, ifelse(iface == 1, T, F)]
                # draw
                lines3d(x = c(iloc.mat[1, ], iloc.mat[1, 1]), y =c(iloc.mat[2, ], iloc.mat[2, 1]), z = c(iloc.mat[3, ], iloc.mat[3, 1]), col="blue",lty=3, lwd=3)
        }
}

lines3d(x = pnt5to3.mat["x", c("k23.bar", "k5", "w4")], y = pnt5to3.mat["y", c("k23.bar", "k5", "w4")], z = pnt5to3.mat["z", c("k23.bar", "k5", "w4")], col="blue",lty=3, lwd=5)

tmp = pnt5to3.mat[, c("5k1", 'k5')]
lines3d(x = tmp[1,], y = tmp[2,], z=tmp[3,], col="blue",lty=3, lwd=5)

tmp = pnt5to3.mat[, c("w4.bar", "5k1", "k23")]
lines3d(x = tmp[1,], y = tmp[2,], z=tmp[3,], col="blue",lty=3, lwd=5)

rgl.snapshot("plots/threeD-five-to-three-node.png")

#----------------------------------------------------#
# add 6->3
#----------------------------------------------------#
points3d(x = pnt6to3.mat[1, ], y = pnt6to3.mat[2, ], z = pnt6to3.mat[3,], col = "darkorange", size = 11)

# polytopes of six-node projection
# non-trial faces, the number of vertexes of the face =< 3
for(i in 1: nrow(six_to_three_hypermat)){
        iface = six_to_three_hypermat[i, ]
        if(sum(iface) < 4){
                iloc.mat = pnt6to3.mat[, ifelse(iface == 1, T, F)]
                # draw
                lines3d(x = c(iloc.mat[1, ], iloc.mat[1, 1]), y =c(iloc.mat[2, ], iloc.mat[2, 1]), z = c(iloc.mat[3, ], iloc.mat[3, 1]), col="darkorange",lty=3, lwd=5)
        }
}

# one special case in six projection, the face is made by 2k3, c6, c6.bar, and 2k3.bar
lines3d(x = pnt6to3.mat["x", c("c6", "Twok3.bar")], y = pnt6to3.mat["y", c("c6", "Twok3.bar")], z = pnt6to3.mat["z", c("c6", "Twok3.bar")], col="darkorange",lty=3, lwd=5)
lines3d(x = pnt6to3.mat["x", c("c6.bar", "Twok3")], y = pnt6to3.mat["y", c("c6.bar", "Twok3")], z = pnt6to3.mat["z", c("c6.bar", "Twok3")], col="darkorange",lty=3, lwd=5)

tmp = pnt6to3.mat[, c("Twok3", "k6", "threeK2.bar")]
lines3d(x = tmp[1,], y = tmp[2,], z=tmp[3,], col="darkorange",lty=3, lwd=5)

tmp = pnt6to3.mat[, c("k6", "k6.bar")]
lines3d(x = tmp[1,], y = tmp[2,], z=tmp[3,], col="darkorange",lty=3, lwd=5)

tmp = pnt6to3.mat[, c("threeK2", "k6.bar", "Twok3.bar")]
lines3d(x = tmp[1,], y = tmp[2,], z=tmp[3,], col="darkorange",lty=3, lwd=5)

rgl.snapshot("plots/threeD-six-to-three-node.png")

#----------------------------------------------------#
# add 7->3
#----------------------------------------------------#  
points3d(x = pnt7to3.mat[1, ], y = pnt7to3.mat[2, ], z = pnt7to3.mat[3,], col = "purple", size = 8)

# connect face with 3 nodes
for(i in 1: nrow(seven_to_three_hypermat)){
        iface = seven_to_three_hypermat[i, ]
        if(sum(iface) < 4){
                iloc.mat = pnt7to3.mat[, ifelse(iface == 1, T, F)]
                # draw
                lines3d(x = c(iloc.mat[1, ], iloc.mat[1, 1]), y =c(iloc.mat[2, ], iloc.mat[2, 1]), z = c(iloc.mat[3, ], iloc.mat[3, 1]), col="purple",lty=3, lwd = 3)
        }
}

# adding untrivial face
# face with four faces
tmp_idx = which(apply(seven_to_three_hypermat, 1, sum) == 4)
seven_to_three_hypermat[tmp_idx,]
for(i in 1:length(tmp_idx)){
        i_face = seven_to_three_hypermat[tmp_idx[i], ]
        print(which(i_face == 1))
}
lines3d(x = pnt7to3.mat[1, c("class175", "class627", "class571")], y = pnt7to3.mat[2, c("class175", "class627", "class571")], z = pnt7to3.mat[3, c("class175", "class627", "class571")],col="purple",lty=3, lwd = 5)

lines3d(x = pnt7to3.mat[1, c("class903", "class968", "class807","class780")], y = pnt7to3.mat[2, c("class903", "class968", "class807","class780")], z = pnt7to3.mat[3, c("class903", "class968", "class807","class780")],col="purple",lty=3, lwd = 5)

lines3d(x = pnt7to3.mat[1, c("class1038", "class1029", "class1017")], y = pnt7to3.mat[2, c("class1038", "class1029", "class1017")], z = pnt7to3.mat[3, c("class1038", "class1029", "class1017")],col="purple",lty=3, lwd = 5)

# face with five faces
tmp_idx = which(apply(seven_to_three_hypermat, 1, sum) == 5)
seven_to_three_hypermat[tmp_idx,]
for(i in 1:length(tmp_idx)){
        i_face = seven_to_three_hypermat[tmp_idx[i], ]
        print(which(i_face == 1))
}
lines3d(x = pnt7to3.mat[1, c("class1029", "class968")], y = pnt7to3.mat[2, c("class1029", "class968")], z = pnt7to3.mat[3, c("class1029", "class968")],col="purple",lty=3, lwd = 5)
lines3d(x = pnt7to3.mat[1, c("class780", "class627")], y = pnt7to3.mat[2, c("class780", "class627")], z = pnt7to3.mat[3, c("class780", "class627")],col="purple",lty=3, lwd = 5)

# face with six nodes
tmp_idx = which(apply(seven_to_three_hypermat, 1, sum) == 6)
seven_to_three_hypermat[tmp_idx,]
for(i in 1:length(tmp_idx)){
        i_face = seven_to_three_hypermat[tmp_idx[i], ]
        print(which(i_face == 1))
}

tmp = pnt7to3.mat[, c("class1020", "class1044", "class861")]
lines3d(x = tmp[1,], y = tmp[2,], z=tmp[3,],col="purple",lty=3, lwd = 5)

tmp = pnt7to3.mat[, c("class1", "class1044")]
lines3d(x = tmp[1,], y = tmp[2,], z=tmp[3,],col="purple",lty=3, lwd = 5)

tmp = pnt7to3.mat[, c("class895","class1", "class51")]
lines3d(x = tmp[1,], y = tmp[2,], z=tmp[3,],col="purple",lty=3, lwd = 5)


rgl.snapshot("plots/threeD-seven-to-three-node.png")
#----------------------------------------------------#
# asym polytopes
#----------------------------------------------------#


open3d(zoom = 1, userMatrix = user_mat)

# resize window

par3d(windowRect = c(200, 200, 800, 800))
plot3d(c(1,0,0,1),c(0,1,0,0),c(0,0,1,0), ylab="%P3_bar",xlab="%P3",zlab="%K3",par3d(FOV=1)) 
## drow 3-simple
lines3d(x = c(1,0,0,1,0,0),y = c(0,1,0,0,0,0),z = c(0,0,1,0,0,1), col="dark grey",lwd=5)
lines3d(x=c(0,0), y = c(0,1) , z = c(0,0), col="dark grey",lwd=5)
points3d(x = c(1,0,0,0),y = c(0,1,0,0),z = c(0,0,1,0), col="red", size = 12)

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

# add the symmetric line
# lines3d(x = c(0.5, 0), y = c(0.5, 0), z = c(0, 0.5), col="green",lty=3)
# add Erdos Renyi line
p_vec = seq(0, 1, by = 0.01)
p_mat = edge_p(p_vec)
curve_cord = cord.mat %*% t(p_mat)

lines3d(x = curve_cord[1,], y = curve_cord[2,], z = curve_cord[3,], col="purple",lty=3)

for(i in 1: ncol(curve_cord)){
        lines3d(x = c(0, curve_cord[1, i], 0), y = c(0, curve_cord[2, i], 0), z = c(0, curve_cord[3, i], 1), col = "purple", lty = 3)
}
rgl.snapshot("plots/threeD-asym.png")
