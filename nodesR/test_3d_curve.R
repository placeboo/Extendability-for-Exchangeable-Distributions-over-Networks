empty.vec = c(0, 0, 0)
p3_bar.vec = c(1, 0, 0)
p3.vec = c(0, 1, 0)
k3.vec = c(0, 0, 1)

cord.mat = cbind(empty.vec, p3_bar.vec, p3.vec,  k3.vec)

edge_p = function(p_vec){
        # one egde probability p, vec
        p_mat = matrix(0, ncol = 4, nrow = length(p_vec))
        for (i in 1: length(p_vec)){
                p_i = p_vec[i]
                p_mat[i, ] = c((1-p_i)^3, 3* p_i*(1-p_i)^2, 3*p_i^2*(1-p_i), p_i^3)
        }
        colnames(p_mat) = c("3k1", "p3.bar", "p3", "k3")
        return(p_mat)
}
no_edge_p = function(p_vec){
        # one egde probability p, vec
        p_mat = matrix(0, ncol = 4, nrow = length(p_vec))
        for (i in 1: length(p_vec)){
                p_i = p_vec[i]
                p_mat[i, ] = c(p_i^3, 3*p_i^2*(1-p_i), 3*p_i*(1-p_i)^2, (1-p_i)^3)
        }
        colnames(p_mat) = c("3k1", "p3.bar", "p3", "k3")
        return(p_mat)
}
p_vec = seq(0, 1, by = 0.01)

p_mat = edge_p(p_vec)
p1_mat = no_edge_p(p_vec)

# change cord
curve_cord = cord.mat %*% t(p_mat)
curve1_cord = cord.mat %*% t(p1_mat)

open3d()

# resize window

par3d(windowRect = c(200, 200, 800, 800))

# plot 3-simplx
plot3d(x = cord.mat[1, ], y = cord.mat[2, ], z = cord.mat[3,], size = 8, add = T)
# add 3-node network extreme points
text3d(cord.mat[, 3], text = "p3", adj = 1.2)
text3d(cord.mat[, 1], text = "3k1", adj = 1.2)
text3d(cord.mat[, 2], text = "p3.bar", adj = 1.2)
text3d(cord.mat[, 4], text = "k3", adj = 1.2)
## drow 3-simple
lines3d(x = c(0, 0, 1, 0, 0), y = c(0, 0, 0, 1, 0), z = c(1, 0, 0, 0, 0), col="black",lty=3)
lines3d(x = c(0, 0), y = c(0, 1), z = c(1, 0), col="black",lty=3)
lines3d(x = c(0, 1), y = c(0, 0), z = c(1, 0), col="black",lty=3)

# add edge
lines3d(x = curve_cord[1,], y = curve_cord[2,], z = curve_cord[3,], col="black",lty=6)
# remove edge
lines3d(x = curve1_cord[1,], y = curve1_cord[2,], z = curve1_cord[3,], col="red",lty=3)
