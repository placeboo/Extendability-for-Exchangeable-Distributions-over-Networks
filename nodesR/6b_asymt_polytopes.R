#library(extrafont)
#loadfonts(device = "win")

empty.vec = c(0, 0, 0)
p3_bar.vec = c(1, 0, 0)
p3.vec = c(0, 1, 0)
k3.vec = c(0, 0, 1)

cord.mat = cbind(p3.vec, empty.vec, p3_bar.vec, k3.vec)
rownames(cord.mat) = c("x", "y", "z")

open3d()

# resize window

par3d(windowRect = c(200, 200, 800, 800))

## drow 3-simple
lines3d(x = c(0, 0, 1, 0, 0), y = c(0, 0, 0, 1, 0), z = c(1, 0, 0, 0, 0), col="black",lty=3)
lines3d(x = c(0, 0), y = c(0, 1), z = c(1, 0), col="black",lty=3)
lines3d(x = c(0, 1), y = c(0, 0), z = c(1, 0), col="black",lty=3)

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

# plot 3-simplx
plot3d(x = cord.mat[1, ], y = cord.mat[2, ], z = cord.mat[3,], size = 8, add = T)
# add 3-node network extreme points
text3d(cord.mat[, 1], text = "P3", adj = 1.2)
text3d(cord.mat[, 2], text = "3K1", adj = 1.2)
text3d(cord.mat[, 3], text = "P3_bar", adj = 1.2)
text3d(cord.mat[, 4], text = "K3", adj = 1.2)

lines3d(x = curve_cord[1,], y = curve_cord[2,], z = curve_cord[3,], col="purple",lty=3)

for(i in 1: ncol(curve_cord)){
        lines3d(x = c(0, curve_cord[1, i], 0), y = c(0, curve_cord[2, i], 0), z = c(0, curve_cord[3, i], 1), col = "purple", lty = 3)
}


#----------------------------------------------#
# in order to calculate the volumn of the asym convex hull
# move ds
#----------------------------------------------#
open3d()

# resize window

par3d(windowRect = c(200, 200, 800, 800))

# plot 3-simplx
plot3d(x = cord.mat[1, ], y = cord.mat[2, ], z = cord.mat[3,], size = 8, add = T)
# # add 3-node network extreme points
# text3d(cord.mat[, 1], text = "p3", adj = 1.2)
# text3d(cord.mat[, 2], text = "3k1", adj = 1.2)
# text3d(cord.mat[, 3], text = "p3.bar", adj = 1.2)
# text3d(cord.mat[, 4], text = "k3", adj = 1.2)

## drow 3-simple
lines3d(x = c(0, 0, 1, 0, 0), y = c(0, 0, 0, 1, 0), z = c(1, 0, 0, 0, 0), col="black",lty=3)
lines3d(x = c(0, 0), y = c(0, 1), z = c(1, 0), col="black",lty=3)
lines3d(x = c(0, 1), y = c(0, 0), z = c(1, 0), col="black",lty=3)



# add the symmetric line
# lines3d(x = c(0.5, 0), y = c(0.5, 0), z = c(0, 0.5), col="green",lty=3)
# add Erdos Renyi line
p_vec = seq(0, 1, by = 0.01)
p_mat = edge_p(p_vec)
curve_cord = cord.mat %*% t(p_mat)

lines3d(x = curve_cord[1,], y = curve_cord[2,], z = curve_cord[3,], col="purple",lty=3)

lines3d(x = c(0, curve_cord[1, 29], 0), y = c(0, curve_cord[2, 29], 0), z = c(0, curve_cord[3, 29], 1), col = "blue", lty = 3)
lines3d(x = c(0, curve_cord[1, 33], 0), y = c(0, curve_cord[2, 33], 0), z = c(0, curve_cord[3, 31], 1), col = "blue", lty = 3)

lines3d(x = c(curve_cord[1, 29], curve_cord[1, 33]), y = c(curve_cord[2, 29], curve_cord[2, 33]), z = c(curve_cord[3, 29], curve_cord[3, 33]), col = "blue", lty = 3)

plot3d(x = c(curve_cord[1, 29], curve_cord[1, 33]), y = c(curve_cord[2, 29], curve_cord[2, 33]), z = c(curve_cord[3, 29], curve_cord[3, 33]), size = 8, add = T)
# mark the location
#text3d(cord.mat[, 4], text = "expression(x_1)=(0,0,1)", adj = 1.2)
# save
browseURL(paste("file://", writeWebGL(dir=file.path("plots/", "webGL"), width=700), sep=""))
#----------------------------------------------#
# calculate the volumn 
#----------------------------------------------#
volumn = function(n){
        #n: partition
        dp = 1/n
        p_seq = seq(0, 1, by = dp)
        
        tri_area = 1/2 * sqrt(9*p_seq^2*(1-p_seq)^2*(1-2*p_seq+2*p_seq^2))
        
        # height of the tetre
        h = (3*p_seq - 3*p_seq^2) / sqrt(2 * p_seq^2 - 2*p_seq + 1) * dp
        
        # volumn
        return(sum(1/3 * tri_area * h))
}
volumn(100)
volumn(1000)

#----------------------------------------------#
# convex ???
# reture the coefficient after convex combination
# n: the number of points, >1
#----------------------------------------------#
asy_posi = function(n){
        # simulate beta
        Beta = matrix(runif(n * 3), ncol = 3)
        # normalize
        Beta = Beta / apply(Beta, 1, sum)
        
        # simulate alpha
        Alpha = runif(n)
        # normalize
        Alpha = Alpha / sum(Alpha)
        
        # simulate probability parameters
        P = runif(n)
        P1 = 3 * P^2 * (1 - P)
        P2 = 3 * P * (1 - P)^2
        P3 = P^3
        
        # the point after covex combination
        p.star = sum(Alpha * Beta[, 3] * P1) / sum(Alpha * Beta[, 3] * (P1 + P2))
        # coef
        rho3 = sum(Alpha * Beta[, 3] * (P1 + P2))^3 / (3 * sum(Alpha * Beta[, 3] * P1 * sum(Alpha * Beta[, 3] * P2)))
        # q1 = 3 * p.star^2 * (1 - p.star)
        # q2 = 3 * p.star * (1 - p.star)^2
        # q3 = p.star^3
        # 
        #rho3 = sum(Alpha * Beta[,3] * P1) / q1
        rho2 = sum(Alpha * Beta[, 3] * P3) + sum(Alpha * Beta[,2]) - sum(Alpha * Beta[, 3] * P1)^2/(3 * sum(Alpha * Beta[, 3] * P2))
        # rho2 = sum(Alpha * Beta[,3] * P3) + sum(Alpha * Beta[,2]) - rho3 * q3
        
        rst = paste("p=",p.star,"\n", "rho2=",rho2, "\n", "rho3=",rho3, "\n", "rho2+rho3=", rho2 + rho3)
        return(rst)
}
set.seed(17)
asy_posi(2)
asy_posi(3)
asy_posi(8)
asy_posi(100)
asy_posi(179)

#---------------------#
# tangent plane
#---------------------#
# some plots
open3d()

# resize window

par3d(windowRect = c(200, 200, 800, 800))

# plot 3-simplx
plot3d(x = cord.mat[1, ], y = cord.mat[2, ], z = cord.mat[3,], size = 8, add = T)
## drow 3-simple
lines3d(x = c(0, 0, 1, 0, 0), y = c(0, 0, 0, 1, 0), z = c(1, 0, 0, 0, 0), col="black",lty=3)
lines3d(x = c(0, 0), y = c(0, 1), z = c(1, 0), col="black",lty=3)
lines3d(x = c(0, 1), y = c(0, 0), z = c(1, 0), col="black",lty=3)
lines3d(x = curve_cord[1,], y = curve_cord[2,], z = curve_cord[3,], col="purple",lty=3)

poi = 80
lines3d(x = c(0, curve_cord[1, poi], 0), y = c(0, curve_cord[2, poi], 0), z = c(0, curve_cord[3, poi], 1), col = "blue", lty = 3)

plot3d(x = c(curve_cord[1, poi]), y = c(curve_cord[2, poi]), z = c(curve_cord[3, poi]), size = 8, add = T)

# mark the location
#text3d(cord.mat[, 4], text = "expression(x_1)=(0,0,1)", adj = 1.2)
# save
browseURL(paste("file://", writeWebGL(dir=file.path("plots/", "webGL"), width=700), sep=""))
# dist_tangent = function(p){
#         q = p/2
# 
#         x = 3 * q * (1-q)^2
#         y = 3 * q^2 * (1-q)
#         z = q^3
# 
#         g = p^2 * x + 2*(p-1) * p * y + 3 * (p-1)^2 * z
# 
#         return(g)
# }
# 
# p.vec = seq(0,1,0.01)
# dist.vec = rep(0, length(p.vec))
# 
# for(i in 1: length(q.vec)){
#         p_i = p.vec[i]
#         dist.vec[i] = dist_tangent(p_i)
# }
# dist.df = data.frame(p = q.vec, dist = dist.vec)
# 
# ggplot(dist.df, aes(x = p, y = dist)) + geom_point()

#-------------------------#
# draw x-y plane
#-------------------------#
open3d()

# resize window

par3d(windowRect = c(200, 200, 800, 800))
# plot 3-simplx
#plot3d(x = cord.mat[1, ], y = cord.mat[2, ], z = cord.mat[3,], size = 8, add = T)
## drow 3-simple
#lines3d(x = c(0, 0, 1, 0, 0), y = c(0, 0, 0, 1, 0), z = c(1, 0, 0, 0, 0), col="black",lty=3)
#lines3d(x = curve_cord[1,], y = curve_cord[2,], z = curve_cord[3,], col="purple",lty=3)
p_vec = seq(0, 1, by = 0.001)
p_mat = edge_p(p_vec)
curve_cord = cord.mat %*% t(p_mat)

lines3d(x = curve_cord[1,], y = curve_cord[2,], z = curve_cord[3,], col="purple",lty=3)

for(i in 1: ncol(curve_cord)){
        lines3d(x = c(0, curve_cord[1, i], 0), y = c(0, curve_cord[2, i], 0), z = c(0, curve_cord[3, i], 1), col = "purple", lty = 3)
}

lines3d(x = c(0,1), y = c(0,0), z = c(0,0), col="black",lty=3)
lines3d(x = c(0,0), y = c(0,1), z = c(0,0), col="black",lty=3)
rgl.viewpoint(theta = 0, phi = 0, fov = 0, zoom = 1, interactive = TRUE)
rgl.postscript("xy_plane.pdf","pdf")
