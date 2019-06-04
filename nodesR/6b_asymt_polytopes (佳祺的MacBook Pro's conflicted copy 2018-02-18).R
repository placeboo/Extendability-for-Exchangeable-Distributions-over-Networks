library(extrafont)
loadfonts(device = "win")

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
text3d(cord.mat[, 4], text = "expression(x_1)=(0,0,1)", adj = 1.2)

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

asy_posi(1709)


#----------------------------------------------#
# convex ???
# use the geomtry method
#----------------------------------------------#
plane_fct = function(p){
        q = (9 * p^2 + sqrt(24*p^4 - 50*p^3 + 36 * p^2 - 10 *p + 1) - 10 * p + 2) / (19 * p^2 - 18*p + 3)
        
        # location on the curve
        x = 3 * q * (1 - q)^2
        y = 3 * q^2 * (1-q)
        z = q^3
        
        # number on f
        f = (1 - 4*p + 3*p^2) * x + (2*p - 3*p^2) * y + p^2 * z - 3 * p* (1-p)^2 * (1 - 4*p + 3 *p^2) - 3*p^3*(1-p)*(2-3*p) - p^5 - 3 * p^5
        
        return(f)
}

p = seq(0, 1, 0.01)
f = rep(0, length(p))

for(i in 1: length(p)){
        p_i = p[i]
        f[i] = plane_fct(p_i) 
}

plane.mat = data.frame(p, f)

ggplot(plane.mat, aes(x = p, y = f)) + geom_line()
