
#----------------------------------------------------------------#
# manage the location of data, prepare for graph
#----------------------------------------------------------------#
empty.vec = c(1, 0, 0)
p3_bar.vec = c(0, 1, 0)
p3.vec = c(0, 0, 1)
k3.vec = c(1, 1, 1) # dual of (0,0,0,1)

cord.mat = cbind(p3.vec, empty.vec, p3_bar.vec, k3.vec)

pnt5to3.mat = cord.mat %*% Lambdas5to3
pnt6to3.mat = cord.mat %*% Lambdas6to3


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
lines3d(x = c(1, 0, 0, 1, 1), y = c(0, 1, 0, 1, 0), z = c(0, 0, 1, 1, 0), col="black",lty=3)
lines3d(x = c(0, 1), y = c(1, 1), z = c(0, 1), col="black",lty=3)
lines3d(x = c(0, 1), y = c(0, 0), z = c(1, 0), col="black",lty=3)

rgl.postscript("3_simplex.pdf","pdf")

points3d(x = pnt5to3.mat[1, ], y = pnt5to3.mat[2, ], z = pnt5to3.mat[3,], col = "blue", size = 8)

points3d(x = pnt6to3.mat[1, ], y = pnt6to3.mat[2, ], z = pnt6to3.mat[3,], col = "red", size = 15)

# add legend
legend3d("topright", legend = c("5-node network to 3-node", "6-node network to 3-node"), pch = 16, col = c("blue", "red"), cex=2, inset=c(0.02))

# what are those points?
for(i in 1: ncol(pnt5to3.mat)){
        text3d(as.numeric(pnt5to3.mat[,i]), texts = colnames(pnt5to3.mat)[i], adj = 1.2)
}

for(i in 1: ncol(pnt6to3.mat)){
        text3d(as.numeric(pnt6to3.mat[,i]), texts = colnames(pnt6to3.mat)[i], adj = -0.2)
}

six_to_three_combn['c6', 'combination'] 
five_to_three_combn[c('c5', 'p2_p3'), 'combination'] # c6 is in the middle of c5 and p2_p3


rgl.postscript("persp3dd.pdf","pdf")


pair.index = cbind(c(1, 3, 5, 7, 8), c(9, 4, 6, 7, 2))
left.mat = pnt5to3.mat[, pair.index[,1]]
right.mat = pnt5to3.mat[, pair.index[,2]]

# middle
midde.mat = 1/2 * (left.mat + right.mat)

v.mat = cbind(rep(0, 5), rep(1,5), t(midde.mat))
v.mat = d2q(v.mat)

redundant(v.mat, representation = "V") 


