rm(list = ls())

set.seed(100)

load(file = "samplk1_undirected.Rdata")
# load the ineqs which define the pace of 5->3
load(file = "five_to_three_hmat.Rdata")

Nsample = 1000000
coef.mat =  matrix(NA, ncol = 4, nrow = Nsample)
Ninsider = 0

for(i in 1: Nsample){
        groups = sample(18, replace = F)
        mat1 = network[groups[1:3], groups[1:3]] 
        mat2 = network[groups[4:6], groups[4:6]] 
        mat3 = network[groups[7:9], groups[7:9]] 
        mat4 = network[groups[10:12], groups[10:12]] 
        mat5 = network[groups[13:15], groups[13:15]] 
        mat6 = network[groups[16:18], groups[16:18]] 
        
        three_pnt_loc = table(c(sum(mat1), sum(mat2), sum(mat3), sum(mat4), sum(mat5), sum(mat6))/2)/6
        three_pnt_loc_name = dimnames(three_pnt_loc)[[1]]
        
        if("0" %in% three_pnt_loc_name){
                three_pnt_loc = three_pnt_loc[-1]
        }
        
        # the length of three_pnt_loc should be three
        if(length(three_pnt_loc) < 3){
                temp.vec = rep(0, 3)
                temp.vec[as.numeric(dimnames(three_pnt_loc)[[1]])] = three_pnt_loc
                three_pnt_loc = temp.vec # update
        }
        
        coef.mat[i, ] = c(1 - sum(three_pnt_loc), three_pnt_loc)
        
        # test whether the point is inside the 5->3 polytope
        inside.vec = ifelse(all(-five_to_three.hmat[, -c(1,2)] %*% three_pnt_loc - five_to_three.hmat[, 2] <= 0), 1, 0)
        if(inside.vec == 1){
                Ninsider = Ninsider + 1
        }
        print(i)
}

Ninsider/Nsample

save(coef.mat, Ninsider, file = paste("application_Rdata/", "Monk_Sample_", Nsample, ".Rdata", sep = ''))
