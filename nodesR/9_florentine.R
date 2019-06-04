rm(list = ls())

set.seed(100)

load(file = "application_Rdata/florentine_business.Rdata")

# load the ineqs which define the pace of 5->3
load(file = "five_to_three_hmat.Rdata")

Nsample = 1000000
coef.mat =  matrix(NA, ncol = 4, nrow = Nsample)
Ninsider = 0

for(i in 1: Nsample){
        groups = sample(16, replace = F)
        mat1 = flobusiness[groups[1:3], groups[1:3]] 
        mat2 = flobusiness[groups[4:6], groups[4:6]] 
        mat3 = flobusiness[groups[7:9], groups[7:9]] 
        mat4 = flobusiness[groups[10:12], groups[10:12]] 
        mat5 = flobusiness[groups[13:15], groups[13:15]] 
        
        
        three_pnt_loc = table(c(sum(mat1), sum(mat2), sum(mat3), sum(mat4), sum(mat5))/2)/5
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

save(coef.mat, Ninsider, file = paste("application/", "florentine_business_Sample_", Nsample, ".Rdata", sep = ''))
