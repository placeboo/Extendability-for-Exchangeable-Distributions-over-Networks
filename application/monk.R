# ##First read in the arguments listed at the command line
# args=(commandArgs(TRUE))
# 
# ##args is now a list of character vectors
# ## First check to see if arguments are passed.
# ## Then cycle through each element of the list and evaluate the expressions.
# if(length(args)==0){
#         print("No arguments supplied.")
#         ##supply default value
# }else{
#         for(i in 1:length(args)){
#                 eval(parse(text=args[[i]]))
#         }
# }


#######
# change Nov 30
# revise from monk 1]
# I realized I was wrong to calculate the effiecient

rm(list = ls())

load(file = "application_Rdata/samplk1_undirected.Rdata")
# load the ineqs which define the pace of 5->3
load(file = "data/five_to_three_hmat.Rdata")

frontGroup = function(input.vec, leng){
        # input.vec: the vector needed to do combinatioin
        # the length of each combination
        smallest = min(input.vec)
        rst = input.vec[input.vec != smallest]
        t(combn(rst, leng-1))
        firstG = cbind(smallest, t(combn(rst, leng-1)))
        return(firstG)
}

vec = 1:18
leng = 3
group1st = frontGroup(vec, leng)

Niter = 0
Ninsider = 0
saveEvery = 1401400 # how often to write grouping results to disk
saveDir = "group2/"
save_counter = 1

# grouping are saved for future
group.mat = matrix(NA, ncol = length(vec), nrow = saveEvery)
# indicator, whether the point inside the 5->3 polytope
inside.vec = rep(NA, saveEvery)
# coeffiecent of three-node network
coef.mat = matrix(NA, ncol = 4, nrow = saveEvery)

for(i in 1: nrow(group1st)){
        G1_i = group1st[i,]
        G1_i_rst = vec[!vec %in% G1_i]
        
        # second group
        G2_i = frontGroup(G1_i_rst, leng)
        for(j in 1: nrow(G2_i)){
                G2_ij = G2_i[j, ]
                G2_ij_rst = G1_i_rst[!G1_i_rst %in% G2_ij]
                
                # third group
                G3_ij = frontGroup(G2_ij_rst, leng)
                for(k in 1: nrow(G3_ij)){
                        G3_ijk = G3_ij[k, ]
                        G3_ijk_rst = G2_ij_rst[!G2_ij_rst %in% G3_ijk]
                        
                        # forth group
                        G4_ijk = frontGroup(G3_ijk_rst, leng)
                        for(l in 1: nrow(G4_ijk)){
                                G4_ijkl = G4_ijk[l, ]
                                G4_ijkl_rst = G3_ijk_rst[!G3_ijk_rst %in% G4_ijkl]
                                
                                # fifth group
                                G5_ijkl = frontGroup(G4_ijkl_rst, leng)
                                for(m in 1: nrow(G5_ijkl)){
                                        Niter = Niter + 1
                                        
                                        G5_ijklm = G5_ijkl[m, ]
                                        # the rest is sixth group
                                        G6_ijklmn = G4_ijkl_rst[!G4_ijkl_rst %in% G5_ijklm]
                                        
                                        # six groups, G1_i, G2_ij, G3_ijk, G4_ijkl, G5_ijklm, G6_ijklmn
                                        grouping = c(G1_i, G2_ij, G3_ijk, G4_ijkl, G5_ijklm, G6_ijklmn)
                                        #-----------------------#
                                        # testing the place 
                                        # sum(mat)/2: 0: empty; 1: p3.bar; 2:p3; 3: k3
                                        #-----------------------#
                                        mat1 = network[G1_i, G1_i]
                                        mat2 = network[G2_ij, G2_ij]
                                        mat3 = network[G3_ijk, G3_ijk]
                                        mat4 = network[G4_ijkl, G4_ijkl]
                                        mat5 = network[G5_ijklm, G5_ijklm]
                                        mat6 = network[G6_ijklmn,G6_ijklmn]
                                        
                                        three_pnt_loc = table(c(sum(mat1), sum(mat2), sum(mat3), sum(mat4), sum(mat5), sum(mat6))/2)/6
                                        three_pnt_loc_name = dimnames(three_pnt_loc)[[1]]
                                        # if there are portion of 0, representing 3k1. Remove it
                                        if("0" %in% three_pnt_loc_name){
                                                three_pnt_loc = three_pnt_loc[-1]
                                        }
                                        
                                        # the length of three_pnt_loc should be three
                                        if(length(three_pnt_loc) < 3){
                                                temp.vec = rep(0, 3)
                                                temp.vec[as.numeric(dimnames(three_pnt_loc)[[1]])] = three_pnt_loc
                                                three_pnt_loc = temp.vec # update
                                        }
                                        
                                        coef.mat[save_counter, ] = c(1 - sum(three_pnt_loc), three_pnt_loc)
                                        
                                        # test whether the point is inside the 5->3 polytope
                                        inside.vec[save_counter] = ifelse(all(-five_to_three.hmat[, -c(1,2)] %*% three_pnt_loc - five_to_three.hmat[, 2] <= 0), 1, 0)
                                        if(inside.vec[save_counter] == 1){
                                                Ninsider = Ninsider + 1
                                        }
                                        
                                        # store in matrix
                                        group.mat[save_counter, ] = grouping
                                        
                                        save_counter = save_counter + 1
                                        
                                        # check whether the grouping is correct
                                        if(!sum(order(grouping - vec))){
                                                stop(paste("Grouping is wrong! ", "i =", i, ";", "j =", j, ";", "k =", k, ";", "l =",l, ";", "m =", m))
                                        }
                                        if(Niter %% saveEvery == 0){
                                                filename1 = paste(saveDir, "group_i",i,"_", Niter, ".Rdata", sep = "")
                                                save(group.mat, file = filename1)
                                                
                                                filename2 = paste(saveDir, "if_insider_i",i,"_", Niter, ".Rdata", sep = "")
                                                save(inside.vec, file = filename2)
                                                
                                                save(Niter, Ninsider, file = paste(saveDir,"looping_i",i,"_", Niter, ".Rdata", sep = ""))
                                                
                                                save(coef.mat, file = paste(saveDir,"ceof_i",i,"_", Niter, ".Rdata", sep = ""))
                                                # reset store_counter
                                                save_counter = 1
                                        }
                                        
                                }
                                
                        }
                }
        }
}



#-----------------------------------------------#
# load back the ceof and number of insider
#-----------------------------------------------#

# number of inside
N = saveEvery * 20
Ninsider_perRound = rep(NA, 7)
Niter_perRound = rep(NA, 7)
for(n in 1:6){
        i = 20 * n
        filename = paste("monks_server/group2/", "looping_i",i,"_", N, ".Rdata", sep = "")
        load(file = filename)
        Ninsider_perRound[n] = Ninsider
        Niter_perRound[n] = Niter
}
# last file
load(file = "monks_server/group2/looping_i136_22422400.Rdata")
Ninsider_perRound[7] = Ninsider
Niter_perRound[7] = Niter

sum(Niter_perRound)
sum(Ninsider_perRound)/sum(Niter_perRound)

# coef 
sum.mat = matrix(NA, ncol = 4)
num = 0 
counter = 1
for(i in 1: 136){ 
        N = counter * saveEvery
        filename =  paste("monks_server/group2/", "ceof_i",i,"_", N, ".Rdata", sep = "")
        load(file = filename)
        sum.tmp = apply(coef.mat, 2, sum)
        sum.mat = rbind(sum.mat, sum.tmp)
        num = num + nrow(coef.mat)
        counter = counter + 1
        print(i)
        if(i %% 20 == 0){
                counter = 1
        }
}
# remove NA
sum.mat = sum.mat[-1, ]

# mean
coef_mean = apply(sum.mat, 2, sum)/sum(Niter_perRound)

all(-five_to_three.hmat[, -c(1,2)] %*% coef_mean[-1] - five_to_three.hmat[, 2] <= 0)
all(-six_to_three.hmat[, -c(1,2)] %*% coef_mean[-1] - six_to_three.hmat[, 2] <= 0)
# variance
sq_sum.mat = matrix(NA, ncol = 4)
num = 0 
counter = 1
for(i in 1: 136){ 
        N = counter * saveEvery
        filename =  paste("monks_server/group2/", "ceof_i",i,"_", N, ".Rdata", sep = "")
        load(file = filename)
        sq_sum.tmp = apply((coef.mat - coef_mean)^2, 2, sum)
        sq_sum.mat = rbind(sq_sum.mat, sq_sum.tmp)
        counter = counter + 1
        print(i)
        if(i %% 20 == 0){
                counter = 1
        }
}

sq_sum2.mat = sq_sum.mat[-1, ]
# sd
coef_sd = sqrt(apply(sq_sum2.mat,2,sum)/sum(Niter_perRound))


# table
ceof.tab = rbind(coef_mean, coef_sd)
colnames(ceof.tab) = c("3K_1", "P3.bar", "P_3", "K_3")
rownames(ceof.tab) = c("Mean", "Std.")

print(xtable(ceof.tab,  caption = "Coefficients of Monk Data Set in Three-node Network", digits = 3, align = rep("c", ncol(ceof.tab)+1)))
