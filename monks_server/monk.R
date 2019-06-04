rm(list = ls())

load(file = "samplk1_undirected.Rdata")
# load the ineqs which define the pace of 5->3
load(file = "five_to_three_hmat.Rdata")

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
saveEvery = 369600 # how often to write grouping results to disk
saveDir = "grouping2/"

# grouping are saved for future
group.mat = matrix(NA, ncol = length(vec), nrow = saveEvery)
# indicator, whether the point inside the 5->3 polytope
inside.vec = rep(NA, saveEvery)
save_counter = 1

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
                                        # what are they
                                        three_pnt_loc_name = dimnames(three_pnt_loc)[[1]]
                                        # if there are portion of 0, representing 3k1. Remove it
                                        if(length(three_pnt_loc_name) == 4){ # with 3k1
                                                three_pnt_loc = three_pnt_loc[-1]
                                        }
                                        # the length of three_pnt_loc should be three
                                        if(length(three_pnt_loc_name) < 3){
                                                temp.vec = rep(0, 3)
                                                temp.vec[as.numeric(three_pnt_loc_name)] = three_pnt_loc
                                                three_pnt_loc = temp.vec # update
                                        }
                                        
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
                                                filename1 = paste(saveDir, "group_results_", Niter, ".Rdata", sep = "")
                                                save(group.mat, file = filename1)
                                                
                                                filename2 = paste(saveDir, "if_insider_results_", Niter, ".Rdata", sep = "")
                                                save(inside.vec, file = filename2)
                                                # reset store_counter
                                                save_counter = 1
                                        }
                                        Niter = Niter + 1
                                }
                                
                        }
                }
        }
        save(Niter, Ninsider, file = "looping_results.Rdata")
}

190590400

# now
190344000

190590400-190344000

246400

91 * 55 * 28*10

1401400/246400

1401400 - 246400 * 5

169400/15400




