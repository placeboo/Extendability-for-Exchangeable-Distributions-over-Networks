rm(list = ls())

load("application_Rdata/samplk.Rdata")

#-------------------------------------------#
# samplk1
# social network is directed graph, which means there is possible (i,j) = 1 but (j,i)=0. Change directed graph to undirected, will treat these case as (i,j)=(j,i)=1
#-------------------------------------------#
summary(samplk1)
network.org = as.sociomatrix(samplk1)

# directed to undirected
network = ceiling((network.org + t(network.org))/2)
save(network, file = "application_Rdata/samplk1_undirected.Rdata")
# 
# combn(6,2)
# 
# groups = function(input.vec, leng){ 
#         # input.vec: the vector needed to do combinatioin
#         # the length of each combination
#         smallest = min(input.vec)
#         rst = input.vec[input.vec != smallest]
#         t(combn(rst, leng-1))
#         firstG = cbind(smallest, t(combn(rst, leng-1)))
#         restG = apply(firstG, 1, function(x) input.vec[!input.vec %in% x])
#         ls = list()
#         ls[[1]] = firstG
#         ls[[2]] = restG
#         return(ls)
# }

frontGroup = function(input.vec, leng){
        # input.vec: the vector needed to do combinatioin
        # the length of each combination
        smallest = min(input.vec)
        rst = input.vec[input.vec != smallest]
        t(combn(rst, leng-1))
        firstG = cbind(smallest, t(combn(rst, leng-1)))
        return(firstG)
}

# # take 1:6 choosing 2-pairs as example
# vec = 1:12
# leng = 3
# combn.mat = matrix(NA, ncol = length(vec))
# # first group
# ls1 = groups(vec, leng)
# group1 = ls1[[1]]
# rst1 = ls1[[2]]
# 
# # second group
# group2nd = list()
# rst2nd = list()
# for(n in 1: ncol(rst1)){
#         vec_n = rst1[, n]
#         group = groups(vec_n, leng)
#         group2nd[[n]] =group[[1]]
#         rst2nd[[n]] = group[[2]]
# }
# 
# # third group
# group3rd = list()
# rst3rd = list()
# for(i in 1: length(group2nd)){
#         rst_i = rst2nd[[i]]
#         l
#         for(j in 1: ncol(rst_i)){
#                 vec_j = rst_i[, j]
#                 group = groups(vec_j, leng)
#                 group3rd[[i]] = group[[1]]
#                 rst2nd[[i]][j] = group[[2]]
#         }
# }
# 


#-------------------------------------------#
# 1. group
# 2. based on group, testing whether it inside of 5->3, or outside; 
#-------------------------------------------#
# load the ineqs which define the pace of 5->3
load(file = "data/five_to_three_hmat.Rdata")

vec = 1:18
leng = 3
group1st = frontGroup(vec, leng)


Niter = 1
Ninsider = 0
saveEvery = 369600 # how often to write grouping results to disk
saveDir = "grouping/"

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
        
}
save(Ninsider, Niter, file = "looping_results.Rdata")

Ninsider_pct = Ninsider / (Niter - 1)

# calculate the coffienct of three-node network for each grouping.
Nfile = (Niter - 1)/saveEvery



#-------------------------------------------#
# another thought about how to estimate 
#-------------------------------------------#
perm1 = 1: 18  

network_tmp = network[perm1, perm1]

Combn_tmp = combn(perm1, 3)

subgraph_three_node.mat = apply(Combn_tmp, 2, function(x) as.vector(network_tmp[x, x]))

Nedges.vec = apply(subgraph_three_node.mat, 2, sum)/2
coef.mat = table(Nedges.vec)/sum(table(Nedges.vec))

rownames(coef.mat) = c("3k1","p3.bar", "p3", "k3")

all(-five_to_three.hmat[, -c(1,2)] %*% coef.mat[-1] - five_to_three.hmat[, 2] <= 0)


load(file = "data/six_to_three_hmat.Rdata") 
all(-six_to_three.hmat[, -c(1,2)] %*% coef.mat[-1] - six_to_three.hmat[, 2] <= 0)

