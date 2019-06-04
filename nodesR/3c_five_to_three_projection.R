rm(list = ls())

load(file = "data/extreme_points_five_node_graph.Rdata")

#------------------------------------------------------------------------------#
# do projection from five nodes down to three nodes
# keep mobius with meanings
# first two cols + "5k1", "k5_e.bar", "p3_2k1", "k3_2k1"
#------------------------------------------------------------------------------#

vmat.proj.q <- cbind(vmat.prob_mobius.q[, c(1,2)], vmat.prob_mobius.q[, c("5k1", "k2_3k1", "p3_2k1", "k3_2k1")])

## remove redundant constrains
vmat.proj.minimal = redundant(vmat.proj.q, representation = "V")
vmat.proj.minimal.q = vmat.proj.minimal$output
# mark the original index kept
orig.index = c(1: nrow(vmat.proj.q))[-vmat.proj.minimal$redundant]
rownames(vmat.proj.minimal.q) = rownames(vmat.proj.q)[orig.index]

colnames(vmat.proj.minimal.q) = c(NA, NA, "3k1", "p3.bar", "p3", "k3")
Hmat = scdd(vmat.proj.minimal.q)$output
# number of ineq
length(which(Hmat[,1] == 0))


## save
save(vmat.proj.minimal.q, file = 'data/five_nodes_to_three_nodes_projection.Rdata')
write.csv(vmat.proj.minimal.q, file = 'data/five_nodes_to_three_nodes_projection.csv')

#------------------------------------------------------------------------------#
# Compare 
#------------------------------------------------------------------------------#

extre_pnt_three = get(load("data/extreme_points_three_node_graph.Rdata"))
five_to_three_extre_pnt = get(load("data/five_nodes_to_three_nodes_projection.Rdata"))
colnames(five_to_three_extre_pnt)
colnames(extre_pnt_three)

## re-orginize 
extre_pnt_three[, 3: ncol(extre_pnt_three)] = extre_pnt_three[, colnames(five_to_three_extre_pnt)[3: ncol(five_to_three_extre_pnt)]]
colnames(extre_pnt_three) = colnames(five_to_three_extre_pnt)

# remove first two columns
extre_pnt_three2 = extre_pnt_three[ , -c(1,2 )]
five_to_three_extre_pnt2 = five_to_three_extre_pnt[, -c(1,2)]
# search the same rows between three-node extreme points and five-node projection


## index.same: 
## 0: there is no same extreme points after five to three projection
## none-zero: the index row of after-projection (five_to_three_extre_pnt) having the same extreme points
index.same = rep(0, nrow(extre_pnt_three2))
for(i in 1: nrow(extre_pnt_three2)){
        temp_extre_three = extre_pnt_three2[i, ]
        for(j in 1: nrow(five_to_three_extre_pnt2)){
                temp_proj = five_to_three_extre_pnt2[j, ]
                if(all(temp_extre_three == temp_proj)){
                        index.same[i] = j
                }
        }
}

## two-row table, the first row: index of three-node extreme points; the second row: the index of after-projection. 
extend.tab = matrix(c(which(index.same!=0), orig.index[index.same[which(index.same!=0)]]), nrow = 2, byrow = T)

extre_pnt_three2[extend.tab[1,], ]
vmat.prob_mobius.q[extend.tab[2,], ]




