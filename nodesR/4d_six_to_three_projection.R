rm(list = ls())

# load extreme points of six-node graph
load(file = "data/extreme_points_six_node_graph.Rdata")

#------------------------------------------------------------------------------#
# do projection from six nodes down to three nodes
# keep mobius with meanings
# first two cols + "k6.bar", "k2_4k1", "p3_3k1", "k3_3k1"
#------------------------------------------------------------------------------#

vmat.proj.q <- cbind(vmat.prob_mobius.q[, c(1,2)], vmat.prob_mobius.q[, c("k6.bar", "k2_4k1", "p3_3k1", "k3_3k1")]) 

## remove redundant constrains
vmat.proj.minimal = redundant(vmat.proj.q, representation = "V")
vmat.proj.minimal.q = vmat.proj.minimal$output

Hmat = scdd(vmat.proj.minimal.q)$output
# number of ineq
length(which(Hmat[,1] == 0))


# mark the original index kept
orig.index = c(1: nrow(vmat.proj.q))[-vmat.proj.minimal$redundant]

colnames(vmat.proj.minimal.q) = c(NA, NA, "3k1", "p3.bar", "p3", "k3")
rownames(vmat.proj.minimal.q) = rownames(vmat.proj.q)[orig.index]

## save
save(vmat.proj.minimal.q, file = 'data/six_nodes_to_three_nodes_projection.Rdata')
write.csv(vmat.proj.minimal.q[, -c(1,2)], "data/six_nodes_to_three_nodes_projection.csv")
#------------------------------------------------------------------------------#
# Compare
#------------------------------------------------------------------------------#
extre_pnt_three = get(load("data/extreme_points_three_node_graph.Rdata"))[, -c(1,2)]
six_to_three_extre_pnt = get(load("data/six_nodes_to_three_nodes_projection.Rdata"))[, -c(1,2)]
colnames(extre_pnt_three)
colnames(six_to_three_extre_pnt)

# search the same rows between four-node extreme points and six-node projection
## re-orginize

extre_pnt_three2 = extre_pnt_three[, colnames(six_to_three_extre_pnt)]

## index.same:
## 0: there is no same extreme points after six to four projection
## none-zero: the index row of after-projection (six_to_three_extre_pnt) having the same extreme points
index.same = rep(0, nrow(extre_pnt_three2))
for(i in 1: nrow(extre_pnt_three2)){
      temp_extre_three = extre_pnt_three2[i, ]
      for(j in 1: nrow(six_to_three_extre_pnt)){
            temp_proj = six_to_three_extre_pnt[j, ]
            if(all(temp_extre_three == temp_proj)){
                  index.same[i] = j
            }
      }
}

## two-row table, the first row: index of four-node extreme points; the second row: the index of after-projection.
extend.tab = matrix(c(which(index.same!=0), orig.index[index.same[which(index.same!=0)]]), nrow = 2, byrow = T)

