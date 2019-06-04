rm(list = ls())

# load extreme points of six-node graph
load(file = "data/extreme_points_six_node_graph.Rdata")

#------------------------------------------------------------------------------#
# do projection from six nodes down to four nodes
# keep mobius with meanings
# first two cols + "k6.bar", "k2_4k1", "p3_3k1", "twoK2_2k1", "p4_2k1", "claw_2k1", "k3_3k1", "c4_2k1", "co-dart_k1", "diamond_2k1", "k4_2k1"
#------------------------------------------------------------------------------#

vmat.proj.q <- cbind(vmat.prob_mobius.q[, c(1,2)], vmat.prob_mobius.q[, c("k6.bar", "k2_4k1", "p3_3k1", "twoK2_2k1", "p4_2k1", "claw_2k1", "k3_3k1", "c4_2k1", "co-dart_k1", "diamond_2k1", "k4_2k1")])

## remove redundant constrains
vmat.proj.minimal = redundant(vmat.proj.q, representation = "V")
vmat.proj.minimal.q = vmat.proj.minimal$output


# mark the original index kept
orig.index = c(1: nrow(vmat.proj.q))[-vmat.proj.minimal$redundant]


colnames(vmat.proj.minimal.q) = c(NA, NA, "4k1", "co-diamond", "co-paw", "c4.bar", "p4", "claw", "co-claw", "c4", "paw", "diamond", "k4")

## save
save(vmat.proj.minimal.q, file = 'data/six_nodes_to_four_nodes_projection.Rdata')
write.csv(vmat.proj.minimal.q[, -c(1,2)], "data/six_nodes_to_four_nodes_projection.csv")
#------------------------------------------------------------------------------#
# Compare
#------------------------------------------------------------------------------#
extre_pnt_four = read.csv("data/extreme_points_four_node_graph.csv", check.names = F, stringsAsFactors = F)[, -1]
six_to_four_extre_pnt = read.csv("data/six_nodes_to_four_nodes_projection.csv", check.names = F, stringsAsFactors = F)[, -1]
colnames(extre_pnt_four)
colnames(six_to_four_extre_pnt)

# search the same rows between four-node extreme points and six-node projection
## re-orginize
extre_pnt_four2 = extre_pnt_four[, colnames(six_to_four_extre_pnt)]

## index.same:
## 0: there is no same extreme points after six to four projection
## none-zero: the index row of after-projection (six_to_four_extre_pnt) having the same extreme points
index.same = rep(0, nrow(extre_pnt_four2))
for(i in 1: nrow(extre_pnt_four2)){
      temp_extre_four = extre_pnt_four2[i, ]
      for(j in 1: nrow(six_to_four_extre_pnt)){
            temp_proj = six_to_four_extre_pnt[j, ]
            if(all(temp_extre_four == temp_proj)){
                  index.same[i] = j
            }
      }
}

## two-row table, the first row: index of four-node extreme points; the second row: the index of after-projection.
extend.tab = matrix(c(which(index.same!=0), orig.index[index.same[which(index.same!=0)]]), nrow = 2, byrow = T)

