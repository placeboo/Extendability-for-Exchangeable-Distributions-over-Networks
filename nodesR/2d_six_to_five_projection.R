rm(list = ls())

load(file = "data/extreme_points_six_node_graph.Rdata")


#------------------------------------------------------------------------------#
# do projection from six nodes down to five nodes
# keep mobius with meanings
# first two cols + "k6.bar", "k2_4k1", "p3_3k1", "twoK2_2k1", "p4_2k1", "claw_2k1", "x197", "k3_3k1", "c4_2k1", "k14_k1", "chair_k1", "co-dart_k1", "p5_k1", "k1_k2_k3", "x198", "w5.bar", "bull_k1", "cricket_k1", "P_k1", "k5_k1", "diamond_2k1", "co-fork_k1", "butterfly_k1", "co-4-fan", "k4_2k1", "dart_k1", "k23_k1","k5-e_k1","claw_k1.bar_k1", "p2_p3.bar_k1", "gem_k1", "w4_k1", "k3_2k1.bar_k1", "p3_2k1.bar_k1"
#------------------------------------------------------------------------------#

vmat.proj.q <- cbind(vmat.prob_mobius.q[, c(1,2)], vmat.prob_mobius.q[, c("k6.bar", "k2_4k1", "p3_3k1", "twoK2_2k1", "p4_2k1", "claw_2k1", "x197", "k3_3k1", "c4_2k1", "k14_k1", "chair_k1", "co-dart_k1", "p5_k1", "k1_k2_k3", "x198", "w5.bar", "bull_k1", "cricket_k1", "P_k1", "k5_k1", "diamond_2k1", "co-fork_k1", "butterfly_k1", "co-4-fan", "k4_2k1", "dart_k1", "k23_k1","k5-e_k1","claw_k1.bar_k1", "p2_p3.bar_k1", "gem_k1", "w4_k1", "k3_2k1.bar_k1", "p3_2k1.bar_k1")])

## remove redundant constrains
vmat.proj.minimal = redundant(vmat.proj.q, representation = "V")
vmat.proj.minimal.q = vmat.proj.minimal$output

rownames(vmat.proj.minimal.q) = rownames(vmat.proj.q)

# three is no redundant constrains after projection!

colnames(vmat.proj.minimal.q) = c(NA, NA, "5k1", "k2_3k1", "p3_2k1", "w4.bar", "co-gem", "claw_k1", "p2_p3", "k3_2k1", "co-butterfly", "k14", "fork", "co-dart", "p5", "k23.bar", "P.bar", "c5", "bull", "cricket", "P", "k5", "co-cricket", "co-fork", "butterfly", "p5.bar", "k14.bar", "dart", "k23","k2_3k1.bar","claw_k1.bar", "p2_p3.bar", "gem", "w4", "k3_2k1.bar", "p3_2k1.bar")

## save
save(vmat.proj.minimal.q, file = 'data/six_nodes_to_five_nodes_projection.Rdata')
write.csv(vmat.proj.minimal.q, "data/six_nodes_to_five_nodes_projection.csv")

#------------------------------------------------------------------------------#
# Compare
#------------------------------------------------------------------------------#
extre_pnt_five = get(load("data/extreme_points_five_node_graph.Rdata"))
six_to_five_extre_pnt = get(load("data/six_nodes_to_five_nodes_projection.Rdata"))

colnames(extre_pnt_five)
colnames(six_to_five_extre_pnt)

# search the same rows between five-node extreme points and six-node projection
## re-orginize
# search the same rows between four-node extreme points and five-node projection
## re-orginize 
extre_pnt_five[, 3: ncol(extre_pnt_five)] = extre_pnt_five[, colnames(six_to_five_extre_pnt)[3: ncol(six_to_five_extre_pnt)]]
colnames(extre_pnt_five) = colnames(six_to_five_extre_pnt)

# remove first two columns
extre_pnt_five2 = extre_pnt_five[, -c(1,2)]
six_to_five_extre_pnt2 = six_to_five_extre_pnt[, -c(1,2)]


## index.same:
## 0: there is no same extreme points after six to five projection
## none-zero: the index row of after-projection (six_to_five_extre_pnt) having the same extreme points
index.same = rep(0, nrow(extre_pnt_five2))
for(i in 1: nrow(extre_pnt_five2)){
      temp_extre_five = extre_pnt_five2[i, ]
      for(j in 1: nrow(six_to_five_extre_pnt2)){
            temp_proj = six_to_five_extre_pnt2[j, ]
            if(all(temp_extre_five == temp_proj)){
                  index.same[i] = j
            }
      }
}

## two-row table, the first row: index of four-node extreme points; the second row: the index of after-projection.
extend.tab = matrix(c(which(index.same!=0), index.same[which(index.same!=0)]), nrow = 2, byrow = T)


# save
five_node_kept = extre_pnt_five2[extend.tab[1,], ]
six_extend = six_to_five_extre_pnt2[extend.tab[2,], ]
write.csv(five_node_kept, "data/five_node_kept_from_six.csv")
write.csv(six_extend, "data/six_node_graphs_extendable.csv")

# graph kept when compared
graphKept = rownames(five_node_kept)

#------------------------------------------------------------------------------#
# convex combination of the left extreme points after projectioin
# 
# any point in the convex hull has a unique representation as a mixture of all the extreme points.
# for example, E: extreme poionts matrix, x: some vector in the convex hull. Then lambda = E^-1 %*% x, once the entry of lambda is zero, which means its coressponding extreme point not involving.
#------------------------------------------------------------------------------#

extre_pnt_five3 = q2d(extre_pnt_five2)
extre_five_name = colnames(extre_pnt_five3)
extre_five_rowname = rownames(extre_pnt_five3)

# find the subgraph-name for each row (extreme point vectors)

six_to_five_extre_pnt3 = q2d(six_to_five_extre_pnt2)

Lambdas = round(solve(t(extre_pnt_five3)) %*% t(six_to_five_extre_pnt3), 2)

# check lambdas
all(Lambdas > 0 || Lambdas == 0) # lambda >= 0
apply(Lambdas, 2, sum) # sum = 1

# make a table to describe the extreme points after projection (five-node to four-node). What are they, how many extreme points used from four-node, how the convex combination of four-node, are there "p4", "k4", "4k1"
combina.vec = rep(0,  ncol(Lambdas))
number_cmbn = rep(0,  ncol(Lambdas))
is_kept = rep(0,  ncol(Lambdas)) # binary, 1: yes, 0: no
graph_kept = rep(0,  ncol(Lambdas)) # what graph 
for(i in 1: ncol(Lambdas)){
        iLambda = Lambdas[, i]
        graphs = extre_five_rowname[iLambda != 0]
        is_kept[i] = ifelse(any(graphs %in% graphKept), 1, 0)
        graph_kept[i] = paste(graphs[graphs %in% graphKept], collapse = ", ")
        vec = paste(iLambda[iLambda != 0], graphs)
        combina.vec[i] = paste(vec, collapse = ", ")
        number_cmbn[i] = length(vec)
}
six_to_five_extre_cmbn = data.frame(number = number_cmbn, combination = combina.vec, isKept = is_kept, graphKept = graph_kept, six_to_five_extre_pnt2, check.names = F)

# save
save(six_to_five_extre_cmbn, file = "data/six_to_five_combination.Rdata")
write.csv(six_to_five_extre_cmbn, file = "data/six_to_five_combination.csv")
