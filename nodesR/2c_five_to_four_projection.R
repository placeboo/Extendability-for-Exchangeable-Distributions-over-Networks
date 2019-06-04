rm(list = ls())

load(file = "data/extreme_points_five_node_graph.Rdata")
#------------------------------------------------------------------------------#
# do projection from five nodes down to four nodes
# keep mobius with meanings
# first two cols + "5k1", "k5_e.bar", "p3_2k1", "w4.bar", "claw_k1", "co_gem", "k3_2k1", "k14.bar","co-butterfly", "co-dart", 'co-cricket'
#------------------------------------------------------------------------------#

vmat.proj.q <- cbind(vmat.prob_mobius.q[, c(1,2)], vmat.prob_mobius.q[, c("5k1", "k2_3k1", "p3_2k1", "w4.bar", "claw_k1", "co-gem", "k3_2k1", "k14.bar","co-butterfly", "co-dart", 'co-cricket')])

## remove redundant constrains
vmat.proj.minimal = redundant(vmat.proj.q, representation = "V")
vmat.proj.minimal.q = vmat.proj.minimal$output

hmat = scdd(vmat.proj.minimal.q)$out
length(which(hmat[, 1]==0))

rownames(vmat.proj.minimal.q) = rownames(vmat.proj.q)

# three is no redundant constrains after projection!

colnames(vmat.proj.minimal.q) = c(NA, NA, "4k1", "co-diamond", "co-paw", "c4.bar", "claw", "p4", "co-claw", "k4","c4", "paw", 'diamond')

hmat = scdd(vmat.proj.minimal.q)$out
length(which(hmat[, 1]==0))



## save
save(vmat.proj.minimal.q, file = 'data/five_nodes_to_four_nodes_projection.Rdata')
write.csv(vmat.proj.minimal.q, "data/five_nodes_to_four_nodes_projection.csv")
#------------------------------------------------------------------------------#
# Compare
#------------------------------------------------------------------------------#
View(vmat.proj.minimal.q)

extre_pnt_four = get(load("data/extreme_points_four_node_graph.Rdata"))

five_to_four_extre_pnt = get(load("data/five_nodes_to_four_nodes_projection.Rdata"))
                           
# search the same rows between four-node extreme points and five-node projection
## re-orginize 
extre_pnt_four[, 3: ncol(extre_pnt_four)] = extre_pnt_four[, colnames(five_to_four_extre_pnt)[3: ncol(five_to_four_extre_pnt)]]
colnames(extre_pnt_four) = colnames(five_to_four_extre_pnt)

# remove first two columns
extre_pnt_four2 = extre_pnt_four[, -c(1,2)]
five_to_four_extre_pnt2 = five_to_four_extre_pnt[, -c(1,2)]
## index.same:
## 0: there is no same extreme points after five to four projection
## none-zero: the index row of after-projection (five_to_four_extre_pnt) having the same extreme points
index.same = rep(0, nrow(extre_pnt_four2))
for(i in 1: nrow(extre_pnt_four2)){
        temp_extre_four = extre_pnt_four2[i, ]
        for(j in 1: nrow(five_to_four_extre_pnt2)){
              temp_proj = five_to_four_extre_pnt2[j, ]
              if(all(temp_extre_four == temp_proj)){
                      index.same[i] = j
              }
        }
}

## two-row table, the first row: index of four-node extreme points; the second row: the index of after-projection.
extend.tab = matrix(c(which(index.same!=0), index.same[which(index.same!=0)]), nrow = 2, byrow = T)

extre_pnt_four2[extend.tab[1,], ]
five_to_four_extre_pnt2[extend.tab[2,], ]

# graph kept
graphKept = rownames(extre_pnt_four2[extend.tab[1,], ])

# #------------------------------------------------------------------------------#
# # convex combination of the left extreme points after projectioin
# #------------------------------------------------------------------------------#
# # combine two vmat togather, since "fixed" extreme points appear in two matrix, remove ones from projection matrix
# two.vmat = as.matrix(rbind(extre_pnt_four, five_to_four_extre_pnt[-extend.tab[2,], ]))
# twoRdt = redundant(two.vmat, representation = "V")
# 
# rm.mat = two.vmat[twoRdt$redundant, ]
# # take one row as an example
# Exrm = rm.mat[4, ]


#------------------------------------------------------------------------------#
# any point in the convex hull has a unique representation as a mixture of all the extreme points.
# for example, E: extreme poionts matrix, x: some vector in the convex hull. Then lambda = E^-1 %*% x, once the entry of lambda is zero, which means its coressponding extreme point not involving.
#------------------------------------------------------------------------------#

extre_pnt_four3 = q2d(extre_pnt_four2)
extre_four_name = colnames(extre_pnt_four3)

extre_pnt_four_rowname = rownames(extre_pnt_four3)

five_to_four_extre_pnt3 = q2d(five_to_four_extre_pnt2)

Lambdas = round(round(solve(t(extre_pnt_four3))) %*% t(five_to_four_extre_pnt3), 3)

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
        graphs = extre_pnt_four_rowname[iLambda != 0]
        is_kept[i] = ifelse(any(graphs %in% graphKept), 1, 0)
        graph_kept[i] = paste(graphs[graphs %in% graphKept], collapse = ", ")
        vec = paste(iLambda[iLambda != 0], graphs)
        combina.vec[i] = paste(vec, collapse = ", ")
        number_cmbn[i] = length(vec)
}
five_to_four_extre_cmbn = data.frame(number = number_cmbn, combination = combina.vec, isKept = is_kept, graphKept = graph_kept, five_to_four_extre_pnt2, check.names = F)

# save
save(five_to_four_extre_cmbn, file = "data/five_to_four_combination.Rdata")
write.csv(five_to_four_extre_cmbn, file = "data/five_to_four_combination.csv")
