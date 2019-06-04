rm(list = ls())

load(file = "application/samplk1_undirected.Rdata")
# load the ineqs which define the pace of 5->3
load(file = "data/five_to_three_hmat.Rdata")

combn_monk = combn(1:18, 3)

graph.vec= rep(NA, ncol(combn_monk))
for(i in 1: ncol(combn_monk)){
        monk3_i = combn_monk[, i]
        graph_i = network[monk3_i, monk3_i]
        graph.vec[i] = sum(graph_i) /2
}

graph_tab = table(graph.vec)
graph_tab = graph_tab/sum(graph_tab)
names(graph_tab) = c("3k1", "p3.bar", "p3", "k3")

