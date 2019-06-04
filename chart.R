rm(list = ls())
#------------------------------------------------------------------------------#
# Make a table, combine all the informatin of graphs
#------------------------------------------------------------------------------#
six_node= read.csv(file = "data/six_nodes_unique_graphs.csv", check.names = F, stringsAsFactors = F, header = T, row.names = 1)
five_node= read.csv(file = "data/five_nodes_unique_graphs.csv", check.names = F, stringsAsFactors = F, header = T, row.names = 1)
four_node= read.csv(file = "data/four_nodes_unique_graphs.csv", check.names = F, stringsAsFactors = F, header = T, row.names = 1)
three_node= read.csv(file = "data/three_nodes_unique_graphs.csv", check.names = F, stringsAsFactors = F, header = T, row.names = 1)

grabGraph =  function(mat){
        graph.vec = c()
        col.name = colnames(mat)
        for(i in 1: nrow(mat)){
                vec = mat[i, ]
                lines = col.name[vec != 0]
                graph.vec = c(graph.vec, paste(lines,collapse=", "))
        }
        return(graph.vec)
}

six.tab = rbind(name =  rownames(six_node), node = rep(6, nrow(six_node)), edge = apply(six_node, 1, sum), graph = grabGraph(six_node))


five.tab = rbind(name =  rownames(five_node), node = rep(5, nrow(five_node)), edge = apply(five_node, 1, sum), graph = grabGraph(five_node))

four.tab = rbind(name =  rownames(four_node), node = rep(4, nrow(four_node)), edge = apply(four_node, 1, sum), graph = grabGraph(four_node))

three.tab = rbind(name =  rownames(three_node), node = rep(4, nrow(three_node)), edge = apply(three_node, 1, sum), graph = grabGraph(three_node))

graphTable = cbind(three.tab, four.tab, five.tab, six.tab)

graphTable = as.data.frame(graphTable)

write.csv(graphTable, "graphTable.csv")
save(graphTable, file = "graphTable.Rdata")
