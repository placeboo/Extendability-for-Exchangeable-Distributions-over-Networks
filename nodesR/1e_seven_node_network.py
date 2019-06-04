
# coding: utf-8

# In[16]:


import numpy as np
from itertools import combinations
import networkx as nx
import matplotlib.pyplot as plt
import csv


# In[2]:


N = 3 # node 
NN = range(1, N+1)
comb =combinations(NN, 2)
edges_ls = [] # edges name
for edge in comb:
    edges_ls.append([edge[0], edge[1]])
edges_arr = np.array(edges_ls)
num_edge = len(edges_ls)
#print(edges_arr)


# In[3]:


count = 1
edgenode_class_dict = dict() # (edge_num, node_num) -> subgraph
class_edge_dict = dict() # subgraph -> edge graph matrix, need to be saved
#class_graph_dict = dict() # subgraph -> graph name
class_Graph_dict = dict() # subgraph -> graph object
values = np.arange(0,2**num_edge)



# In[4]:


for value in values:
    graph_i = np.zeros(num_edge, dtype=int)
    incidence_mat = np.zeros((N,N), dtype=int)
    #graph_i = [0] * num_edge
    tmp = [int(num) for num in bin(value)[2:]] # change bin to int arrary
    graph_i[num_edge-len(tmp):] = tmp # [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1]
    #edges_ls[np.where(graph_i==1)]
    #print(graph_i)
    edges = edges_arr[np.where(graph_i == 1)] # what are edges for the graph_i
    edge_num = edges.shape[0]
    node_num = np.unique(edges).shape[0]
    #print(edges)
    incidence_mat[edges[:,0]-1, edges[:,1]-1] = 1
    incidence_mat[edges[:,1]-1, edges[:,0]-1] = 1
    #print(incidence_mat)
    G_tmp = nx.Graph(incidence_mat)
    tmp_class = "class"+str(count)
    if((edge_num, node_num) in edgenode_class_dict.keys()):
        classes = edgenode_class_dict[(edge_num, node_num)] # List
        #print(classes)
        not_in_count = 0
        for class_j in classes: # check whether they are isomorphism.
            if(nx.is_isomorphic(G_tmp, class_Graph_dict[class_j])): # if isomorphism
                #print(class_edge_dict)
                class_edge_dict[class_j] = np.vstack([class_edge_dict[class_j] , graph_i])
                break
            else: 
                not_in_count += 1
        if(not_in_count == len(classes)): # there is no iso
            edgenode_class_dict[(edge_num, node_num)].append(tmp_class)
            class_edge_dict[tmp_class] = graph_i 
            #class_graph_dict[tmp_class] = [value]
            class_Graph_dict[tmp_class] = G_tmp
            count += 1         
    else: # add new keys to all dictionary
        edgenode_class_dict[(edge_num,node_num)] = [tmp_class]
        class_edge_dict[tmp_class] = graph_i 
        #class_graph_dict[tmp_class] = [value]
        class_Graph_dict[tmp_class] = G_tmp
        count += 1


# In[5]:


print(edgenode_class_dict)
print(class_edge_dict) 


# In[41]:


edge_graph_arr = np.empty((0,num_edge), int)
class_name = []
class_edge_dict.keys()
for keys in class_edge_dict.keys():
    class_tmp_arr = class_edge_dict[keys]
    #print(class_arr)
    class_name.extend([keys] * int(class_tmp_arr.size/num_edge))
    edge_graph_arr = np.vstack([edge_graph_arr, class_tmp_arr])
    

np.savetxt("edge_graph_arr.csv", edge_graph_arr, fmt='%i',delimiter=",")
#np.savetxt("edge_num_arr.csv", edge_num_arr, delimiter=",")

with open("class_name.csv", 'w') as myfile:
    for nm in class_name:
        myfile.write(nm + '\n')
        
egde_node_num_list = list(edgenode_class_dict.keys())
with open("egde_node_num.csv", 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for nm in egde_node_num_list:
        wr.writerow(nm)

    
#np.savetxt("file_name.csv", data1, delimiter=",", fmt='%s', header=header)
#np.savetxt("class_name.csv", class_name,delimiter=",")
#np.savetxt("edge_num_arr.csv", edge_num_arr,delimiter=",")


# In[38]:


values
edges_arr
print(edge_graph_arr)
print(class_name)
np.array(class_name, dtype= str)
edge_num_list

