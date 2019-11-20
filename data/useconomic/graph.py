import numpy as np
import networkx as nx
import torch as th

with open('1997.txt', 'r') as f:
    l = [[float(num) for num in line.split(' ')[:-1]] for line in f]
mat=np.matrix(l)
mat.resize((15, 15))
#print(mat.shape)
G=nx.from_numpy_matrix(mat, create_using=nx.DiGraph)
