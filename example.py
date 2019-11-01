import csv
import pandas as pd 
import collections
from numpy import matrix, rank
from numpy.linalg import matrix_rank
import networkx as nx

#import matplotlib.pyplot as plt
import numpy
import sys
import heapq
#from GF import GF
import sympy 
from scipy.linalg import lu
from minBasis import minBasis
from PPH import persHomo

filename='NeuronConnect.xls'
neurontypename='neuronname.xls'
cs=pd.read_excel(filename)
G=nx.DiGraph()
print(cs.at[0, 'Neuron 1'])
for i in range(cs.shape[0]):
	if cs.at[i, 'Type']=='S' or cs.at[i, 'Type']=='Sp':
		G.add_edge(cs.at[i, 'Neuron 1'], cs.at[i, 'Neuron 2'], weight=cs.at[i, 'Nbr'])
print(max(G.edges))
g=persHomo(G)
print("done1")
g.perHom(50)
print('# of non-tree edges: ', len(g.NonTreeEdge))
print('# of edges: ', len(g.curG.edges))
minB=minBasis(g)
minB.ComputeAnnotation()
minB.ComputeMinimalBasis('celegans1')