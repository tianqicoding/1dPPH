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
import timeit

filename='NeuronConnect.xls'
neurontypename='neuronname.xls'
cs=pd.read_excel(filename)
G=nx.DiGraph()
print(cs.at[0, 'Neuron 1'])
for i in range(cs.shape[0]):
	if cs.at[i, 'Type']=='S' or cs.at[i, 'Type']=='Sp':
		G.add_edge(cs.at[i, 'Neuron 1'], cs.at[i, 'Neuron 2'], weight=cs.at[i, 'Nbr'])
print(len(G.edges), len(G.nodes))
# g=persHomo(G)
# print("done1")
# start = timeit.default_timer()
# g.perHom(50)
# stop = timeit.default_timer()
# g.PrintPair()
# print('Time: ', stop - start)
# print('# of non-tree edges: ', len(g.NonTreeEdge))
# print('# of edges: ', len(g.curG.edges))
# minB=minBasis(g)
# minB.ComputeAnnotation()
# minB.ComputeMinimalBasis('celegans1')

# nn=50
# G=nx.DiGraph()
# G.add_nodes_from(list(range(nn)))
# E=[]
# for i in range(nn-1):
# 	#E.append((i, i+1))
# 	G.add_edge(i, i+1, weight=1)
# G.add_edge(nn-1, 0, weight=1)
# g=persHomo(G)
# start = timeit.default_timer()
# g.perHom(2)
# #print(g.Root)
# stop = timeit.default_timer()
# g.PrintPair()
# print('Time: ', stop - start)

# G=nx.DiGraph()
# G=nx.read_weighted_edgelist("cora", create_using = nx.DiGraph())
# print(len(G.nodes), len(G.edges))
# G=nx.DiGraph()
# G=nx.read_weighted_edgelist("citeseer", create_using = nx.DiGraph())
# print(len(G.nodes), len(G.edges))
# g=persHomo(G)
# start = timeit.default_timer()
# g.perHom(2)
# stop = timeit.default_timer()
# print('Time: ', stop - start)
# print(len(g.TreeEdge), len(g.curG.edges), len(g.Bnd))
# g.PrintPair()
# g.PrintCycle()