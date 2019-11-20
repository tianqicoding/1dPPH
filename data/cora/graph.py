import numpy as np
import dgl
import networkx as nx
import torch as th
from dgl.data import citation_graph as citegrh
from scipy import sparse
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import dgl
import networkx as nx
import torch as th

g=nx.DiGraph()

citationname='cora.cites'
citationfile=open(citationname, 'r')
citationstring=citationfile.read()
citationl=citationstring.split('\n')
citationl=citationl[:-1]
for i in range(len(citationl)):
	citationl[i]=citationl[i].split('\t')
	if citationl[i][0]==citationl[i][1]:
		continue
	print(citationl[i])
	g.add_edge(citationl[i][0], citationl[i][1], weight=1)
	
	
nx.write_weighted_edgelist(g, "../../cora", comments='#', delimiter=' ', encoding='utf-8')

print(len(g.edges), len(g.nodes))


# ##############input is a list of graphs
# data=citegrh.load_citeseer()
# print(len(data.graph.edges))
# GMulti=dgl.DGLGraph(data.graph)
# print(data.features.shape)
# features=sparse.csr_matrix(data.features)
# similarities = cosine_similarity(features)
# print(similarities)
# edgeset1=set(data.graph.edges)
# print(len(edgeset1))
# #print(np.count_nonzero(similarities))
# EdgeCount=0
# edgeset2=set()
# for i in range(similarities.shape[0]):
# 	for j in range(i+1,similarities.shape[0]):
# 		if similarities[i, j]>0.2:
# 			EdgeCount+=1
# 			edgeset2.add((i, j))
# 			GMulti.add_edge(i, j)

# print(GMulti)
# print(EdgeCount)



#print(features.toarray().shape)
# for i in range(features.shape[0]):
# 	for j in range(i+1, features.shape[0]):

# filename='cora.content'
# file=open(filename, 'r')
# string=file.read()
# l=string.split('\n')
# NodeDic={}

# for i in range(len(l)):
# 	l[i]=l[i].split('\t')
# 	NodeDic[l[i][0]]=len(NodeDic)
# print(len(NodeDic))


# citationname='cora.cites'
# citationfile=open(citationname, 'r')
# citationstring=citationfile.read()
# citationl=citationstring.split('\n')
# citationl=citationl[:-1]
# for i in range(len(citationl)):
# 	citationl[i]=citationl[i].split('\t')
	
# 	#print(citationl[i])
	
# 	if citationl[i][0] not in NodeDic:
# 		#NodeDic[citationl[i][0]]=len(NodeDic)
# 		print(citationl[i][0])
# 	if citationl[i][1] not in NodeDic:
# 		#NodeDic[citationl[i][1]]=len(NodeDic)
# 		print(citationl[i][1])


# NumberOfGraphs=2
# GMulti=dgl.DGLGraph()
# graphMat=[]
# NumberOfNodes=len(NodeDic)
# GMulti.add_nodes(NumberOfNodes)
# edgeset=set()
# print(NumberOfNodes)
# for i in range(len(citationl)):
	
# 	u=citationl[i][0]
# 	v=citationl[i][1]
# 	GMulti.add_edge(NodeDic[u], NodeDic[v])
# 	edgeset.add((NodeDic[u], NodeDic[v]))
# print(GMulti)

#GMulti.edata['att'] = th.zeros((len(edgeset), 2))




	
# for i in range(len(graphMat)):
# 	for e in graphMat[i].edges:
# 		if e not in edgeset:
# 			GMulti.add_edge(e[0], e[1])
# 			edgeset.add(e)
# 			#GMulti.edges[e[0], e[1]].data['att']=th.tensor(th.zeros(NumberOfGraphs))
# 		#GMulti.edges[e].data['att'][i]=1
# GMulti.edata['att'] = th.zeros((len(edgeset), NumberOfGraphs))
# for i in range(len(graphMat)):
# 	for e in graphMat[i].edges:
# 		a=GMulti.edges[e[0], e[1]].data['att']
# 		a[0][i]=1
# 		#print(a)
# 		GMulti.edges[e[0], e[1]].data['att']=a

#print(GMulti.edata)



# # GList=[]
# # for file in filelist:
# # 	mat=np.loadtxt(file)
# # 	#print(mat.shape)
# # 	nxg=nx.from_numpy_matrix(mat)
# # 	#print(G.nodes)
# # 	g = dgl.DGLGraph()
# # 	g.from_networkx(nxg)
# # 	print(g)
# # 	GList.append(g)


# ###################a multi graph
