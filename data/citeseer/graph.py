import numpy as np
import networkx as nx
import torch as th

##############input is a list of graphs
g=nx.DiGraph()

citationname='citeseer.cites'
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
	
	
nx.write_weighted_edgelist(g, "../../citeseer", comments='#', delimiter=' ', encoding='utf-8')

print(g.edges)

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
