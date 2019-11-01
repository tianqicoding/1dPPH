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
import timeit
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
class Tree:
	def __init__(self):
		self.tree=collections.defaultdict(list)
		self.cycle={}
class Graph:
	def __init__(self, V, E): #V is the vertex set, E is the edge set (u,v) u->v
		self.G=nx.DiGraph()
		for i in range(len(V)):
			self.G.add_node(V[i])
		self.G.add_edges_from(E)


class PathHomo:
	def __init__(self, G):
		self.nontreeEdge=[]
		self.G=nx.DiGraph(G)
		
		self.Boundary=[]
		self.Quad=[]
		self.TreeEdge=set()
		self.Bnd=None
		self.H=[]
		self.nontreeEdgeId={}		
		self.Hlist=[]
		self.Tree={}
		self.TreeRoot=[]
		self.MinHomoBasis=[]
		self.HomoBasis=[]
		self.Support=[]
		self.Candidate={}
		self.HomoBasisId={}
		self.cycles=[]
		self.dim=0
	def ComputeCycle(self):#g is the graph
		#compute spanning trees, return to a set of nontree edges
		s=set(self.G.nodes)
		se=set(self.G.edges)
		
		#print(se)
		
		#vid='ADAL'#
		#s.remove(vid)
		vid=s.pop()
		
		###Compute Tree
		while True:
			self.TreeRoot.append(vid)
			q=collections.deque([vid])

			self.Tree[vid]=(vid, 0, 0)
			V=0
			
			while q:
				v=q.popleft()
				_,level,_a=self.Tree[v]
				#print(level)
				for w in self.G.successors(v):
					if w in s:
						se.remove((v,w))
						s.remove(w)
						q.append(w)
						V+=1
						self.TreeEdge.add((v,w))
						self.Tree[w]=(v, level+1, 1) ####v->w


				for w in self.G.predecessors(v):
					if w in s:
						#print()
						#print(w,v,eid, se, 'in')
						s.remove(w)
						se.remove((w,v))
						q.append(w)
						V+=1
						self.TreeEdge.add((w,v))
						self.Tree[w]=(v, level+1, -1)#####w->v
			if s:
				vid=s.pop()
			else:
				break
		print("tree edge", len(self.TreeEdge))
		#for every nontree edge, take the cycle
		for i,e in enumerate(se):
			self.nontreeEdge.append(e)
			self.nontreeEdgeId[e]=i
		print("nontreeedge", len(self.nontreeEdge))
		#print(se)
	def BiEdge(self):
		d=len(self.nontreeEdge)
		E=set(self.G.edges)
		for e in self.G.edges:
			u,v=e[0], e[1]

			if (v, u) in E:
				bi=[0 for i in range(d)]
				if e not in self.TreeEdge:
					bi[self.nontreeEdgeId[(u, v)]]=1
				if (v, u) not in self.TreeEdge:
					bi[self.nontreeEdgeId[(v, u)]]=-1
				self.Boundary.append(bi)
		print(matrix_rank(self.Boundary))
	def Triangle(self):
		d=len(self.nontreeEdge)
		for u in self.G.nodes:
			su=set(self.G.successors(u))
			for v in su:
				for w in self.G.successors(v):
					if w in su:
						tri=[0 for i in range(d)]
						#u->v, u->w, v->w
						if (u,v) not in self.TreeEdge:
							tri[self.nontreeEdgeId[(u,v)]]=1
						
						if (u,w) not in self.TreeEdge:
							tri[self.nontreeEdgeId[(u,w)]]= -1
							
						if (v,w) not in self.TreeEdge:
							tri[self.nontreeEdgeId[(v,w)]]=1
						self.Boundary.append(tri)
		print(matrix_rank(self.Boundary))
	def Quadrangle(self):
		d=len(self.nontreeEdge)
		for vi in self.G.nodes:
			U=collections.defaultdict(list)
			for u in self.G.successors(vi):
				for w in self.G.successors(u):
					if w!=vi:
						U[w].append(u)
			for w in U:
				if len(U[w])>=2:
					for i in range(len(U[w])-1):
						quad=[0 for i in range(d)]
						#vi->U[w][i]->w vi->U[w][i+1]->w
						if (vi, U[w][i]) not in self.TreeEdge:
							quad[self.nontreeEdgeId[(vi, U[w][i])]]=1
						if (U[w][i], w) not in self.TreeEdge:
							quad[self.nontreeEdgeId[(U[w][i], w)]]= 1

						if (vi, U[w][i+1]) not in self.TreeEdge:
							quad[self.nontreeEdgeId[(vi, U[w][i+1])]]=-1
						if (U[w][i+1], w) not in self.TreeEdge:
							quad[self.nontreeEdgeId[(U[w][i+1], w)]]=-1
						self.Boundary.append(quad)
		print(matrix_rank(self.Boundary))
	def ComputeBasis(self):
		self.ComputeCycle()
		A=[]
		d=len(self.nontreeEdge)
		self.Boundary=[]
		#print(len(self.TreeEdge), len(self.nontreeEdge), len(self.G.edges))
		self.BiEdge()
		
		self.Quadrangle()
		
		self.Triangle()
		
		#print(d)
		
		#print(matrix_rank(A))
		#p, l, u = lu(matrix(A).T)
		#print(matrix_rank(l))
		#print(matrix_rank(l))
		self.Bnd=matrix(self.Boundary)
		
		#if self.Bnd.size==0:
		#	print("dimension of the boundary group is : 0", "dimension of the homology group is :", d)
		#else:
		#	dim=matrix_rank(self.Bnd)
		#	print("dimension of the boundary group is :", dim, "dimension of the homology group is :", d-dim)
	

	

	def Homology(self):
		
		#print(self.adjcentTree)
		#print(self.TreeEdge)
		#print(self.nontreeEdge)
		#print(len(self.TreeEdge), len(self.nontreeEdge))
		#print(self.Boundary)
		self.ComputeBasis()
		#print(self.Boundary)
		H=[]
		A=self.Bnd

		m=A.shape[0]
		n=A.shape[1]
		I=set(range(m))
		pivot=numpy.array([None for i in range(m)])
		#print('pivit', pivot)
		A=numpy.float_(A)
		t=0
		tmp=0
		while I:
			i=I.pop()
			
			j=None
			#print(A.shape)
			aaaa=numpy.where(A[i]!=0)
			#print(aaaa)
			if aaaa[1].shape[0]>0:
				j=aaaa[1][0]


			# for k in range(A.shape[1]):
			# 	if A[i,k]!=0:
			# 		j=k
			# 		break
			
			if j!=None:
				tmp+=1
				pivot[i]=j
				#print(111)
				for k in range(m):
					if k!=i and A[k,j]!=0:
						A[k]=A[k]-A[i]*1.0/A[i,j]*A[k,j]

				
		#print(tmp)
		#print(matrix_rank(A))
		cnt=0
		ss=set(range(A.shape[1]))
		for x in pivot:
			if x!=None:
				cnt+=1
				ss.remove(x)
		print(cnt, len(ss))
		#ss=numpy.where(pivot==None)[0]
		#print('ss', ss.shape)
		
		cnt=0
		
		self.HomoBasis=[[0 for i in range(len(ss))] for i in range(len(self.nontreeEdge))]
		self.HomoBasis=matrix(self.HomoBasis)
		#print(self.HomoBasis.shape)
		self.HomoBasis=numpy.float_(self.HomoBasis)
		for item in ss:
			H.append(self.nontreeEdge[item])
			self.HomoBasisId[self.nontreeEdge[item]]=cnt
			self.HomoBasis[item, cnt]=1
			#print(item, cnt)
			cnt+=1
		aaaa=[]
		for i in range(pivot.shape[0]):
			if pivot[i]:
				aaaa.append((pivot[i], i))
		aaaa.sort(key=lambda x:x[1], reverse=True)
		self.dim=len(ss)
		for t, i in aaaa:
			#t=pivot[i]

			if t!=None:
				
				for j in range(A.shape[1]):
					if j!=t and A[i,j]!=0:
						self.HomoBasis[t]+=-A[i,j]*1.0/A[i, t]*self.HomoBasis[j]
						#print(A[i,j], A[i, t], -A[i,j]*1.0/A[i, t])

		
		
		#print(self.HomoBasis)

		
		
		
		

		
		# 
		# attrs=collections.defaultdict(dict)

		# for i, H in enumerate(self.Hlist):
		# 	for e in H:
		# 		attrs[e]['attr'+str(i)]=1
		
		# nx.set_edge_attributes(self.G, attrs)
		# #print(i, self.MinHomoBasis)
		# file=open("testfile1.txt","w")
		# for H in self.Hlist:
		# 	for e in H:
		# 		file.write('(')
		# 		file.write(e[0])
		# 		file.write('\t')
		# 		file.write(e[1])
		# 		file.write(')')
		# 		file.write('\t')
		# 	file.write('\n')
		# file.close()
		
		

		
		

	def GraphPlot(self):
		return
		#nx.write_gml(self.G,'g1.gml')
		#nx.write_graphml(self.G,'g1.xml')
		# fig = plt.figure()
		# xrandom=[1]

		# pos = nx.spring_layout(self.G)
	
		# nx.draw_networkx_nodes(self.G, pos, node_size=15, node_color='b')
		# nx.draw_networkx_edges(self.G, pos)
		# nx.draw_networkx_edges(self.G, pos, edgelist=self.Hlist[0], edge_color='r', width=2, alpha=0.5)

		# plt.show()
	def ComputeCandidate(self):
		d=len(self.HomoBasisId)
		print("homodimension", d)
		for i,vi in enumerate(self.G.nodes):
			self.Candidate[vi]=Tree()
			s=set(self.G.nodes)
			s.remove(vi)
			q=collections.deque([vi])
			se=set(self.G.edges)
			lev={}

			lev[vi]=(vi,0,0)
			while q:
				v=q.popleft()
				_, level, _=lev[v]
				for w in self.G.successors(v):
					if w in s:
						self.Candidate[vi].tree[v].append(w)
						se.remove((v,w))
						s.remove(w)
						q.append(w)
						lev[w]=(v,level+1, 1)
				for w in self.G.predecessors(v):
					if w in s:
						self.Candidate[vi].tree[v].append(w)
						s.remove(w)
						se.remove((w,v))
						q.append(w)
						lev[w]=(v,level+1, -1)
			for edge in se:
				edgelist=[0.0 for i in range(d)]
				edgelist=numpy.array(edgelist)
				v1=edge[0]
				v2=edge[1]

				_, d1, _=lev[v1]
				_, d2, _=lev[v2]
				
				edgecycle=[]
				if edge in self.nontreeEdgeId:
					#print(self.HomoBasis)
					edgelist=edgelist+self.HomoBasis[self.nontreeEdgeId[edge]]
				edgecycle.append(edge)
				w=1
				#print(d1, d2)
				while d1>d2:
					parent, d1, _in=lev[v1]
					#print(d1)
					if _in==-1:
						edgecycle.append((v1, parent))
						if (v1,parent) in self.nontreeEdgeId: ##v1->parent
							edgelist=edgelist-self.HomoBasis[self.nontreeEdgeId[(v1,parent)]]
					elif _in==1:
						edgecycle.append((parent, v1))
						if (parent, v1) in self.nontreeEdgeId:
							edgelist=edgelist+self.HomoBasis[self.nontreeEdgeId[(parent, v1)]]
						
					v1=parent
					d1-=1
					w+=1
				while d1<d2:
					parent, d2, _in=lev[v2]
					#print(d1)
					if _in==-1:
						if (v2,parent) in self.nontreeEdgeId: ##v1->parent
							edgelist=edgelist+self.HomoBasis[self.nontreeEdgeId[(v2,parent)]]
						edgecycle.append((v2, parent))
					elif _in==1:
						if (parent, v2) in self.nontreeEdgeId:
							edgelist=edgelist-self.HomoBasis[self.nontreeEdgeId[(parent, v2)]]
						edgecycle.append((parent, v2))
					v2=parent
					d2-=1
					w+=1
				#print(d1, d2, edgelist)

				while v1!=v2 and d1>0:
					parent1, d1, _=lev[v1]

					if _==-1:
						edgecycle.append((v1, parent1))
						if (v1,parent1) in self.nontreeEdgeId: ##v1->parent
							edgelist=edgelist-self.HomoBasis[self.nontreeEdgeId[(v1,parent1)]]
						
					elif _==1:
						if (parent1, v1) in self.nontreeEdgeId:
							edgelist=edgelist+self.HomoBasis[self.nontreeEdgeId[(parent1, v1)]]
						edgecycle.append((parent1, v1))
					w+=1
					parent2, d2, _=lev[v2]
					if _==-1:
						if (v2,parent2) in self.nontreeEdgeId: ##v2->parent
							edgelist=edgelist+self.HomoBasis[self.nontreeEdgeId[(v2,parent2)]]
						edgecycle.append((v2, parent2))
					elif _==1: 
						if (parent2, v2) in self.nontreeEdgeId:
							edgelist=edgelist-self.HomoBasis[self.nontreeEdgeId[(parent2, v2)]]
						edgecycle.append((parent2, v2))
					w+=1
					v1,v2=parent1,parent2
				if len(edgecycle)==2 and edgelist.any():
					print(1)
				self.cycles.append((w,edgelist, edgecycle))
				#print(edgecycle)

						


	def ComputeMinimalBasis(self, filename):
		self.ComputeCandidate()
		dim_homo=len(self.HomoBasisId)
		
		self.cycles.sort(key=lambda x:x[0])
		KK=[a for _,a, _ in self.cycles]
		print(matrix_rank(KK))
		
		
		i=0
		A=[]
		B=[]
		while i<len(self.cycles):
			if self.cycles[i][1].any():
				A.append(self.cycles[i][1].tolist()[0])
				B.append(self.cycles[i][2])
			i+=1
		#A=self.cycles[i][1]
		a=matrix_rank(A)
		print(a)
		# while i<len(self.cycles) and a<dim_homo:
		# 	w, cyclist, cycedge=self.cycles[i]
			
		# 	if not cyclist.any():
		# 		print('continue')
		# 		continue
		# 	#print(A)
		# 	A=numpy.append(A, cyclist, axis=0)
		# 	print(A)
		# 	rk=matrix_rank(A)
		# 	print(rk)
		# 	if rk>a:
		# 		a=rk
		# 		print(a)
		# 		self.MinHomoBasis.append(cycedge)
		# 	else:
		# 		A=numpy.delete(A, A.shape[0]-1, 0)
		# 	i+=1
		# #print(A)
		A=matrix(A)
		m=A.shape[0]
		pivot=[None for i in range(m)]
		A=numpy.float_(A)
		t=0
		tmp=0
		for i in range(m):
			j=None
			#print(A.shape)
			for k in range(A.shape[1]):
				if A[i,k]!=0:
					j=k
					break
			
			if j!=None:
				tmp+=1
				pivot[i]=j
				#print(111)
				for k in range(m):
					if k!=i and A[k,j]!=0:
						A[k]=A[k]-A[i]*1.0/A[i,j]*A[k,j]
		for i in range(len(pivot)):
			if pivot[i]!=None:
				self.MinHomoBasis.append(B[i])
		attrs=collections.defaultdict(dict)

		for i, H in enumerate(self.MinHomoBasis):
			for e in H:
				attrs[e]['attr'+str(i)]=1
		
		nx.set_edge_attributes(self.G, attrs)
		#print(i, self.MinHomoBasis)
		file=open(filename,"w")
		#print(self.MinHomoBasis)
		for H in self.MinHomoBasis:
			for e in H:
				file.write('(')
				file.write(str(e[0]))
				file.write('\t')
				file.write(str(e[1]))
				file.write(')')
				file.write('\t')
			file.write('\n')
		file.close()
		# #print(A)
		# return A

	

#GraphPlot()
# def GraphPlot(filename):
# 	lines = []
# 	with open(filename) as file:
# 		for line in file:
# 			line = line.strip('\n').split('\t')
# 			l=[]
# 			i=0
# 			while i<len(line)-1:
# 				l.append((line[i], line[i+1]))
# 				i+=2
# 			lines.append(l)
# 		t=len(lines)
# 		if t>10:
# 			t=10
# 		pos = nx.spring_layout(self.G)
# 		for i in range(1,t+1):
# 			plt.subplot(2, 3, i)
#     		#nx.draw(self.G.G)
# 			nx.draw_networkx_nodes(self.G, pos, node_size=15)
# 			nx.draw_networkx_edges(self.G, pos)
# 			nx.draw_networkx_edges(self.G, pos, edgelist=self.Hlist[i-1], edge_color='r', width=2, alpha=0.5)

# 		plt.show()


# g=PathHomo(['a','b','c','d','e','f','g'],[('a','b'), ('b','d'), ('b','c'), ('d', 'e'), ('d', 'a'), ('c', 'e'), ('c', 'a'), ('e', 'f'), ('f', 'c'), ('c', 'g'), ('g', 'f')])
# g.ComputeCycle()
# g.Triangle()
# g.Quadrangle()
# g.ComputeBasis()
# g.Homology()
# g.GraphPlot()
# for _ in range(10):
# 	filename='NeuronConnect.xls'
# 	neurontypename='neuronname.xls'
# 	cs=pd.read_excel(filename)
# 	neurontype=pd.read_excel(neurontypename)
# 	attrs=collections.defaultdict(dict)
# 	for i in range(len(neurontype)):
# 		row=neurontype.iloc[i]
		
# 		at1=row['motor']
# 		at2=row['inter']
# 		at3=row['sensory']
# 		at4=row['unknown']
# 		if at4==1:
# 			attrs[row['Neuron']]['type']=0
# 		elif at1==1 and not at2==1 and not at3==1:
# 			attrs[row['Neuron']]['type']=1
# 		elif not at1==1 and at2==1 and not at3==1:
# 			attrs[row['Neuron']]['type']=2
# 		elif not at1==1 and not at2==1 and at3==1:
# 			attrs[row['Neuron']]['type']=3
# 		elif at1==1 and at2==1 and not at3==1:
# 			attrs[row['Neuron']]['type']=4
# 		elif at1==1 and not at2==1 and at3==1:
# 			attrs[row['Neuron']]['type']=5
# 		elif not at1==1 and at2==1 and at3==1:
# 			attrs[row['Neuron']]['type']=6
# 		else:
# 			attrs[row['Neuron']]['type']=7

		

# 	V=set()
# 	E=set()
# 	for i in range(len(cs)):
# 		row=cs.iloc[i]
# 		edgetype=row['Type']
# 		v1=row['Neuron 1']
# 		v2=row['Neuron 2']
# 		#print(edgetype)
# 		if edgetype=='S' or edgetype=='Sp':# or edgetype=='R' or edgetype=='Rp':
# 			V.add(v1)
# 			V.add(v2)
# 			#a=(v2, v1)
# 			#if edgetype=='S' or edgetype=='Sp':
# 			a=(v1, v2)
# 			E.add(a)
# 	G=nx.DiGraph()
# 	for i in V:
# 		G.add_node(i)
# 	nx.set_node_attributes(G, attrs)
# 	G.add_edges_from(E)

# 	print(len(V),len(E))


# 	g=PathHomo(G)
# 	g.ComputeBasis()
# 	g=None
# 	G=None
# start = timeit.default_timer()
# g.Homology()
# stop = timeit.default_timer()
# print('Time: ', stop - start)

# g.ComputeMinimalBasis("minbasis_neuron.txt")

# #
# g.GraphPlot()
# l=[]
# def GraphPlot(l):
# 	#nx.write_gml(self.G,'g1.gml')
# 	#nx.write_graphml(self.G,'g1.xml')
# 	#fig = plt.figure()
# 	xrandom=[0 for i in l]
# 	x = numpy.array([0,1])
# 	plt.scatter(xrandom, l, marker='^', color='blue')
# 	plt.scatter([1], [17], marker='o', color='red')
# 	plt.xticks(x, ['random graph', 'C.elegance'])
# 	plt.show()
# def Detect(h):
# 	if len(h)!=4:
# 		return False
# 	a=set([x for (x, _) in h])
# 	if len(a)!=4:
# 		return False
# 	b=set([x for (_, x) in h])
# 	if len(b)!=4:
# 		return False
# 	return True


# cnt_random=0
# l=[]
# for i in range(1000):
# 	RG=nx.generators.random_graphs.erdos_renyi_graph(279, 0.028, directed=True)
# 	rg=PathHomo(RG)
# 	rg.Homology()
	

# 	#rg.ComputeMinimalBasis(str(i)+'.txt')
# 	#for h in rg.MinHomoBasis:
# 	#	if Detect(h):
# 	#		cnt_random+=1
# 	#		break
# 	l.append(rg.dim)
# 	print(i, 'cnt', rg.dim)
# l=numpy.array(l, dtype=numpy.int32)

# fileobj = open('random.csv', mode='wb')
# l.tofile(fileobj)
# fileobj.close()
# print(l)


# GraphPlot(l)

#ef GraphPlot(l):
# g.Homology()
# g.ComputeMinimalBasis()
# #
# g.GraphPlot()
# filename='industry.csv'
# with open(filename) as csvfile:
#     readCSV = csv.reader(csvfile, delimiter=',')
#     dates = []
#     colors = []
#     cnt=0
#     A=[]
#     V=[]
#     E=[]
#     vid=0
#     #print(readCSV)
#     for row in readCSV:
#     	if cnt==0:
#     		nn=len(row)-4
#     		A=[[0 for i in range(nn)] for i in range(nn)]
#     		for i in range(2,len(row)-2):
#     			V.append(row[i])
#     		#print(V)
#     	elif cnt>1:
#     		for i in range(2, len(row)-2):
#     			if vid!=i-2:
#     				#E.append((V[vid], V[i-2], float(row[i])))
#     				A[vid][i-2]=float(row[i])
#     		vid+=1
#     	cnt+=1
#     for i in range(nn):
#     	sm=0
#     	for j in range(nn):
#     		sm+=A[j][i]
#     	for j in range(nn):
#     		A[j][i]=A[j][i]/sm
#     for i in range(nn):
#     	for j in range(nn):
#     		if i!=j:
#     			A[i][j]=1-A[i][j]
#     for i in range(nn):
#     	for j in range(nn):
#     		if 0<A[i][j]<=0.9:
#     			E.append((V[i], V[j], A[i][j]))
#     	#print(1)
#     print(V)
#     G=nx.DiGraph()
#     G.add_nodes_from(V)
#     for e in E:
#     	G.add_edge(e[0], e[1], weight=e[2])
#     g=PathHomo(G)

# nn=100
# G=nx.DiGraph()
# G.add_nodes_from(list(range(nn)))
# E=[]
# for i in range(nn-1):
# 	#E.append((i, i+1))
# 	G.add_edge(i, i+1, weight=1)
# G.add_edge(nn-1, 0, weight=1)
# g=PathHomo(G)
# start = timeit.default_timer()
# g.Homology()
# stop = timeit.default_timer()
# print('Time: ', stop - start)

G=nx.DiGraph()
G=nx.read_weighted_edgelist("citeseer")
g=PathHomo(G)
g.ComputeBasis()
start = timeit.default_timer()
#g.Homology()
stop = timeit.default_timer()
print('Time: ', stop - start)
print(len(g.TreeEdge), len(g.Bnd))

#     g.Homology()
    
#     g.ComputeMinimalBasis()

#     g.GraphPlot()
