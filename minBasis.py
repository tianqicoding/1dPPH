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
from PPH import persHomo
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


class minBasis:
	def __init__(self, persG):
		self.persG=persG
		self.Candidate={}
		self.G=persG.G
		self.HomoBasisId={}
		self.cycles=[]
		self.MinHomoBasis=[]
		self.HomoBasis=numpy.float_(matrix([[0 for i in range(len(persG.HomoEdgeId))] for i in range(len(persG.NonTreeEdge))]))
	def ComputeAnnotation(self):
		cnt=0
		for edge in self.persG.HomoEdgeId:
			#print(edge)
			#.append(self.NonTreeEdge[item])
			edgeId=self.persG.EdgeId[edge]
			#print(edge, edgeId)
			self.HomoBasisId[edge]=cnt
			self.HomoBasis[edgeId, cnt]=1
			#print(item, cnt)
			cnt+=1

		A=numpy.float_(matrix([self.persG.Bnd[i][1].tolist()[0] for i in range(len(self.persG.Bnd))]))
		print(cnt, "number of basis")
		pivot=[self.persG.Bnd[i][0]for i in range(len(self.persG.Bnd))]
		aaaa=[]
		for i in range(len(pivot)):
			aaaa.append((pivot[i], i))
		aaaa.sort(key=lambda x:x[1], reverse=True)
		for t, i in aaaa:
			for j in range(A.shape[1]):
				if j!=t and A[i,j]!=0:
					self.HomoBasis[t]+=-A[i,j]*1.0/A[i, t]*self.HomoBasis[j]
		print(matrix_rank(self.HomoBasis))

	def ComputeCandidate(self):
		persG=self.persG
		#print(type(persG.G))
		d=len(self.HomoBasisId)
		print("homodimension", d)
		for i,vi in enumerate(persG.G.nodes):
			self.Candidate[vi]=Tree()
			treeNodes=set()
			s=set(persG.G.nodes)
			s.remove(vi)
			q=collections.deque([vi])

			se=set(persG.G.edges)
			lev={}
			treeNodes.add(vi)
			lev[vi]=(vi,0,0)
			while q:
				v=q.popleft()

				_, level, _=lev[v]
				for w in persG.G.successors(v):
					if w in s:
						self.Candidate[vi].tree[v].append(w)
						se.remove((v,w))
						treeNodes.add(w)
						s.remove(w)
						q.append(w)
						lev[w]=(v,level+1, 1)
				for w in persG.G.predecessors(v):
					if w in s:
						self.Candidate[vi].tree[v].append(w)
						s.remove(w)
						treeNodes.add(w)
						se.remove((w,v))
						q.append(w)
						lev[w]=(v,level+1, -1)
			for edge in se:
				if edge[0] not in treeNodes or edge[1] not in treeNodes:
					continue
				edgelist=[0.0 for i in range(d)]
				edgelist=numpy.array(edgelist)
				v1=edge[0]
				v2=edge[1]

				_, d1, _=lev[v1]
				_, d2, _=lev[v2]
				
				edgecycle=[]
				if edge in persG.EdgeId:
					#print(self.HomoBasis)
					edgelist=edgelist+self.HomoBasis[persG.EdgeId[edge]]
				edgecycle.append(edge)
				w=1
				#print(d1, d2)
				while d1>d2:
					parent, d1, _in=lev[v1]
					#print(d1)
					if _in==-1:
						edgecycle.append((v1, parent))
						if (v1,parent) in persG.EdgeId: ##v1->parent
							edgelist=edgelist-self.HomoBasis[persG.EdgeId[(v1,parent)]]
					elif _in==1:
						edgecycle.append((parent, v1))
						if (parent, v1) in persG.EdgeId:
							edgelist=edgelist+self.HomoBasis[persG.EdgeId[(parent, v1)]]
						
					v1=parent
					d1-=1
					w+=1
				while d1<d2:
					parent, d2, _in=lev[v2]
					#print(d1)
					if _in==-1:
						if (v2,parent) in persG.EdgeId: ##v1->parent
							edgelist=edgelist+self.HomoBasis[persG.EdgeId[(v2,parent)]]
						edgecycle.append((v2, parent))
					elif _in==1:
						if (parent, v2) in persG.EdgeId:
							edgelist=edgelist-self.HomoBasis[persG.EdgeId[(parent, v2)]]
						edgecycle.append((parent, v2))
					v2=parent
					d2-=1
					w+=1
				#print(d1, d2, edgelist)

				while v1!=v2 and d1>0:
					parent1, d1, _=lev[v1]

					if _==-1:
						edgecycle.append((v1, parent1))
						if (v1,parent1) in persG.EdgeId: ##v1->parent
							edgelist=edgelist-self.HomoBasis[persG.EdgeId[(v1,parent1)]]
						
					elif _==1:
						if (parent1, v1) in persG.EdgeId:
							edgelist=edgelist+self.HomoBasis[persG.EdgeId[(parent1, v1)]]
						edgecycle.append((parent1, v1))
					w+=1
					parent2, d2, _=lev[v2]
					if _==-1:
						if (v2,parent2) in persG.EdgeId: ##v2->parent
							edgelist=edgelist+self.HomoBasis[persG.EdgeId[(v2,parent2)]]
						edgecycle.append((v2, parent2))
					elif _==1: 
						if (parent2, v2) in persG.EdgeId:
							edgelist=edgelist-self.HomoBasis[persG.EdgeId[(parent2, v2)]]
						edgecycle.append((parent2, v2))
					w+=1
					v1,v2=parent1,parent2
				# if len(edgecycle)==2 and edgelist.any():
				# 	print(1)
				self.cycles.append((w,edgelist, edgecycle))
				#print(edgecycle)	

			print(i)	


	def ComputeMinimalBasis(self, filename):
		self.ComputeCandidate()
		persG=self.persG
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
