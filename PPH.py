#input: graph G with weights, time t
#all nodes in, edges sorted from smallest to largest, only pick <=t
import pandas as pd
import csv
import networkx as nx
import collections
from numpy import matrix, rank
import numpy
from numpy.linalg import matrix_rank
import timeit
import gudhi
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
class persHomo:
	def __init__(self, G):
		self.G=G # weighted, directed
		#print(type(self.G))
		self.Root={}
		for v in G.nodes:
			self.Root[v]=None
		self.parent={}
		for v in G.nodes:
			self.parent[v]=None
		self.generators=[]
		self.TreeEdge=set() #tree edge set
		self.EdgeId={} # for every edge, give one id
		self.NonTreeEdge=[] #non-tree edge list
		self.EdgeCycle={} 
		self.HomoCycle={}
		self.HomoEdgeId=[]
		self.Bnd=[] #bnd list, every bnd is a vector
		self.pair=[] #pair list
		self.curG=nx.DiGraph() #current graph
		self.curG.add_nodes_from(self.G.nodes) 
		self.ns_nt_pair={} #u->x->v
		self.s_pair={} #u<-x->v
		self.t_pair={} #u->x<-v
		self.treepath={}
		self.generatorCycles=[]
	def find(self,u):
		if u not in self.Root:
			return u
		return self.find(self.Root[u])

	def findRoot(self, u):
		return self.Root[u]

	def PrintPair(self):
		print('pair', self.pair)
		print('number of cycles:', len(self.EdgeCycle))
		print('number of nuntreeedge:', len(self.EdgeId))
		

	def PrintCycle(self):
		print(self.HomoCycle)
	def perHom(self, t): #the function to compute persistent path homology, t is the time
		self.Root={}
		self.TreeEdge=set()
		self.EdgeId={}
		self.NonTreeEdge=[]
		self.EdgeCycle={}
		self.HomoCycle={}
		self.Bnd=[]
		self.generators=[]
		self.pair=[]
		self.curG=nx.DiGraph()
		self.curG.add_nodes_from(self.G.nodes)
		self.ns_nt_pair={} #u->x->v
		self.generatorCycles=[]
		self.s_pair=collections.defaultdict(list) #u<-x->v
		self.t_pair=collections.defaultdict(list) #u->x<-v
		# t is the time
		edgeSet=collections.defaultdict(list)
		bnd=[]
		ws=set()
		#sort edges
		for e in self.G.edges:
			w=self.G.edges[e]['weight']
			if w<=t:
				edgeSet[w].append(e)
				ws.add(w)
		#edgeSet.sort(key=lambda x:self.G.edges[x]['weight'])
		wsl=[w for w in ws]
		wsl.sort()
		cnt=0
		for weight in wsl:
			
			last_cnt=cnt
			edges=edgeSet[weight]
			s=[]
			tmpcycle=set()
			bi=[]
			tri=[]
			quad=[]
			anothercnt=0
			for e in edges:
				self.curG.add_edge(e[0], e[1])
				#print(e)
				u=e[0]
				v=e[1]
				ru=self.find(u)
				rv=self.find(v)
				if ru!=rv: #this edge is a tree edge, they will not form a cycle
					#anothercnt+=1
					
					self.Root[rv]=ru
					self.Root[v]=ru
					# if not self.parent[v]:
					# 	self.parent[v]=u
					#print(anothercnt, u, v, ru, rv, self.Root[rv], self.Root[v])
					self.TreeEdge.add(e)
					
				else:#nontree edge case-------which means a cycle
					if ru!=rv:
						print('not equal')
					
					#lu=self.findrootlist(u)
					#lv=self.findrootlist(v)
					self.EdgeId[(u,v)]=cnt
					self.NonTreeEdge.append((u, v))
					tmpcycle.add(cnt)
					cnt+=1
					
					#bi-gon
					if (v, u) in self.curG.edges():
						tmp=[0 for i in self.EdgeId]
						tmp[cnt-1]=1
						if (v,u) in self.EdgeId:
							tmp[self.EdgeId[(v, u)]]=-1
						s.append(tmp)
						bi.append(tmp)
					#triangle
					if (u, v) in self.ns_nt_pair:
						time, e1, e2=self.ns_nt_pair[(u, v)] #u->x->v
						a=[]
						tmp=[0 for i in self.EdgeId]
						tmp[cnt-1]=1
						if e1 in self.EdgeId:
							tmp[self.EdgeId[e1]]=-1
						if e2 in self.EdgeId:
							tmp[self.EdgeId[e2]]=-1
						
						s.append(tmp)
						tri.append(tmp)

					if (u, v) in self.s_pair:
						#print(self.s_pair[(u, v)])
						for (time, e1, e2) in self.s_pair[(u, v)]: #u<-x->v
							tmp=[0 for i in self.EdgeId]
							tmp[cnt-1]=1
							if e1 in self.EdgeId:
								tmp[self.EdgeId[e1]]=1
							if e2 in self.EdgeId:
								tmp[self.EdgeId[e2]]=-1
						
							s.append(tmp)
							tri.append(tmp)

					if (u, v) in self.t_pair:
						for (time, e1, e2) in self.t_pair[(u, v)]: #u->x<-v
							tmp=[0 for i in self.EdgeId]
							tmp[cnt-1]=1
							if e1 in self.EdgeId:
								tmp[self.EdgeId[e1]]=-1
							if e2 in self.EdgeId:
								tmp[self.EdgeId[e2]]=1
							s.append(tmp)
							tri.append(tmp)
					#####Quad

					for w in self.curG.successors(v): #u->x->w
						if w==u or (u, w) not in self.ns_nt_pair:
							continue
						time, e1, e2=self.ns_nt_pair[(u, w)]
						tmp=[0 for i in self.EdgeId]
						tmp[cnt-1]=1
						if (v, w) in self.EdgeId:
							tmp[self.EdgeId[(v, w)]]=1
						if e1 in self.EdgeId:
							tmp[self.EdgeId[e1]]=-1
						if e2 in self.EdgeId:
							tmp[self.EdgeId[e2]]=-1
						s.append(tmp)
						quad.append(tmp)
						# self.ns_nt_pair[(u,w)]=(weight, (u, v), (v, w))

					for w in self.curG.predecessors(u): #w->x->v
						if w==v or (w, v) not in self.ns_nt_pair:
							continue
						time, e1, e2=self.ns_nt_pair[(w, v)]
						tmp=[0 for i in self.EdgeId]
						tmp[cnt-1]=1
						if (w, u) in self.EdgeId:
							tmp[self.EdgeId[(w, u)]]=1
						if e1 in self.EdgeId:
							tmp[self.EdgeId[e1]]=-1
						if e2 in self.EdgeId:
							tmp[self.EdgeId[e2]]=-1
						s.append(tmp)
						quad.append(tmp)
						# self.ns_nt_pair[(w, v)]=(weight, (w, u), (u, v))
					# for w in self.curG.predecessors(v):
					# 	if w!=u:
					# 		self.t_pair[(u, w)]=(weight, (u, v), (w, v))
					# 		self.t_pair[(w, u)]=(weight, (w, v), (u, v))
					# for w in self.curG.successors(u):
					# 	if w!=v:
					# 		self.s_pair[(w, v)]=(weight, (u, w), (u, v))
					# 		self.s_pair[(v, w)]=(weight, (u, v), (u, w))


				for w in self.curG.predecessors(v): #u->v<-w
					if w!=u:
						self.t_pair[(u, w)].append((weight, (u, v), (w, v)))
						self.t_pair[(w, u)].append((weight, (w, v), (u, v)))
				for w in self.curG.successors(u):#w<-u->v
					if w!=v:
						self.s_pair[(w, v)].append((weight, (u, w), (u, v)))
						self.s_pair[(v, w)].append((weight, (u, v), (u, w)))
				for w in self.curG.predecessors(u): #w->u->v
					self.ns_nt_pair[(w, v)]=(weight, (w, u), (u, v))
				for w in self.curG.successors(v):
					self.ns_nt_pair[(u, w)]=(weight, (u, v), (v, w))

			if not s:
				for e in tmpcycle:
					self.EdgeCycle[self.NonTreeEdge[e]]=(weight, 0)#####0 means nontrivial cycle has not been paired
				continue
			#self.EdgeCycle[(u,v)]=(weight, 2) #####2 means trivial cycle
			for i in range(len(s)):
				tmp=[0 for j in range(cnt)]
				tmp[:len(s[i])]=s[i]
				s[i]=tmp
			for i in range(len(bi)):
				tmp=[0 for j in range(cnt)]
				tmp[:len(bi[i])]=bi[i]
				bi[i]=tmp
			for i in range(len(tri)):
				tmp=[0 for j in range(cnt)]
				tmp[:len(tri[i])]=tri[i]
				tri[i]=tmp
			for i in range(len(quad)):
				tmp=[0 for j in range(cnt)]
				tmp[:len(quad[i])]=quad[i]
				quad[i]=tmp
				#print(s[i])
			A=matrix(s, numpy.float32)
			m=A.shape[0]

			for i in range(m):
				bnd.append(A[i])
			
			for i, (tpivot, tind) in enumerate(self.Bnd):
				#print(tpivot, tind, cnt)
				tind=tind.tolist()[0]
				tmp=[0 for j in range(cnt)]
				tmp[:len(tind)]=tind
				#tind=numpy.resize(tind, (1, cnt))
				#=tind=numpy.pad(tind[0], (0,cnt-tind.shape[1]), 'constant')
				#print(tind, 1)
				self.Bnd[i]=(tpivot, matrix(tmp))
				#print(self.Bnd)
			#print(A)
			I=set(range(m))
			pivot=[None for i in range(m)]
			#A=numpy.float_(A)
		

			for i in range(m):
				for tpivot, tind in self.Bnd:
					if A[i, tpivot]!=0:
						A[i]=A[i]-tind*1.0/tind[0,tpivot]*A[i, tpivot]
					
			
			while I:
				i=I.pop()
				
				j=None
				
				for k in range(A.shape[1]-1, -1, -1):
					if A[i,k]!=0:
						j=k
						break
				
				if j!=None:
					pivot[i]=j
					#print(111)
					for k in range(m):
						if k!=i and A[k,j]!=0:
							A[k]=A[k]-A[i]*1.0/A[i,j]*A[k,j]
						


					for tpivot, tind in self.Bnd:	
						if tind[0, j]!=0:
							tind=tind-A[i]*1.0/A[i, j]*tind[0, j]
						
			
			for i in range(len(pivot)):	###every cycle die at weight
				if pivot[i]!=None and pivot[i]>=last_cnt:##meaning the cycle born now and die now
					edge=self.NonTreeEdge[pivot[i]]

					self.EdgeCycle[edge]=(weight, 2)
					tmpcycle.remove(pivot[i])
					self.Bnd.append((pivot[i], A[i]))
				elif pivot[i]!=None and pivot[i]<last_cnt:
					edge=self.NonTreeEdge[pivot[i]]
					birthtime, status=self.EdgeCycle[edge]

					deathtime=weight
					#tmpcycle.remove(pivot[i])
					self.pair.append((birthtime, deathtime))
					self.EdgeCycle[edge]=(birthtime, 1)
					self.Bnd.append((pivot[i], A[i]))
					self.generators.append(A[i])
					# tt=[]
					# for a,b in self.Bnd:
					# 	#print(b)
					# 	tt.append(b.tolist()[0])
					# #print(tt)
					# q1=matrix_rank(matrix(tt))
					# q2=len(self.pair)
					# print(q1, q2)
					# if q1!=q2:
					# 	return
			for i in tmpcycle:
				self.EdgeCycle[self.NonTreeEdge[i]]=(weight, 0)

		homocount=0
		homo1=0
		homo2=0
		homo0=0
		for edge in self.EdgeCycle:
			birthtime,status=self.EdgeCycle[edge]
			if status==0:
				homo0+=1
			elif status==1:

				homo1+=1
			else:
				homo2+=1

		for edge in self.EdgeCycle:

			birthtime,status=self.EdgeCycle[edge]
			if status==0:
				self.pair.append((birthtime, t))
				homocount+=1
				self.HomoEdgeId.append(edge)
				self.EdgeCycle[edge]=(birthtime, 1)

		print("homocount", homocount)
		print("homo0:", homo0)
		print("homo1:", homo1)
		print("homo2:", homo2)

	def findCycle(self, filename):#####canbe done only when there is only one component
		treeNodes=set(self.G.nodes)
		lev={}
		l=list(self.TreeEdge)
		lp=[0 for i in range(len(l))]
		while treeNodes:
			root=treeNodes.pop()
			#print(root)
			lev[root]=(root,0,0, root)
			q=collections.deque()
			q.append(root)
			while q:
				v=q.popleft()
				_, level, _, r=lev[v]
				for i in range(len(l)):
					if lp[i]==1:
						continue
					

					e=l[i]
					#print(e, v)
					if e[0]==v:
						lp[i]=1
						q.append(e[1])
						lev[e[1]]=(v,level+1, 1, root)
						treeNodes.remove(e[1])
					elif e[1]==v:
						lp[i]=1
						q.append(e[0])
						lev[e[0]]=(v,level+1, -1, root)
						treeNodes.remove(e[0])
			#print(len(treeNodes))
		print(lev)
		for v in self.generators:
			edgecycle=set()
			for index in range(v.shape[1]):
				if v[0,index]==0:
					continue
				
				edge=self.NonTreeEdge[index]
				
				
				v1=edge[0]
				v2=edge[1]
				if v1==v2:
					print('wrong')
				_, d1, _, r1=lev[v1]
				_, d2, _, r2=lev[v2]
				if r1!=r2:
					print('wrong')
				
				if edge not in edgecycle:
					edgecycle.add(edge)
				else:
					edgecycle.remove(edge)
				
				while d1>d2:
					parent, d1, _in, r1=lev[v1]
					#print(d1)
					if v1==parent:
						print(v1, parent1, edgecycle)
					if _in==-1:
						if (v1, parent) not in edgecycle:
							edgecycle.add((v1, parent))
						else:
							edgecycle.remove((v1, parent))
						
					elif _in==1:
						if (parent, v1) not in edgecycle:
							edgecycle.add((parent, v1))
						else:
							edgecycle.remove((parent, v1))
						
					v1=parent
					d1-=1
				while d1<d2:
					parent, d2, _in, r=lev[v2]
					#print(d1)
					if v2==parent:
						print(v2, parent2, edgecycle)
					if _in==-1:
						if (v2, parent) not in edgecycle:
							edgecycle.add((v2, parent))
						else:
							edgecycle.remove((v2, parent))
						
					elif _in==1:
						if (parent, v2) not in edgecycle:
							edgecycle.add((parent, v2))
						else:
							edgecycle.remove((parent, v2))
					
					v2=parent
					d2-=1

				while v1!=v2 and d1>0:
					parent1, d1, _, r=lev[v1]
					if v1==parent1:
						print(v1, parent1, edgecycle)
					if _in==-1:

						if (v1, parent1) not in edgecycle:
							edgecycle.add((v1, parent1))
						else:
							edgecycle.remove((v1, parent1))
						
					elif _in==1:
						if (parent1, v1) not in edgecycle:
							edgecycle.add((parent1, v1))
						else:
							edgecycle.remove((parent1, v1))
					
					parent2, d2, _, r=lev[v2]
					if v2==parent2:
						print(v2, parent2, edgecycle)
					if _in==-1:
						if (v2, parent2) not in edgecycle:
							edgecycle.add((v2, parent2))
						else:
							edgecycle.remove((v2, parent2))
						
					elif _in==1:
						if (parent2, v2) not in edgecycle:
							edgecycle.add((parent2, v2))
						else:
							edgecycle.remove((parent2, v2))
					v1,v2=parent1,parent2
				#print(edgecycle)
			self.generatorCycles.append(edgecycle)
		file=open(filename,"w")
			#print(self.MinHomoBasis)
		for cycle in self.generatorCycles:
			for e in cycle:
				file.write('(')
				file.write(str(e[0]))
				file.write('\t')
				file.write(str(e[1]))
				file.write(')')
				file.write('\t')
			file.write('\n')
		file.close()


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
# print(g.Root)
# stop = timeit.default_timer()
# g.PrintPair()
# print('Time: ', stop - start)

# G=nx.DiGraph()
# G=nx.read_weighted_edgelist("citeseer", create_using = nx.DiGraph())
# print(len(G.edges))
# g=persHomo(G)
# start = timeit.default_timer()
# g.perHom(2)
# stop = timeit.default_timer()
# print('Time: ', stop - start)
# print(len(g.TreeEdge), len(g.curG.edges), len(g.Bnd))
# g.PrintPair()
# g.PrintCycle()