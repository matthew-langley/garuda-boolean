# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 16:57:24 2017

@author: Matthew
"""

import networkx as nx
import numpy as np

# Make a networkx representation of the state transition graph
G = nx.read_gml('test/ayako-boolean-psc-jun2017-output-2iL.gml')

# CLear any existing membership annotations
for n in G.nodes():
    if 'membership' in G.node[n]:
        del G.node[n]['membership']

# Determine the SCCs and steady states in the graph
sccGenerator = nx.strongly_connected_component_subgraphs(G)
scc = None
for x in sccGenerator:
    if len(x.nodes()) > 1:
        scc = x
        
### CALCULATE NODE PROBABILITIES ###

for node in scc:
    TotalSum = sum([float(scc[s][t]['weight']) for s,t in scc.edges_iter(node)])
    #print "node,total passing in the SCC%s= %s,%s" %(SCCnum,node,TotalSum)
    for [s,t] in list(scc.edges_iter(node)): #weight = (times of passing)/(total passing)
        scc[s][t]['internalweight'] = float(scc[s][t]['weight'])/float(TotalSum)

N = len(scc.nodes())
left = np.ones([N,N], dtype=float)
np.fill_diagonal(left, 0.0)
nodes = scc.nodes()
for i in range(N):
    node_t = nodes[i]
    node_s_list = list(scc.predecessors(node_t))
    for node_s in node_s_list:
        j = nodes.index(node_s)
        left[i,j] += scc[node_s][node_t]['internalweight']

right = np.array(np.ones((len(scc),1)))
nodeProbabilities = np.linalg.solve(left,right)