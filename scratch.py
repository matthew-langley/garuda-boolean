# 
# Copyright 2017 Matthew Langley. All rights reserved.
#


"""
Scratchpad for testing out new snippets of code.
"""

import networkx as nx
import numpy as np

def calc_node_probabilities(scc):
    """
    Calculates the steady-state probability of reaching each node in the SCC, 
    by simulating a chain of random walks along the SCC graph. This is a less
    computationally-intensive method than solving the associated eigenvalue
    problem, Ax = x. Returns a dictionary of the probabilities, keyed by node
    name. Note that this method also assign the probability of each node as 
    an attribute to the node itself, i.e. scc.nodes()[i]['probability']
    """
    
    # Determine the weight of each edge in the SCC, based on how many 
    # times it was passed during the simulation. The edge weights are 
    # normalized, such that the sum of the edge weights leaving each 
    # node is 1.
    # Note: Attempt to speed up this part by converting scc.edges()
    # into a set.
    # for (s,t) in scc.edges():
    #     scc[s][t]['weight'] = edgeDict[(s,t)]
    for node in scc.nodes():
        totalSum = 0.0
        for (s,t) in scc.edges(node):
            totalSum += scc[s][t]['weight']
        for (s,t) in scc.edges(node):
            scc[s][t]['weight'] /= totalSum

    # Construct the adjacency matrix for the SCC graph.
    # Note: row = target, column = source
    adjacencyMatrix = []
    for t in scc.nodes():
        adjacencyRow = []
        for s in scc.nodes():
            if s in scc.predecessors(t):
                adjacencyRow.append(scc[s][t]['weight'])
            else:
                adjacencyRow.append(0)
        adjacencyMatrix.append(adjacencyRow)

    # Let A be the adjacency matrix that was just constructed.
    # Applying A onto any vector, x, simulates a iteration of a walk
    # on the graph, such that Ax is the probability of arriving at
    # each node after one step from the initial conditions. It
    # follows that (A^n)x is the probability of each node after n
    # steps. We are looking for the steady-state probabilities
    # (i.e. the 'x' for which Ax = x). Rather than solving the
    # eigenvalue problem, which is computationally intensive for
    # large matrices, we instead find the steady-state numerically
    # by calculating (A^n)x for increasing values of n. We stop when
    # the value of (A^n)x equals (A^(n-1))x to 4 decimal places.
    A = np.matrix(adjacencyMatrix)
    x = np.matrix([[1.0]] * len(adjacencyMatrix)) / len(adjacencyMatrix)
    P_old = (A) * x
    P_new = (A**2) * x
    n = 0
    while (n < 1000 and not(np.array_equal(np.round(P_old, decimals=4), np.round(P_new, decimals=4)))):
        P_old = P_new
        P_new = A * P_old
        n += 1
    nodeProbabilities = P_new.tolist()
    for i in range(len(scc.nodes())):
        nodeID = scc.nodes()[i]
        scc.node[nodeID]['probability'] = nodeProbabilities[i][0]
    return scc, A, nodeProbabilities
    
def ayako_scc(i):
    
    equationList = []
    count_target = 0
    for node_t in i:

        #print 'equation for ----%s' %(node_t)

        equation = []
        count_source = 0
        ### (+1.0) at the end of equation -> summation restriction
        for node_s in i:  #create equation for target node for every possible pred.(if not: insert 0)

            if node_s in set(i.predecessors(node_t)):                    
                if count_target == count_source: # "= node_t" -> move from right term to left term
                    equation.append(i[node_s][node_t]['weight'])#i[node_s][node_t]-1.0+1.0
                else:
                    equation.append(i[node_s][node_t]['weight']+1.0)
            else:
                if count_target == count_source: # "= node_t" -> move from right term to left term
                    equation.append(0.0) #(-1.0+1.0)
                else:
                    equation.append(1.0) #(0.0+1.0)

            count_source += 1

        equationList.append(equation)
        count_target += 1

    ##### solve simultaneous equations ==== Array of Node Probability for (i.nodes())####
    left = np.array(equationList)
    right = np.array(np.ones((len(i),1)))
    #print sum(numpy.linalg.solve(left,right)) -> must be 1.0

    nodeprob_array = np.linalg.solve(left,right)
    return i, left, nodeprob_array
    

if __name__ == '__main__':
    
    G = nx.read_gml('test/test_graph_simple.gml')
    sccList = list(nx.strongly_connected_component_subgraphs(G, copy=True))
    
    ML_scc, ML_adjacencyMatrix, ML_nodeProbabilities = calc_node_probabilities(sccList[0])
    AY_scc, AY_adjacencyMatrix, AY_nodeProbabilities = ayako_scc(sccList[0])