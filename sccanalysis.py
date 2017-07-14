# 
# Copyright 2017 Matthew Langley. All rights reserved.
#


"""
Analysis of state transition graphs arising from Boolean network simulations, 
including identification of steady states and strongly connected components
"""

import argparse
import networkx as nx
import numpy as np
import os
import pandas as pd


def expression_profile(scc):
    """
    Returns the (average) gene expression profile of an attractor, either a 
    steady state or an SCC. In the case of a steady state, the gene expression 
    profile is simply the vector of all gene values for that state. In the 
    case of SCCs, the expression profile is the weighted average of gene 
    expression values over all states (nodes) in the SCC. The weighting 
    corresponds to each state's probability in the stationary distribution of 
    the SCC (i.e. the limit of probability after an infinite-step random walk 
    on the SCC). Notably considering only the internal edges to the SCC are 
    considered; edges which lead to exit from the SCC are excluded.
    
    Args:
        scc (nx.DiGraph)
            The SCC/steady state subgraph, in which each entry in its 
            corresponding dictionary corresponds to a gene and its expression 
            value
    
    Returns:
        (pd.Series)
            Gene expression profile
    """
    
    # If steady state, convert dictionary to Series and return
    if len(scc.nodes()) == 1:
        n = scc.nodes()[0]
        return pd.Series(scc.node[n], dtype=float, name=scc.graph['name'])
        
    # Else calculate the node probabilities of SCC and do the weighted average
    else:
        
        # Start by normalizing the edges in the SCC
        # Needed because we're calculating an *internal* SCC node probability,
        # considering only edges within SCC. The edge weights on the argument 
        # 'scc' are weighted by edges external to the SCC as well. Thus, the
        # edges need to be recalculated.
        for node in scc:
            TotalSum = 0
            for [s,t] in list(scc.edges_iter(node)): #list of edges of "starting from the node"
                TotalSum = TotalSum + float(scc[s][t]['weight'])
            #print "node,total passing in the SCC%s= %s,%s" %(SCCnum,node,TotalSum)

            for [s,t] in list(scc.edges_iter(node)): #weight = (times of passing)/(total passing)
                scc[s][t]['weight'] = float(scc[s][t]['weight'])/float(TotalSum)        
        
        equationList = []
        count_target = 0
        for node_t in scc:
            equation = []
            count_source = 0
            ### (+1.0) at the end of equation -> summation restriction
            for node_s in scc:  #create equation for target node for every possible pred.(if not: insert 0)
                if node_s in set(scc.predecessors(node_t)):                    
                    if count_target == count_source: # "= node_t" -> move from right term to left term
                        equation.append(scc[node_s][node_t]['weight'])#i[node_s][node_t]-1.0+1.0
                    else:
                        equation.append(scc[node_s][node_t]['weight']+1.0)
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
        right = np.array(np.ones((len(scc),1)))
        nodeProbabilities = np.linalg.solve(left,right)
        
        # Calculate the probability-weighted average
        nodeProfiles = pd.DataFrame([pd.Series(vals, name=n, dtype=float) for n, vals in scc.nodes(data=True)]).transpose()
        averageExpr = np.dot(nodeProfiles.as_matrix(), nodeProbabilities)
        return pd.Series(averageExpr.flatten(), index=nodeProfiles.index, name=scc.graph['name'], dtype=float)
        

if __name__ == '__main__':
    """
    Code runs only if module is executed as mainline. Supports command line 
    arguments.
    """
    
    # Configure command line arguments
    description = ('Find strongly connected components (SCCs) and steady '
                   'states within a state transition graph, such as those '
                   'produced by Boolean network simulation, and calculates '
                   'their expression profile')
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('graph',
                        help='State transition graph in GML format')
    parser.add_argument('--output',
                        help='Output file destination for SCC expression profiles')
    parser.add_argument('--writeSubgraphs',
                        dest='writeSubgraphs',
                        action='store_true',
                        default=False,
                        help='Write GML representations of identified SCCs and SSs to file')
    parser.add_argument('--annotateGraph',
                        dest='annotateGraph',
                        action='store_true',
                        default=False,
                        help='Annotate nodes in input graph with SCC/SS membership')
                           
    # Switch between listening to command-line arguments or hard-coded 
    # arguments depending on whether running in IDE or from cmd.
    if any([name.startswith('SPYDER') for name in os.environ]):
        # myArgs = 'test/test_output.gml --output test/test_output_scc_profiles.csv --annotateGraph'
        # myArgs = 'test/test_graph_simple.gml --output test/test_graph_simple_scc_profiles.csv --annotateGraph'
        # myArgs = 'test/ES_2iL+B-A_output.gml --output test/ES_2iL+B-A_output_scc_profiles.csv --annotateGraph'
        myArgs = 'test/test_graph_simple.gml --output Output/test_graph_simple_scc_profiles.csv --annotateGraph'
        args = parser.parse_args(myArgs.split())
    else:
        args = parser.parse_args()
        
    # Make a networkx representation of the state transition graph
    G = nx.read_gml(args.graph)
    
    # CLear any existing membership annotations
    for n in G.nodes():
        if 'membership' in G.node[n]:
            del G.node[n]['membership']
    
    # Determine the SCCs and steady states in the graph
    sccGenerator = nx.strongly_connected_component_subgraphs(G)
    sccList = []
    ssList = []
    for scc in sccGenerator:
        if len(scc.nodes()) > 1:
            scc.graph['name'] = 'SCC' + str(len(sccList) + 1)
            sccList.append(scc)
        else:
            # NB: Not all SCCs of size = 1 are steady states! Transition states
            # count as SCCs too in the purest graph-theory sense.
            # http://www.geeksforgeeks.org/strongly-connected-components/
            # So need to check if SCC consists of only one node with only one
            # (self-edge). i.e. set of 
            if scc.edges() != []:
                sourceNode = scc.edges()[0][0]
                if set([sourceNode]) == set(G.successors(sourceNode)):
                    scc.graph['name'] = 'SS' + str(len(ssList) + 1)
                    ssList.append(scc)
    print 'Found %i strongly connected component(s) and %i steady states' %(len(sccList), len(ssList))
    
    # Calculate the expression profiles of each SCC and SS
    expressionProfiles = pd.DataFrame([expression_profile(s) for s in ssList + sccList])
    if args.output is None:
        args.output = '../Output/' + args.graph.rsplit('.', 1)[0] + '_scc_expression.csv'
    expressionProfiles.to_csv(args.output)
    
    # Write GML for each SCC and SS if desired
    if args.writeSubgraphs is True:
        for scc in sccList + ssList:
            if '\\' in args.graph:
                args.graph = args.graph.replace('\\', '/')
            if '/' in args.graph:
                graphName = args.graph.rsplit('/', 1)[1].rsplit('.', 1)[0]
            else:
                graphName = args.graph.rsplit('.', 1)[0]
            gmlOutputFile = args.output.rsplit('/', 1)[0] + '/' + graphName + '_' + scc.graph['name'] + '.gml'
            nx.write_gml(scc, gmlOutputFile)
    
    # Annotate original graph file with SCC and SS membership if desired, 
    # and pass along chain regardless
    if args.annotateGraph is True:
        for scc in sccList + ssList:
            for n in scc.nodes():
                G.node[n]['membership'] = scc.graph['name']
    passthroughGraphFile = args.output.rsplit('/', 1)[0] + '/' + args.graph.rsplit('/', 1)[-1]
    nx.write_gml(G, passthroughGraphFile, stringizer=str)
