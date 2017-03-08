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
        return pd.Series(scc.nodes[0], dtype=float)
        
    # Else calculate the node probabilities of SCC and do the weighted average
    else:
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
                           
    # Switch between listening to command-line arguments or hard-coded 
    # arguments depending on whether running in IDE or from cmd.
    if any([name.startswith('SPYDER') for name in os.environ]):
        myArgs = 'test/test_graph_simple.gml --output test/test_graph_scc_expression.csv --writeSubgraphs'
        args = parser.parse_args(myArgs.split())
    else:
        args = parser.parse_args()
        
    # Make a networkx representation of the state transition graph
    G = nx.read_gml(args.graph)
    
    # Determine the SCCs and steady states in the graph
    sccGenerator = nx.strongly_connected_component_subgraphs(G)
    sccList = []
    ssList = []
    for scc in sccGenerator:
        if len(scc.nodes()) == 1:
            scc.graph['name'] = 'SS' + str(len(ssList) + 1)
            ssList.append(scc)
        elif len(scc.nodes()) > 1:
            scc.graph['name'] = 'SCC' + str(len(sccList) + 1)
            sccList.append(scc)
    print 'Found %i strongly connected component(s) and %i steady states' %(len(sccList), len(ssList))
    
    # Calculate the expression profiles of each SCC and SS
    expressionProfiles = pd.DataFrame([expression_profile(s) for s in ssList + sccList])
    if args.output is None:
        args.output = '../Output/' + args.graph.rsplit('.', 1)[0] + '_scc_expression.csv'
    expressionProfiles.to_csv(args.output)
    
    # Write GML for each SCC and SS if desired
    if args.writeSubgraphs is True:
        for scc in sccList + ssList:
            gmlOutputFile = args.output.rsplit('/', 1)[0] + '/' + args.graph.rsplit('/', 1)[1].rsplit('.', 1)[0] + '_' + scc.graph['name'] + '.gml'
            nx.write_gml(scc, gmlOutputFile)
    
    