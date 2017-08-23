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


def read_parameters(args, paramFile):
    """
    Reads command line arguments specifiying parameters for SCC analysis 
    from a specified file
    
    Args:
        args (argparse.Namespace)
            Current arguments object from argparse
        paramFile (str)
            Filepath to parameters file
    
    Returns:
        (argparse.Namespace)
            Updated arguments object, with arguments added as specified in 
            paramFile
    """
    
    with open(paramFile, 'r') as f:
        for line in f:
            # Ignore any comments and blank lines
            if len(line.strip()) > 0 and (not line.startswith('#')):
                param, value = [x.strip() for x in line.split('=')]
                if param == 'minSize':
                    args.minSize = int(value)
                elif param == 'minSustainability':
                    args.minSustainability = float(value)
                elif param == 'writeSubgraphs':
                    if value in ['True', 'true', 'T', 't', '1', 'Yes', 'Y', 'y']:
                        args.writeSubgraphs = True
                    else:
                        args.writeSubgraphs = False
                elif param == 'annotateGraph':
                    if value in ['True', 'true', 'T', 't', '1', 'Yes', 'Y', 'y']:
                        args.annotateGraph = True
                    else:
                        args.annotateGraph = False
    
    return args


def calculate_node_probabilities(scc):
    """
    Calculates the stationary distribution of the SCC (i.e. the limit of 
    probability after an infinite-step random walk on the SCC) and returns 
    an array of the probability of each node. Notably, only the internal 
    edges to the SCC are considered; edges which lead to exit from the SCC are 
    excluded.
    
    Args:
        scc (nx.DiGraph)
            The SCC/steady state subgraph
            
    Returns:
        (np.array)
            Array of node probabilities in stationary distribution of SCC
    """
    
    # Check if SCC is a steady state, and return dictionary {node: 1.0} if yes
    if len(scc.nodes()) == 1:
        return np.array([1.0])
        
    else:
        
        # Start by normalizing the edges in the SCC, and save as 
        # 'internal_weight' attribute.
        # Needed because we're calculating an *internal* SCC node probability,
        # considering only edges within SCC. The edge weights on the argument 
        # 'scc' are weighted by edges external to the SCC as well. Thus, the
        # edges need to be recalculated.
        for node in scc:
            TotalSum = sum([float(scc[s][t]['weight']) for s,t in scc.edges_iter(node)])
            #print "node,total passing in the SCC%s= %s,%s" %(SCCnum,node,TotalSum)
            for [s,t] in list(scc.edges_iter(node)): #weight = (times of passing)/(total passing)
                scc[s][t]['internalweight'] = float(scc[s][t]['weight'])/float(TotalSum)
        
        # Solve simultaneous equations based on array of node probabilities
        # (+1.0) at the end of equation -> summation restriction
        # Basic calculation:
        # if node_s in scc.predecessors(node_t):                    
        #   if node_s == node_t -> move from right term to left term
        #     equation.append(scc[node_s][node_t]['internalweight'])#i[node_s][node_t]-1.0+1.0
        #   else:
        #     equation.append(scc[node_s][node_t]['internalweight']+1.0)
        # else:
        #   if node_s == node_t -> move from right term to left term
        #     equation.append(0.0) #(-1.0+1.0)
        #   else:
        #     equation.append(1.0) #(0.0+1.0)
        # equationList.append(equation)
        # Code section below optimizes this calculation by pre-initializing
        # the matrix 'left' with {1.0 if off-diagonal, 0.0 if on-diagonal},
        # then adding edge weights.
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
        return nodeProbabilities


def calculate_expression_profile(scc, nodeProbabilities=None):
    """
    Returns the (average) gene expression profile of an attractor, either a 
    steady state or an SCC. In the case of a steady state, the gene expression 
    profile is simply the vector of all gene values for that state. In the 
    case of SCCs, the expression profile is the probability-weighted average 
    of gene expression values over all states (nodes) in the SCC.
    NOTE: By default, the expression profile is not weighted by the 
    sustainability score of the SCC (i.e. closed SCC evaluation).
    
    Args:
        scc (nx.DiGraph)
            The SCC/steady state subgraph, in which each entry in its 
            corresponding dictionary corresponds to a gene and its expression 
            value
        nodeProbabilities (np.array)
            The stationary probabilities of each node in the SCC, in order of 
            scc.nodes(). If none provided, the function will internally call 
            'calculate_node_probabilities(scc)' to calculate it.
    
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
        
        # Calculate node probabilities if None provided as argument.
        if nodeProbabilities is None:
            nodeProbabilities = calculate_node_probabilities(scc)
        
        # Calculate the probability-weighted average
        nodeProfiles = pd.DataFrame([pd.Series(vals, name=n, dtype=float) for n, vals in scc.nodes(data=True)]).transpose()
        averageExpr = np.dot(nodeProfiles.as_matrix(), nodeProbabilities)
        return pd.Series(averageExpr.flatten(), index=nodeProfiles.index, name=scc.graph['name'], dtype=float)
        

def calculate_sustainability(G, scc, nodeProbabilities=None):
    """
    Calculates the sustainability score of an SCC.
    
    Args:
        g (nx.DiGraph)
            The complete state transition graph (needed to calculate outgoing 
            trajectories)
        scc (nx.DiGraph)
            The SCC/steady state subgraph
        nodeProbabilities (np.array)
            The stationary probabilities of each node in the SCC, in order of 
            scc.nodes(). If none provided, the function will internally call 
            'calculate_node_probabilities(scc)' to calculate it.
    
    Returns:
        (float)
            Sustainability score
    """

    # Calculate node probabilities if None provided as argument.
    if nodeProbabilities is None:
        nodeProbabilities = calculate_node_probabilities(scc)
    
    sccSet = set(scc)
    totaloutweight = 0.0
    for i in range(len(scc.nodes())):
        s = scc.nodes()[i]
        edges = {t: v['weight'] for (t,v) in G[s].items()}
        pass_count = sum(edges.values())
        internal_count = sum(v for t,v in edges.items() if t in sccSet)
        outweight = float(nodeProbabilities[i]) * float(internal_count) / float(pass_count)
        totaloutweight += outweight
    
    # Return total outweight, accounting for possible float rounding errors 
    # that would lead to slightly higher sustainability values than 1.
    return min([totaloutweight, 1.0])
        

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
    parser.add_argument('parameters',
                        nargs='?',
                        help='(optional) File containing simulation parameters described below')
    parser.add_argument('--profilesOutput',
                        help='Output file destination for SCC expression profiles')
    parser.add_argument('--metricsOutput',
                        help='Output file destination for SCC metrics, such as size and sustainability')
    parser.add_argument('--minSize',
                        dest='minSize',
                        type=int,
                        default=0,
                        help='Minimum size of SCCs to be considered, in terms of numbers of profiles')
    parser.add_argument('--minSustainability',
                        dest='minSustainability',
                        type=float,
                        default=0.0,
                        help='Minimum sustainability score of SCCs to be considered')
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
        # myArgs = 'test/ES_2iL+B-A_output.gml --output test/ES_2iL+B-A_output_scc_profiles.csv --annotateGraph'
        myArgs = 'test/test_graph_simple.gml --profilesOutput test/test_graph_simple_scc_profiles.csv --metricsOutput test/test_graph_simple_scc_metrics.csv --annotateGraph --writeSubgraphs'
        # myArgs = 'test/2017-06-30-LS/ayako-boolean-psc-jun2017-boolean-functions_output.gml --output test/2017-06-30-LS/ayako-boolean-psc-jun2017-boolean-functions_output_scc_profiles.csv'        
        # myArgs = 'test/ayako-boolean-psc-jun2017-output-2iL.gml test/test_scc_params.txt'
        args = parser.parse_args(myArgs.split())
    else:
        args = parser.parse_args()
        
    # If a parameter file is specified, read in parameters from there.
    if args.parameters is not None:
        args = read_parameters(args, args.parameters)
        
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
            # Only keep SCCs larger than the minimum size (user-specified)
            if len(scc.nodes()) > args.minSize:
                # Calculate node probabilities and sustainability
                nodeProbabilities = calculate_node_probabilities(scc)
                sustainability = calculate_sustainability(G, scc, nodeProbabilities)
                # Only keep SCCs with sustainability larger than minimum 
                # threshold (user-specified)
                if sustainability > args.minSustainability:
                    scc.graph['name'] = 'SCC' + str(len(sccList) + 1)
                    scc.graph['sustainability'] = sustainability
                    expressionProfile = calculate_expression_profile(scc, nodeProbabilities)
                    scc.graph['expression'] = expressionProfile
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
                    scc.graph['sustainability'] = 1.0  # By definition
                    expressionProfile = calculate_expression_profile(scc)
                    scc.graph['expression'] = expressionProfile
                    ssList.append(scc)
    
    # Prepare to output to file
    graphName = os.path.basename(args.graph).rsplit('.', 1)[0]
    if args.profilesOutput is None:
        args.profilesOutput = '../Output/' + graphName + '_scc_expression.csv'
    if args.metricsOutput is None:
        args.metricsOutput = '../Output/' + graphName + '_scc_metrics.csv'   
    dirName = os.path.dirname(args.profilesOutput).replace('\\', '/')
    
    # Output expression profiles and metrics of each SCC and SS to file
    if len(sccList) > 0:
        # Calculate average expression profile over all SCCs in the model.
        # Weight the contribution of each SCC based on their sustainability and 
        # their size.
        # NOTE: Steady states excluded from this calculation
        modelAverageProfile = pd.Series(sum([scc.graph['expression'] * scc.graph['sustainability'] * len(scc.nodes()) for scc in sccList]) / sum([len(scc.nodes()) for scc in sccList]),
                                        name='Weighted Average')
        sccProfiles = pd.concat([scc.graph['expression'] for scc in sccList + ssList], axis=1)
        expressionProfiles = pd.concat([sccProfiles, modelAverageProfile], axis=1).transpose()
        expressionProfiles.to_csv(args.profilesOutput)
        # Output metrics of each SCC and SS to file
        metrics = pd.DataFrame([pd.Series({'Size': len(scc.nodes()), 
                                           'Edges': len(scc.edges()), 
                                           'Sustainability': scc.graph['sustainability']}, 
                                name=scc.graph['name']) for scc in sccList + ssList])
        metrics.to_csv(args.metricsOutput)
    else:
        # Write blank csv files if no SCCs or steady states
        with open(args.profilesOutput, 'w') as f:
            f.write('No SCCs or steady states matching criteria were found.')
        with open(args.metricsOutput, 'w') as f:
            f.write('No SCCs or steady states matching criteria were found.')
    
    # Write GML for each SCC and SS if desired
    if args.writeSubgraphs is True:
        for scc in sccList + ssList:
            del scc.graph['expression'] # Not supported by GML format
            gmlOutputFile = dirName + '/' + graphName + '_' + scc.graph['name'] + '.gml'
            nx.write_gml(scc, gmlOutputFile)
    
    # Annotate original graph file with SCC and SS membership if desired, 
    # and pass along chain regardless
    if args.annotateGraph is True:
        for scc in sccList + ssList:
            for n in scc.nodes():
                G.node[n]['membership'] = scc.graph['name']
    passthroughGraphFile = dirName + '/' + graphName + '.gml'
    nx.write_gml(G, passthroughGraphFile, stringizer=str)
    