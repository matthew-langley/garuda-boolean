# 
# Copyright 2017 Matthew Langley. All rights reserved.
#


"""
Visualization of state transition graphs arising from Boolean network 
simulation.
"""

import argparse
import os
import pandas as pd
import requests
import seaborn as sns
from py2cytoscape.data import BASE_URL, HEADERS
from py2cytoscape.data.cyrest_client import CyRestClient
from py2cytoscape.data.util_network import NetworkUtil as util
from sklearn.decomposition import PCA as sklearnPCA
import traceback
import time
# Necessary workaround for pyinstaller since igraph loads this dynamically...
import igraph.vendor.texttable
import sklearn.neighbors.typedefs


def read_parameters(args, paramFile):
    """
    Reads command line arguments specifiying parameters for Boolean simulation 
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
                if value != '':
                    if param == 'layout':
                        args.layout = str(value)
                    elif param == 'refstates':
                        args.refstates = str(value)
                    elif param == 'coloursubgraphs':
                        if value in ['True', 'true', 'T', 't', '1', 'Yes', 'Y', 'y']:
                            args.coloursubgraphs = True
                        else:
                            args.coloursubgraphs = False
                    elif param == 'clearsession':
                        if value in ['True', 'true', 'T', 't', '1', 'Yes', 'Y', 'y']:
                            args.clearsession = True
                        else:
                            args.clearsession = False
                    elif param == 'outputLayout':
                        args.outputLayout = str(value)
                    elif param == 'outputImage':
                        args.outputImage = str(value)
    
    return args


def manual_rest(command, httpType='get', headers=None):
    """
    Performs a manual request to Cytoscape via the REST API. Allows access to 
    other parts of the REST API that aren't implemented in py2cytoscape Python 
    package.
    
    Args:
        command (str)
            Text of request command following the BASE (i.e. after '/v1/')
        httpType (str)
            Type of http request (choices: get, put, post, delete)
    Return:
        response (requests.models.Response)
            Response object from 'requests' package
    """
    
    myHeaders = HEADERS.copy()
    if headers is not None:
        myHeaders.update(headers)
        
    if httpType == 'get':
        return requests.get(BASE_URL + command, headers=myHeaders)
    elif httpType == 'put':
        return requests.put(BASE_URL + command, headers=myHeaders)
    elif httpType == 'post':
        return requests.post(BASE_URL + command, headers=myHeaders)
    elif httpType == 'delete':
        return requests.delete(BASE_URL + command, headers=myHeaders)
    else:
        raise ValueError('Unsupported httpType: %s' %(httpType))


def pca_locations(graph, xscale=1.0, yscale=1.0):
    """
    Computes a graph layout based on the principal components of all profiles
    
    Args:
        graph (py2cytoscape.data.cynetwork.CyNetwork)
            Cytoscape network object to layout
        xscale, yscale (float)
            Factor to multiply post-PCA (x,y) locations by to increase space 
            occupied by layout
    Return:
        (list of 3-tuples)
            (SUID, x, y) for each node in graph
    """
    
    # Read the graph's node table from Cytoscape
    nodeTable = graph.get_node_table()
    
    # Format the data frame to only genes shared across all nodes
    df = nodeTable.copy()
    df.index = df['name']
    columnsToRemove = ['shared name', 'name', 'selected', 'isExcludedFromPaths', 'datasource', 'membership']
    for column in columnsToRemove:
        df = df.drop(column, axis=1)
    df = df.dropna(axis=1, how='any')
    
    pca = sklearnPCA().fit(df.as_matrix())
    dfPCA = pd.DataFrame(pca.transform(df.as_matrix()), index=df.index)
    
    # Make the locations list
    idmap = util.name2suid(graph)
    locations = []
    for k in dfPCA.index:
        locations.append([int(idmap[k]),
                          dfPCA.loc[k, 0] * xscale,
                          dfPCA.loc[k, 1] * yscale])
    
    return locations
    

def reference_pca_locations(graph, xscale=1.0, yscale=1.0):
    """
    Computes a graph layout based on the principal components of the reference 
    profiles.
    
    Args:
        graph (py2cytoscape.data.cynetwork.CyNetwork)
            Cytoscape network object to layout
        xscale, yscale (float)
            Factor to multiply post-PCA (x,y) locations by to increase space 
            occupied by layout
    Return:
        (list of 3-tuples)
            (SUID, x, y) for each node in graph
    """
    
    # Read the graph's node table from Cytoscape
    nodeTable = graph.get_node_table()
    
    # Format the data frame to only genes shared across all nodes
    df = nodeTable.copy()
    df.index = df['name']
    columnsToRemove = ['shared name', 'name', 'selected', 'isExcludedFromPaths', 'membership']
    for column in columnsToRemove:
        df = df.drop(column, axis=1)
    df = df.dropna(axis=1, how='any')
    
    referenceData = df[df['datasource'] == 'reference']
    referenceData = referenceData.drop('datasource', axis=1)
    simulationData = df[df['datasource'] == 'simulation']
    simulationData = simulationData.drop('datasource', axis=1)
    
    pca = sklearnPCA().fit(referenceData.as_matrix())
    referencePCA = pd.DataFrame(pca.transform(referenceData.as_matrix()),
                                index=referenceData.index)
    simulationPCA = pd.DataFrame(pca.transform(simulationData.as_matrix()),
                                 index=simulationData.index)
    
    # Make the locations list
    idmap = util.name2suid(graph)
    locations = []
    for k in simulationPCA.index:
        locations.append([int(idmap[k]),
                          simulationPCA.loc[k, 0] * xscale,
                          simulationPCA.loc[k, 1] * yscale])
    for k in referencePCA.index:
        locations.append([int(idmap[k]),
                          referencePCA.loc[k, 0] * xscale,
                          referencePCA.loc[k, 1] * yscale])
    
    return locations
    

if __name__ == '__main__':
    
    # Connect to Cytoscape instance via REST
    try:
        cy = CyRestClient()
    except requests.ConnectionError:
        tb = traceback.format_exc()
        print tb        
        raise requests.ConnectionError('Cytoscape not started. Please start Cytoscape before running script.')
    
    # Configure layout choices
    originalLayouts = cy.layout.get_all()
    availableLayouts = list(originalLayouts)
    availableLayouts.append('PCA')
    availableLayouts.append('referencePCA')
    
    # Configure the command line interface
    parser = argparse.ArgumentParser(description='Visualize resulting state transition graph of Boolean network simulation',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('graph',
                        help='State transition graph to visualize')
    parser.add_argument('parameters',
                        nargs='?',
                        help='(optional) File containing visualization parameters described below')
    parser.add_argument('--layout',
                        dest='layout',
                        type=str,
                        default='circular',
                        choices=availableLayouts,
                        help='Graph layout method')
    parser.add_argument('--refstates',
                        dest='refstates',
                        type=str,
                        help='Reference states from experimental data')
    parser.add_argument('--coloursubgraphs',
                        dest='coloursubgraphs',
                        action='store_true',
                        default=False,
                        help='Assign node colours based on SCCs, steady states, etc.')
    parser.add_argument('--clearsession',
                        dest='clearsession',
                        action='store_true',
                        default=False,
                        help='Clear existing Cytoscape sessions and start with new session')
    parser.add_argument('--outputLayout',
                        dest='outputLayout',
                        default=None,
                        help='Output graph with node layout positions to this file')
    parser.add_argument('--outputImage',
                        dest='outputImage',
                        default=None,
                        help='Output image of graph with layout and styles to this file')
    
    # Switch between listening to command-line arguments or hard-coded 
    # arguments depending on whether running in IDE or from cmd.
    if any([name.startswith('SPYDER') for name in os.environ]):
        # myArgs = 'test/test_output.gml --layout referencePCA --refstates test/test_refstates.csv --coloursubgraphs --clearsession'
        # myArgs = 'test/ES_2iL+B-A_output.gml --layout PCA --coloursubgraphs --clearsession'
        """        
        myArgs = ('test/test_graph_simple.gml --layout grid --coloursubgraphs '
                  '--clearsession --outputLayout test/test_graph_simple_layout.json '
                  '--outputImage test/test_graph_simple_layout.png')
        """
        
        myArgs = ('test/test_output.gml --layout referencePCA '
                  '--refstates test/test_refstates.csv --coloursubgraphs '
                  '--clearsession --outputLayout test/test_ouput_layout.json '
                  '--outputImage test/test_output_layout.png')
        
        # myArgs = 'test/test_graph_simple.gml test/test_graph_simple_params.txt'
        args = parser.parse_args(myArgs.split())
    else:
        args = parser.parse_args()
    
    # If a parameter file is specified, read in parameters from there.
    if args.parameters is not None:
        args = read_parameters(args, args.parameters)    
    
    # Delete old session and start a new session if desired.
    if args.clearsession is True:
        cy.session.delete()
    
    # Load state transition graph in Cytoscape
    #
    # This didn't work if there are spaces in the absolute filepath...!
    # Changing the source code of py2cytoscape to fix this
    cyGraph = cy.network.create_from(args.graph, collection='Boolean')
    cyGraphTable = cyGraph.get_node_table()
    cyGraphTable['datasource'] = 'simulation'
    cyGraph.update_node_table(cyGraphTable, data_key_col='name')
        
    # Load reference profiles into Cytoscape
    if args.refstates is not None:
        refTable = pd.read_csv(args.refstates, index_col=0)
        refTable['datasource'] = 'reference'
        refTable['membership'] = 'reference'
        cyGraph.add_nodes(refTable.index.tolist())
        cyGraph.update_node_table(refTable)
    
    # Apply styles
    booleanStyle = cy.style.create('Boolean')
    booleanStyle.update_defaults({
        'NODE_FILL_COLOR': '#9E9E9E',
        'NODE_WIDTH': '20',
        'NODE_HEIGHT': '20',
        'NODE_SHAPE': 'ELLIPSE',
        'NODE_BORDER_PAINT': '#424242',
        'NODE_BORDER_WIDTH': '1.0',
        'NODE_LABEL_COLOR': '#424242',
        'NODE_LABEL_FONT_SIZE': '8',
        'NODE_SIZE': '20',
    })
    booleanStyle.create_discrete_mapping(column='datasource',
                                         vp='NODE_HEIGHT',
                                         mappings={
                                             'reference': '40',
                                             'simulation': '20'
                                         })
    booleanStyle.create_discrete_mapping(column='datasource',
                                         vp='NODE_WIDTH',
                                         mappings={
                                             'reference': '40',
                                             'simulation': '20'
                                         })
    booleanStyle.create_discrete_mapping(column='datasource',
                                         vp='NODE_SIZE',
                                         mappings={
                                             'reference': '40',
                                             'simulation': '20'
                                         })
    booleanStyle.create_discrete_mapping(column='datasource',
                                         vp='NODE_LABEL_COLOR',
                                         mappings={
                                             'reference': '#000000',
                                             'simulation': '#424242'
                                         })
    booleanStyle.create_discrete_mapping(column='datasource',
                                         vp='NODE_LABEL_FONT_SIZE',
                                         mappings={
                                             'reference': '12',
                                             'simulation': '8'
                                         })
    booleanStyle.create_passthrough_mapping(column='name',
                                            vp='NODE_LABEL')
                                            
    # Apply colour based on SCC / SS
    if args.coloursubgraphs is True:
        cyGraphTable = cyGraph.get_node_table()
        hasAttractors = False
        try:
            nSCC = cyGraphTable[cyGraphTable['membership'].str.contains('SCC')].shape[0]
            nSS = cyGraphTable[cyGraphTable['membership'].str.contains('SS')].shape[0]
            hasAttractors = True
        except ValueError:
            hasAttractors = False
        
        if hasAttractors is True:
            sccColors = sns.color_palette('husl', nSCC).as_hex()
            ssColors = sns.color_palette('husl', nSS).as_hex()
            colorMappings = {}
            i = 0
            j = 0
            for label in cyGraphTable['membership'].unique().tolist():
                if label not in colorMappings:
                    if 'SCC' in label:
                        colorMappings[label] = sccColors[i]
                        i += 1
                    elif 'SS' in label:
                        colorMappings[label] = ssColors[j]
                        j += 1
            colorMappings['reference'] = '#9575CD'
            booleanStyle.create_discrete_mapping(column='membership',
                                                 vp='NODE_FILL_COLOR',
                                                 mappings=colorMappings)
    
    # Apply colour by reference vs. simulation
    booleanStyle.create_discrete_mapping(column='datasource',
                                     vp='NODE_FILL_COLOR',
                                     mappings={
                                         'reference': '#9575CD',
                                         'simulation': '#4DD0E1'
                                     })
    
    # Push styles to Cytoscape
    cy.style.apply(style=booleanStyle, network=cyGraph)
    
    
    
    # Apply layout
    if args.layout in originalLayouts:
        cy.layout.apply(args.layout, cyGraph)
        cy.layout.fit(cyGraph)
    elif args.layout == 'referencePCA':
        locations = reference_pca_locations(cyGraph, xscale=200.0, yscale=-200.0)
        cy.layout.apply_from_presets(cyGraph, positions=locations)
        cy.layout.fit(cyGraph)
    elif (args.layout == 'PCA') or (args.layout == 'pca'):
        locations = pca_locations(cyGraph, xscale=5000.0, yscale=5000.0)
        cy.layout.apply_from_presets(cyGraph, positions=locations)
        cy.layout.fit(cyGraph)
    else:
        raise ValueError('Unknown layout method: %s' %(args.layout))
        
    
    # Output layout image to file
    if args.outputImage is not None:
        imageExt = args.outputImage.rsplit('.')[-1]        
        acceptableExt = set(['cx', 'pdf', 'png', 'svg'])
        if imageExt in acceptableExt:
            # Wait for a bit, otherwise the image that gets returned is not complete
            time.sleep(5)
            networkId = str(cyGraph.get_id())
            command = 'networks/' + networkId + '/views/first.' + imageExt
            res = manual_rest(command, httpType='get')
            with open(args.outputImage, 'wb') as f:
                f.write(res.content)
        else:
            raise ValueError('Unsupported image format: %s' %(imageExt))
            
    # Output graph layout to file as JSON
    if args.outputLayout is not None:
        layoutExt = args.outputLayout.rsplit('.')[-1]
        acceptableLayoutExt = set(['json', 'cyjs'])
        if layoutExt in acceptableLayoutExt:
            json = cyGraph.get_first_view()
            jsonString = str(json).replace("u'", "'").replace("'", '"').replace('False', "false").replace('True', "true")
            with open(args.outputLayout, 'wb') as f:
                f.write(jsonString)
        else:
            raise ValueError('Unsupported graph layout format: %s' %(layoutExt))
    