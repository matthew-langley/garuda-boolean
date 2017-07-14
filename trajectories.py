"""
Trajectory finding in Boolean state transition graph, as well as clustering and
thematic analysis of those clusters (in terms of average gene expression)
"""


import argparse


def read_parameters(args, paramFile):
    """
    Reads command line arguments specifiying parameters for Boolean 
    trajectory analysis from a specified file
    
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
                    if param == 'source':
                        args.source = int(value)
                    elif param == 'target':
                        args.target = int(value)
                    elif param == 'clusters':
                        args.clusters = int(value)
                    elif param == 'linkage':
                        args.linkage = str(linkage)
                    elif param == 'distance':
                        args.outputImage = str(distance)
    
    return args


if __name__ == '__main__':
    """
    Code runs only if module is executed as mainline. Supports command line 
    arguments.
    """
    
    # Configure command line arguments
    description = ('Find trajectories between states in a Boolean network '
                   'simulation, and optionally cluster those trajectories '
                   'and discover themes in average gene expression')
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('graph',
                        help='State transition graph in GML format')
    parser.add_argument('parameters',
                        nargs='?',
                        help='(optional) File containing simulation parameters described below')
    parser.add_argument('--source',
                        type=int,
                        dest='source',
                        help='ID of starting state or SCC')
    parser.add_argument('--target',
                        type=int,
                        dest='target',
                        help='ID of final state or SCC')
    parser.add_argument('--clusters',
                        dest='clusters',
                        type=int,
                        default=0,
                        help=('Analyze the top (N) clusters of trajectories. '
                              'By default, do not perform clustering.'))
    parser.add_argument('--linkage',
                        dest='linkage',
                        type=str,
                        default='average',
                        help='Linkage metric used for clustering')
    parser.add_argument('--distance',
                        dest='distance',
                        type=str,
                        default='frechet',
                        help='Distance metric used for clustering')
    parser.add_argument('--refprofiles',
                        dest='refprofiles',
                        help='')
    parser.add_argument('--outputTrajectories',
                        help='Output file destination for trajectory list')
    parser.add_argument('--outputTrajectoryGraph',
                        help='Output file destination for trajectory graph')
    parser.add_argument('--outputClusters',
                        help='Output file destination for cluster mapping')
    parser.add_argument('--outputClusterTimecourses',
                        help='Output file destination for cluster timecourse average data')
    parser.add_argument('--outputClusterGraphs',
                        help='Output file destination for cluster timecourse line plots')
                        
    # If a parameter file is specified, read in parameters from there.
    if args.parameters is not None:
        args = read_parameters(args, args.parameters)