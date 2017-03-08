# 
# Copyright 2017 Matthew Langley. All rights reserved.
#


"""
Simulation of Boolean networks
"""

import argparse
import boolean2 as bn
import itertools as it
import multiprocessing as mp
import multiprocessing.forking as forking
import os
import re
import sys


def fix_multiprocessing_for_exe():
    """
    Code snippet to make multiprocessing work properly on Windows when frozen  
    into a .exe using pyinstaller.
    [1] https://github.com/pyinstaller/pyinstaller/wiki/Recipe-Multiprocessing
    """
    mp.freeze_support()
    if sys.platform.startswith('win'):
        # First define a modified version of Popen.
        class _Popen(forking.Popen):
            def __init__(self, *args, **kw):
                if hasattr(sys, 'frozen'):
                    # We have to set original _MEIPASS2 value from sys._MEIPASS
                    # to get --onefile mode working.
                    os.putenv('_MEIPASS2', sys._MEIPASS)
                try:
                    super(_Popen, self).__init__(*args, **kw)
                finally:
                    if hasattr(sys, 'frozen'):
                        # On some platforms (e.g. AIX) 'os.unsetenv()' is not
                        # available. In those cases we cannot delete the variable
                        # but only set it to the empty string. The bootloader
                        # can handle this case.
                        if hasattr(os, 'unsetenv'):
                            os.unsetenv('_MEIPASS2')
                        else:
                            os.putenv('_MEIPASS2', '')
                            
        # Second override 'Popen' class with our modified version.
        forking.Popen = _Popen
  

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
                if param == 'runs':
                    args.runs = int(value)
                elif param == 'steps':
                    args.steps = int(value)
                elif param == 'mode':
                    if value == 'sync':
                        args.mode = 'sync'
                    else:
                        args.mode = 'async'
                elif param == 'writeStateDict':
                    if value in ['True', 'true', 'T', 't', '1', 'Yes', 'Y', 'y']:
                        args.writeStateDict = True
                    else:
                        args.writeStateDict = False
                elif param == 'writeEdgeDict':
                    if value in ['True', 'true', 'T', 't', '1', 'Yes', 'Y', 'y']:
                        args.writeEdgeDict = True
                    else:
                        args.writeEdgeDict = False
    
    return args
  

def simulate(modelText, mode='async', steps=100, runs=100, nprocesses=None):
    """
    Simulates a Boolean model over multiple instances. Each instance (or run) 
    consists of a specified number of discrete time steps in which some 
    (async) or all (sync) of the genes in the model are updated based on their 
    associated Boolean update function. Since these Boolean simulations can be 
    time-intensive, the method allows for run instances to be performed on 
    parallel processes.
    
    Args:
        modelText (str)
            The Boolean model to simulate, either as a filepath or the
            complete definition as text.
        mode (str)
            The update strategy used for simulation. Either 'async' or 'sync'. 
            Defaults to 'async'.
        steps (int)
            The number of update time steps performed in each run. Defaults to 
            100.
        runs (int)
            The number of run instances to simulate. Defaults to 100.
        nprocesses (int)
            The number of processes used for the parallel worker pool. Set to 
            1 to run without multiprocessing. If unspecifed, the number of 
            worker processes used will be equal to the number of CPU cores on
            the target machine.
            
    Returns:
        (list)
            A list of all of trajectories produced in the batch of simulation
            instances. Each trajectory in the list is itself a list of the 
            Boolean states traversed by that run instance of the simulation, 
            up to the specified number of steps or the end of any limit cycle 
            that was encountered.
    """
    
    # If nprocesses not specified, set to the number of CPU cores on target
    # machine.
    if nprocesses is None:
        nprocesses = mp.cpu_count()
        
    # Run the simulations by calling into 'simulate_one_run' and retrieve the 
    # results. If mulitple processes are to be used, set up a multiprocessing 
    # Pool of workers to perform the simulations.
    if nprocesses > 1:
        pool = mp.Pool(processes=nprocesses)
        handles = [pool.apply_async(simulate_one_run, 
                   args=(modelText, mode, steps)) for i in range(runs)]
        pool.close()
        results = [r.get() for r in handles]
    else:
        results = [simulate_one_run(modelText, mode, steps) for i in range(runs)]
        
    # Return the list of boolean2.state.State trajectories traversed by each
    # run instance.
    return results


def simulate_one_run(modelText, mode, steps):
    """
    Simulates a single run instance of a Boolean model, i.e. only one 
    trajectory from one initial condition. The simulation consists of a 
    specified number of discrete time steps in which some (async) or all 
    (sync) of the genes in the model are updated based on their associated 
    Boolean update function.
    
    Args:
        modelText (str)
            The Boolean model to simulate, either as a filepath or the
            complete definition as text.
        mode (str)
            The update strategy used for simulation. Either 'async' or 'sync'. 
            Defaults to 'async'.
        steps (int)
            The number of update time steps performed in each run. Defaults to 
            100.
            
    Returns:
        (list)
            The trajectory of Boolean states traversed by the simulation, up 
            to the specified number of steps or the end of any limit cycle 
            that was encountered.
    """
    
    # Setup and run model.
    model = bn.Model(modelText, mode=mode)
    model.initialize()
    model.iterate(steps)
    
    # Determine length of path to return based on limit cycle detection.
    # Need to return (index + size + 1) instead of (index + size) so that the
    # returning edge from the last state in cycle to the first one is included 
    # in the edge list.
    index, size = model.detect_cycles()
    if index + size == 0 or index + size >= steps:
        pathLen = steps
    else:
        pathLen = index + size + 1
    
    return model.states[:pathLen]


def assign_ids(stateList):
    """
    Given a list of Boolean State objects from a simulation, assigns a 
    globally consistent fingerprint to that State object as a new attribute,
    'id'. This is necessary because we can't trust how BooleanNet handles 
    state.fp() when using multiprocessing, since it relies on global variables 
    for the State class that are not shared across processes.
    
    Args:
        stateList (list)
            A list of all Boolean states to which IDs should be assigned. 
            States with the same binary profile will be assigned identical IDs.
    """
    
    # We can't trust how BooleanNet handles state.fp() when using 
    # multiprocessing, since it relies on global variables for the State class 
    # that are not shared across processes. Therefore, we need to implement 
    # our own version of fingerprinting that is calculated after all simulation
    # runs have completed.
    # stateDict = {s.fp(): s for s in stateList}
    
    myMapper = {}
    myCounter = 0
    for state in stateList:
        # Use hashed value of string representation of State for indexing in 
        # the dictionary. Need to do this instead of the State object itself 
        # because even if two State objects represent the same binary profile,
        # they are still different *objects*.
        # i.e. Even if stateA.bin() == stateB.bin();
        #      myMapper[stateA] != myMapper[stateB]
        key = hash(str(state))
        if key not in myMapper:
            myMapper[key] = myCounter
            myCounter += 1
        state.id = myMapper[key]


def make_state_dictionary(stateList):
    """
    Given a list of Boolean states, returns a dictionary in which the keys are 
    an ID (fingerprint) unique to that state, and the associated values are 
    the Boolean state. Note that IDs for a state are consistent across runs 
    completed at the same time; i.e. '101' in run #1 and '101' in run #2 will 
    have the same ID. However, consistent IDs are not guaranteed across 
    experiments; i.e. if I run a model today, and that same model tomorrow (new
    instance of Python), the IDs will be different.
    
    Args:
        stateList (list)
            A list of Boolean states to include in the dictionary.
            
    Returns:
        (dict)
            A dictionary in which the keys are an ID (fingerprint) unique to 
            that state, and the associated values are the Boolean state.
    """
    
    # We can't trust how BooleanNet handles state.fp() when using 
    # multiprocessing, since it relies on global variables for the State class 
    # that are not shared across processes. Therefore, we must use our own
    # State.id attribute which is created after all simulations are completed 
    # and control returns to the main (parent) process.
    
    stateDict = {s.id: s for s in stateList}    
    return stateDict


def make_edge_dictionary(trajectoryList):
    """
    Given a list of Boolean state trajectories, returns a dictionary in which 
    the keys are a tuple of (sourceID, targetID) for each edge in the 
    trajectories, and the values are the number of times that edge was 
    encountered in the trajectory list.
    
    Args:
        trajectoryList (list)
            A list of Boolean state trajectories.
            
    Returns:
        (dict)
            A dictionary in which the keys are a tuple of (sourceID, targetID) 
            for each edge in the trajectories, and the values are the number 
            of times that edge was encountered in the trajectory list.
    """
    
    edgeDict = {}
    for trajectory in trajectoryList:
        for source, target in it.izip(trajectory, trajectory[1:]):
            s = source.id
            t = target.id
            if (s,t) in edgeDict.keys():
                edgeDict[(s,t)] += 1
            else:
                edgeDict[(s,t)] = 1
    return edgeDict
    

def write_state_dictionary(stateDict, outputFile):
    """
    Writes a simulation's state dictionary to file.
    
    Args:
        stateDict (dict)
            The dictionary to be written to file.
        outputFile (str)
            The destination filepath.
    """
    
    # Start by writing the header row.
    genes = sorted(stateDict.values()[0].keys())
    genes.remove('id')
    outputText = 'ID,' + ','.join(genes) + '\n'
    
    # For each state, write the state ID and value for each gene (as 1 or 0) 
    # as a new row.
    for k, v in sorted(stateDict.items()):
        outputText += str(k) + ','
        outputText += ','.join([str(int(v[i])) for i in genes])
        outputText += '\n'
    
    # Write output to file.
    with open(outputFile, 'w') as f:
        f.write(outputText)
        

def write_edge_dictionary(edgeDict, outputFile):
    """
    Writes a simulation's edge dictionary to file.
    
    Args:
        edgeDict (dict)
            The dictionary to be written to file.
        outputFile (str)
            The destination filepath.
    """
    
    # Start by writing the header row.
    outputText = 'Source,Target,Frequency\n'
    
    # For each edge, write the source state ID, target state ID, and 
    # frequency as a new row.
    for (s, t), f in sorted(edgeDict.items()):
        outputText += ','.join([str(s), str(t), str(f)])
        outputText += '\n'
        
    # Write output to file.
    with open(outputFile, 'w') as f:
        f.write(outputText)
        
        
def write_gml(stateDict, edgeDict, outputFile, graphName):
    """
    Writes the results of the simulation (i.e. the directed graph representing 
    the simulation's state transition graph) to file in the Graph Modeling 
    Language, GML. Details of the GML format available at:
    [1] http://www.fim.uni-passau.de/fileadmin/files/lehrstuhl/brandenburg/projekte/gml/gml-technical-report.pdf
    [2] https://en.wikipedia.org/wiki/Graph_Modelling_Language
    
    Note: Due to limitations of the GML spec, any non alpha-numeric characters 
    will be ommitted from the gene names
    
    Args:
        stateDict (dict)
            The state dictionary for the simulation
        edgeDict (dict)
            The edge dictionary for the simulation
        outputFile (str)
            The destination filepath
    """
    
    genes = sorted(stateDict.values()[0].keys())
    genes.remove('id')
    
    lines = []
    lines.append('graph [')
    lines.append('        directed 1')
    lines.append('        name "%s"' %(graphName))
    for k,v in stateDict.items():
        lines.append('        node [')
        lines.append('               id %i' %(k))
        lines.append('               label "%i"' %(k))
        for gene in genes:
            gene0 = re.sub(r'[^A-Za-z0-9#]+', "", gene)
            lines.append('               %s %i' %(gene0, int(v[gene])))
        lines.append('        ]')
    for (s,t), f in edgeDict.items():
        lines.append('        edge [')
        lines.append('               source %i' %(s))
        lines.append('               target %i' %(t))
        lines.append('               weight %i' %(f))
        lines.append('        ]')
    lines.append(']')
    
    with open(outputFile, 'w') as f:
        f.write('\n'.join(lines))


if __name__ == '__main__':
    """
    Code runs only if module is executed as mainline. Supports command line 
    arguments.
    """
    
    # Configure multiprocessing to work properly in frozen .exe version
    fix_multiprocessing_for_exe()
    
    # Configure command line arguments
    parser = argparse.ArgumentParser(description='Simulate a Boolean network.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('model',
                        help='File describing the Boolean network to be simulated')
    parser.add_argument('parameters',
                        nargs='?',
                        help='(optional) File containing simulation parameters described below')
    parser.add_argument('--output',
                        help='Output file destination')
    parser.add_argument('--runs',
                        dest='runs',
                        type=int,
                        default=1000,
                        help='Number of simulation runs to perform')
    parser.add_argument('--steps',
                        dest='steps',
                        type=int,
                        default=100,
                        help='Number of Boolean update steps per simulation')
    parser.add_argument('--processes',
                        dest='nprocesses',
                        type=int,
                        default=1,
                        help='Number of processes for parallel computation')
    modeGroup = parser.add_mutually_exclusive_group()
    modeGroup.add_argument('--sync',
                           dest='mode',
                           action='store_const',
                           const='sync',
                           help='All genes are updated at each step')
    modeGroup.add_argument('--async',
                           dest='mode',
                           action='store_const',
                           const='async',
                           help='Random subset of genes updated at each step')
    parser.add_argument('--writeStateDict',
                         dest='writeStateDict',
                         action='store_true',
                         default=False,
                         help='Write resulting state dictionary to file')
    parser.add_argument('--writeEdgeDict',
                         dest='writeEdgeDict',
                         action='store_true',
                         default=False,
                         help='Write resulting edge dictionary to file')
                           
    # Switch between listening to command-line arguments or hard-coded 
    # arguments depending on whether running in IDE or from cmd.
    if any([name.startswith('SPYDER') for name in os.environ]):
        myArgs = 'test/test_model.txt --output test/test_output.gml --runs 50 --steps 10 --processes 2 --async --writeStateDict --writeEdgeDict'
        args = parser.parse_args(myArgs.split())
    else:
        args = parser.parse_args()
    
    # If a parameter file is specified, read in parameters from there.
    if args.parameters is not None:
        args = read_parameters(args, args.parameters)
    
    # Default to asyncrhonous simulation mode.
    if args.mode is None:
        args.mode = 'async'
    
    # Make sure the specified model file actually exists.
    if not os.path.isfile(args.model):
        parser.error('The specified model file does not exist: %s' %(args.model))
    
    # Simulate the Boolean network
    trajectories = simulate(modelText=args.model,
                            mode=args.mode,
                            steps=args.steps,
                            runs=args.runs,
                            nprocesses=args.nprocesses)
                            
    # Aggregate the results of all the runs into dictionaries, ensuring that 
    # expression states that are common between runs are given the same ID, 
    # and that all the edges that are common between runs are counted.
    assign_ids(list(it.chain(*trajectories)))
    stateDict = make_state_dictionary(list(it.chain(*trajectories)))
    edgeDict = make_edge_dictionary(trajectories)
    
    # Write the results (state transition graph) to file.
    if args.output is None:
        args.output = '../Output/' + args.model.rsplit('.', 1)[0] + '_output.gml'
    graphName = args.model.rsplit('/', 1)[1].rsplit('.', 1)[0]
    write_gml(stateDict, edgeDict, args.output, graphName)
    
    # If requested, write state and/or edge dictionaries to file.
    if args.writeStateDict is True:
        stateDictFile = args.output.rsplit('.', 1)[0] + '_nodes.csv'
        write_state_dictionary(stateDict, stateDictFile)
    if args.writeEdgeDict is True:
        edgeDictFile = args.output.rsplit('.', 1)[0] + '_edges.csv'
        write_edge_dictionary(edgeDict, edgeDictFile)
    
    