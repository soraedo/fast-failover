import sys
from typing import List, Any, Union

import networkx as nx
import numpy as np
import random
import time
import math
from objective_function_experiments import *

import outerplanar_graph as outerplanar
import routing
DEBUG = True

# Data structure containing the algorithms under scrutiny. Each entry contains a name and a pair of algorithms.
#
# The first algorithm is used for any precomputation to produce data structures later needed for routing on the graph passed along in args. If the precomputation
# fails, the algorithm must return -1. Examples for precomputation algorithms can be found in arborescences.py
#
# The second algorithm decides how to forward a packet from source s to destination d, despite the link failures fails using data structures from precomputation
# Examples for precomputation algorithms can be found in routing.py
#
# In this example we compare the RightHandRule on outerplanar networks and EDP on the original topology.
# You can add more algorithms to this data structure to compare the performance of additional algorithms.
algos = {'Outerplanar': [outerplanar.createOuterplanar, outerplanar.RightHandRule],
         'EdgeDisjointPaths': [outerplanar.original, outerplanar.EdgeDisjointPaths]}

# run one experiment with graph g
# out denotes file handle to write results to
# seed is used for pseudorandom number generation in this run
# returns a score for the performance:
#       if precomputation fails : 0
#       if d becomes disconnected from the rest of the graph : -1
#       otherwise score is the mean of all success scores (closest to 1 is best)
def one_experiment(g, seed, out, algo, f_list, first_string):
    [precomputation_algo, routing_algo] = algo[:2]
    if DEBUG: print('experiment for ', algo[0])

    # precomputation
    random.seed(seed)
    t = time.time()
    precomputation = precomputation_algo(g)
    pt = time.time() - t
    if precomputation == -1:  # error...
        out.write(', %f, %f, %f, %f, %f, %f\n' %
                  (float('inf'), float('inf'), float('inf'), 0, 0, pt))
        print('precomputation didnt work', algo[:2], g.graph['seed'], seed)
        sys.exit()
    stat = routing.Statistic(routing_algo, str(routing_algo), g=None)
    success_ratios = []

    for f_num in f_list:
        out.write(first_string)
        if attack == "CLUSTER":
            g.graph['fails'] = targeted_attacks_against_clusters(g,f_num)
        cc_size = len(set(connected_component_nodes_with_d_after_failures(g,g.graph['fails'][:f_num],g.graph['root'])))
        stat.reset(g.nodes())
        random.seed(seed)
        t_start = time.time()

        # simulate the routing of the packet
        if -1 == routing.SimulateGraph(g, True, [stat], f_num, samplesize, dest = g.graph['root'], precomputation=precomputation, targeted=False):
            print('simulate graph returns -1', algo[1], len(g.edges()), f_num )
            sys.exit()
        t_end = time.time()
        num_experiments = stat.succ + stat.fails
        if num_experiments < min(len(g.nodes())-1,samplesize):
            print(num_experiments, samplesize, cc_size, algo[1], "sample size doesn't match number of experiments") #TODO investigate
            sys.exit()
        rt = (t_end - t_start)/max(num_experiments,1)

        # success can be from 0 to 1 per rep with 1 representing 100% success.
        success_ratio = stat.succ/max(num_experiments,1)
        if stat.succ > cc_size -1:
            print('more success %i, than cc size %i, investigate seed'  % (algo[1], (stat.succ, cc_size, seed)))
            sys.exit()
        # write results
        if stat.succ > 0:
            if DEBUG: print('success', stat.succ, algo[0])
            # stretch, load, hops, success, cc_size/n, routing time, precomputation time
            out.write(', %i, %i, %i, %i, %f, %i, %f, %f\n' %
                      (f_num, np.max(stat.stretch), stat.load, np.max(stat.hops),
                       success_ratio, cc_size, rt, pt))
            out.flush()
        else:
            if DEBUG:
                print('no success_ratio', algo[0], 'seed', g.graph['seed'],'expected ratio <=', cc_size/n )
            out.write(', %i, %f, %f, %f, %f, %i, %f, %f\n' %
                      (f_num, float('inf'), float('inf'), float('inf'), float(0), cc_size, float(rt), float(pt)))
            
        # cc_size == 1 means d is unreachable and isolated from the rest of the graph.
        if cc_size > 1:
            success_ratios.append(success_ratio)
        else:
            success_ratios.append(-1)  
    return success_ratios

'''
# run experiments with AS graphs
# out denotes file handle to write results to
# seed is used for pseudorandom number generation in this run
# rep denotes the number of repetitions in the shuffle for loop
def run_AS(out=None, seed=0, rep=5):
    for i in range(4, 5):
        generate_trimmed_AS(i)
    files = glob.glob('./infocom_graphs/AS*.csv')
    original_params = [n, rep, k, samplesize, f_num, seed, name]
    for x in files:
        random.seed(seed)
        kk = int(x[-5:-4])
        g = nx.read_edgelist(x).to_directed()
        g.graph['k'] = kk
        nn = len(g.nodes())
        mm = len(g.edges())
        ss = min(int(nn / 2), samplesize)
        fn = min(int(mm / 4), f_num)
        fails = random.sample(list(g.edges()), fn)
        g.graph['fails'] = fails
        set_infocom_parameters([nn, rep, kk, ss, fn, seed, name + "AS-"])
        shuffle_and_run(g, out, seed, rep, x)
        set_infocom_parameters(original_params)
'''
        
# run experiments with zoo graphs
# out denotes file handle to write results to
# seed is used for pseudorandom number generation in this run
# rep denotes the number of repetitions in the shuffle for loop
def run_zoo(out=None, seed=0, rep=5, min_connectivity=1):
    original_params = [n, rep, k, samplesize, f_num, seed, name]
    if DEBUG:
        print('n_before, n_after, m_after, connectivity, degree')
    # 263 graphs
    for i in range(262):
        random.seed(seed)
        g = read_zoo(i, min_connectivity) # this is a directed graph
        if g is None:
            continue
        g = g.to_undirected()
        nn = len(g.nodes())
        # only run it on graphs with more than 20 and less than 100 nodes
        if nn < 20 or nn > 50 or (short and len(g.nodes()) > 2500):
            continue
        kk = nx.edge_connectivity(g)
        mm = len(g.edges())
        ss = min(int(nn / 1.1), samplesize)
        fn = min(int(mm / 4), f_num)
        f_list = [int(math.ceil(0.01*f_num*mm))] # denotes the percentage of edges in the graph that fail. f_num is given in percentage
        print(i, f_list, 'failure number list', time.asctime(time.localtime(time.time())))
        set_infocom_parameters([nn, rep, kk, ss, -1, seed, name + "zoo-"])
        results = shuffle_and_run(g, out, seed, rep, str(i), f_list)
        set_infocom_parameters(original_params)
        for fn in f_list:
            print(fn, "failures")
            print('#connected destination component of size 1:', results[fn]['small_cc']/len(algos), 'repetitions', rep)
            if results[fn]['small_cc']/len(algos) == rep:
                continue
            print('(mean(cc_size)-1)/(n-1):', (np.mean(results[fn]['cc_size'])-1)/(n-1))
            print('mean success ratio, min, std')
            for (algoname, algo) in algos.items():
                scores = results[fn]['scores'][algoname]
                print('%.4E, %.4E, %.4E : %s' % (np.mean(scores), np.max(scores), np.std(scores), algoname))
                sys.stdout.flush()
        if short and i > 5:
            break

# shuffle root nodes and run algorithm
def shuffle_and_run(g, out, seed, rep, x, f_list):
    random.seed(seed)
    results = {f : {'scores':{algoname:[] for algoname in algos.keys()}, 'small_cc':0, 'cc_size':[]}  for f in f_list}
    nodes = list(g.nodes())
    for count in range(rep):
        g.graph['root'] = nodes[count%len(nodes)]
        for (algoname, algo) in algos.items():
            # graph, size, connectivity, algorithm, index,
            first_string = '%s, %i, %i, %s, %i' % (x, len(g.nodes()), g.graph['k'], algoname, count)
            scores = one_experiment(g, seed + rep, out, algo, f_list, first_string)
            print(scores)
            for i in range(len(f_list)):
                f = f_list[i]
                if scores[i] > -1:
                    algos[algoname] += [scores[i]]
                    c = [len(set(connected_component_nodes_with_d_after_failures(g,g.graph['fails'][:f],g.graph['root'])))]
                    results[f]['scores'][algoname] += [scores[i]]
                    results[f]['cc_size'] += c
                else:
                    results[f]['small_cc'] += 1
    return results

# return nodes in connected component with node d (after failures have been removed)
def connected_component_nodes_with_d_after_failures(g, failures, d):
    G =g.to_undirected()
    G.remove_edges_from(failures)
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
    for i in range(len(Gcc)):
        if d in Gcc[i]:
            return list(Gcc[i])

'''
# run experiments with d-regular graphs
# out denotes file handle to write results to
# seed is used for pseudorandom number generation in this run
# rep denotes the number of repetitions in the secondary for loop
def run_regular(out=None, seed=0, rep=5):
    ss = min(int(n / 2), samplesize)
    fn = min(int(n * k / 4), f_num)
    set_infocom_parameters([n, rep, k, ss, fn, seed, name + "regular-"])
    write_graphs()
    f_list = [f_num]
    for i in range(rep):
        random.seed(seed + i)
        g = read_graph(i)
        random.seed(seed + i)
        for (algoname, algo) in algos.items():
            # graph, size, connectivity, algorithm, index,
            first_string = '%s, %i, %i, %s, %i' % ("regular", n, k, algoname, i)
            algos[algoname] += [one_experiment(g, seed + i, out, algo, f_list, first_string)]
'''

# Custom targeted link failure model
# Return the selected links incident to nodes with a non-zero value of the clustering coefficient
def targeted_attacks_against_clusters(g, f_num):
    candidate_links_to_fail = list()
    links_to_fail = list()
    clustering_coefficients = nx.clustering(g)
    for ( v, cc ) in clustering_coefficients.items():
        if cc == 0.0:
            continue
        neighbors = nx.neighbors( g, v )
        for w in neighbors:
            if not ( v, w ) in candidate_links_to_fail and not ( w, v ) in candidate_links_to_fail:
                candidate_links_to_fail.append( ( v, w ) )
    # Select up to f_num bi-directional links that should be disabled
    if len( candidate_links_to_fail ) > f_num:
        links_to_fail = random.sample( candidate_links_to_fail, f_num )
    else:
        links_to_fail.extend( candidate_links_to_fail )
    # Append the opposite arcs to the list (we assume failures affect links in both directions)
    for ( v, w ) in links_to_fail:
        if not ( w, v ) in candidate_links_to_fail:
            links_to_fail.append( ( w, v ) )
    return links_to_fail

# start file to capture results
def start_file(filename):
    out = open(filename + ".txt", 'w')
    out.write(
        "graph, size, connectivity, algorithm, repetition, failures, " +
        "stretch, load, hops, success, cc_size," +
        "routing_time, precomputation_time\n")
    out.write(
        "#" + str(time.asctime(time.localtime(time.time()))) + "\n")
    return out

# set global parameters in this file and in objective_function_experiments
def set_infocom_parameters(params):
    global seed, n, rep, k, samplesize, name, f_num
    [n, rep, k, samplesize, f_num, seed, name] = params
    set_objective_parameters(params)

# print global parameters in this file and in routing_stats
def print_infocom_parameters():
    print(n, rep, k, samplesize, f_num, seed, name)

# run experiments
# seed is used for pseudorandom number generation in this run
# switch determines which experiments are run
def experiments(switch="all", seed=0, rep=100):
    #if switch in ["regular", "all"]:
    #    out = start_file("results/outputREG_seed45_rep" + str(n) + "_" + str(f_num) + "f")
    #    run_regular(out=out, seed=seed, rep=rep)
    #    out.close()

    if switch in ["zoo", "all"]:
        out = start_file("results/output_45seed_" + str(f_num) + "f" + "_" + "50-100n_" + attack)
        run_zoo(out=out, seed=seed, rep=rep, min_connectivity=1)
        out.close()

    #if switch in ["AS"]:
    #    out = start_file("results/outputAS_seed45_rep" + str(n) + "_" + str(f_num) + "f")
    #    run_AS(out=out, seed=seed, rep=rep)
    #    out.close()

    print(attack)
    for (algoname, algo) in algos.items():
        print('%.5E %s' % (np.mean(algo[2:]), algoname))
    print("\closest to 1.0 is better")



if __name__ == "__main__":
    #default values
    f_num = 10 #number of failed links
    n = 100 # number of nodes
    k = 5 #base connectivity
    samplesize = 20 #number of sources to route a packet to destination
    rep = 100 #number of experiments
    switch = 'zoo' #which experiments to run with same parameters
    seed = 45 #random seed
    name = "output-" #result files start with this name
    attack = 'RANDOM' #how edge failures are chosen. Can be RANDOM or CLUSTER
    short = False  #if true only small zoo graphs < 25 nodes are run
    start = time.time()
    print(time.asctime(time.localtime(start)))
    if len(sys.argv) > 1:
        switch = sys.argv[1]
    if len(sys.argv) > 2:
        seed = int(sys.argv[2])
    if len(sys.argv) > 3:
        rep = int(sys.argv[3])
    if len(sys.argv) > 4:
        n = int(sys.argv[4])
    if len(sys.argv) > 5:
        samplesize = int(sys.argv[5])
    if len(sys.argv) > 6:
        f_num = int(sys.argv[6]) #as a percentage of failed edges
    if len(sys.argv) > 7:
        attack = sys.argv[7] #RANDOM, EDGECUT, CLUSTER
    if len(sys.argv) > 8:
        short = sys.argv[8] == 'True' #True or False
    print(time.asctime(time.localtime(start)), 'attack', attack)
    random.seed(seed)
    set_parameters([n, rep, k, samplesize, f_num, seed, "benchmark-"])
    experiments(seed=seed, rep=rep)
    end = time.time()
    print("time elapsed", end - start)
    print("start time", time.asctime(time.localtime(start)))
    print("end time", time.asctime(time.localtime(end)))

    #example call: python nancy_benchmark_template.py zoo 45 100 100 20 10 RANDOM False