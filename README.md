This repository is a fork from [https://github.com/inet-tub/fast-failover](https://github.com/inet-tub/fast-failover) and contains the source code for my bachelor thesis: "Improving the Resilience of Fast Failover Routing: Planar Graphs". 
I am indebted to [Klaus-Tycho Foerster](https://ktfoerster.github.io/), [Andrzej Kamisinski](https://home.agh.edu.pl/~andrzejk/), [Yvonne-Anne Pignolet](http://yvonneanne.pignolet.ch/), [Stefan Schmid](https://www.inet.tu-berlin.de/menue/people/profs0/stefan/), [Gilles Tredan](https://homepages.laas.fr/gtredan/) for publishing their simulation framework source code publicly. 

## Overview

* benchmark_graphs: directory to be filled with network topologies used in the experiments (prefilled with Topology Zoo dataset)
* results: directory to which csv and other output files are written

* arborescence.py: arborescence decomposition and helper algorithms
* routing_stats.py: routing algorithms, simulation and statistic framework
* objective_function_experiments.py: objective functions, independence and SRLG experiments
* srds2019_experiments.py: experiments for SRDS 2019 / TSDC 2022 paper
* dsn2019_experiments.py: experiments for DSN 2019 paper
* infocom2019_experiments.py: experiments for INFOCOM 2019 paper
* infocom2021_experiments.py: experiments for INFOCOM 2021 paper
* benchmark_template.py: template to compare algorithms
* nancy_benchmark_template.py: experiments for my bachelor thesis: "Improving the Resilience of Fast Failover Routing: Planar Graphs"
* outerplanar_graph.py: my implementations of an Outerplanar Routing Scheme and an edge disjoint path routing algorithm for comparison

To run the experiments, execute the corresponding python file:
```
python nancy_benchmark_template.py
```
With additional arguments the experiments can be customised (see main function of the python file). E.g., 
```
python nancy_benchmark_template.py zoo 45 100 100 20 10 RANDOM False
```
executes 100 repetitions of the Topology Zoo experiments with seed 45.
