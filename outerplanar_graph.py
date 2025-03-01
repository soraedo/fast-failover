import networkx as nx
import matplotlib.pyplot as plt
DEBUG = False


# Implementation of an Outerplanar Routing Scheme.
# An algorithm was created to create outerplanar subgraphs out of any network graph based on the paper "An improved algorithm for finding maximum outerplanar
# subgraphs" from 2024 by Gruia Călinescu, Hemanshu Kaul, Bahareh Kudarzi, as well as its outerplanar embedding.
# The routing algorithm is the skipping Right-Hand Rule mentioned in https://ktfoerster.github.io//paper/2021-apocs.pdf
# To use it in benchmark_template.py, import outerplanar_graph and extend algos with 'Outerplanar': [outerplanar.createOuterplanar, outerplanar.RightHandRule]
# The first entry of the list is the precomputation algorithm, the second the actual routing algorithm
# Implemented by Nancy Mey Ching Hou

'''
PRECOMPUTATION ALGORITHM
'''
# function which returns the embedding of an outerplanar subgraph of graph g
# outerplanar subgraph creation based on the paper "An improved algorithm for finding maximum outerplanar subgraphs" (2024) by Gruia Călinescu, Hemanshu Kaul, Bahareh Kudarzi
# given undirected graph g
def createOuterplanar(g):
    if DEBUG:
        print(f'Topology zoo edges: {g.edges}')
    # copy g
    G = nx.Graph()
    G.add_nodes_from(g)
    G.add_edges_from(g.edges)

    # convert graph labels to integers, remove self-loops, and make sure graph is undirected
    mapping = dict()
    for i in range(101):
        mapping[str(i)] = i
    G = nx.relabel_nodes(G, mapping, copy = True)
    G.remove_edges_from(nx.selfloop_edges(G))
    G.to_undirected()
    if DEBUG:
        print(f"G {G.edges}")

    # Step 1
    G1 = find_triangular_cactus(G)
    if DEBUG:
        print(f"G1 {G1.edges}")
    # Step 2
    G2 = find_square_cactus(G, G1)
    if DEBUG:
        print(f"G2 {G2.edges}")
    # Step 3
    G3 = find_edges(G, G2)
    if DEBUG:
        print(f"G3 {G3.edges}")
    # Step 4
    G4 = greedy_edges(G, G3)
    if DEBUG:
        print(f"G4 {G4.edges}")
    # G4 should be outerplanar. If not, exit with error -1.
    if not is_outerplanar(G4):
        print(" Output subgraph was not outerplanar")
        return -1
    
    # now we create a planar embedding of G3
    # Step 5
    biconnected = make_biconnected(G4)
    # Step 6
    path = generate_outer_face(biconnected)
    # Step 7
    return convert_to_outerplanar_embedding(G4, path)

# dummy function which returns the original topology. Used for the EDP Routing Scheme, which routes on the original graph.
def original(g):
    return g

# TODO: createOuterplanar() function that repeats the outerplanar subgraph creation steps, each time on the leftover edges, to create multiple graphs.

'''
ROUTING ALGORITHMS
'''
# source s
# destination d
# link failure set fails 
# g is an outerplanar embedding, when graph is outerplanar. Otherwise the original graph.
# returns result as boolean (True means routing failed, False means it succeeded) and hop count as int. 
# Switches and detour_edges are not relevant to our algorithm, but will be returned empty to fit in the benchmark_template structure.
def RightHandRule(s, d, fails, g):
    # Assert that graph is outerplanar
    if not is_outerplanar(g):
        print("Embedding is not outerplanar")
        return False, hops, 0, []
    
    print(f"s {s}, d {d}, failures {fails}")
    if DEBUG:
        print(f"Input graph: {g.edges}")
    # convert types so it works with template
    s = int(s) 
    d = int(d)
    source = s # save source node, because it has more specific forwarding rules
    p = source # previous node, to keep track of the incoming link
    hops = 0
    path = [source]

    if source not in g.nodes:
        print("Failure: invalid or isolated source")
        return True, hops, 0, []

    # generate clockwise list of neighbors of source. To be used later for specific source forwarding rules. Skip neighbors if edge is in fails
    neighbor = g.neighbors_cw_order(source)
    source_cw_neighbors = [] # list holding neighbors
    for n in neighbor:
        if((s, n) not in fails and (n, s) not in fails):
            source_cw_neighbors.append(n)
    if DEBUG:
        print(f"Starting node neighbors {source_cw_neighbors}")

    # if all outgoing links from source are failures, source is disconnected from destination and it returns a failure.
    if not source_cw_neighbors:
        print("Failure: source is disconnected from graph")
        return True, hops, 0, []
    
    # if source is destination, we are done
    if(s==d):
        return False, hops, 0, []
    while (s != d):
        if DEBUG:
            print(f"s {s} p {p}")
        # specific source forwarding rules
        if (s == source):
            # if previous node is also source, it means the routing algorithm just started. Route to first ccw neighbor that is not in fails
            if (p == source):
                if DEBUG:
                    print("p == source identified")
                source_ccw_neighbors = list(reversed(source_cw_neighbors))
                i = 0
                s = source_ccw_neighbors[i]
                while ((p,s) in fails or (s,p) in fails):
                    s = source_ccw_neighbors[(i+1)%len(source_ccw_neighbors)]
                hops += 1
                path.append(s)
                if DEBUG:
                    print(f"Starting: {path}")
                continue
            # if previous node is first in the cw_neighbors (last in ccw), it means we looped back to source without finding the destination. 
            elif (p == source_cw_neighbors[0]):
                # Otherwise, destination was not found.
                print(f"Failure path: {path}")
                print("Failure: d was not found")
                return True, hops, 0, []
            
        # route to the first ccw edge from the incoming link that is not in fails
        ccw_node = g.edges[s,p]['ccw']
        if ((s,ccw_node) in fails or (ccw_node,s) in fails):
            neighbor = g.neighbors_cw_order(s)
            cw_neighbors = []
            for n in neighbor:
                cw_neighbors.append(n)
            ccw_neighbors = list(reversed(cw_neighbors))
            while ((s,ccw_node) in fails or (ccw_node,s) in fails):
                if DEBUG:
                    print(f"{s} to {ccw_node} in fails")
                ccw_node = ccw_neighbors[(ccw_neighbors.index(ccw_node)+1) % len(ccw_neighbors)]
                if DEBUG:
                    print(f"try {s} to {ccw_node}") 
        p = s
        s = ccw_node
        hops += 1
        path.append(s)
        if DEBUG:
            print(path)

    # Success
    print(f"Success path: {path}")
    return False, hops, 0, []

# routing algorithm that generates edge disjoint paths between source and target and sorts them based on their length
# if a link on any path fails, packets bounce back to the source and another path is used for routing, until the destination is reached
# source s
# destination d
# link failure set fails 
# graph g
# returns result as boolean (True means routing failed, False means it succeeded) and hop count as int. 
# Switches and detour_edges are not relevant to our algorithm, but will be returned empty to fit in the benchmark_template structure.
def EdgeDisjointPaths(s, d, fails, g):
    # generate edge disjoint paths and sort them ascending
    paths = list(nx.edge_disjoint_paths(g, s, d))
    paths.sort(key=len)
    if DEBUG:
        print(f"EDPs: {paths}")

    hops = 0
    # attempt to route through each path until one is successful
    for path in paths:
        if DEBUG:
           print(f"Considering path: {path}")
        back = 0
        p = s
        for n in path:
            # ignore the first node, since it is the source
            if (n == s):
               continue
            # if the path encounters a failure, try the next path
            if ((p,n) in fails or (n,p) in fails):
               # add the amount of hops it would take to go back to the source
               hops += back
               break
            hops += 1
            back += 1
            # success if destination is found
            if (n==d):
                return False, hops, 0, []
            # otherwise current node becomes previous node
            p = n
    # else failure, if destination is not found
    return True, hops, 0, []

'''
HELPER METHODS METHOD
'''

# Step 1.
# Starting with a spanning graph with no edges, repeatedly (as long as possible) find a triangle T whose vertices are in different components of subgraph,
# and add the edges of T to the subgraph
def find_triangular_cactus(graph):
    # Create an empty spanning graph G[E]
    subgraph = nx.Graph()
    subgraph.add_nodes_from(graph.nodes)

    # Continue finding triangles as long as possible
    T = [cycle for cycle in nx.cycle_basis(graph) if len(cycle) == 3]
    for triangle in T:
        # Check if the triangle's vertices are in different components of G[E]
        components = list(nx.connected_components(subgraph))
        component_map = {node: i for i, comp in enumerate(components) for node in comp}

        # Count how many unique components the triangle's vertices belong to
        unique_components = len(set(component_map.get(node, -1) for node in triangle))

        # If all vertices are in different components (or not in any component yet), add the triangle
        if unique_components == 3:
            subgraph.add_edges_from([(triangle[i], triangle[(i + 1) % 3]) for i in range(3)])

    return subgraph

# Step 2.
# Starting with a triangular cactus_graph, repeatedly (as long as possible) find a square S in graph whose vertices are all in different components of 
# the cactus_graph, and add the edges of S to the cactus_graph. Graph is the original graph and cactus_graph is the result of find_triangular_cactus
def find_square_cactus(graph, cactus_graph):
    # Copy the cactus_graph
    subgraph = nx.Graph()
    subgraph.add_nodes_from(cactus_graph)
    subgraph.add_edges_from(cactus_graph.edges)

    # Continue finding squares as long as possible
    S = [cycle for cycle in nx.cycle_basis(graph) if len(cycle) == 4]
    for square in S:
        # Check if the square's vertices are in different components of G[E]
        components = list(nx.connected_components(subgraph))
        component_map = {node: i for i, comp in enumerate(components) for node in comp}

        # Count how many unique components the square's vertices belong to
        unique_components = len(set(component_map.get(node, -1) for node in square))

        # If all vertices are in different components (or not in any component yet), add the square
        if unique_components == 4:
            subgraph.add_edges_from([(square[i], square[(i + 1) % 4]) for i in range(4)])

    return subgraph

# Step 3.
# Repeatedly (as long as possible) find an edge e in graph whose endpoints are in different components of cactus_graph, and add e to cactus_graph. 
def find_edges(graph, cactus_graph):
    # Copy the triangular square cactus subgraph
    subgraph = nx.Graph()
    subgraph.add_nodes_from(cactus_graph)
    subgraph.add_edges_from(cactus_graph.edges)

    # For each edge
    e = graph.edges
    for edge in e:
        # Check if the edge's vertices are in different components of G[E]
        components = list(nx.connected_components(subgraph))
        component_map = {node: i for i, comp in enumerate(components) for node in comp}

        # Count how many unique components the edge's vertices belong to
        unique_components = len(set(component_map.get(node, -1) for node in edge))

        # If all vertices are in different components (or not in any component yet), add the edge
        if unique_components == 2:
            subgraph.add_edge(*edge)

    return subgraph

# Step 4. My own addition to augment the outerplanar subgraph.
# Repeatedly (as long as possible) find an edge e in graph and add it to the subgraph only if the end result stays outerplanar. 
def greedy_edges(graph, outerplanar):
    # Copy the outerplanar subgraph
    subgraph = nx.Graph()
    subgraph.add_nodes_from(outerplanar)
    subgraph.add_edges_from(outerplanar.edges)

    # Greedily add edges as long as the result stays outerplanar
    e = graph.edges
    for u,v in e:
        if not subgraph.has_edge(u, v):
            subgraph.add_edge(u, v)
            if is_outerplanar(subgraph):
                continue
            subgraph.remove_edge(u, v)
    return subgraph

# Augment outerplanar graph to biconnected outerplanar graph using BICONNECT(G) from Augmenting Outerplanar Graphs (1996) by Goos Kant
# graph is an outerplanar graph
# returns a biconnected outerplanar graph
def make_biconnected(graph):
    graph = graph.copy()

    # Find cut vertices
    cut_vertices = list(nx.articulation_points(graph))

    for v in cut_vertices:
        neighbors = list(graph.neighbors(v))

        # Partition neighbors into blocks based on biconnected components
        subgraph = nx.Graph(graph)
        subgraph.remove_node(v)
        blocks = []
        for component in nx.connected_components(subgraph):
            block_neighbors = [n for n in neighbors if n in component]
            if block_neighbors:
                blocks.append(block_neighbors)

        # Add edges between blocks to make the graph biconnected
        for j in range(len(blocks) - 1):
            u = blocks[j][-1]
            w = blocks[j + 1][0]
            if not graph.has_edge(u, w):
                graph.add_edge(u, w)

                # If the edge breaks outerplanarity, revert and try other combinations
                if not is_outerplanar(graph):
                    graph.remove_edge(u, w)
                    for u_alt in blocks[j]:
                        for w_alt in blocks[j + 1]:
                            if not graph.has_edge(u_alt, w_alt):
                                graph.add_edge(u_alt, w_alt)
                                if is_outerplanar(graph):
                                    break
                                graph.remove_edge(u_alt, w_alt)

    return graph

# function which returns a list of nodes (sorted in order of appearance), which form the outer face of the graph
# graph is a biconnected outerplanar graph
def generate_outer_face(graph):
    print("calculating simple_cycles")
    cycles = list(nx.simple_cycles(graph))
    print("finished calculating simple_cycles")
    #maxLength = max( len(l) for l in cycles )
    #path = list(l for l in cycles if len(l) == maxLength)
    path = max(cycles, key=len)
    if DEBUG:
        print(f"largest cycle: {path}")
    return path

# return planar embedding of an outerplanar graph using an ordered node list.
# graph is an outerplanar graph
# path is a list of nodes (sorted in order of appearance), which form the outer face of the graph
def convert_to_outerplanar_embedding(graph, path):
    # Step 1: Verify that the graph is outerplanar
    if not is_outerplanar(graph):
        raise ValueError("The graph is not outerplanar.")

    # Step 2: Create a PlanarEmbedding object
    planar_embedding = nx.PlanarEmbedding()
    
    if DEBUG:
        print(f"Convert to outerplanar embedding input graph edges: {graph.edges}")

    # Step 3: Assign the cyclic order of edges around each node based on positions
    for node in graph.nodes:
        if(node in nx.isolates(graph)):
            continue
        neighbors = list(graph.neighbors(node))
        if DEBUG:
            print(f"Neighbors of {node}: {neighbors}")
        neighbors.append(node)
        
        # Sort neighbors based on their order in the outer face, starting from the first node cw of the current node
        neighbors = sorted(neighbors, key=lambda n: path.index(n))
        neighbors = neighbors[neighbors.index(node):] + neighbors[:neighbors.index(node)]
        neighbors.remove(node)

        # Add edges in sorted order to the embedding
        neighbors = list(reversed(neighbors))
        planar_embedding.add_half_edge(node, neighbors[0])
        for i in range(1, len(neighbors)):
            u = neighbors[i-1]
            v = neighbors[i] 
            planar_embedding.add_half_edge(node, v, ccw = u)
    return planar_embedding

# function to test if a graph is outerplanar. Based on the fact that a graph G is outerplanar iff K+G (a new vertex K is joined with all vertices of G) is planar.
def is_outerplanar(graph):
    # Copy graph
    G = nx.Graph()
    G.add_nodes_from(graph)
    G.add_edges_from(graph.edges)

    # Create a new vertex and join it to all vertices of G
    G.add_node("k")
    for n in G:
        if(n!="k"): 
            G.add_edge("k",n)

    return nx.is_planar(G)

#Test the Outerplanar Routing Scheme using a Zoo Topology graph. 
#Some suitable graphs with 20<|V|<50 and not already outerplanar: Geant2009, Renater2010, SwitchL3

if DEBUG:
    file_path = "./benchmark_graphs/Surfnet.graphml" 
    G = nx.Graph(nx.read_graphml(file_path))

    # Check if graph is planar. If yes, if it is outerplanar
    print(f"Is G outerplanar? {is_outerplanar(G)}")

    embedding = createOuterplanar(G)
    print(f"Is output outerplanar: {is_outerplanar(embedding)}")

    # Simulate experiment
    s = 26
    d = 0
    fails = [(3, 49), (16, 17), (21, 28), (42, 46), (34, 35), (30, 31), (11, 23), (47, 48), (14, 45), (24, 27)]

    # Original graph without the edges
    g = G.copy(as_view=False)
    g = nx.convert_node_labels_to_integers(g)
    g.remove_edges_from(nx.selfloop_edges(g))
    G = g.to_undirected()
    G.remove_edges_from(fails)

    # Draw visualizations
    plt.subplot(212)
    nx.draw(g, with_labels=True, node_color='lightblue', edge_color='black')
    plt.title('Original Graph')

    plt.subplot(221)
    nx.draw_planar(embedding, with_labels=True, node_color='lightgreen', edge_color='red')
    plt.title('Outerplanar subgraph')

    plt.subplot(222)
    nx.draw(G, with_labels=True, node_color='lightgreen', edge_color='red')
    plt.title('Link failures removed')
    
    plt.show()
    
    # Test the RightHandRule
    result = RightHandRule(s,d,fails,embedding)
    print(f"RightHandRule: {result[0]}, Hops: {result[1]}")

    # Test EdgeDisjointPaths
    result2 = EdgeDisjointPaths(s,d,fails,embedding)
    print(f"EDG: {result2[0]}, Hops: {result2[1]}")

    # Test for connectivity
    dist = nx.shortest_path_length(G, target=d)
    print(f"shortest paths: {dist}")
    print(f"Is s in dist? {s in dist}")