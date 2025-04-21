import networkx as nx
from networkx.algorithms import approximation as ap
import matplotlib.pyplot as plt
DEBUG = True


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
    #nx.draw_planar(biconnected, with_labels=True, node_color='lightgreen', edge_color='red')
    #plt.show()
    
    # Step 6
    #path = generate_outer_face(biconnected)
    

    #ChatGPT Experimentations
    #G5, outer_face = augment_to_maximal_outerplanar(biconnected)
    #outer_face = outerplanar_embedding_from_dual(G5)
    #G_augmented, embedding_aug, outer_face = get_embedding_after_triangulation(biconnected)
    #return convert_to_outerplanar_embedding_from_augmented(G4, embedding_aug, set(G4.edges()))
    
    # Step 7
    return convert_to_outerplanar_embedding(G4, outer_face)

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

# Experimentations for finding an algorithm that computers the outer face of an outerplanar subgraph with less time complexity than nx.simple_cycles
'''
# Experiment: Build dual tree of maximal outerplanar graph. Traverse the dual to order the faces. Reconstruct outer face based on the order from the last step.
# Results are incomplete. Errors in the outerplanar building stage.

# function to build a dual tree from a maximal outerplanar graph.
def build_dual_tree(G):
    _, embedding = nx.check_planarity(G)
    visited = set()
    face_map = {}
    face_id = 0
    face_nodes = []

    # Extract faces and map shared edges
    for u in embedding:
        for v in embedding[u]:
            if (u, v) not in visited:
                face = tuple(embedding.traverse_face(u, v))
                face_nodes.append(face)
                face_map[face_id] = face
                for i in range(len(face)):
                    a, b = face[i], face[(i + 1) % len(face)]
                    visited.add((a, b))
                face_id += 1

    # Identify outer face (largest)
    outer_face_id = max(face_map, key=lambda fid: len(face_map[fid]))

    # Build dual graph (connect faces via shared edges)
    dual = nx.Graph()
    for i in face_map:
        if i == outer_face_id:
            continue  # skip outer face
        dual.add_node(i)

    # Build adjacency
    edge_to_faces = {}
    for fid, face in face_map.items():
        if fid == outer_face_id:
            continue
        for i in range(len(face)):
            a, b = sorted((face[i], face[(i + 1) % len(face)]))
            edge_to_faces.setdefault((a, b), []).append(fid)

    for face_list in edge_to_faces.values():
        if len(face_list) == 2:
            f1, f2 = face_list
            dual.add_edge(f1, f2)

    return dual, face_map, outer_face_id

# function to build the outerplanar embedding using the dual of a maximal outerplanar graph.
def outerplanar_embedding_from_dual(G):
    dual, face_map, outer_face_id = build_dual_tree(G)

    # Start from any leaf node in the dual tree
    start = [node for node in dual.nodes if dual.degree[node] == 1][0]
    visited_faces = set()
    stack = [(start, None)]
    face_order = []

    while stack:
        curr, parent = stack.pop()
        if curr in visited_faces:
            continue
        visited_faces.add(curr)
        face_order.append(curr)
        for neighbor in dual.neighbors(curr):
            if neighbor != parent:
                stack.append((neighbor, curr))

    # Reconstruct vertex ordering from face sequence
    outer_cycle = []
    added = set()
    for fid in face_order:
        face = face_map[fid]
        for v in face:
            if v not in added:
                outer_cycle.append(v)
                added.add(v)

    return outer_cycle

# Experiment 2: Based on the paper "Finding Hamiltonian cycles in certain planar graphs" by Robert J. Cimikowski from 1990. Only works on inner triangulations: 2-connected 
# planar graphs, where every interior face is a triangle (maximal outerplanar graps fall into this category).

# Step 1: Assign nodes to levels using BFS
def assign_levels(graph):
    levels = {}
    start_node = list(graph.nodes)[0]
    queue = [(start_node, 0)]
    visited = set()
    
    while queue:
        node, level = queue.pop(0)
        if node not in visited:
            visited.add(node)
            levels[node] = level
            for neighbor in graph.neighbors(node):
                if neighbor not in visited:
                    queue.append((neighbor, level + 1))
    
    return levels

# Step 2: Find ramps, which are edges that unify adjacent levels
def find_ramps(graph, levels):
    ramps = []
    for u, v in graph.edges():
        if abs(levels[u] - levels[v]) == 1:
            ramps.append((u, v))
    return ramps

# Step 3: Construct the Hamiltonian Cycle by joining all levels with ramps.
def construct_hamiltonian_cycle(graph):
    levels = assign_levels(graph)
    ramps = find_ramps(graph, levels)
    
    if not ramps:
        return None
    
    path = []
    level_nodes = sorted(levels.keys(), key=lambda x: levels[x])
    
    for i in range(len(level_nodes) - 1):
        path.append(level_nodes[i])
    
    return path
'''
    
#Test the Outerplanar Routing Scheme using a Zoo Topology graph. 
#Some suitable graphs with 20<|V|<50 and not already outerplanar: Geant2009, Renater2010, SwitchL3, HostwayInternational

if DEBUG:
    file_path = "./benchmark_graphs/Renater2010.graphml" 
    g = nx.Graph(nx.read_graphml(file_path))
    
    # Prepare graph
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

    # Find an outerplanar subgraph
    G1 = find_triangular_cactus(G)
    G2 = find_square_cactus(G, G1)
    G3 = find_edges(G, G2)
    G4 = greedy_edges(G, G3)

    # Augment to biconnected outerplanar
    biconnected = make_biconnected(G4)

    # Augment to maximal outerplanar from biconnected
    maximal, alpha = nx.complete_to_chordal_graph(biconnected)
    
    is_planar, embedding = nx.check_planarity(maximal)

    # Experiments with TSP. Does not return the correct outer face.
    #path = ap.traveling_salesman_problem(biconnected, weight='weight', nodes=set(biconnected.nodes), cycle=True)
    #print(path)
    #unique_path = []
    #for a in path:
    #    if a not in unique_path:
    #        unique_path.append(a)
    #print(unique_path)
    #embedding = convert_to_outerplanar_embedding(G4, unique_path)
    #path2 = generate_outer_face(biconnected)
    #print(path2)
    
    #is_planar, embedding2 = nx.check_planarity(maximal)

    plt.subplot(212)
    nx.draw_circular(embedding, with_labels=True)
    plt.show()

    #embedding = createOuterplanar(G)
    #print(f"Is output outerplanar: {is_outerplanar(embedding)}")
    
    #print(f"Is output a correct planar embedding: {embedding.check_structure()}")

    #nx.draw_planar(embedding, with_labels=True)
    #plt.show()

    '''
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
    print(f"Is s in dist? {s in dist}")'
    '''

