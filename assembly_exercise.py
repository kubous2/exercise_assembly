# ========== Exercise 4 - functions copies
def in_out_degree(graph):
    in_degree = {vertex: 0 for vertex in graph.keys()}
    out_degree = {vertex: len(neighbors) for vertex, neighbors in graph.items()}

    for vertex, neighbors in graph.items():
        for neighbor in neighbors:
            if neighbor not in in_degree:
                in_degree[neighbor] = 0
            in_degree[neighbor] += 1

    return in_degree, out_degree


def dfs(graph, start):
    """
    Depth-first search WITHOUT recursion.
    Uses an explicit stack (a list).
    """
    visited = []
    visited_set = set()
    stack = [start]

    while stack:
        vertex = stack.pop()

        if vertex not in visited_set:
            visited.append(vertex)
            visited_set.add(vertex)

            # To visit neighbors in alphabetical order
            for neighbor in sorted(graph.get(vertex, []), reverse=True):
                if neighbor not in visited_set:
                    stack.append(neighbor)

    return visited


def connected_components(graph):
    components = []
    visited_all = set()

    for vertex in graph:
        if vertex not in visited_all:
            component_list = dfs(graph, vertex)
            component_set = set(component_list)
            components.append(component_set)
            visited_all.update(component_set)

    return components


def is_weakly_connected(graph):
    if not graph:
        return True

    shadow = {vertex: [] for vertex in graph}
    for u, neighbors in graph.items():
        for v in neighbors:
            if v not in shadow:
                shadow[v] = []
            shadow[u].append(v)
            shadow[v].append(u)

    components = connected_components(shadow)
    return len(components) <= 1


def degrees(graph):
    # Helper for undirected branch of eulerian_type
    return {vertex: len(neighbors) for vertex, neighbors in graph.items()}


def eulerian_type(graph, directed=False):
    if not directed:
        components = connected_components(graph)
        degs = degrees(graph)
        active_components = []

        for comp in components:
            has_edges = False
            for v in comp:
                if degs[v] > 0:
                    has_edges = True
                    break
            if has_edges:
                active_components.append(comp)

        if len(active_components) > 1:
            return "none"

        odd_count = sum(1 for degree in degs.values() if degree % 2 != 0)

        if odd_count == 0:
            return "circuit"
        elif odd_count == 2:
            return "path"
        else:
            return "none"

    else:
        shadow = {vertex: [] for vertex in graph}
        for u, neighbors in graph.items():
            for v in neighbors:
                if v not in shadow:
                    shadow[v] = []
                shadow[u].append(v)
                shadow[v].append(u)

        weak_components = connected_components(shadow)
        in_deg, out_deg = in_out_degree(graph)

        active_weak_components = []
        for comp in weak_components:
            has_edges = False
            for v in comp:
                if in_deg.get(v, 0) + out_deg.get(v, 0) > 0:
                    has_edges = True
                    break
            if has_edges:
                active_weak_components.append(comp)

        if len(active_weak_components) > 1:
            return "none"

        start_nodes, end_nodes = 0, 0
        all_vertices = set(graph.keys()).union(in_deg.keys()).union(out_deg.keys())

        for vertex in all_vertices:
            excess = out_deg.get(vertex, 0) - in_deg.get(vertex, 0)
            if excess == 1:
                start_nodes += 1
            elif excess == -1:
                end_nodes += 1
            elif excess != 0:
                return "none"

        if start_nodes == 0 and end_nodes == 0:
            return "circuit"
        elif start_nodes == 1 and end_nodes == 1:
            return "path"
        else:
            return "none"



# =================TASK 2 - DE BRUIJN GRAPH=============

def kmer_composition(text, k):
    """
    Returns every substring of length k of text, in the order they occur.
    """
    kmers = []
    for i in range(len(text) - k + 1):
        kmers.append(text[i:i + k])
    return kmers


def debruijn_graph(kmers):
    """
    Builds a directed de Bruijn graph from a list of k-mers.
    Every k-mer is one edge from its prefix to its suffix.
    """
    graph = {}

    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]

        # Ensure both prefix and suffix are initialized in the dictionary
        if prefix not in graph:
            graph[prefix] = []
        if suffix not in graph:
            graph[suffix] = []

        # Add the directed edge
        graph[prefix].append(suffix)

    return graph


def spell_path(nodes):
    """
    Spells the string from a walk in a de Bruijn graph.
    Starts with the whole first node, then adds the last letter of every later node.
    """
    if not nodes:
        return ""

    result = nodes[0]
    for node in nodes[1:]:
        result += node[-1]

    return result


# ==========TASK 3 - HIERHOLZER'S ALGORITHM =======================
def hierholzer(graph, start, verbose=False):
    """
    Finds an Eulerian path using Hierholzer's algorithm with strict rules:
    - Always picks the alphabetically first unused edge.
    - Splices new walks at the first node in the current walk that has unused edges.
    """
    # Create a tracking dictionary for unused edges.
    # Sorting in reverse alphabetical order allows us to use .pop() to efficiently
    # retrieve and remove the alphabetically FIRST edge.

    unused = {}
    for node, neighbors in graph.items():
        # sorting backwards
        sorted_neighbors = sorted(neighbors, reverse=True)
        unused[node] = sorted_neighbors

    def get_walk(curr):
        w = [curr]
        while unused.get(curr):
            nxt = unused[curr].pop()
            w.append(nxt)
            curr = nxt
        return w

    # Step 1: Initial walk
    walk = get_walk(start)
    if verbose:
        print(f"walk 1 (from the start node): {walk}")

    # Steps 2-4: Repeat until all edges are used
    while any(unused.values()):
        # Step 3: Find the first node in the current walk that has unused edges
        restart_idx = -1
        restart_node = None
        for i, node in enumerate(walk):
            if unused.get(node):
                restart_idx = i
                restart_node = node
                break

        # Walk from the restart node
        new_walk = get_walk(restart_node)
        if verbose:
            print(f"Restart at: {restart_node}")
            print(f"Walk: {new_walk}")

        # Step 4: Splice the new walk into the exact occurrence of the restart node
        walk = walk[:restart_idx] + new_walk + walk[restart_idx + 1:]
        if verbose:
            print(f"Walk after splicing: {walk}")

    return walk


def eulerian_path(graph, verbose=False):
    """
    Determines the correct start node and returns the Eulerian path.
    """
    in_deg, out_deg = in_out_degree(graph)

    start_node = None
    all_vertices = set(graph.keys()).union(in_deg.keys()).union(out_deg.keys())

    # Check for the node with out - in = 1
    for v in all_vertices:
        if out_deg.get(v, 0) - in_deg.get(v, 0) == 1:
            start_node = v
            break

    # If no such node exists, all nodes are balanced
    # use alphabetically first node
    if start_node is None:
        start_node = min(all_vertices)

    return hierholzer(graph, start_node, verbose)

if __name__ == "__main__":
    from collections import Counter

    # ==========================TASK 2  ==============================
    print ("============ TASK 2 ================")
    # Test kmer_composition
    text = "GATTACA"
    kmers_4 = kmer_composition(text, 4)
    print("kmer_composition:", kmers_4)  #['GATT', 'ATTA', 'TTAC', 'TACA']

    #Test debruijn_graph with Toy A 3-mers
    toy_A_reads = ["ACT", "AGA", "CAG", "CCA", "CTG", "GAA", "GAC", "GGT", "GTC", "TCC", "TGA"]
    db_graph = debruijn_graph(toy_A_reads)
    print(db_graph)

    edges_count = sum(len(neighbors) for neighbors in db_graph.values())
    print(f"Edges count (should be {len(toy_A_reads)}):", edges_count)

    print("Eulerian type for Toy A:", eulerian_type(db_graph, directed=True))  # 'path'

    # Test spell_path
    path = ["GAT", "ATT", "TTA"]
    print("spell_path:", spell_path(path))  #'GATTA'
    print("\n")

    # ===================== TASK 3 ===================
    print("============ TASK 3 ================")
    # Toy B (k-mers length 3)
    toy_B_reads = ["ACT", "AGC", "ATT", "CCG", "CGA", "CGT", "CTG", "GAC", "GAG", "GAT", "GTC", "TCC", "TCG", "TGA",
                   "TTC"]
    toy_B_graph = debruijn_graph(toy_B_reads)

    print("Testing eulerian_path on Toy B (verbose=True):")
    final_walk = eulerian_path(toy_B_graph, verbose=True)

    # K-mer composition and Eulere path
    print("\nSpelled genome:", spell_path(final_walk))


