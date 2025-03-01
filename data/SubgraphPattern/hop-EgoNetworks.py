import random
import sys
from collections import defaultdict

def read_edgelist(filename):
    """ Reads an edge list file and returns an adjacency list. """
    graph = defaultdict(set)

    with open(filename, 'r') as file:
        for line in file:
            node1, node2 = map(int, line.strip().split())
            graph[node1].add(node2)
            graph[node2].add(node1)  # Since it's an undirected graph

    return graph

def hop_network(graph, start_node, max_hops=10):
    ego_net = defaultdict(set)
    return_nets = []
    visited_nodes = set([start_node])
    num_edges = 0
    hops = 0
    distance = {start_node: 0}
    cur_dist = 0

    queue = [start_node]
    while queue and hops < max_hops:
        node = queue.pop(0)
        if distance[node] != cur_dist:
            return_nets.append([dict(ego_net), num_edges])
            cur_dist = distance[node]

        if distance[node] >= max_hops:
            break

        for neighbor in graph[node]:
            if neighbor not in visited_nodes:
                visited_nodes.add(neighbor)
                queue.append(neighbor)
                distance[neighbor] = distance[node]+1

            # Add edge to ego network
            if neighbor not in ego_net[node]:
                num_edges += 1
            ego_net[node].add(neighbor)
            ego_net[neighbor].add(node)



    return return_nets

def select_start_nodes(graph, num_nodes=10):
    """ Selects 10 random starting nodes from the graph. """
    return random.sample(list(graph.keys()), min(num_nodes, len(graph)))



if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <edge_list_file> <dataset_name>")
        sys.exit(1)

    filename = sys.argv[1]
    graph = read_edgelist(filename)

    start_nodes = select_start_nodes(graph, num_nodes=2)
    print("Selected start nodes:", start_nodes)


    for i, start_node in enumerate(start_nodes):
        hop_networks = hop_network(graph, start_node, max_hops=1000)
        for j, [ego_net, num_edges] in enumerate(hop_networks):
            if len(ego_net.keys()) < 200:
                continue
            if j > 0:
                if abs(len(ego_net.keys()) - len(hop_networks[j-1][0].keys())) < 50:
                    continue
            print(f"\nhop Network {i+1} (Starting Node: {start_node}), n={len(ego_net.keys())}, m={num_edges}:")
            with open(str(sys.argv[2])+"_hop_"+str(start_node)+"_"+str(len(ego_net.keys()))+"_"+str(num_edges)+".edgelist", 'w') as pattern_file:
                lines = [f"{node} {neigh}" for node, neighbors in sorted(ego_net.items()) for neigh in neighbors if node < neigh]
                pattern_file.write("\n".join(lines))
