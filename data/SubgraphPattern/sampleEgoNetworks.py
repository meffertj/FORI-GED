import random
import sys
from collections import defaultdict
import copy

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
        # if distance[node] != cur_dist:
        #     return_nets.append([copy.deepcopy(ego_net), num_edges])
        #     cur_dist = distance[node]

        if distance[node] >= max_hops:
            return_nets.append([copy.deepcopy(ego_net), num_edges])
            cur_dist = distance[node]
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

def ego_network(graph, start_node, max_nodes=10):
    ego_net = defaultdict(set)
    visited_nodes = set([start_node])
    num_edges = 0

    queue = [start_node]
    while queue and len(visited_nodes) < max_nodes:
        node = queue.pop(0)

        for neighbor in graph[node]:
            if len(visited_nodes) >= max_nodes:
                break

            if neighbor not in visited_nodes:
                visited_nodes.add(neighbor)
                queue.append(neighbor)

            if neighbor not in ego_net[node]:
                num_edges += 1
            ego_net[node].add(neighbor)
            ego_net[neighbor].add(node)


    return ego_net, num_edges



def select_start_nodes(graph, num_nodes=10):
    """ Selects 10 random starting nodes from the graph. """
    return random.sample(list(graph.keys()), min(num_nodes, len(graph)))



# if __name__ == "__main__":
#     if len(sys.argv) != 3:
#         print("Usage: python script.py <edge_list_file> <dataset_name>")
#         sys.exit(1)
#
#     filename = sys.argv[1]
#     graph = read_edgelist(filename)
#
#     start_nodes = select_start_nodes(graph, num_nodes=5)
#     print("Selected start nodes:", start_nodes)
#
#     for n in [5,10,15,20]:
#         for i, start_node in enumerate(start_nodes):
#             ego_net = hop_network(graph, start_node, max_hops=n)
#             if len(ego_net.keys()) != n:
#                 change = random.sample(list(graph.keys()), 1)
#                 for j, entry in enumerate(start_nodes):
#                     if entry == start_node:
#                         start_nodes[j] = change[0]
#                         break
#                 ego_net, m = ego_network(graph, change[0], max_nodes=n)
#                 print(f"\n############################### swapped (Starting Node: {start_node}), n={n}, m={m}:")
#             with open(str(sys.argv[2])+"_"+str(start_node)+"_"+str(n)+"_"+str(m)+".edgelist", 'w') as pattern_file:
#                 lines = [f"{node} {neigh}" for node, neighbors in sorted(ego_net.items()) for neigh in neighbors if node < neigh]
#                 pattern_file.write("\n".join(lines))

if __name__ == "__main__":
    

    filename = "largeGraphs/Pubmed-Diabetes_19717_44327.edgelist"
    graph = read_edgelist(filename)

    start_nodes = select_start_nodes(graph, num_nodes=2)
    print("Selected start nodes:", start_nodes)


    for i, start_node in enumerate(start_nodes):
        hop_networks = hop_network(graph, start_node, max_hops=5)
        for j, [ego_net, num_edges] in enumerate(hop_networks):
            #if 30 < len(ego_net.keys()): continue
            if len(ego_net.keys()) < 5: continue
            if j > 0:
                if abs(len(ego_net.keys()) - len(hop_networks[j-1][0].keys())) < 5:
                    continue
            print(f"\nhop Network {i+1} (Starting Node: {start_node}), n={len(ego_net.keys())}, m={num_edges}:")
            with open("large_pubmed_hop_"+str(start_node)+"_"+str(len(ego_net.keys()))+"_"+str(num_edges)+".edgelist", 'w') as pattern_file:
                lines = [f"{node} {neigh}" for node, neighbors in sorted(ego_net.items()) for neigh in neighbors if node < neigh]
                pattern_file.write("\n".join(lines))