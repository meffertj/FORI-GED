import xml.etree.ElementTree as ET
import os

def parse_gxl(file_path):
    tree = ET.parse(file_path)
    root = tree.getroot()
    graph = root.find("graph")

    node_label_map = {}
    label_counter = 0
    node_map = {}
    edge_label_map = {}
    edge_label_counter = 0

    nodes = []
    edges = []

    # Process nodes
    for i, node in enumerate(graph.findall("node")):
        node_id = node.get("id")
        label = node.find("attr/string").text

        if label not in node_label_map:
            node_label_map[label] = label_counter
            label_counter += 1

        node_map[node_id] = i  # Assign zero-based index to node
        nodes.append(f"{i} {node_label_map[label]}")

    # Process edges
    for edge in graph.findall("edge"):
        from_id = node_map[edge.get("from")]
        to_id = node_map[edge.get("to")]
        edge_label = edge.find("attr/int").text

        if edge_label not in edge_label_map:
            edge_label_map[edge_label] = edge_label_counter
            edge_label_counter += 1

        edges.append(f"{from_id} {to_id} {edge_label_map[edge_label]}")

    return nodes, edges


if __name__ == "__main__":
    for file in os.listdir():
        if file.endswith(".gxl"):
            nodes, edges = parse_gxl(file)
            for node in nodes:
                print(node)
            for edge in edges:
                print(edge)