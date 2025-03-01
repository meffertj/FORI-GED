import random
import sys
from collections import defaultdict
import xml.etree.ElementTree as ET
import os
import re

def read_edgelist(filename):
    """ Reads an edge list file and returns an adjacency list. """
    graph = defaultdict(set)

    with open(filename, 'r') as file:
        for line in file:
            node1, node2 = map(int, line.strip().split())
            graph[node1].add(node2)
            graph[node2].add(node1)  # Since it's an undirected graph

    return graph



def graph_to_gxl(graph, output_file):
    """ Converts a graph into a GXL file with random node and edge labels. """
    gxl = ET.Element("gxl")
    graph_elem = ET.SubElement(gxl, "graph", id="ego_network", edgeids="true", edgemode="undirected")

    node_mapping = {node: idx+1 for idx, node in enumerate(graph.keys())}
    node_labels = {node: random.choice(["C", "O", "H", "N"]) for node in graph.keys()}

    # Create nodes
    for node, label in node_labels.items():
        mapped_id = node_mapping[node]
        node_elem = ET.SubElement(graph_elem, "node", id=str(mapped_id))
        attr_elem = ET.SubElement(node_elem, "attr", name="chem")
        ET.SubElement(attr_elem, "string").text = label

    # Create edges
    for node, neighbors in graph.items():
        for neighbor in neighbors:
            if node < neighbor:  # Avoid duplicate edges
                edge_elem = ET.SubElement(graph_elem, "edge", attrib={"from": str(node_mapping[node]), "to": str(node_mapping[neighbor])})
                attr_elem = ET.SubElement(edge_elem, "attr", name="valence")
                ET.SubElement(attr_elem, "int").text = str(random.choice([1, 2]))  # Random edge label

    # Write to file with proper formatting
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("<?xml version=\"1.0\"?>\n")
        f.write("<gxl>\n")
        f.write("<graph id=\"ego_network\" edgeids=\"true\" edgemode=\"undirected\">\n")

        for node, label in node_labels.items():
            mapped_id = node_mapping[node]
            f.write(f"<node id=\"{mapped_id}\"><attr name=\"chem\"><string>{label}</string></attr></node>\n")

        for node, neighbors in graph.items():
            for neighbor in neighbors:
                if node < neighbor:  # Avoid duplicate edges
                    label = random.choice([1, 2])
                    f.write(f"<edge from=\"{node_mapping[node]}\" to=\"{node_mapping[neighbor]}\"><attr name=\"valence\"><int>{label}</int></attr></edge>\n")

        f.write("</graph>\n")
        f.write("</gxl>\n")


if __name__ == "__main__":
    """ converts all edge list files from the largeGraphs directory to randomly labeled gxl graphs with the same label alphabet as mutagenicity """
    pattern = re.compile(r"(.*)_(\d+)_(\d+)_(\d+)\.edgelist")
    for file in os.listdir("./largeGraphs/"):
        match = pattern.match(os.path.basename(file))
        if match:
            graph = read_edgelist("./largeGraphs/"+file)
            edges = 0
            for i, [node, neighbors] in enumerate(graph.items()):
                for neigh in neighbors:
                    if node < neigh:
                        edges += 1
            prefix, num1, num2, num3 = match.groups()
            num1, num2, num3 = int(num1), int(num2), int(num3)

            new_path = file
            if (num2, num3) != (len(graph.keys()), edges):
                new_name = f"{prefix}_{num1}_{len(graph.keys())}_{edges}.edgelist"
                new_path = os.path.join(os.path.dirname(file), new_name)
                os.rename("./largeGraphs/"+file, "./largeGraphs/"+new_path)
                print(f"Renamed: {file} -> {new_path}")

            graph_to_gxl(graph, "./largeGraphs/"+new_path.replace(".edgelist", ".gxl"))
        elif file in ["Pubmed-Diabetes_19717_44327.edgelist", "cora_edges_2708_5278.edgelist"]:
            graph_to_gxl(graph, "./largeGraphs/"+file.replace(".edgelist", ".gxl"))

