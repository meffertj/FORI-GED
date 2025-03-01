import os
import sys

def transform_tab_to_edgelist(input_file, output_file):
    paper_to_index = {}
    index_counter = 0
    edges = set()

    # script should be in the same directory as the data files
    script_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(script_dir, input_filename)
    output_file = os.path.join(script_dir, output_filename)

    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        return

    with open(input_file, 'r') as file:
        next(file)
        next(file) # first two lines in tab format are comments

        for line in file:
            parts = line.strip().split()
            if len(parts) < 3:
                continue  # Skip lines not fitting the format

            paper1_id = parts[1].split(":")[1]
            paper2_id = parts[3].split(":")[1]

            if paper1_id not in paper_to_index: # paperIDs become zero-based nodeIDs
                paper_to_index[paper1_id] = index_counter
                index_counter += 1
            if paper2_id not in paper_to_index:
                paper_to_index[paper2_id] = index_counter
                index_counter += 1

            # edges are stored from lower nodeID to higher nodeID, no parallel edges
            node1 = paper_to_index[paper1_id]
            node2 = paper_to_index[paper2_id]
            edge = (min(node1, node2), max(node1, node2))
            edges.add(edge)

    # write to edgelist format
    with open(output_file, 'w') as out_file:
        for node1, node2 in sorted(edges):
            out_file.write(f"{node1} {node2}\n")

    print(f"Edgelist saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python transformTabToEdgeslist.py <input_file.tab> <output_file.txt>")
        sys.exit(1)

    input_filename = sys.argv[1]  # .tab file
    output_filename = sys.argv[2]  # edgelist filename

    transform_tab_to_edgelist(input_filename, output_filename)