import os
import sys
import csv

def transform_csv_to_edgelist(input_filename, output_filename):
    paper_to_index = {}
    index_counter = 0
    edges = set()

    # script and datafiles should be in the same directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(script_dir, input_filename)
    output_file = os.path.join(script_dir, output_filename)

    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        return

    with open(input_file, 'r', newline='') as file:
        reader = csv.reader(file)

        next(reader, None) # skip header
        for row in reader:
            if len(row) != 4:
                print(f"Skipping malformed line: {row}")
                continue

            try:
                source_id = row[1].strip()
                target_id = row[2].strip()

                if source_id not in paper_to_index:
                    paper_to_index[source_id] = index_counter
                    index_counter += 1
                if target_id not in paper_to_index:
                    paper_to_index[target_id] = index_counter
                    index_counter += 1


                # does not add parallel edges
                node1 = paper_to_index[source_id]
                node2 = paper_to_index[target_id]
                edge = (min(node1, node2), max(node1, node2))
                edges.add(edge)

            except Exception as e:
                print(f"Error processing row {row}: {e}")
                continue

    # Write edges to output file in "Node1 Node2" format
    with open(output_file, 'w') as out_file:
        for node1, node2 in sorted(edges):
            out_file.write(f"{node1} {node2}\n")

    print(f"Edgelist saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python transformCsvToEdgelist.py <input_file.csv> <output_file.txt>")
        sys.exit(1)

    input_filename = sys.argv[1]  # Input CSV file
    output_filename = sys.argv[2]  # Output edge list file

    transform_csv_to_edgelist(input_filename, output_filename)
