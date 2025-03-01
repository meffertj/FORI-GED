import re
import os
import sys

def parse_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            line = line.strip()
            if "OPTIMAL" in line:
                outfile.write("* OPTIMAL\n")
            elif "TIMELIMIT" in line:
                outfile.write("* SUBOPTIMAL\n")
            else:
                # Match the structure for node and edge lines
                match = re.match(r'(node|edge)_(sub|del|ins)\[(\d+)(?:\]\[(\d+))?\]', line)
                if match:
                    element_type, action, first_num, second_num = match.groups()
                    e_or_n = 'n' if element_type == 'node' else 'e'
                    s_d_i = 's' if action == 'sub' else ('d' if action == 'del' else 'i')
                    if second_num:
                        # For sub, there are two node ids
                        outfile.write(f"{e_or_n} {s_d_i} {first_num} {second_num}\n")
                    else:
                        # For del/ins, only one node id
                        outfile.write(f"{e_or_n} {s_d_i} {first_num}\n")


def rename_file(input_filename, data_name):
    match = re.match(r"scip"+data_name+"_(.*?)_(\d+)_(\d+)_.*_solution\.txt", input_filename)
    if not match:
        raise ValueError("Filename format is incorrect")

    affine, num1, num2 = match.groups()
    base_name = f"{data_name}_{num1}_{num2}_{affine}"

    # different seeds could return different solution files
    integer = 1
    while os.path.exists(f"{base_name}_{integer}.gedsol"):
        integer += 1

    return f"{base_name}_{integer}.gedsol"
    
    
def rename_files_in_directory():
	for filename in os.listdir():
		new_filename = filename
		if 'affineF2' in filename:
			new_filename = filename.replace('affineF2', 'F2+')
		elif 'affineF1' in filename:
			new_filename = filename.replace('affineF1', 'F1+')
		rename_stub(filename, new_filename)
            
            
def rename_stub(filename, new_filename):
    if filename == new_filename:
    	return
    else:
    	old_path = os.path.join(os.getcwd(), filename)
    	new_path = os.path.join(os.getcwd(), new_filename)
    	os.rename(old_path, new_path)


# def filter_files_by_muta(directory, muta90vldbj):
# 	valid_pairs = set()
# 	for i in range(len(muta90vldbj)):
# 		for j in range(i, len(muta90vldbj)):
# 			num1, num2 = muta90vldbj[i], muta90vldbj[j]
# 			valid_pairs.add(f"{num1}_{num2}")
# 			valid_pairs.add(f"{num2}_{num1}")
#
# 	for filename in os.listdir(directory):
# 		if filename.endswith("_solution.txt"):
# 			# Check if file is a muta90 solution or an aids solution with wrong filename
# 			if not any(pair in filename for pair in valid_pairs):
# 				file_path = os.path.join(directory, filename)
# 				os.remove(file_path)
# 				print(f"Deleted: {file_path}")
#
#
# muta90vldbj = ['2490', '643', '267', '292', '800', '4084', '3970', '4204', '3691', '3755']
# filter_files_by_muta(os.getcwd(), muta90vldbj)

print("Usage: python3 convert_solution.py <dataset name> <path/to/solutionFiles/dir> \n[dataset name: aids, muta] as in the solution filename \"scipaids_..\" \n[path: if none is passed will try to convert _solution.txt files in current working directory]\n")

data_name = ""
dir_path = os.getcwd()
if len(sys.argv) < 2:
    print("No dataset name passed, exiting...")
    exit(1)
elif len(sys.argv) < 3:
    data_name = sys.argv[1]
    print("Trying to convert _solution.txt files in current working directory")
else:
    data_name = sys.argv[1]
    dir_path = sys.argv[2]
    print("Trying to convert _solution.txt files in ", dir_path)

os.chdir(dir_path)
rename_files_in_directory()

for file in os.listdir():
	if file.endswith("_solution.txt"):
		parse_file(file, rename_file(file, data_name))
