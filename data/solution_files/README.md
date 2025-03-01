.gedsol solution files have the following syntax:
lines with * are comments
<n/e> (depeding on node or edge operation) <s/d/i> (depending on deletion insertion substitution) <node_1> (<node_2> only when substitution) 


to convert the raw scip solution files to .gedsol place "convert_solution.py" into the dir containing the _solution.txt files and run it with python3 convert_solution.py
