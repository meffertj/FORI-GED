gi=0
for g in 395 167 524 74 442 174 309 563 195 43
do
hi=0
for h in 395 167 524 74 442 174 309 563 195 43 
do
if [ $hi -gt $gi ]
then
./gurobi_protein -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/Protein/Protein_41_50
./gurobi_protein -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/Protein/Protein_41_50
./gurobi_protein -f F1+ -g ${g} -h ${h}  -t 8 -l 600 -w F1Plus/Protein/Protein_41_50
./gurobi_protein -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/Protein/Protein_41_50
fi
((hi++))
done
((gi++))
done
