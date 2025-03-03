gi=0
for g in 88 289 148 139 10 580 277 247 427 230
do
hi=0
for h in 88 289 148 139 10 580 277 247 427 230 
do
if [ $hi -gt $gi ]
then
./gurobi_protein -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/Protein/Protein_31_40
./gurobi_protein -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/Protein/Protein_31_40
./gurobi_protein -f F1+ -g ${g} -h ${h}  -t 8 -l 600 -w F1Plus/Protein/Protein_31_40
./gurobi_protein -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/Protein/Protein_31_40
fi
((hi++))
done
((gi++))
done
